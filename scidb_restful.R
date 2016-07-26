# This demo walks through an entire usage workflow of SciDB's RESTful API (Shim)
# Specifically, the script demonstrates the following:
# 1. connecting to SciDB (server) from your client (R in this case) via Shim
# 2. executing any query and saving the result of that query in a binary file (on the server)
# 3. reading that binary file over HTTP into client (R memory)
# 4. also we show how to convert the binary payload into an R dataframe

# NOTES:
# - the code for this demo is adapted from the code and documentation of `query_to_df()` function 
#         of SciDBR (https://github.com/Paradigm4/SciDBR/blob/master/R/utility.R).
# - in #2 above, we run a query and ask SciDB to save the result into a binary files. One can also 
#         execute queries on SciDB that do not generate any results. For examples of such usage of Shim, 
#         refer the Shim documentation 
#         https://htmlpreview.github.io/?https://raw.github.com/Paradigm4/shim/master/wwwroot/help.html
# - #4 above (converting SciDB output to R dataframes) is just one way of interpreting the SciDB output. 
#         Other users might decide to use the binary data in a completely different way

### PREP ### (Run only once)
# # Connect to one R session and execute the following commands
# # The store(apply(..)) query creates a temporary array `temp_demo` for demo purposes
# host = "40.121.92.238"  #options("scidb.default_shim_host")[[1]]
# port = 8080             #options("scidb.default_shim_port")[[1]]
# library(scidb)
# scidbconnect(host, port)   # or just `scidbconnect()` if you are running this demo on the same machine as the Shim server
# iquery("store(apply(build(<val1:double>[i=0:4,5,0, j=0:3,2,0], random()), val2, i*10+j), temp_demo)")


rm(list=ls())
### A CONVENIENCE FUNCTION ###
# This is useful for converting SciDB queries to URL friendly strings

# Some versions of RCurl seem to contain a broken URLencode function, so a version is added here
# From documentation of URLencode (https://stat.ethz.ch/R-manual/R-devel/library/utils/html/URLencode.html): 
# "Characters in a URL other than the English alphanumeric characters and - _ . ~ should be encoded 
# as % plus a two-digit hexadecimal representation, and any single-byte character can be so encoded. 
# (Multi-byte characters are encoded byte-by-byte.) The standard refers to this as ‘percent-encoding’."
oldURLencode = function (URL, reserved = FALSE) 
{
  OK = paste0("[^-ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz0123456789$_.+!*'(),", 
              if (!reserved) 
                ";/?:@=&", "]")
  x = strsplit(URL, "")[[1L]]
  z = grep(OK, x)
  if (length(z)) {
    y = sapply(x[z], function(x) paste0("%", as.character(charToRaw(x)), 
                                        collapse = ""))
    x[z] = y
  }
  paste(x, collapse = "")
}

### THE DEMO ###
host = "40.121.92.238"  #options("scidb.default_shim_host")[[1]]
port = 8080             #options("scidb.default_shim_port")[[1]]
aflstr="between(temp_demo,0,0,2,2)"    # a SELECT operation
TYPES=c("double", "int64")
NULLABILITY=c(FALSE, FALSE)
attribs=c("val1", "val2")
dims=c("i","j")
DEBUG=TRUE



aflstr = sprintf("project(apply(%s,%s),%s,%s)", 
                 aflstr, 
                 paste(dims, dims, sep=",", collapse=","), 
                 paste(dims, collapse=",") , 
                 paste(attribs, collapse=",")
)

# Some basic string cleaning on the AFL string so that it can be passed to httr::GET e.g. SPACE is replaced with %20
aflstr = oldURLencode(aflstr)

# Formulate the shim url
shimurl=paste("http://", host, ":", port, sep = "")

# Connect to SciDB
# /new_session
id = httr::GET(paste(shimurl, "/new_session", sep=""))
session = gsub("([\r\n])", "", rawToChar(id$content)) # the gsub is added to remove some trailing characters if present 

# AFL query and save format
TYPES=c(rep("int64", length(dims)), TYPES)
NULLABILITY=c(rep(FALSE, length(dims)), NULLABILITY)
savestr=sprintf("(%s)", 
                paste(TYPES, 
                      ifelse(NULLABILITY,
                             "+null", ""), 
                      sep="", collapse = ",")
)

# Now formulate and execute the query
# /execute_query
query = sprintf("%s/execute_query?id=%s&query=%s&save=%s",
                shimurl, session, aflstr, savestr)
resp = httr::GET(query)

# View the status of the result
if(DEBUG) { cat("save status: "); cat(sprintf("%d\n", resp$status_code)) }

# Get the data returned by the previous save
# /read_bytes
resp = httr::GET(sprintf("%s/read_bytes?id=%s&n=0", shimurl, session))
if (DEBUG) {
  cat("read status: "); cat(sprintf("%d\n", resp$status_code))
  print(resp)
}

# Release the session
httr::GET(sprintf("%s/release_session?id=%s", shimurl, session))

#BEGIN: Copied from scidbr/R/internal.R >> scidb_unpack_to_dataframe()
len = length(resp$content)
p = 0 
ans = c()

# BEGIN: Hard coded for tests
buffer = 100000L
# END: HARD CODED

cnames = c(dims, attribs, "lines", "p")  # we are unpacking to a SciDB array, ignore dims
n = length(dims) + length(attribs)

while(p < len)
{
  dt2 = proc.time()
  tmp   = .Call("scidb_parse", as.integer(buffer), TYPES, NULLABILITY, resp$content, as.double(p), PACKAGE="scidb")
  names(tmp) = cnames
  lines = tmp[[n+1]]
  p_old = p
  p     = tmp[[n+2]]
  if(DEBUG) cat("  R buffer ", p, "/", len, " bytes parsing time", (proc.time() - dt2)[3], "\n")
  dt2 = proc.time()
  if(lines > 0)
  {
    if("binary" %in% TYPES)
    {
      if(DEBUG) cat("  R rbind/df assembly time", (proc.time() - dt2)[3], "\n")
      return(lapply(1:n, function(j) tmp[[j]][1:lines]))
    }
    len_out = length(tmp[[1]])
    if(lines < len_out) tmp = lapply(tmp[1:n], function(x) x[1:lines])
    # adaptively re-estimate a buffer size
    avg_bytes_per_line = ceiling((p - p_old) / lines)
    buffer = min(getOption("scidb.buffer_size"), ceiling(1.3 * (len - p) / avg_bytes_per_line)) # Engineering factors
    # Assemble the data frame
    if(is.null(ans)) ans = data.table::data.table(data.frame(tmp[1:n], stringsAsFactors=FALSE))
    else ans = rbind(ans, data.table::data.table(data.frame(tmp[1:n], stringsAsFactors=FALSE)))
  }
  if(DEBUG) cat("  R rbind/df assembly time", (proc.time() - dt2)[3], "\n")
}
ans = as.data.frame(ans)
print(ans)
#END: Copied from scidbr/R/internal.R >> scidb_unpack_to_dataframe()