# BEGIN_COPYRIGHT
# 
# Copyright © 2016 Paradigm4, Inc.
# This script is used in conjunction with the Community Edition of SciDB.
# SciDB is free software: you can redistribute it and/or modify it under the terms of the Affero General Public License, version 3, as published by the Free Software Foundation.
#
# END_COPYRIGHT

#Welcome! If you are new to RStudio, note you can run a line of code by clicking on it and then pressing CTRL+R or Command+R. 
#This is a set of some demo workflows for the 1000 Genomes dataset (http://www.1000genomes.org)
#!!! !!! Run these lines first to set up the environment and connect to the database:
library(scidb)
library('threejs')
library('ggplot2')
source('/home/scidb/ami_functions.R') #Will output "creating a generic function for 'image'... that is normal
scidbconnect("localhost")

#Here are some statements in a function, but do run them one after another to get a feel for how the SciDB-R package works:
demo = function()
{
  #Create a reference to an array called KG_VARIANT. 
  #This array stores all the variations in the 1000 Genomes dataset (chromosomes 7 and 8)
  #When we run this line, no data is moved - only a description of the array is retrieved and stored:
  KG_VARIANT    = scidb("KG_VARIANT")
  #Now examine the array schema. Here we chose chromosome/start/end as dimensions. There is another dimension
  #called alternate_id for multiple variants that occur on the same position, i.e. A->T or A->G:
  str(KG_VARIANT)
  #Check out some counts and stats:
  summarize(KG_VARIANT)
  #Finally, check out the first few entries:
  head(KG_VARIANT)
  
  #Now the GENOTYPE array had all the variant calls in the 1000 Genomes: all the variants in KG_VARIANT,
  #called with a true/false on both alleles across 2504 samples. This array is larger! 
  KG_GENOTYPE   = scidb("KG_GENOTYPE")
  str(KG_GENOTYPE)
  summarize(KG_GENOTYPE)
  head(KG_GENOTYPE, n=20)  
  #Be a little careful with larger arrays; some whole-array statistics may take some time ~ 10s of minutes.
  #If stuck, hitting the red "stop" button at the top of the console will properly cancel the query
  
  #Let's run our first query. Note: SciDB excels at slicing by dimension, so this kind of selection
  #is very efficient. Given the region 55Mbp - 56Mbp on Chromosome 7 (we store them 0-based), 
  #return all the variants that overlap this region:
  selection = subset(KG_VARIANT, chromosome_id==6 && end >= 55000000 && start <= 56000000)
  #How many variants in this region?
  count(selection)
  #Look at a few:
  head(selection)
  
  #Note: "selection" is actually a reference to a SciDB query! These stored queries are similar to 
  #array references. The query is not evaluated until the user asks for data (or part of it).
  #Check out the SciDB AFL query string like so:
  selection@name
  
  #Now let's make it interesting!
  #For all the non-reference variants calls:
  gt_selection = subset(KG_GENOTYPE, allele_1==TRUE || allele_2==TRUE)
  #Join them with the "selection" variants as above
  gt_selection = merge(gt_selection, selection)
  #Using another array to filter for South Asian population only
  KG_POPULATION = scidb("KG_POPULATION")
  str(KG_POPULATION)
  gt_selection = merge(gt_selection, subset(KG_POPULATION, population=="SAS"))
  #For those variant calls, return just these 5 fields:
  gt_selection = project(gt_selection, c("id","reference","alternate","allele_1","allele_2"))
  #And store the result into a new temporary array in SciDB:
  gt_selection = scidbeval(gt_selection, temp=TRUE)
  #Check it out:
  gt_selection@name
  count(gt_selection)
  head(gt_selection) 
  
  #Now do a count of non-reference calls, per variant:
  aggregated = grouped_aggregate(gt_selection, FUN="count(*)", by=list("start","reference", "alternate", "id"))
  #Sort and return the top 100 most frequent:
  aggregated = sort(aggregated, attributes="count", decreasing=TRUE)
  aggregated = subset(aggregated, n<100)
  top_100_SAS_vars = iquery(aggregated, return=TRUE)
  head(top_100_SAS_vars) 
  
  ###
  ### Here's another example of a PCA visualization workflow.
  
  #Earlier, we selected 1 million variants randomly and made this dense matrix:
  #1 million variants x 2504 samples:
  matrix = scidb("KG_VAR_CENTERED")
  str(matrix)
  head(matrix)
  #We can regrid large matrices down to display-able size for visualization:
  image(matrix)  
  #The code to make a matrix like this is in the PCA function below.
  
  #And we computed the sample-by-sample covariance matrix:
  covariance = scidb("KG_VAR_CV")
  str(covariance)
  image(covariance)
  
  #And we also computed the left singular vectors. Let's plot the first 3 in 3D:
  svded = scidb("KG_VAR_SVD")
  str(svded)
  #Download just the 3 vectors into R and make a matrix out of them:
  svd_top = df2xyvm(iqdf(subset(svded, i<=2), n=Inf))
  #Do kmeans clustering of these vectors in R now:
  clustering = kmeans(svd_top, 5, nstart=50)
  #Convert the kmeans cluster assignments to colors
  color=gsub("[0-9]","",palette()[clustering$cluster+1])
  #Plot. You can mouse over the plot and spin it around:
  scatterplot3js(svd_top, size=0.2, color=color, renderer="canvas")  
  #Upload the clusters into SciDB and join them against the population array, then download back:
  sample_clustering = merge(as.scidb(clustering$cluster), scidb("KG_POPULATION"), by=1)[]
  sample_clustering$super_population = factor(sample_clustering$population)
  sample_clustering$val = factor(color)
  #Now examine clusters versus actual populations - not 100% but pretty good:
  table(sample_clustering$val, sample_clustering$super_population)
}

#This concludes the quick intro. The workloads presented below get a little more advanced. If you found the demo
#interesting so far, we of P4 would love to hear more about your interest and your use cases! Chat with us via:
#lifesciences@paradigm4.com
#http://forum.paradigm4.com
#http://www.paradigm4.com/

#Read on for more workflows:
#1) compute the Hardy-Weinberg equilibrium for the variants in a particular location
#2) perform a range join of variants and genes: find the genes that overlap the most variants
#3) compute the PCA of samples using K randomly chosen variants (i.e. create the matrix we've plotted above)

#Create links to a few of the arrays:
KG_CHROMOSOME = scidb("KG_CHROMOSOME")
KG_VARIANT    = scidb("KG_VARIANT")
KG_GENOTYPE   = scidb("KG_GENOTYPE")
KG_SAMPLE     = scidb("KG_SAMPLE")
KG_POPULATION = scidb("KG_POPULATION")
GENE = scidb("GENE_37")

######
#Hardy-Weinberg Equilibrium Workflow:
#Given a genomic region, extract all the variants that overlap it. For those variants, compute the number of
#Homozygous-reference,
#Heterozygous and
#Homozygous-variant calls
#Then use the Hardy-Weinberg principle to estimate the expected number of heterozygous calls; compare with the actual.
#The variants out of equlibrium may be due to technical measurement problems, or, in fact, interesting evolutionary reasons.
#ETA: ~5 minutes for the default arguments
hardy_weinberg_equilibrium = function(chromosome=7, start_coord = 10000000, end_coord = 19999999)
{
  print("Computing initial counts")
  selection = transform(KG_GENOTYPE, aa="iif(allele_1 = TRUE and allele_2 = TRUE, 1, 0)",
                        Aa="iif((allele_1 = FALSE and allele_2 = TRUE) or (allele_1 = TRUE and allele_2 = FALSE), 1, 0)",
                        AA="iif(allele_1 = FALSE and allele_2 = FALSE, 1, 0)")
  chrm = chromosome-1;
  selection = subset(selection, chromosome_id==chrm && end>=start_coord  && start<=end_coord)
  counts = grouped_aggregate(selection, FUN="sum(AA) as AA, sum(Aa) as Aa, sum(aa) as aa", 
                             by=list("chromosome_id", "start", "end", "alternate_id"))
  counts = scidbeval(counts, temp=TRUE)
  print("Computed counts")
  counts = redimension(counts, sprintf("<aa:int64 null, Aa:int64 null, AA:int64 null>%s", scidb:::build_dim_schema(KG_VARIANT)))
  counts = scidbeval(counts, temp=TRUE)
  print("Redimensioned counts to match KG_VARIANT")
  stats = transform(counts, N    = "AA+Aa+aa")
  stats = transform(stats,  p    = "double(2*AA + Aa)/(2*N)",
                    q    = "double(2*aa + Aa)/(2*N)")
  stats = transform(stats,  E_AA = "N*p*p",
                    E_Aa = "2*N*p*q",
                    E_aa = "N*q*q")
  stats = transform(stats,  chi_2 = "pow(AA-E_AA,2)/E_AA + pow(Aa-E_Aa,2)/E_Aa + pow(aa-E_aa,2)/E_aa")
  stats = scidbeval(stats, temp=TRUE)
  print("Computed per-variant stats")
  
  #Plot a sampling of 10K variants
  download = transform(stats, r="random()")
  download = sort(download, attributes="r")
  download = head(download, n=10000)
  print(qplot(as(download$E_Aa,"vector"), as(download$Aa,"vector"), xlab="E_Aa", ylab="Aa"))
  
  #Examine some variants with the most deviation
  most_deviant = merge(KG_VARIANT, stats)
  most_deviant = sort(unpack(most_deviant), attributes="chi_2", decreasing=TRUE)
  most_deviant = project(most_deviant, c("chromosome_id", "start", "reference", "alternate",
                                          "aa", "Aa", "AA", "E_aa", "E_Aa", "E_AA", "chi_2"))
  subset(most_deviant, n<20)[]
}

######
#Range join workflow:
#Locate variants that overlap with genes, and count the number of overlapping variants for each gene.
#Output the top 20 genes with the most variants, plus plot distribution in first 200.
#This represents a range join, which has many different use cases in genomics. For example,
#we can use the same procedure to compare two sets of variants and reference ranges, i.e. 
#give me variants in A that are fully enclosed by reference ranges in B - to detect de novo variants.
#ETA: ~1 minute with default arguments
#Note: range_join is a pretty sophisticated routine in ~/scidb/ami_functions.R. We're considering dropping
#it into a C++ SciDB operator in the future.
rank_genes_by_variants = function(top_n = 200, overlaps = TRUE)
{
  res = range_join(KG_VARIANT, 
                   index_lookup(GENE, KG_CHROMOSOME, attr="chromosome", new_attr="chromosome_id"), 
                   left_attributes="qual", 
                   right_attributes="gene", 
                   overlaps = overlaps)
  res = grouped_aggregate(res, FUN="count(*) as num_variants", by="gene")
  res = sort(res, attributes="num_variants", decreasing=TRUE)
  res = subset(res, n<top_n)[]
  print(qplot(x=res$n, y=res$num_variants))
  head(res, n=20)
}

######
#PCA workflow:
#Select all variants whose af is greater than min_af and less than max_af; if variant_limit is
#specified, then select up to that many variants, chosen randomly. 
#<value> [sample_id, variant_number] where value is 2 if the variant is present in both alleles,
#1 if the variant is present in one allele or 0 otherwise. Then compute the left singular vectors
#of this matrix and store the result in a new SciDB array.
#ETA: ~25 minutes for the default arguments:
#~3 min to select and assemble the matrix, ~20 min to compute the SVD
#This is the community edition of SciDB, not using the optimized MKL BLAS. In the enterprise
#edition, the SVD phase uses optimized MKL BLAS and so completes ~5x faster!
#The EE truncated tsvd() can be used to return an approximation even sooner.
compute_pca = function(min_af = 0.1, max_af = 0.9, variant_limit=100000, result_name="KG_VAR_SVD_2")
{
  t1=proc.time();
  selected_variants = project(unpack(subset(KG_VARIANT, af>min_af && af<max_af)), c("chromosome_id", "start", "end", "alternate_id"))
  selected_variants = transform(selected_variants, randomizer="random()")
  selected_variants = sort(selected_variants, attributes="randomizer")
  selected_variants = rename_fields(selected_variants, dims="variant_number")
  selected_variants = subset(selected_variants, variant_number<variant_limit)
  selected_variants = scidbeval(selected_variants, temp=TRUE)
  
  num_variants = count(selected_variants)
  num_samples = count(KG_SAMPLE)
  print(sprintf("Found %i variants that fit the criteria; Extracting [%i x %i] matrix", 
                num_variants, num_variants, num_samples))
  print(proc.time()-t1);
  
  t1=proc.time();
  redim_guide = redimension(selected_variants, sprintf("<variant_number:int64> %s", scidb:::build_dim_schema(KG_VARIANT)))
  
  redim_matrix = subset(KG_GENOTYPE, allele_1 == TRUE || allele_2 == TRUE)
  redim_matrix = transform(redim_matrix, v="double(iif(allele_1 and allele_2, 2.0, 1.0))")
  redim_matrix = merge(redim_matrix, redim_guide)
  redim_matrix = redimension(redim_matrix, sprintf("<v:double> [sample_id = 0:%i,512,0, variant_number=0:%i,512,0]", 
                                                   num_samples-1, num_variants-1))
  redim_matrix = merge(redim_matrix, build(0, redim_matrix), merge=TRUE)
  redim_matrix = scidbeval(redim_matrix, temp=TRUE)
  print(proc.time()-t1);
  
  t1=proc.time();
  print("Centering matrix")
  centered = project(transform(merge(redim_matrix, aggregate(redim_matrix, FUN=mean, by="variant_number"), by="variant_number"), c="v-v_avg"), "c")
  centered = scidbeval(centered, temp=TRUE)
  print(proc.time()-t1);
  
  #Covariance if you like
  #Y = gemm(centered, transpose(centered))
  #CV = bind(Y, "cv", sprintf("gemm/%i", num_samples))
  #CV = project(CV, "cv")
  #CV= scidbeval(CV)
  #image(CV)
  
  t1=proc.time();
  print("Computing SVD")
  svded = gesvd(centered, 'left')
  svded = scidbeval(svded, name = result_name, gc=0)
  print(proc.time()-t1);
  
  print(sprintf("Result stored as %s", result_name))
  plot(x= subset(svded,i==0)[]$u, y=subset(svded, i==1)[]$u)
}
