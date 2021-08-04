library(Matrix)
library(RCTD)
library(doParallel)
library(ggplot2)
library(xlsx)
library(tidyverse)
library(dplyr)
library(DEGLAM)
library(GSA)
# Data is obtained from the following article.
# https://science.sciencemag.org/content/362/6416/eaau5324/tab-figures-data
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248 merfish data
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113576 reference data
# normalization data was obtained by emailing the authors of the article
# Directory setup
pwd = getwd()
datadir <- paste0(pwd,'/data/moffitt','/')
resultsdir <- paste0(pwd,'/ResultsMerfish','/')

# Load in spatialRNA data and create spatialRNA object
# colnames are Cell_ID, Animal_ID, Animal_sex, Behavior, Bregma, Centroid_X, Centroid_Y, Neuron_cluster_ID, All_the_other_genes
# 161 genes in the columns here
merfish_data = read.csv(paste0(datadir,'Moffitt_and_Bambah-Mukku_et_al_merfish_all_cells.csv'))
merfish_barcodes = merfish_data[,1]
normalization_data = read.csv(paste0(datadir,'volume_and_batch_correction.csv'))
normalization_data = select(normalization_data, abs_volume, batch_correction)
rownames(normalization_data) = merfish_barcodes

# Create SpatialRNA object using coords, counts, and nUMI
coords = select(merfish_data,Cell_ID,Centroid_X,Centroid_Y)
rownames(coords) = coords[,1]; coords[,1] = NULL; colnames(coords) = c("xcoord","ycoord");

counts = select(merfish_data, -Cell_ID, -Animal_ID, -Animal_sex, -Behavior, -Bregma, -Centroid_X, -Centroid_Y, -Neuron_cluster_ID, -Cell_class)
counts = counts[1:140] # Throw out genes at and after Adcyap1
rownames(counts) = merfish_barcodes
counts[is.na(counts)] <- 0

# Filter coords, counts, and normalization_data so they only contain data where
# animal ID = 1 and bregma = .01.
# This represents data from one slice of tissue
# nUMI is calculated later from counts
relevant_barcodes = merfish_data$Cell_ID[merfish_data$Animal_ID==1 & merfish_data$Bregma==.01]
coords = filter(coords, rownames(coords) %in% relevant_barcodes)
counts = filter(counts, rownames(counts) %in% relevant_barcodes)
normalization_data = filter(normalization_data, rownames(normalization_data) %in% relevant_barcodes)

# Reverse the normalization done on expression values to get their raw integer counts.
# Counts is barcodes x gene, original counts were  (counts / volume) * batch_correction value.
counts = counts / normalization_data[,2] # reverse batch correction
counts = counts * normalization_data[,1] # reverse division by cell volume

# Verify most of the values in counts are now basically integers
close_integer = function(val){
	if(val==0){return(FALSE)} # don't care really about the 0's
	eps = 1e-2
	remainder = val %% 1
	if (remainder + eps > 1 | remainder - eps < 0)
		return(TRUE)
	print(val)
	return(FALSE)
}
non_zero_count = sum(counts > 0)
close_integers_matrix = apply(counts, MARGIN=c(1,2) ,FUN=close_integer)
total_close_integers = sum(close_integers_matrix)
close_integer_ratio = total_close_integers / non_zero_count # should be close to 1

# Make the counts values integers
counts = as.matrix(round(counts))
mode(counts) = "integer"
# counts = data.frame(counts) # If counts needs to be in dataframe form.

counts = t(counts)

nUMI <- colSums(counts)

# Clear massive variables not going to be used again. Anything >1GB
# rm(merfish_data) 

puck = SpatialRNA(coords,counts,nUMI) 

# Examine spatialRNA object (optional)
print(dim(counts)) # 140 genes x 6111 cells
hist(log(nUMI,2))
puck_barcodes <- colnames(puck@counts)
plot_puck_continuous(puck, puck_barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), 
                     title ='plot of nUMI') 

# Load in reference data and then create Reference Object
# 27998x31299 genes x barcodes
reference_matrix_data = readMM(paste0(datadir,'GSE113576_matrix.mtx'))
reference_barcodes_data = read.table(file=paste0(datadir,'GSE113576_barcodes.tsv'))
reference_genes_data = read.table(file=paste0(datadir,'GSE113576_genes.tsv'))
cell_type_data = read.xlsx(file=paste0(datadir,'aau5324_Moffitt_Table-S1.xlsx'),sheetIndex =1 ,header=TRUE)
colnames(cell_type_data) = c("barcodes","sex","replicate_number","cell_type","non_neuronal_cluster","neuronal_cluster")
cell_type_data = cell_type_data[-1,]; # Set actual column names since file is loaded oddly

# Create Reference Object using counts, and cell_types
counts = as.matrix(reference_matrix_data)
mode(counts) = "integer" # convert matrix data from numeric to integers
reference_barcodes = reference_barcodes_data[,1]
colnames(counts) = reference_barcodes
reference_genes = reference_genes_data[,2]
rownames(counts) = reference_genes

# Create nUMI before filtering out any data from counts
nUMI = colSums(counts)
names(nUMI) <- colnames(counts)

# There are duplicated genes in counts. Remove them.
duplicate_rows_boolean = duplicated(rownames(counts))
duplicate_row_indices = which(duplicate_rows_boolean)
duplicate_gene_names = unique(rownames(counts)[duplicate_row_indices])
# Can't convert to dataframe and filter out rows because conversion creates unique rownames.
counts = counts[!(rownames(counts) %in% duplicate_gene_names), ]

cell_types = cell_type_data$cell_type
names(cell_types) = cell_type_data$barcodes
cell_types <- as.factor(cell_types)

# Clear massive variables no longer being used, anything >1GB
# rm(reference_matrix_data)

# 27877 genes x 31299 cells
reference <- Reference(counts, cell_types, nUMI)

# Create and run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 8)
# Save/Load RCTD object
# saveRDS(myRCTD,file.path(resultsdir,'preRCTD.rds'))
myRCTD = readRDS(paste0(resultsdir,'preRCTD.rds')) # 53 DE genes found from 140

myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
# Save/Load postRCTD object
# saveRDS(myRCTD,file.path(resultsdir,'postRCTD.rds'))
myRCTD = readRDS(paste0(resultsdir,'postRCTD.rds'))

# Create explanatory.variable
cell_types = myRCTD@cell_type_info$info[[2]]
barcodes <- colnames(myRCTD@spatialRNA@counts)
explanatory.variable = c(rep(0,length(barcodes)))
names(explanatory.variable) = barcodes
x_coords = myRCTD@spatialRNA@coords[,1]
midline_x = (max(x_coords)+min(x_coords)) / 2
for(i in 1:length(barcodes)) {
  if(i%%50==0) { print(i) } # Track progress
  barcode = barcodes[i]
  curr_x = myRCTD@spatialRNA@coords[barcode,1]
  distance = abs(curr_x - midline_x)
  explanatory.variable[i]=distance
}
explanatory.variable = explanatory.variable / max(explanatory.variable)
hist(explanatory.variable) # Observe the distribution of explanatory.variable
plot_puck_continuous(myRCTD@spatialRNA, colnames(myRCTD@spatialRNA@counts), explanatory.variable, # ylimit = c(0,round(quantile(myRCTD@spatialRNA@nUMI,0.9))), 
                     title ='plot of explanatory variable') 


# run DE
myRCTD <- run.de.single(myRCTD, explanatory.variable) # outputs plots into pwd() + "/de_plots/"

# Save/Load and plot resulting RCTDde
saveRDS(myRCTD,file.path(resultsdir,'myRCTDde.rds'))
myRCTD = readRDS(paste0(resultsdir,'myRCTDde.rds'))
make_all_de_plots(myRCTD, resultsdir)
# Extra stuff
# Histogram of n.iter used by DEGLAM -- myRCTD@de_results$gene_fits$n.iter
hist(myRCTD@de_results$gene_fits$n.iter)
# Check for convergence -- look at myRCTD@de_results$gene_fits$con_val and myRCTD@de_results$gene_fits$precision_val
con_val = myRCTD@de_results$gene_fits$con_val
proportion_genes_converged = sum(con_val)/length(con_val)

# Precision vals should be low
precision_val = myRCTD@de_results$gene_fits$precision_val
hist(precision_val)

# Make sure the reject cells were discarded
dim(myRCTD@internal_vars_de$my_beta)
dim(myRCTD@spatialRNA@coords)
length(myRCTD@internal_vars_de$all_barc)

########################################################################
# Dylans Gene ontoogy code, adapted to analysis here
# Determine what cell type to look at
source("~/DEGLAM/R/de_utils.R", echo=TRUE)
cell_aggregate = names(aggregate_cell_types(myRCTD))
res_gene_list_cell_types = names(de_results$res_gene_list)
cell_type_options = intersect(cell_aggregate,res_gene_list_cell_types)
sorted_options = agg_types[order(-agg_types)][cell_type_options]
# Astrocytes              Endothelial               Excitatory Immature oligodendrocyte 
# 800                      321                      763                      186 
# Inhibitory   Mature oligodendrocyte 
# 580                      448 
# Note: abandoning gene enrichment analysis for merfish because 140 genes is too few and we find no pathways of interest to look at

de_results = myRCTD@de_results
# Get a list of the overexpressed and under expressed genes in CAF cells
over_genes <- tolower(rownames(de_results$res_gene_list$Inhibitory[de_results$res_gene_list$Inhibitory$log_fc > 0,]))
under_genes <- tolower(rownames(de_results$res_gene_list$Inhibitory[de_results$res_gene_list$Inhibitory$log_fc < 0,]))

# Load hallmark gene sets and change the genes to all lowercase so it's easier to work with.
gene_sets = GSA.read.gmt(file.path(datadir,'hallmark_genesets.gmt'))
gene_set_names = gene_sets$geneset.names
gene_set_descriptions = gene_sets$geneset.descriptions
gene_sets = gene_sets$genesets
names(gene_sets)=gene_set_names
gene_sets = lapply(gene_sets, tolower)
n_sets = length(gene_sets)

# Optional: Check for intersections between hallmark sets.
for(i in 1:n_sets) {
  for(j in 1:50){
    if(i!=j){
      print(intersect(gene_sets[[i]], gene_sets[[j]]))
    }
  }
}
# Conclusion, there's definitely intersections between hallmark sets

# Optional: How many of the gene set genes do we have spatial counts information on
genes_list <- tolower(rownames(myRCTD@originalSpatialRNA@counts))
coverage = c()
for(i in 1:n_sets){
  coverage[i] = length(intersect(gene_sets[[i]], genes_list))/length(gene_sets[[i]])
}
hist(coverage)

# Optional: Print the percent of genes in over and under genes that are in the hallmark gene sets
for(i in 1:n_sets) {
  print(paste0(i,"-----"))
  print(length(intersect(gene_sets[[i]], over_genes))/length(gene_sets[[i]]))
  print(length(intersect(gene_sets[[i]], under_genes))/length(gene_sets[[i]]))
}

# Use binomial tests to see if the ratio of genes over/under expressed in our data is significantly different
# when compared to the ratio of over/under expressed genes in the hallmark gene sets.
# If the ratio is significantly different, we can assume that gene sets pathway is significant somehow in
# how the cells are differentially expressed.
# Given a random gene in the differentially expressed genes, it should have a p_avg chance of being over expressed.
# So for n_o + n_u random genes, we expect p_avg amount of them to be overexpressed but we actually have n_o overexpressed.
# If we find too many more overexpressed it implies our overexpressed genes seem to fit this gene set well as in this gene set is enriched.
# two.sided makes the opposite true as well, so it also tests for underexpression
p_vals = numeric(n_sets)
p_avg = length(over_genes) / (length(over_genes) + length(under_genes))
n_o_vals = numeric(n_sets)
n_u_vals = numeric(n_sets)

max_n = 0; # Used later when generating a chart to know how many rows it needs to have, rows contain the over/under expressed genes
for(i in 1:n_sets) {
  print(i)
  n_o = length(intersect(gene_sets[[i]], over_genes))
  n_u = length(intersect(gene_sets[[i]], under_genes))
  n_o_vals[i] = n_o
  n_u_vals[i] = n_u
  total = n_o+n_u
  # arguments are heads, flips, rate. it tells you the chance that you got n_e when the fair rate was p_avg
  if(total>0)
    p_vals[i] <- binom.test(n_o, total, p_avg, alternative="two.sided")$p.value
  else
    p_vals[i] = 1.0
  max_n = max(max_n,n_o,n_u)
}

# Thorough analysis of pvalues taking indo account false discovery rate
# https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini%E2%80%93Hochberg_procedure
fdr = 0.1 
N_sig <- max(which(p_vals[order(p_vals)] < fdr * 1:n_sets / n_sets)) # which pvalues are less than .1 scaled from 0 t .1???
sig_list <- order(p_vals)[1:N_sig]
sig_pathways = gene_set_names[sig_list]
sig_list_df <- data.frame(sig_list, p_vals[sig_list], n_o_vals[sig_list], n_u_vals[sig_list])
colnames(sig_list_df) <- c('index', 'p_val', 'n_over', 'n_under')
sig_list_df$name <- gene_set_names[sig_list]
sig_list_df$desc <- gene_set_descriptions[sig_list]
gene_list_df <- data.frame()
# 2 columns per gene set
df <- matrix(nrow = max_n, ncol = 2*length(sig_list_df$name))
colnames(df) <- c(rbind(paste("UnderExpressed_" ,sig_list_df$name,sep=""),paste("OverExpressed_" ,sig_list_df$name,sep="")))
for(i in 1:N_sig) {
  gl <- intersect(gene_sets[[sig_list[i]]], under_genes)
  if(length(gl) > 0)
    df[1:length(gl),2*i-1] <- gl
  gl <- intersect(gene_sets[[sig_list[i]]], over_genes)
  if(length(gl) > 0)
    df[1:length(gl),2*i] <- gl # i+N_sig because the over expressed gene columns are shifted over by a constant from their under expressed column counterpart
}
# TODO: create a new folder for the gene enrichment analysis results
write.csv(df,file = file.path(resultsdir, 'gene_pathways.csv'))
write.csv(sig_list_df, file = file.path(resultsdir, 'sig_list_df.csv'))

# Looking at the list
pathways = substring(tolower(gene_set_names[sig_list]),10,)






