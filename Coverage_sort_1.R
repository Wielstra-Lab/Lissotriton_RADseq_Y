####
#Basic script to produce a table of the coverage of a group of RADseqed samples for certain loci

args_in <- commandArgs(trailingOnly = TRUE)

popmap_file <- args_in[1]
coverage_directory <- args_in[2]
BLAST_file <- args_in[3]
output_file <- args_in[4]

#No need to to load lines beyond the highest numbered loci we use

loci_table <- read.table(BLAST_file)
loci_list <- loci_table[,2]
loci_num <- max(loci_list)

print(loci_list)

files_list <- list.files(coverage_directory, full.names = TRUE)

sample_list_table <- read.table(popmap_file)
sample_list <- sample_list_table[,1]

print(sample_list)

#create empty coverage table of appropriate length

coverage_table <- data.frame(matrix(nrow = loci_num, ncol = 0))

#find the appropriate coverage file for each sample, open it as its own table and then copy the coverage values as a new column in the main table

for(sample in sample_list) {
  print(sample)
  sample_index <- grep(sample, files_list)
  sample_file <- files_list[sample_index]
  sample_table <- read.table(sample_file)
  #sample_table <- read.table(sample_file, nrows = loci_num)
  coverage_table[,sample] <- sample_table[1:loci_num,2]
} 

#filter the main coverage table for the loci we are interested in and then output it

filtered_cov_table <- coverage_table[loci_list,]

#set any NA to 0
filtered_cov_table[is.na(filtered_cov_table)] <- 0

rownames(filtered_cov_table) <- loci_table$V1

write.table(filtered_cov_table, output_file, quote = FALSE, sep = "\t")