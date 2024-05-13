# Script designed to import coverage data for RAD markers into R and filter for possible male markers
# Script caculates median coverage for both sexes and removes any markers with less than 20 in males and more than 0 in females

### Read arguments ### 

args_in <- commandArgs(trailingOnly = TRUE)

popmap_file <- args_in[1]
coverage_directory <- args_in[2]
output_name <- args_in[3]

### Load lists of coverage files for samples ###

popmap_table <- read.table(popmap_file)
male_sample_table <- subset(popmap_table, V2 == "Male")
female_sample_table <- subset(popmap_table, V2 == "Female")

male_list <- male_sample_table$V1
female_list <- female_sample_table$V1
files_list <- list.files(coverage_directory, full.names = TRUE)

list_all <- c(male_list,female_list)

### note how many samples of each sex ###

num_males <- length(male_list)
num_females <- length(female_list)

print(c("number of males: ", num_males))
print(c("number of females: ", num_females))

### Get the required length for the main table (this is the same as the highest marker name in the dataset - easiest to get this with unix shell commands) ###

max_vals <- c()

for(sample in list_all) {
  sample_index <- grep(sample, files_list)
  sample_file <- files_list[sample_index]
  max_vals <- c(max_vals, as.integer(system(paste("tail -1", sample_file, "| cut -f 1"), intern = TRUE)))
}

loci_num <- max(max_vals)

print(paste("Highest numbered locus:", loci_num))

### create main table ### 

coverage_table <- data.frame(matrix(nrow = loci_num, ncol = 0))

### Load coverage data from male samples ###

for(sample in male_list) {
  sample_index <- grep(sample, files_list)
  sample_file <- files_list[sample_index]
  sample_table <- read.table(sample_file)
  coverage_table[,sample] <- sample_table[1:loci_num,2]
  print(c("imported ",sample))
} 

### Load coverage data from female samples ###

for(sample in female_list) {
  sample_index <- grep(sample, files_list)
  sample_file <- files_list[sample_index]
  sample_table <- read.table(sample_file)
  coverage_table[,sample] <- sample_table[1:loci_num,2]
  print(c("imported ",sample))
} 

print("finished importing samples")
print(c("number of samples imported: ", ncol(coverage_table)))

coverage_table[is.na(coverage_table)] <- 0 #change any missing data to 0 

print("calculating medians")

### split coverage data into male and female tables, and calculate a median for each row in each sex ###

male_table <- coverage_table[,1:num_males]
female_table <- coverage_table[,(num_males + 1):(num_males + num_females)]

male_table$median = apply(male_table, 1, median, na.rm=T)

print("finished calculating male medians")

female_table$median = apply(female_table, 1, median, na.rm=T)

print("finished calculating female medians")

### add the male and female median columns to the main table ###

coverage_table[,"male_median"] <- male_table[,num_males + 1]
coverage_table[,"female_median"] <- female_table[,num_females + 1]

### initial rough subset ###

print("subsetting table: step 1")

filtered_coverage_table <- subset.data.frame(coverage_table, male_median > 20 & female_median < 1)

print("number of remaining loci")
print(nrow(filtered_coverage_table))

### subset for candiate Y markers ###

print("subsetting table: step 2")

total_loci <- nrow(filtered_coverage_table)

fem_max <- num_females - 3
male_min <- num_males - 3

for(row in 1:total_loci) {
    filtered_coverage_table[row,"male_hits"] <- num_males - sum(filtered_coverage_table[row,1:num_males] == 0)
    filtered_coverage_table[row,"female_zeroes"] <- sum(filtered_coverage_table[row,(num_males + 1):(num_males + num_females)] == 0)
    }

markers_missing_females <- subset.data.frame(filtered_coverage_table, female_zeroes >= fem_max & male_hits >= male_min)

print("number of remaining loci")
print(nrow(markers_missing_females))

for(samp in 1:(num_males + num_females)) {
    samp_hits <- sum(markers_missing_females[,samp] != 0)
    print(colnames(markers_missing_females)[samp])
    print(samp_hits)
    }

### export ###

output_table_file_name <- paste(output_name, "_table.txt", sep = "")
output_table_list_name <- paste(output_name, "_markers.txt", sep = "")

ordered_table <- markers_missing_females[ with(markers_missing_females, order(female_zeroes, male_hits, male_median, decreasing = TRUE)),]

write.table(ordered_table, output_table_file_name, quote = FALSE, sep = "\t")

write(rownames(ordered_table), output_table_list_name)


