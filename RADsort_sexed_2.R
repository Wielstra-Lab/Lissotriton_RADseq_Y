#Script to filter candidate markers according to the number of paralogs found in a blast search

args_in <- commandArgs(trailingOnly = TRUE)

# Intially load the relavant files, the input table of makers, the blast results and also determine what the output should be called

marker_table_file <- args_in[1]
blast_results <- args_in[2]
output_code <- args_in[3]

marker_table <- read.table(marker_table_file)
blast_table <- read.table(blast_results)

long_output <- paste(output_code, "all.txt", sep = "_")
short_output <- paste(output_code, "top_10.txt", sep = "_")

# Now count the number of times each marker name occurs in the blast results

marker_table[,"Paralogs"] <- 0

for(i in 1:nrow(marker_table)) {
  marker <- row.names(marker_table)[i]
  x <- sum(blast_table$V1 == marker)
  marker_table[marker,"Paralogs"] <- x
}

# Rank the table by first absence in females, then fewest paralogs, then presence in males and finally median coverage in males

ordered_1 <- marker_table[ with(marker_table, order(female_zeroes,-Paralogs,male_hits,male_median,decreasing = TRUE)),]

# export both the full ranked table and a summary of the top 10

num_cols <- ncol(marker_table)

Top_10_summary <- ordered_1[1:10,(num_cols - 4):num_cols]

write.table(ordered_1, long_output, quote = FALSE, sep = "\t")

write.table(Top_10_summary, short_output, quote = FALSE, sep = "\t")