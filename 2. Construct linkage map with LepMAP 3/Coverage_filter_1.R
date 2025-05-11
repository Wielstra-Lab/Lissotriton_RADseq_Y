
args_in <- commandArgs(trailingOnly = TRUE)

input_file <- args_in[1]
output_file <-  args_in[2]
min_mean_cov <- as.numeric(args_in[3])
mother_ID <- args_in[4]

if(length(args_in) < 3) {

  mother_ID <- NA 
  
} else mother_ID <- args_in[4]

input_table <- read.table(input_file)

N_samples <- ncol(input_table)

input_table[,"mean_cov"] <- rowMeans(input_table)

input_table[,"keep"] <- TRUE

print(mother_ID)

print(input_table[,mother_ID])

for(i in 1:nrow(input_table)){
  
  if(as.numeric(input_table[i,"mean_cov"]) < min_mean_cov) {
    input_table[i,"keep"] <- FALSE
    print(paste("Ditching", rownames(input_table)[i], "due to low cov:", input_table[i,"mean_cov"]))
  }
  
  if( is.na(mother_ID) == FALSE) {
    if( input_table[i,mother_ID] > 0) {
      input_table[i,"keep"] <- FALSE
      print(paste("Ditching", rownames(input_table)[i], "due to reads in mother:", input_table[i,mother_ID]))
    }
  }
}

filtered_table_1 <- subset(input_table, keep == TRUE)
filterd_table_2 <- filtered_table_1[,1:N_samples]

write.table(filterd_table_2, output_file, quote = FALSE, sep = "\t")  
  