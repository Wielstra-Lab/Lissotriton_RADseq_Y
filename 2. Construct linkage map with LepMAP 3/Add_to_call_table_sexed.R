
#script takes a table of coverage for a linkage map family and converts coverage/presence absence data, into psuedoSNPs in the parent.call format used by LepMAP 3
#Designed for candidate Y markers
#input file must be formated with sample names as header (in the same order as PED file/parent.call file) and marker designations as row_names
#markers added by this scripted can be noted by a prefix

args_in <- commandArgs(trailingOnly = TRUE)

input_file <- args_in[1]
output_file <-  args_in[2]
prefix <- args_in[3]

Y_marker_cov_table <- read.table(input_file)

XX_code <- "1.0 0 0 0 0 0 0 0 0 0"
XY_code <- "0 1.0 0 0 0 0 0 0 0 0"


Y_marker_format <- Y_marker_cov_table[,1:2]
colnames(Y_marker_format) <- c("CHR","POS")

Y_marker_calls <- cbind(Y_marker_format,Y_marker_cov_table)

for(i in 1:nrow(Y_marker_calls)){
  Y_marker_calls[i,1] <- paste(prefix,row.names(Y_marker_calls)[i],sep = "")
  Y_marker_calls[i,2] <- 100
  for(j in 3:ncol(Y_marker_calls)){
    if(as.numeric(Y_marker_calls[i,j]) == 0){
      Y_marker_calls[i,j] <- XX_code
    }
    else{Y_marker_calls[i,j] <- XY_code}
  }
}
write.table(Y_marker_calls,output_file,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE) 
