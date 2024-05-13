
### Script to filter BLAST results from sorted BLAST outputs
### Produces table of best hits, where that hit is 5 orders of magnitude more specific than the next best hit
### (or where alll similar scoring hits are tightly clusterd in the same region of genome as the best hit ~10 Mbp)

# Get the input and output file names:

args_in <- commandArgs(trailingOnly = TRUE)

input_table_file <- args_in[1]
output_file <- args_in[2]

# Read the input BLAST table, add a column to keep track of the hits to keep

master_table <- read.table(input_table_file)

master_table[,"Keep"] <- FALSE

# Initialise variables

gene_id <- 0
gene_top_score <- 0
gene_top_row <- 0
gene_chr <- 0
gene_loc <- 0
accepting <- FALSE

# Main loop: two branches, one for the best hit for each query (which will be encountered first, as the blast is sorted), one for all secondary hits
    
for(i in 1:nrow(master_table)){
    
    #Branch for best hits, used if the gene_id from the new row does not match the current gene of interest
    if (master_table[i,1] != gene_id) {
    
        print(master_table[i,1])
    
        #If the current gene of interest is still a candidate (after we've processed all its hits) confirm it for keeping
        if (accepting == TRUE) {
        
            master_table[gene_top_row,"Keep"] <- TRUE
            print(paste("Keeping (only copy/all entries co-located)", gene_id))
            }
    
        #Now switch to the new row and update all variables
        gene_id <- master_table[i,1]
        gene_top_score <- master_table[i,11]
        gene_top_row <- i
        gene_chr <- master_table[i,2]
        gene_loc <- master_table[i,9]
        accepting <- FALSE
        
        #If the new gene's best hit is significant enough award it "candidate" status 
        if (gene_top_score < 1.0e-25) {accepting <- TRUE}
        
    }
    
    #Branch for supplementary hits, skip this if the best hit for gene is no longer a candidate for keeping (if accepting is false)   
    else if (accepting == TRUE) {
    
        #If the supplementary hit is more than 5 orders of magnitude less significant than the best hit, set the best hit to KEEP (and close candidacy) 
        if (gene_top_score * 100000 < master_table[i,11]) {
        
            accepting <- FALSE
            master_table[gene_top_row,"Keep"] <- TRUE
            print(paste("Keeping (highest value by far)", gene_id))
            }
            
        #If the supplementary hit is not on the same chromosome, close the candidacy (without setting the best hit to KEEP)
        else if (master_table[i,2] != gene_chr) {accepting <- FALSE}
        
        #If the supplementary hit is not close to the best hit locus, close the candidacy (without setting the best hit to KEEP) 
        else if (abs(master_table[i,9] - gene_loc) > 10000000) {accepting <- FALSE}
    }
}

# Write the results to the given output file

final_table <- subset(master_table, Keep==TRUE)
write.table(final_table, output_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
            
        
            
    
    