
### Script for automatically trimming loosely bound outlying markers from linkage map, and re-ordering map by group length

args_in <- commandArgs(trailingOnly = TRUE)

Raw_map_file <- args_in[1]

output_file <- args_in[2]

if(length(args_in) < 3) {

  priority_group <- NA #if this is an integer, that group will always become group 1
  
} else priority_group <- as.integer(args_in[3])

if(length(args_in) < 4) {

  keep_weight <- 1 #default keep weight is 1
  
} else keep_weight <- as.numeric(args_in[4])

if(is.na(priority_group)){
  print("No group prioritised")
} else print(paste("Make raw group ", priority_group ," group 1 in final map"))

print(paste("Keep weight: ", keep_weight))

# This final argument allows for keeping the bp info (within the target sequences) of the SNPs in the linkage map

if(length(args_in) < 5) {
  keep_bp_info <- FALSE
  
} else if(args_in[5] == "keep") {

  keep_bp_info <- TRUE
  print("keep bp info")
  
} else keep_bp_info <- FALSE
  

print("")

# The number of markers to be cut is multiplied by the keep weight, and compared to the 
# Length in cm to be cut, if the length is longer, and at least 5, then the cut is made

Raw_map <- read.table(Raw_map_file)

Trimmed_map <- Raw_map[FALSE,]

nGroups <- max(Raw_map[,4])

Trimmed_group_lenghts <- data.frame(1:nGroups)
Trimmed_group_lenghts[,2] <- NA
Trimmed_group_lenghts[,3] <- 0

if(is.na(priority_group) == FALSE) {Trimmed_group_lenghts[priority_group,3] <- 1}

### Main trimming loop ###

for(Raw_group_id in 1:nGroups){
  
  # First: subset the group from the combined map

  Raw_group_map <- subset(Raw_map, V4 == Raw_group_id)
  
  Raw_group_marker_num <- nrow(Raw_group_map)
  Raw_group_length_cm <- Raw_group_map[Raw_group_marker_num,3]
  
  print(paste("markers in raw group ",Raw_group_id,": ",Raw_group_marker_num))
  print(paste("length of raw group ",Raw_group_id,": ",Raw_group_length_cm))
  
  # Set default cut points (that keep the entire group), in case nothing is found
  
  Keep_below <- 0
  Keep_above <- nrow(Raw_group_map) + 1
  
  upper_length_removed <- 0
  lower_length_removed <- 0
  
  # Find upper cut point
  
  for(i in 1:(nrow(Raw_group_map) - 1)){
    
    distance_to_next_marker <- Raw_group_map[i + 1,3] - Raw_group_map[i,3]
    
    if (distance_to_next_marker > i * keep_weight & distance_to_next_marker > 5) {
      Keep_below <- i
      upper_length_removed <- Raw_group_map[i + 1,3]
      }
  }
  
  print(paste("cutting ", Keep_below, " markers and ", upper_length_removed, "cM above ", Raw_group_map[(Keep_below + 1),1]))
  
  # Find lower cut point
  
  for(i in nrow(Raw_group_map):2){
    
    distance_to_previous_marker <- Raw_group_map[i, 3] - Raw_group_map[i - 1, 3] 
    
    if (distance_to_previous_marker > (nrow(Raw_group_map) - i) * keep_weight & distance_to_previous_marker > 5) {
      Keep_above <- i
      lower_length_removed <- Raw_group_length_cm - Raw_group_map[i - 1,3]
    }
  }
  
  print(paste("cutting ", (Raw_group_marker_num - Keep_above) + 1, " markers and ", lower_length_removed, "cM below ", Raw_group_map[(Keep_above - 1),1]))
  
  # Subset by cut points to create trimmed map, and then reset co-ordinates so first marker is at 0
  
  Trimmed_group_map <- Raw_group_map[(Keep_below + 1):(Keep_above - 1),]
  
  New_zero_value <- Trimmed_group_map[1,3]
  
  Trimmed_group_map[,3] <- Trimmed_group_map[,3] - New_zero_value
  
  Trimmed_group_marker_num <- nrow(Trimmed_group_map)
  
  print(paste("markers in trimmed group ",Raw_group_id,": ", Trimmed_group_marker_num))
  print(paste("length of trimmed group ",Raw_group_id,": ", Trimmed_group_map[Trimmed_group_marker_num,3]))
  
  # Add trimmed group to overall map and note trimmed group length on table
  
  Trimmed_group_lenghts[Raw_group_id,2] <- Trimmed_group_map[Trimmed_group_marker_num,3]
  
  Trimmed_map <- rbind(Trimmed_map,Trimmed_group_map)
  
  print("")
}

print(paste("Markers in complete final map: ", nrow(Trimmed_map)))
print(paste("Markers removed: ", nrow(Raw_map) - nrow(Trimmed_map)))

### Order by length ###

# sort table of lengths by trimmed group length (and priority)

Sorted_group_lengths <- Trimmed_group_lenghts[order(-Trimmed_group_lenghts$V3,-Trimmed_group_lenghts$V2),]

Reordered_map <- Raw_map[FALSE,]

# subset each group in length order from trimmed table, change group-ID in subset and add subset to re-ordered table

for(length_rank in 1:nrow(Sorted_group_lengths)){
  
  old_ID <- Sorted_group_lengths[length_rank,1]
  old_ID_table <- subset(Trimmed_map, V4 == old_ID)
  old_ID_table[,4] <- length_rank
  
  Reordered_map <- rbind(Reordered_map,old_ID_table)
}

#unless directed to keep delete the second column of the map (which states the bp of the marker where the SNP is) and write the output to the specified file

if(keep_bp_info == FALSE) {Reordered_map$V2 <- NULL}

write.table(Reordered_map, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, na = "")

