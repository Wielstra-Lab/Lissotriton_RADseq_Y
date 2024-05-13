

# first import the primary and secondary map

args_in <- commandArgs(trailingOnly = TRUE)

Primary_map_file <- args_in[1]
Secondary_map_file <- args_in[2]
output_file_name <- args_in[3]

Primary_map_table <- read.table(Primary_map_file)
Secondary_map_table <- read.table(Secondary_map_file)

colnames(Primary_map_table) <- c("marker","position","group")
colnames(Secondary_map_table) <- c("marker","position","group")

# make chart of group lengths of second map

print("make chart")

nGroups_sec <- max(Secondary_map_table$group)

Secondary_groups <- data.frame(1:nGroups_sec)
Secondary_groups[,"Length"] <- NA

for(group in 1:nGroups_sec) {
  
  group_end <- findInterval(group,Secondary_map_table$group)
  Secondary_groups[group,"Length"] <- as.numeric(Secondary_map_table[group_end,2])
}

# match the markers in each group

print("match")

Primary_map_table[,"secondary_group"] <- 0
Primary_map_table[,"secondary_position"] <- 0

Secondary_row_nums <- match(Primary_map_table$marker,Secondary_map_table$marker)

Primary_map_table[,"secondary_group"] <- Secondary_map_table[Secondary_row_nums, "group"]
Primary_map_table[,"secondary_position"] <- Secondary_map_table[Secondary_row_nums, "position"]

#print(Primary_map_table)

# evaluate matching groups and orrientation

print("evaluate")

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

Secondary_groups[,"Primary_group"] <- NA
Secondary_groups[,"Invert"] <- FALSE

for(sec_group in 1:nGroups_sec){
  
  print(paste("sec group:", sec_group))
  
  group_table <- subset(Primary_map_table, secondary_group == sec_group)
  
  print(nrow(group_table))
  
  Primary_group <- getmode(group_table$group)
  
  print(paste("prime group:", Primary_group))
  
  Secondary_groups[sec_group,"Primary_group"] <- Primary_group 
  
  filtered_group_table <- subset(group_table, group == Primary_group)
  
  coefficient <- cor.test(filtered_group_table$position, filtered_group_table$secondary_position)
  if(coefficient$estimate < 0){ Secondary_groups[group,"Invert"] <- TRUE}
}

print(Secondary_groups)

# rebuild secondary map to match primary flipping groups as required

print("build and flip")

nGroups_prime <- max(Primary_map_table$group)

Export_table <- Secondary_map_table[0,]

print(nrow(Export_table))

for(prime_group in 1:nGroups_prime) {
  
  sec_group <- match(prime_group, Secondary_groups$Primary_group)
  
  print(sec_group)
  
  sec_group_table <- subset(Secondary_map_table, group == sec_group)
  
  print(nrow(sec_group_table))
  
  if (Secondary_groups[sec_group,"Invert"] == TRUE) {
    group_length <- max(sec_group_table$position)
    rev_group <- sec_group_table
    rev_group[,"position"] <- sec_group_table[,"position"] * -1 + group_length
    sec_group_table[nrow(sec_group_table):1,] <- rev_group
    }
  Export_table <- rbind(Export_table,sec_group_table)
}

print(nrow(Export_table))

write.table(Export_table, output_file_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, na = "")