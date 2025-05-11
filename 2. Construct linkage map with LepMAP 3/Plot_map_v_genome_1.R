
library(ggplot2)
library(grid)
library(gridtext)

args_in <- commandArgs(trailingOnly = TRUE)

Linkage_map_file <- args_in[1]
Blast_result_file <- args_in[2]
Genome_struct_file <- args_in[3]
output_name <- args_in[4]

# Import the linkage map, filtered blast results and genome structure

Linkage_map_table <- read.table(Linkage_map_file)
Blast_result_table <- read.table(Blast_result_file)
Genome_struct_table <- read.table(Genome_struct_file)

# Get the start positions of the B chromosome segments

Genome_struct_table[,"start_pos"] <- 0

for(i in 1:nrow(Genome_struct_table)) {
  if( Genome_struct_table[i,4] == 2) {Genome_struct_table[i,5] <- Genome_struct_table[i - 1,2] }
}

# Reconstruct blast coordinates

Blast_result_table[,"Raw_scaff"] <- match(Blast_result_table$V2,Genome_struct_table$V1)

for(i in 1:nrow(Blast_result_table)) {
  
  raw_scaff_no <- Blast_result_table[i,"Raw_scaff"]
  
  Blast_result_table[i,9] <- Blast_result_table[i,9] + Genome_struct_table[raw_scaff_no,5]
  Blast_result_table[i,10] <- Blast_result_table[i,10] + Genome_struct_table[raw_scaff_no,5]
  Blast_result_table[i,2] <- Genome_struct_table[raw_scaff_no,3]
  
}

# Then add blast data to map

Blast_result_table[,"Map_no"] <- match(Blast_result_table$V1,Linkage_map_table$V1)

Linkage_map_table[,"Genome_Chr"] <- NA
Linkage_map_table[,"Chr_pos"] <- NA

for(i in 1:nrow(Blast_result_table)) {
  
  Map_no <- Blast_result_table[i,"Map_no"]
  
  if(Map_no %in% row.names(Linkage_map_table)){
    Linkage_map_table[Map_no,"Genome_Chr"] <- Blast_result_table[i,2]
    Linkage_map_table[Map_no,"Chr_pos"] <- Blast_result_table[i,9]
  }
}


# Extract map structure

nGroups <- max(Linkage_map_table$V3)

Group_lenghts <- data.frame(1:nGroups)
Group_lenghts[,"Length"] <- NA

for(group in 1:nGroups) {
  
  group_end <- findInterval(group,Linkage_map_table$V3)
  Group_lenghts[group,"Length"] <- as.numeric(Linkage_map_table[group_end,2])
}

Group_lenghts[,"Start_pos"] <- 0

for(group in 2:nGroups) {
  Group_lenghts[group,"Start_pos"] <- Group_lenghts[group - 1,2] + Group_lenghts[group - 1,3]
}

# Filter down to only blasted loci

Hits_table <- na.omit(Linkage_map_table)

# Process genome structure

Genome_struct_table[,"extra_length"] <- 0

for(i in 2:nrow(Genome_struct_table)) {
  if( Genome_struct_table[i,4] == 2) {Genome_struct_table[i - 1,"extra_length"] <- Genome_struct_table[i,2] }
}

Genome_struct_table[,"total_length"] <- Genome_struct_table$V2 + Genome_struct_table$extra_length

Chr_lengths <- subset(Genome_struct_table, V4 == 1, select = c("V3","total_length"))

# Group stats

stats_table <- data.frame(matrix(NA, nrow = nGroups, ncol = nrow(Chr_lengths)))
colnames(stats_table) <- Chr_lengths$V3

for(i in 1:nGroups){
  
  group_table <- subset(Hits_table, V3==i)
  
  for(j in 1:nrow(Chr_lengths)){
    
    stats_table[i,j] <- sum(group_table$Genome_Chr == Chr_lengths$V3[j], na.rm=T)
  }
}

stats_table[,"Total_hits"] <- rowSums(stats_table)
stats_table[,"Best_match"] <- NA
stats_table[,"recipricol"] <- NA

for(i in 1:nGroups){
  
  best_match_index <- which.max(stats_table[i,1:nrow(Chr_lengths)])
  stats_table[i,"Best_match"] <- colnames(stats_table)[best_match_index]
  stats_table[i,"recipricol"] <- which.max(stats_table[,best_match_index])
  stats_table[i,"Best_match_hits"] <- stats_table[i,best_match_index]
}

stats_table[,"correlation"] <- NA
stats_table[,"inverted"] <- FALSE

for(i in 1:nGroups){
  
  group_table <- subset(Hits_table, V3==i & Genome_Chr == stats_table[i,"Best_match"])
  
  coefficient <- cor.test(group_table$V2, group_table$Chr_pos)
  stats_table[i,"correlation"] <- coefficient$estimate
  
  if(coefficient$estimate < 0) { stats_table[i,"inverted"] <- TRUE}
}

# Flip inverted linkage groups

for(i in 1:nrow(Hits_table)){
  
  group_id <- Hits_table[i,3]
  
  if(stats_table[group_id,"inverted"] == TRUE) {
    
    group_len <- Group_lenghts[group_id,2]
    Hits_table[i,2] <- (Hits_table[i,2] * - 1) + group_len
    
  }
  
}

# Re-order genome structure and calculate start positions

Chr_lengths[,"recipricol"] <- 0

for(i in 1:nrow(Chr_lengths)) {
  
  chrom <- Chr_lengths[i,1]
  Chr_lengths[i,"recipricol"] <- which.max(stats_table[,chrom])
}

Chr_lengths <- Chr_lengths[order(Chr_lengths$recipricol),]

Chr_lengths[,"Start_pos"] <- 0

for(chr in 2:nrow(Chr_lengths)) {
  
  Chr_lengths[chr,"Start_pos"] <- Chr_lengths[chr - 1,2] + Chr_lengths[chr -1,4]
}

# Get X co-ordinates

Hits_table[,"X_coord"] <- NA

for(i in 1:nrow(Hits_table)) {
  
  marker_group <- as.numeric(Hits_table[i,3])
  Hits_table[i,"X_coord"] <- Hits_table[i,2] + Group_lenghts[marker_group,"Start_pos"]
  
}

# Get Y co-ordinates

Hits_table[,"CHR_no"] <- match(Hits_table$Genome_Chr,Chr_lengths$V3)

Hits_table[,"Y_coord"] <- NA

Hits_table[,"Colour"] <- "black"

Hits_table[,"Size"] <- 1

for(i in 1:nrow(Hits_table)) {
  
  marker_chr <- as.numeric(Hits_table[i,"CHR_no"])
  Hits_table[i,"Y_coord"] <- Hits_table[i,5] + Chr_lengths[marker_chr,"Start_pos"]
  if(grepl("Y",Hits_table[i,"V1"])==TRUE) {Hits_table[i,"Colour"] <- "red"}
  if(grepl("Y",Hits_table[i,"V1"])==TRUE) {Hits_table[i,"Size"] <- 3}
  
}

Hits_table <- Hits_table[order(Hits_table$Size),]

# Make plot

Pick_scale <- function(total_length,fraction = 0.25) {
  optimal_scale_length <- total_length * fraction
  scale_mag <- 10^(floor(log10(total_length)) - 1)
  Scale_factors <- c(1,2,5,10)
  Scale_options <- Scale_factors * scale_mag
  Scale_results <- abs(Scale_options - optimal_scale_length)
  scale_best <- Scale_options[which.min(Scale_results)]
  return(scale_best)
}

Group_lenghts[,"Mid_point"] <- Group_lenghts$Start_pos + (0.5 * Group_lenghts$Length)
Chr_lengths[,"Mid_point"] <- Chr_lengths$Start_pos + (0.5 * Chr_lengths$total_length)

Group_font <- 20
Axis_font <- 20

X_text <- "_Lissotriton vulgaris_ Linkage Group"
X_length <- max(Hits_table$X_coord)
X_axis_offset <- -2.0e+09
X_scale <- Pick_scale(X_length)
X_scale_text <- paste(X_scale,"cM")
x_text_offset <- 0.75
x_midpoint <- X_length / 2
  
y_text <- "<i>Plurodeles  waltl</i>  Chromosome"
y_axis_offset <- -50
y_scale <- 5e+09
y_scale_text <- "5 Gbp"
y_text_offset <- 0.75
y_midpoint <- max(Hits_table$Y_coord) / 2

axis_par <- gpar(fontsize = Axis_font, fontfamily = "Helvetica")

final_plot <- ggplot(Hits_table, aes(x=X_coord, y=Y_coord)) + geom_point(colour = Hits_table$Colour, size = Hits_table$Size) + 
  scale_x_continuous(breaks = Group_lenghts$Mid_point, minor_breaks = Group_lenghts$Start_pos, expand = c(0,0), labels=Group_lenghts$X1.nGroups) +
  scale_y_continuous(breaks = Chr_lengths$Mid_point, minor_breaks = Chr_lengths$Start_pos, expand = c(0,0), labels=Chr_lengths$V3) +
  theme(panel.grid.minor = element_line(colour="black", linewidth=0.5), panel.grid.major = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(color= "black", size = Group_font)) +
  theme(axis.text.y = element_text(color= "black", size = Group_font)) +
  
  annotation_custom(grid::linesGrob(gp = gpar(lwd = unit(2, "pt"))), xmin = 0, xmax = X_scale, ymin = X_axis_offset, ymax = X_axis_offset)  +
  annotation_custom(richtext_grob(gp = axis_par, X_scale_text, vjust = -0.5), xmin = 0, xmax = X_scale, ymin = X_axis_offset, ymax = X_axis_offset) +
  annotation_custom(richtext_grob(gp = axis_par, X_text, vjust = -0.5), xmin = x_midpoint, xmax = x_midpoint, ymin = X_axis_offset, ymax = X_axis_offset) +
  
  annotation_custom(grid::linesGrob(gp = gpar(lwd = unit(2, "pt"))), xmin = y_axis_offset, xmax = y_axis_offset, ymin = 0, ymax = y_scale)  +
  annotation_custom(richtext_grob(gp = axis_par, y_scale_text, rot = 90, vjust = -0.5), xmin = y_axis_offset, xmax = y_axis_offset, ymin = 0, ymax = y_scale) +
  annotation_custom(richtext_grob(gp = axis_par, y_text, rot = 90, vjust = -0.5), xmin = y_axis_offset, xmax = y_axis_offset, ymin = y_midpoint, ymax = y_midpoint) +
  
  coord_cartesian(clip = 'off') +
  theme(plot.margin=unit(c(0.25,0.25,1,1), 'in')) +
  labs(title=output_name) +
  theme(plot.title = element_text(size = 25))
  
ggsave(paste(output_name,".png", sep = ""), plot = final_plot, width = 14, height = 12, units = "in")

write.table(stats_table, paste(output_name,"_stats.txt", sep = ""), sep = "\t", quote = FALSE)

