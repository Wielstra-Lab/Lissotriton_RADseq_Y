library(ggplot2)
library(grid)
library(gridExtra)
library(gridtext)
library(ggtext)

# First import the map, and lists of paternal markers and SNPs

args_in <- commandArgs(trailingOnly = TRUE)

Linkage_map_file <- args_in[1]
paternal_marker_file <- args_in[2]
paternal_snp_file <- args_in[3]
output_name <- args_in[4]

Bin_size <- 2

Linkage_map_table <- read.table(Linkage_map_file)
paternal_marker_table <- read.table(paternal_marker_file)
paternal_snp_table <- read.table(paternal_snp_file)

# Score each map marker for paternal status

Linkage_map_table[,"paternal"] <- Linkage_map_table$V1 %in% paternal_marker_table$V1

# Score each map marker for number of paternal SNPs

Linkage_map_table[,"pat_SNPs"] <- 0

num_snps <- table(paternal_snp_table$V1)

for(i in 1:nrow(Linkage_map_table)) {Linkage_map_table[i,"pat_SNPs"] <- num_snps[as.character(Linkage_map_table[i,"V1"])]}

Linkage_map_table$pat_SNPs[is.na(Linkage_map_table$pat_SNPs)] <- 0

Linkage_map_table[,"Y_PAV"] <- grepl("Y", Linkage_map_table$V1)

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

Group_lenghts[,"N_markers"] <- 0
Group_lenghts[,"Pat_markers"] <- 0
Group_lenghts[,"Pat_SNPs"] <- 0
Group_lenghts[,"Y_markers"] <- 0

# Make bins

Group_lenghts[,"Bins"] <- 1 + (Group_lenghts$Length %/% Bin_size)

Master_Bin_table <- data.frame(matrix(NA, nrow = 0, ncol = 8))
colnames(Master_Bin_table) <- c("Group","Start","End","X_pos","N_markers","N_paternal","N_SNPs","Y_markers")

for(i in 1:nrow(Group_lenghts)) {
  
  Bin_table <- data.frame(matrix(NA, nrow = Group_lenghts[i,"Bins"], ncol = 8))
  colnames(Bin_table) <- c("Group","Start","End","X_pos","N_markers","N_paternal","N_SNPs","Y_markers")
  
  Bin_table$Group <- i
  Bin_table$End <- 1:nrow(Bin_table) * Bin_size
  Bin_table$Start <- Bin_table$End - Bin_size
  Bin_table$X_pos <- (0.5 * Bin_size) + Bin_table$Start + Group_lenghts[i,"Start_pos"]
  
  Group_subset <- subset(Linkage_map_table, V3 == i)
  
  Group_lenghts[i,"N_markers"] <- nrow(Group_subset)
  Group_lenghts[i,"Pat_markers"] <- sum(Group_subset$paternal)
  Group_lenghts[i,"Pat_SNPs"] <- sum(Group_subset$pat_SNPs)
  Group_lenghts[i,"Y_markers"] <- sum(Group_subset$Y_PAV)
  
  for(bin in 1:nrow(Bin_table)) {
    
    Bin_subset <- subset(Group_subset, V2 >= Bin_table[bin,"Start"] & V2 < Bin_table[bin,"End"])
    Bin_table[bin,"N_markers"] <- nrow(Bin_subset) - sum(Bin_subset$Y_PAV)
    Bin_table[bin,"N_paternal"] <- sum(Bin_subset$paternal)
    Bin_table[bin,"N_SNPs"] <- sum(Bin_subset$pat_SNPs)
    Bin_table[bin,"Y_markers"] <- sum(Bin_subset$Y_PAV)
  }
  
  Master_Bin_table <- rbind(Master_Bin_table,Bin_table)
}

# Score bins, and compute stats

Pat_marker_freq <- sum(Linkage_map_table$paternal) / nrow(Linkage_map_table)

Pat_SNP_freq <- sum(Linkage_map_table$pat_SNPs) / nrow(Linkage_map_table)

Master_Bin_table[,"Pat_marker_ratio"] <- Master_Bin_table$N_paternal / Master_Bin_table$N_markers

Master_Bin_table[,"Pat_marker_prob"] <- pbinom(Master_Bin_table$N_paternal -1, size=Master_Bin_table$N_markers, prob=Pat_marker_freq, lower.tail = FALSE)

Master_Bin_table[,"Pat_marker_log"] <- -1 * log10(Master_Bin_table$Pat_marker_prob)

Master_Bin_table[,"Pat_SNP_density"] <- Master_Bin_table$N_SNPs / Master_Bin_table$N_markers

#Master_Bin_table[,"Pat_SNP_prob"] <- pbinom((Master_Bin_table$N_SNPs / 5) - 1, size=Master_Bin_table$N_markers, prob=(Pat_SNP_freq / 5), lower.tail = FALSE)

Master_Bin_table[,"Pat_SNP_prob"] <- ppois(Master_Bin_table$N_SNPs, lambda=(Master_Bin_table$N_markers * Pat_SNP_freq), lower=FALSE)

Master_Bin_table[,"Pat_SNP_log"] <- -1 * log10(Master_Bin_table$Pat_SNP_prob)

Master_Bin_table[sapply(Master_Bin_table, is.infinite)] <- NA

Master_Bin_table[,"Bin_colour"] <- "black"

Master_Bin_table$Bin_colour[Master_Bin_table$Group %% 2 == 0] <- "dark grey"

Master_Bin_table$Bin_colour[Master_Bin_table$Y_markers > 0] <- "red"

Master_Bin_table[,"dot_size"] <- 2

Master_Bin_table$dot_size[Master_Bin_table$Y_markers > 0] <- 4

Master_Bin_table <- Master_Bin_table[order(Master_Bin_table$dot_size),]

Y_subset <- subset(Master_Bin_table, Y_markers > 0)
max_Y_marker_log <- max(Y_subset$Pat_marker_log)
max_Y_SNP_log <- max(Y_subset$Pat_SNP_log)

max_Y_marker_prob <- min(Y_subset$Pat_marker_prob)
max_Y_SNP_prob <- min(Y_subset$Pat_SNP_prob)

print(paste("Total bins:", nrow(Master_Bin_table)))
print(paste("Bins with higher max marker log than highest Y bin:", sum(Master_Bin_table$Pat_marker_log > max_Y_marker_log)))
print(paste("Bins with higher max SNP log than highest Y bin:", sum(Master_Bin_table$Pat_SNP_log > max_Y_SNP_log)))

# Plot



Group_lenghts[,"Mid_point"] <- Group_lenghts$Start_pos + (0.5 * Group_lenghts$Length)

text_size <- 15

Plot_A <- ggplot(Master_Bin_table, aes(x=X_pos, y=Pat_marker_log)) + 
  geom_hline(yintercept = max_Y_marker_log, linetype="dashed") +
  geom_point(colour = Master_Bin_table$Bin_colour, size = Master_Bin_table$dot_size) +
  scale_x_continuous(breaks = Group_lenghts$Mid_point, minor_breaks = Group_lenghts$Start_pos,
                     labels=Group_lenghts$X1.nGroups, expand = c(0.01,0)) +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA)) +
  ggtitle("Fig 6A") +
  labs(x = "<i>Lissotriton vulgaris</i> Linkage Group") +
  labs(y = "Paternal Specific Marker -log<sub>10</sub>(<i>P</i>)") +
  theme(axis.title.y = ggtext::element_markdown(size = text_size)) +
  theme(axis.title.x = ggtext::element_markdown(size = text_size)) +
  theme(axis.text.x = element_text(size = text_size - 2),
        axis.text.y = element_text(size = text_size - 2),
        plot.title = element_text(size = text_size + 2))

Plot_B <- ggplot(Master_Bin_table, aes(x=X_pos, y=Pat_SNP_log)) + 
  geom_hline(yintercept = max_Y_SNP_log, linetype="dashed") +
  geom_point(colour = Master_Bin_table$Bin_colour, size = Master_Bin_table$dot_size) +
  scale_x_continuous(breaks = Group_lenghts$Mid_point, minor_breaks = Group_lenghts$Start_pos,
                     labels=Group_lenghts$X1.nGroups, expand = c(0.01,0)) +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA)) +
  ggtitle("Fig 6B") +
  labs(x = "<i>Lissotriton vulgaris</i> Linkage Group") +
  labs(y = "Paternal Specific SNP -log<sub>10</sub>(<i>P</i>)") +
  theme(axis.title.y = ggtext::element_markdown(size = text_size)) +
  theme(axis.title.x = ggtext::element_markdown(size = text_size)) +
  theme(axis.text.x = element_text(size = text_size - 2),
        axis.text.y = element_text(size = text_size - 2),
        plot.title = element_text(size = text_size + 2))
  

ggsave(output_name, arrangeGrob(Plot_A, Plot_B), width = 14, height = 8, units = "in")
  
