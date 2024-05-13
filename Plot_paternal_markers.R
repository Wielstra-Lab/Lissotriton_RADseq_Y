
# Plot number of paternal markers and SNPs per linkage group

library(ggplot2)
library(ggrepel)
library(gridExtra)

# First import the map, and lists of paternal markers and SNPs

args_in <- commandArgs(trailingOnly = TRUE)

Linkage_map_file <- args_in[1]
paternal_marker_file <- args_in[2]
paternal_snp_file <- args_in[3]
output_name <- args_in[4]

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

Groups <- data.frame(1:nGroups)

Groups[,"N_markers"] <- 0
Groups[,"Pat_markers"] <- 0
Groups[,"Pat_SNPs"] <- 0
Groups[,"Y_markers"] <- 0
Groups[,"Colour"] <- "black"

for(i in 1:nGroups) {
  
  Group_subset <- subset(Linkage_map_table, V3 == i)
  
  Groups[i,"N_markers"] <- nrow(Group_subset)
  Groups[i,"Pat_markers"] <- sum(Group_subset$paternal)
  Groups[i,"Pat_SNPs"] <- sum(Group_subset$pat_SNPs)
  Groups[i,"Y_markers"] <- sum(Group_subset$Y_PAV)
  
  if(Groups[i,"Y_markers"] > 0) {Groups[i,"Colour"] <- "red"}
}

# Make plots

Xlim <- max(Groups$N_markers) * 1.05
Ylim1 <- max(Groups$Pat_markers) * 1.05
Ylim2 <- max(Groups$Pat_SNPs) * 1.05

text_size <- 15
marker_label_size <- 7
  
Plot_A <- ggplot(Groups, aes(x=N_markers, y=Pat_markers, label = X1.nGroups)) + 
  geom_smooth(method=lm, formula=y~x-1, se=FALSE, fullrange = TRUE, linetype = "dashed", color = "grey") +
  geom_point(colour = Groups$Colour, size = 4) +
  geom_text_repel(vjust=1.5, size = marker_label_size) +
  scale_x_continuous(expand = c(0,0), limits = c(0,Xlim)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,Ylim1)) +
  ggtitle("Fig 5A") +
  labs(x = "Total Markers") +
  labs(y = "Paternal Specific Markers") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
        axis.title.x = element_text(size = text_size), axis.text.x = element_text(size = text_size - 2),
        axis.title.y = element_text(size = text_size), axis.text.y = element_text(size = text_size - 2),
        plot.title = element_text(size = text_size + 2))

Plot_B <- ggplot(Groups, aes(x=N_markers, y=Pat_SNPs, label = X1.nGroups)) + 
  geom_smooth(method=lm, formula=y~x-1, se=FALSE, fullrange = TRUE, linetype = "dashed", color = "grey") +
  geom_point(colour = Groups$Colour, size = 4) +
  geom_text_repel(vjust=1.5, size = marker_label_size) +
  scale_x_continuous(expand = c(0,0), limits = c(0,Xlim)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,Ylim2)) +
  ggtitle("Fig 5B") +
  labs(x = "Total Markers") +
  labs(y = "Paternal Specific SNPs") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA),
        axis.title.x = element_text(size = text_size), axis.text.x = element_text(size = text_size - 2),
        axis.title.y = element_text(size = text_size), axis.text.y = element_text(size = text_size - 2),
        plot.title = element_text(size = text_size + 2))

ggsave(output_name, arrangeGrob(Plot_A, Plot_B, ncol=2), width = 14, height = 6, units = "in")

sum(Groups$Pat_markers)   
sum(Groups$Pat_SNPs)
