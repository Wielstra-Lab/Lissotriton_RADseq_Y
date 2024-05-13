
# Script that ingests a linkage map (with paternal and maternal positions) and plots marker density

library(ggplot2)
library(ggtext)

args_in <- commandArgs(trailingOnly = TRUE)

Paternal_map_file <- args_in[1]
Maternal_map_file <- args_in[2]
output_name <- args_in[3]

custom_smooth <- function(s_vector, s_mag = 5){
  
  s_len <- length(s_vector)
  
  return_vec <- c()
  
  # First calculate weights, and store as vector
  
  weight_calc_points <- (s_mag - 1) / 2
  weight_calc_dist <- 2 / s_mag
  mid_point <- weight_calc_points + 1
  
  raw_weights <- c((pnorm(weight_calc_dist) - 0.5) * 2)
  
  for(i in 1:weight_calc_points) {
    new_weight <- pnorm((2 * i + 1) * weight_calc_dist) - pnorm((2 * i - 1) * weight_calc_dist)
    raw_weights <- c(raw_weights,new_weight)
  }
  
  for_weights <- rev(raw_weights[2:length(raw_weights)])
  
  all_weights <- c(for_weights,raw_weights)
  
  # Then for each point in vector to be smoothed, would the weights vector over/undershoot?
  
  for(i in 1:s_len) {
    
    undershoot <- mid_point - i
    overshoot <- (mid_point + i) - (s_len + 1)
    
    if(undershoot > 0) {
      
      weight_start <- undershoot + 1
      weight_end <- s_mag
      tab_start <- (i - weight_calc_points) + undershoot
      tab_end <- i + weight_calc_points
    }
    
    else if(overshoot > 0) {
      
      weight_start <- 1
      weight_end <- s_mag - overshoot
      tab_start <- i - weight_calc_points
      tab_end <- (i + weight_calc_points) - overshoot
    }
    
    else {
      
      weight_start <- 1
      weight_end <- s_mag
      tab_start <- i - weight_calc_points
      tab_end <- i + weight_calc_points
    }
    
    vec_weights <- all_weights[weight_start : weight_end]
    vec_vals <- s_vector[tab_start : tab_end]
    smoothed_val <- sum(vec_vals * vec_weights) / sum(vec_weights)
    return_vec <- c(return_vec,smoothed_val)
  }
  
  return(return_vec)
}

# load map files

Male_table <- read.table(Paternal_map_file, col.names = c("Marker", "Position", "Group"))
Female_table <- read.table(Maternal_map_file, col.names = c("Marker", "Position", "Group"))
Male_table[,"Sex"] <- "Paternal"
Female_table[,"Sex"] <- "Maternal"
colnames(Male_table) <- c("Marker", "Position", "Group", "Sex")
colnames(Female_table) <- c("Marker", "Position", "Group", "Sex")

Sexed_table <- rbind(Female_table,Male_table)

# subset and bin


nGroups <- max(Sexed_table$Group)
Bin_size <- 2

Master_Bin_table <- data.frame(matrix(NA, nrow = 0, ncol = 8))
colnames(Master_Bin_table) <- c("Group","Sex","Start","End","N_markers","Smooth_markers","Y_linked","Y_smooth")

padding_bins <- 5

for (sex in c("Paternal","Maternal")) {
  
  for (group in 1:nGroups) {
    
    Group_subset <- subset(Sexed_table, Sex == sex & Group == group)
    
    Group_length <- Group_subset[nrow(Group_subset),"Position"]
    
    nBins <- 1 + Group_length %/% Bin_size + padding_bins * 2
    
    print("step 0")
    
    Bin_table <- Master_Bin_table[1:nBins,]
    
    Bin_table$Group <- group
    Bin_table$End <- (1:nBins * Bin_size) - padding_bins * Bin_size
    Bin_table$Start <- Bin_table$End - Bin_size
    Bin_table$Sex <- sex
    
    print("step 1")
    
    for(bin in 1:nrow(Bin_table)) {
      
      Bin_subset <- subset(Group_subset, Position >= Bin_table[bin,"Start"] & Position < Bin_table[bin,"End"])
      Bin_table[bin,"N_markers"] <- nrow(Bin_subset)
      Bin_table[bin,"Y_linked"] <- sum(grepl("Y",Bin_subset$Marker))
    }
    
    print("step 2")
    
    Bin_table$Smooth_markers <- custom_smooth(Bin_table$N_markers, 11) / Bin_size
    
    Bin_table$Y_smooth <- custom_smooth(Bin_table$Y_linked, 7) / Bin_size
    
    Master_Bin_table <- rbind(Master_Bin_table,Bin_table)
  }
  
}

Master_Bin_table[,"Y_linked"] <- as.numeric(Master_Bin_table$Y_linked)

Master_Bin_table[,"Sex"] <- factor(Master_Bin_table$Sex, levels = c("Paternal","Maternal"))

text_size <- 15

plot <- ggplot(Master_Bin_table) + 
  geom_area(aes(x=Start, y=Smooth_markers), linewidth = 1) + 
  geom_area(aes(x=Start, y=Y_smooth), linewidth = 1, fill = "red") +
  facet_grid(cols = vars(Master_Bin_table$Group), rows = vars(Master_Bin_table$Sex), scales = "free", space = "free") +
  scale_x_continuous(breaks = c(0,100),expand=c(0.1,0)) +
  scale_y_continuous(breaks = c(0,10,20,30,40,50), minor_breaks = c(-0.01), expand=c(0.01,0)) + 
  labs(x = "Position in Linkage Group (cM from origin)") +
  labs(y = "Marker Density (markers/cM)") +
  labs(title='Fig 4') +
  labs(subtitle='_Lissotriton vulgaris_ Linkage Group') +
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_line(colour="grey", linewidth=1),
        panel.spacing.x = unit(0, "null"),
        panel.background = element_blank(),
        strip.background = element_rect(fill="white"),
        strip.text.y.right = element_text(angle = 90, size = text_size),
        strip.text.x.top = element_text(size = text_size),
        plot.title.position = "plot",
        plot.title = element_text(size = text_size + 2),
        plot.subtitle = ggtext::element_markdown(hjust = 0.5, size = text_size + 2),
        axis.title.x = element_text(size = text_size),
        axis.title.y = element_text(size = text_size),
        axis.text.y = element_text(size = text_size - 2))

ggsave(output_name, plot,  width = 14, height = 6, units = "in")

max(Master_Bin_table$N_markers)

