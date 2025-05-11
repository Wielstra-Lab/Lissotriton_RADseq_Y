library(grid)
library(gridtext)

args_in <- commandArgs(trailingOnly = TRUE)

input_file <- args_in[1]
special_marker_file <- args_in[2]
output_file <- args_in[3]
title_text <- args_in[4]

Load_map <- function(map_file, secondary_map_file = NA, special_markers = NA, default_colour = "black") {
  
  input_table <- read.table(map_file, col.names = c("Marker","Pos","Group"))
  
  if(is.na(secondary_map_file) == FALSE) {
    
    input_table_2 <- read.table(secondary_map_file, col.names = c("Marker","Pos","Group"))
    
    input_table[,"Group_2"] <- NA
    input_table[,"Pos_2"] <- NA
    
    for (i in 1:nrow(input_table)) {
      input_table[i,"Group_2"] <- input_table_2[match(input_table[i,"Marker"], input_table_2$Marker), "Group"]
    }
    
    for (i in 1:nrow(input_table)) {
      input_table[i,"Pos_2"] <- input_table_2[match(input_table[i,"Marker"], input_table_2$Marker), "Pos"]
    }
  }
  
  input_table["Colour"] <- default_colour
  
  input_table[,"Special"] <- FALSE
  
  if(is.na(special_markers) == FALSE) {
    
    input_table[,"Special"] <- input_table$Marker %in% special_markers$Marker
    
    for(i in 1:nrow(input_table)) {
      if(input_table[i,"Special"] == TRUE){
        input_table[i,"Colour"] <- special_markers[match(input_table[i,"Marker"],special_markers$Marker),"Colour"]
      }
    }
  }
  
  return(input_table)
}

plot_group <- function(group_table,x_origin = 0, y_origin = 1, Gname = "G", scale_factor = 200) {
  
  # Determine dimensions
  
  line_weight <- 2
  
  linkage_group_length <- max(group_table$Pos) / scale_factor
  linkage_group_width <- 0.03
  
  x_left <- x_origin
  x_right <- x_origin + linkage_group_width
  y_top <- y_origin
  y_bottom <- y_top - linkage_group_length
  
  x_mid <- x_left + (linkage_group_width / 2)
  y_text <- y_top + linkage_group_width + 0.03
  
  # Draw perimeter
  
  left_edge <- segmentsGrob(x_left,y_top,x_left,y_bottom)
  right_edge <- segmentsGrob(x_right,y_top,x_right,y_bottom)
  top_cap <- curveGrob(x_left,y_top,x_right,y_top, curvature=-arcCurvature(180),square=F,ncp=10)
  bottom_cap <- curveGrob(x_left,y_bottom,x_right,y_bottom, curvature=arcCurvature(180),square=F,ncp=10)
  
  perimeter <- gTree(children = gList(left_edge,right_edge,top_cap,bottom_cap), gp = gpar(lwd = line_weight))
  
  # Add group title
  
  group_title <- textGrob(Gname, x_mid, y_text) 
  
  # Add marker bars
  
  group_table[,"Plot_y"] <- y_top - (group_table$Pos / scale_factor)
  
  reg_table <- subset(group_table, Special == FALSE)
  
  reg_bars <- gTree()
  
  for(i in 1:nrow(reg_table)) {
    
    y_pos <- reg_table[i,"Plot_y"]
    bar <- segmentsGrob(x_left,y_pos,x_right,y_pos, gp = gpar(col = reg_table[i,"Colour"]))
    reg_bars <- addGrob(reg_bars, bar)
  }
  
  spec_bars <- gTree()
  
  if(sum(group_table$Special == TRUE) > 0){
    
    spec_table <- subset(group_table, Special == TRUE)
    
    for(i in 1:nrow(spec_table)) {
      
      y_pos <- spec_table[i,"Plot_y"]
      bar <- segmentsGrob(x_left,y_pos,x_right,y_pos, gp = gpar(col = spec_table[i,"Colour"], lwd = 2))
      spec_bars <- addGrob(spec_bars, bar)
    }
  }
  
  group_grob <- gTree(children = gList(perimeter,group_title,reg_bars,spec_bars))
  return(group_grob)
}

plot_scale <- function(x_position, y_top, y_length, tick_space = 20, scale_factor = 200, reverse = FALSE) {
  
  y_bottom <- y_top - y_length / scale_factor
  
  y_mid <- y_top - 0.5 * y_length / scale_factor
  
  # make the backbone
  
  grid.segments(x_position, y_top, x_position, y_bottom)
  
  # make the tick marks
  
  major_tick_length <- 0.02
  
  if(reverse == TRUE){major_tick_length <- major_tick_length * -1}
  
  major_tick_end <- x_position + major_tick_length 
  
  major_tick_seq <- seq(0,y_length,tick_space)
  
  major_tick_pos <- y_top - (major_tick_seq / scale_factor) 
  
  for(i in 1:length(major_tick_seq)){
    grid.segments(x_position, major_tick_pos[i], major_tick_end, major_tick_pos[i])
  }
  
  minor_tick_space <- tick_space / 4
  
  minor_tick_length <- major_tick_length / 2
  
  minor_tick_end <- x_position + minor_tick_length 
  
  minor_tick_seq <- seq(0,y_length,minor_tick_space)
  
  minor_tick_pos <- y_top - (minor_tick_seq / scale_factor) 
  
  for(i in 1:length(minor_tick_seq)){
    grid.segments(x_position, minor_tick_pos[i], minor_tick_end, minor_tick_pos[i])
  }
  
  # draw the labels
  
  text_just <- "left"
  
  if(reverse == TRUE){text_just <- "right"}
  
  for(i in 1:length(major_tick_seq)){
    grid.text(major_tick_seq[i], major_tick_end + major_tick_length * 0.5, major_tick_pos[i], just = text_just, gp = gpar(fontsize = 8))
  }
  
  # draw the title:
  
  Title_text <- "Length in cM"
  grid.text(Title_text, x_position + major_tick_length * -0.5, y_mid, rot = 90, just = "bottom")
  
}

Y_markers <- read.table(special_marker_file, col.names = c("Marker","Colour"))

Map_table <- Load_map(input_file, special_markers = Y_markers)

pdf(output_file)

grid.newpage()

plot_scale(0.05, 0.8, 140)

for(i in 1:12){
  
  group_i <- Map_table[Map_table$Group == i,]
  
  grid.draw(plot_group(group_i, (i * 0.075) + 0.05, 0.8, i))
}

Title <- richtext_grob(title_text, 0.5, 0.92)
grid.draw(Title)

dev.off()


