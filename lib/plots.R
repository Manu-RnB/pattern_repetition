# ------ Theme definition
apa7 <- theme(panel.grid = element_blank(), # to delete the grid
              panel.background = element_rect(fill = "#FFFFFF"), # to add white background to panel
              plot.background = element_rect(fill = "#FFFFFF"), # to add white background to plot
              axis.text.y = element_text(color = "black"),
              axis.text.x = element_text(color = "black"), # to have black text
              axis.line = element_line(color ="black",
                                       linewidth = 0.231), # to add a black line
              axis.ticks = element_line(color ="black",
                                        linewidth = 0.231), # to add black ticks
              axis.title.y = element_text(face = "bold", # to have bold title
                                          margin = margin(t = 0, r = 15, b = 0, l = 0)), # to increase space between title and text label
              axis.title.x = element_text(face = "bold",
                                          margin = margin(t = 15, r = 0, b = 0, l = 0)),
              plot.subtitle = element_text(hjust = 0.5), # to have centered subtitle
              legend.position = "top", # to modify legend position
              rect = element_blank() # to delete background color
) 

# To change police size
theme_set(
  theme_classic(base_size = 20)
)