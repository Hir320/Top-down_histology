rm(list = ls())
library(ggplot2)
library(cowplot)
library(rlang)
# Define the directory where you want to search for the files
setwd('/home/rstudio/')


# Function ----------------------------------------------------------------

make_zoom_plots <- function(
    df,
    xvar, yvar,
    main_xlim    = c(0, 250),
    main_ylim    = c(-0.1, 1),
    main_xbreaks = seq(0, 300, by = 50),
    main_ybreaks = seq(-0.2, 1,   by = 0.2),
    inset_xlim,
    inset_ylim,
    inset_xbreaks,
    inset_ybreaks,
    inset_pos    = c(x = 0.45, y = 0.45, width = 0.45, height = 0.45),
    rect_args    = list(fill = NA, color = "black", linetype = "dashed", size =  0.5/2.13),
    inset_theme  = theme(
      axis.title      = element_blank(),
      axis.text       = element_text(size = 7),
      axis.line                = element_line(linewidth = 0.5/2.13),
      axis.ticks               = element_line(linewidth = 0.5/2.13),
      plot.background          = element_rect(fill = "transparent", color = NA),
      panel.background         = element_rect(fill = "transparent")
    )
) {
  # 1) Main plot with proper aes()
  p_main <- ggplot(df, aes(
    x = .data[[xvar]],
    y = .data[[yvar]]
  )) +
    geom_line(linewidth = 0.5/2.13) +
    #geom_point() +
    scale_x_continuous(
      limits = main_xlim,
      breaks = main_xbreaks,
      expand = expansion(mult = c(0, 0.1))
    ) +
    scale_y_continuous(
      limits = main_ylim,
      breaks = main_ybreaks
    ) +
    theme_classic() +
    theme(
      axis.title        = element_blank(),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.ticks.length = unit(0.5, "mm"),
      axis.line                = element_line(linewidth = 0.5/2.13),
      axis.ticks               = element_line(linewidth = 0.5/2.13),
      plot.background          = element_rect(fill = "transparent", color = NA),
      panel.background         = element_rect(fill = "transparent")
    ) +
    #labs(x = "Distance (µm)", y = "Normalized Autocorrelation") +
    geom_hline(yintercept = 0, linewidth =  0.5/2.13)
  
  # 2) Add the rectangle showing inset region
  p_main_rect <- p_main +
    annotate("rect",
             xmin     = inset_xlim[1], xmax = inset_xlim[2],
             ymin     = inset_ylim[1], ymax = inset_ylim[2],
             fill     = rect_args$fill,
             color    = rect_args$color,
             linetype = rect_args$linetype,
             size     = rect_args$size
    )
  
  # 3) Build the inset panel
  p_inset <- p_main +
    scale_x_continuous(limits = inset_xlim, breaks = inset_xbreaks) +
    scale_y_continuous(limits = inset_ylim, breaks = inset_ybreaks) +
    inset_theme
  
  # 4) Return both
  list(
    main_plot_rect = p_main_rect,
    inset_plot     = p_inset,
    inset_position = inset_pos
  )
}



# Load your data ----------------------------------------------------------
# 1) Load your data
df_RSC <- read_csv("./Data/Histology/Retinotopic_quantification/Clustering_individual/csv/118_RSC.csv")
df_ACC <- read_csv("./Data/Histology/Retinotopic_quantification/Clustering_individual/csv/118_ACC.csv")
df_LGN <- read_csv("./Data/Histology/Retinotopic_quantification/Clustering_individual/csv/118_LGN.csv")


# Plot settings -----------------------------------------------------------
res_RSC <- make_zoom_plots(
  df_RSC,
  xvar         = "Radius_[microns]",
  yvar         = "Normalized_Integrated_Intensity",
  main_xlim    = c(0, 250),
  main_ylim    = c(-0.1, 1),
  main_xbreaks = seq(0, 300, 100),
  main_ybreaks = seq(-0.5, 1,   0.5),
  inset_xlim   = c(40, 230),
  inset_ylim   = c(-0.05, 0.04),
  inset_xbreaks= seq(0, 300, 100),
  inset_ybreaks= seq(-0.1, 0.1, 0.02),
  inset_pos    = c(x = 0.25, y = 0.35, width = 0.75, height = 0.55)
)

RSC_plot <- ggdraw() +
  draw_plot(res_RSC$main_plot_rect) +
  draw_plot(res_RSC$inset_plot,
            x      = res$inset_position["x"],
            y      = res$inset_position["y"],
            width  = res$inset_position["width"],
            height = res$inset_position["height"])



res_ACC <- make_zoom_plots(
  df_ACC,
  xvar         = "Radius_[microns]",
  yvar         = "Normalized_Integrated_Intensity",
  main_xlim    = c(0, 250),
  main_ylim    = c(-0.1, 1),
  main_xbreaks = seq(0, 300, 100),
  main_ybreaks = seq(-0.5, 1,   0.5),
  inset_xlim   = c(15, 205),
  inset_ylim   = c(-0.05, 0.04),
  inset_xbreaks= seq(0, 300, 100),
  inset_ybreaks= seq(-0.1, 0.1, 0.02),
  inset_pos    = c(x = 0.25, y = 0.35, width = 0.75, height = 0.55)
)

ACC_plot <- ggdraw() +
  draw_plot(res_ACC$main_plot_rect) +
  draw_plot(res_ACC$inset_plot,
            x      = res$inset_position["x"],
            y      = res$inset_position["y"],
            width  = res$inset_position["width"],
            height = res$inset_position["height"])



res_LGN <- make_zoom_plots(
  df_LGN,
  xvar         = "Radius_[microns]",
  yvar         = "Normalized_Integrated_Intensity",
  main_xlim    = c(0, 250),
  main_ylim    = c(-0.1, 1),
  main_xbreaks = seq(0, 300, 100),
  main_ybreaks = seq(-0.5, 1,   0.5),
  inset_xlim   = c(20, 210),
  inset_ylim   = c(-0.05, 0.04),
  inset_xbreaks= seq(0, 300, 100),
  inset_ybreaks= seq(-0.1, 0.1, 0.02),
  inset_pos    = c(x = 0.25, y = 0.35, width = 0.75, height = 0.55)
)

LGN_plot <- ggdraw() +
  draw_plot(res_LGN$main_plot_rect) +
  draw_plot(res_LGN$inset_plot,
            x      = res$inset_position["x"],
            y      = res$inset_position["y"],
            width  = res$inset_position["width"],
            height = res$inset_position["height"])



# Show plots --------------------------------------------------------------------


plot(RSC_plot)
plot(ACC_plot)
plot(LGN_plot)


# Save plots --------------------------------------------------------------

ggsave("./Data/Histology/Retinotopic_quantification/Clustering_individual/118_RSC_radialprofile.svg",
       RSC_plot,
       width = 38, height = 32, units = "mm")
ggsave("./Data/Histology/Retinotopic_quantification/Clustering_individual/118_ACC_radialprofile.svg",
       ACC_plot,
       width = 38, height = 32, units = "mm")
ggsave("./Data/Histology/Retinotopic_quantification/Clustering_individual/118_LGN_radialprofile.svg",
       LGN_plot,
       width = 38, height = 32, units = "mm")
