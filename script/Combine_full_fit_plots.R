library(ggplot2)
library(cowplot)

# Run from the `script` folder.
sym_env <- new.env(parent = globalenv())
sys.source(file.path("Symptomatic", "Model_Efficacy_Nab_Symptomatic.R"), envir = sym_env)

asym_env <- new.env(parent = globalenv())
sys.source(file.path("Asymptomatic", "Model_Efficacy_Nab_Asymptomatic.R"), envir = asym_env)

panel_a <- sym_env$Full_fit_sym +
  theme(plot.margin = margin(t = 14, r = 8, b = 6, l = 24))

panel_b <- asym_env$Full_fit_asym +
  theme(plot.margin = margin(t = 14, r = 24, b = 6, l = 8))

Combined_full_fit_sym_asym <- cowplot::plot_grid(
  panel_a,
  NULL,
  panel_b,
  labels = c("A.", "", "B."),
  label_size = 18,
  label_fontface = "bold",
  label_x = c(0.004, 0, 0.004),
  label_y = c(0.99, 0.99, 0.99),
  hjust = 0,
  vjust = 1,
  ncol = 3,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 0.08, 1)
)

output_dir <- "Output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(output_dir, "Full_fit_sym_asym_combined.pdf"),
  plot = Combined_full_fit_sym_asym,
  width = 11,
  height = 5.1,
  units = "in"
)

ggsave(
  filename = file.path(output_dir, "Full_fit_sym_asym_combined.png"),
  plot = Combined_full_fit_sym_asym,
  width = 11,
  height = 5.1,
  units = "in",
  dpi = 300
)

message("Saved combined full-fit plots to: ", output_dir)
