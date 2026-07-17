qmd_path <- qmd_path <- "C:/Users/sergasvi/Documents/Y2025.Wolbachia.CI/reports/02_Hatch_rate_analysis_batch8_9.qmd"

stopifnot(file.exists(qmd_path))

quarto::quarto_render(
  input = qmd_path,
  output_format = "html",
  execute_dir = dirname(qmd_path),
  quiet = FALSE
)
