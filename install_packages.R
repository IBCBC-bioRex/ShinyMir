# install_packages.R
pkgs <- c(
  # UI framework
  "shiny", "shinydashboard", "shinyjs", "shinycssloaders",
  # Database
  "DBI", "duckdb",
  # Data manipulation
  "dplyr", "tidyr", "purrr", "glue", "scales",
  # Tables & charts
  "DT", "plotly", "visNetwork", "heatmaply", "igraph",
  # Visuals & report
  "RColorBrewer", "htmlwidgets", "htmltools", "kableExtra",
  # File I/O
  "readxl"
)

missing <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(missing) == 0) {
  message("All packages already installed.")
} else {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)
  message("Done.")
}
