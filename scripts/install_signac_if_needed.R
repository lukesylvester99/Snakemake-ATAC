#!/usr/bin/env Rscript
if (!requireNamespace("Signac", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("timoast/signac", upgrade = "never", dependencies = FALSE)
}
