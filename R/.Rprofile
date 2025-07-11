# if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
#   if (requireNamespace("httpgd", quietly = TRUE)) {
#     options(vsc.plot = FALSE)
#     options(device = function(...) {
#       httpgd::hgd(silent = TRUE)
#       .vsc.browser(httpgd::hgd_url(history = FALSE), viewer = "Beside")
#     })
#   }
# }

if (interactive() && Sys.getenv("TERM_PROGRAM") == "vscode") {
  if ("httpgd" %in% .packages(all.available = TRUE)) {
    options(vsc.plot = FALSE)
    options(device = function(...) {
      httpgd::hgd(silent = TRUE)
      .vsc.browser(httpgd::hgd_url(), viewer = "Beside")
    })
  }
}

source("renv/activate.R")

# Load required packages
library(data.table)
library(MASS)  # For negative binomial
library(pROC)  # For ROC curves
library(dplyr) # For pipe operators
library(ggplot2)
library(patchwork)