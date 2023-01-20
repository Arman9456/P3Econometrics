#' @description Function that sets up the project

Init <- function() {
  # Load packages
  packages <- c("matlib", "expm", "tidyverse")
  sapply(packages, Package_fctn)
  # Source functions
  Scripts <- paste0("R/", list.files(path = paste0(getwd(), "/R")))
  sapply(Scripts, source)
}

#' @description Function that loads required packages or installs them if necessary
#' @param string of a package name

Package_fctn <- function(pckg) {
  if (!require(pckg, character.only = TRUE)) {
    cat("Installing required packages \n")
    install.packages(pckg, dep = TRUE, quiet = T)
  }
  require(pckg, character.only = TRUE, quietly = T)
}
