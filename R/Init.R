#' @description Function that sets up the project

Init <- function() {
  # Load packages
  packages <- c("tseries", "tidyverse", "Rcpp", "RcppArmadillo", "numDeriv", "asbio", "matrixcalc")
  sapply(packages, Package_fctn)
  # Source functions
  RScripts <- paste0("R/", list.files(path = paste0(getwd(), "/R")))
  CScripts <- paste0("src/", list.files(path = paste0(getwd(), "/src")))
  sapply(RScripts, source)
  sapply(CScripts, Rcpp::sourceCpp)
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
