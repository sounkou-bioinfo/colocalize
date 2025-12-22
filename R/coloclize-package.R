## usethis namespace: start
#' @import coloc
#' @import seqminer
#' @import reticulate
#' @import rix
#' @import vroom
#' @import dplyr
#' @import ggplot2
#' @import data.table
#' @import rmarkdown
#' @import S7
#' @import vcfppR
#' @import RBCFLib
## usethis namespace: end
NULL


#' magma executable path
#' @return path to magma executable
#' @export
magma_exec <- function() {
    system.file("bin", "magma", "magma", package = "coloclize", mustWork = TRUE)
}

#' make pops python package venv
#' @param  environment_name name of the environment to create
#' @param environment_path path to create the environment in
#' @return list with environment name and path
#' @export
make_pops_venv <- function(environment_name = "coloclize_pops_env",
                           environment_path = tempdir()) {
    requirement_file <- system.file("python", "pops", "requirements.txt", package = "coloclize", mustWork = TRUE)
    venv_dir <- file.path(environment_path, environment_name)
    reticulate::virtualenv_create(
        envname = environment_name,
        python = NULL,
        path = environment_path
    )
    reticulate::virtualenv_install(
        envname = environment_name,
        packages = requirement_file,
        ignore_installed = TRUE,
        path = environment_path
    )
    message(paste0("Virtual environment created at: ", venv_dir))
    invisible(list(
        environment_name = environment_name,
        environment_path = venv_dir
    ))
}

#' make FLAMES python package venv
#' @param  environment_name name of the environment to create
#' @param environment_path path to create the environment in
#' @return list with environment name and path
#' @export
make_flames_venv <- function(environment_name = "coloclize_flames_env",
                             environment_path = tempdir()) {
    requirement_file <- system.file("python", "FLAMES", "requirements.txt", package = "coloclize", mustWork = TRUE)
    venv_dir <- file.path(environment_path, environment_name)
    reticulate::virtualenv_create(
        envname = environment_name,
        python = NULL,
        path = environment_path
    )
    reticulate::virtualenv_install(
        envname = environment_name,
        packages = requirement_file,
        ignore_installed = TRUE,
        path = environment_path
    )
    message(paste0("Virtual environment created at: ", venv_dir))
    invisible(list(
        environment_name = environment_name,
        environment_path = venv_dir
    ))
}
