#' Compute Sobol indices from a finished PEcAn run
#'
#' Loads standardized ensemble output for a Sobol run, computes first-order and
#' total-order Sobol indices with `sensobol`, and saves both the Sobol design
#' metadata and the computed indices using PEcAn-style filenames.
#'
#' First-order indices quantify the share of output variance attributable to a
#' parameter alone, while total-order indices summarize the full contribution of
#' that parameter including its interactions with other parameters. PEcAn uses
#' `sensobol` for these variance-based estimators; readers wanting estimator
#' details or broader background should refer to Puy et al. (2022).
#'
#' @param outdir PEcAn run output directory containing `ensemble.output.*.Rdata`
#'   files.
#' @param sobol_obj object produced by
#'   `PEcAn.uncertainty::generate_joint_ensemble_design(..., sobol = TRUE)`.
#' @param var Variable name to summarize (default `"GPP"`).
#' @param stat_fun Summary statistic applied to `var`. Retained for backwards
#'   compatibility; standardized ensemble outputs are already scalar summaries.
#'
#' @return A tibble of Sobol first-order and total-order indices with attached
#'   factor metadata.
#' @references Puy, A., Lo Piano, S., Saltelli, A., and Levin, S. A. (2022).
#'   sensobol: An R Package to Compute Variance-Based Sensitivity Indices.
#'   Journal of Statistical Software, 102(5), 1-37.
#'   \doi{10.18637/jss.v102.i05}
#' @export
compute_sobol_indices <- function(outdir,
                                  sobol_obj,
                                  var = "GPP",
                                  stat_fun = mean) {
  if (is.null(sobol_obj$backend) || sobol_obj$backend != "sensobol") {
    PEcAn.logger::logger.error(
      "compute_sobol_indices expects a sensobol design object returned by ",
      "generate_joint_ensemble_design(..., sobol = TRUE)"
    )
  }

  output_files <- list.files(
    outdir,
    pattern = "^ensemble\\.output\\..*\\.Rdata$",
    full.names = TRUE
  )
  if (length(output_files) == 0) {
    PEcAn.logger::logger.error("No ensemble.output.*.Rdata files found in ", outdir)
  }

  output_var <- vapply(
    strsplit(basename(output_files), "\\."),
    function(x) if (length(x) >= 4) x[[4]] else NA_character_,
    character(1)
  )
  matched_files <- output_files[output_var == var]

  if (length(matched_files) == 0) {
    PEcAn.logger::logger.error(
      "No standardized ensemble output found for variable '", var,
      "' in ", outdir
    )
  }
  if (length(matched_files) > 1) {
    PEcAn.logger::logger.error(
      "Multiple standardized ensemble outputs found for variable '", var,
      "' in ", outdir, ". Please keep only one matching file."
    )
  }

  output_file <- matched_files[[1]]
  output_env <- new.env(parent = emptyenv())
  load(output_file, envir = output_env)
  if (is.null(output_env$ensemble.output)) {
    PEcAn.logger::logger.error(
      "Object `ensemble.output` missing from standardized output file ",
      output_file
    )
  }

  y <- as.numeric(unlist(output_env$ensemble.output, use.names = FALSE))
  expected_length <- sobol_obj$N * (length(sobol_obj$params) + 2L)
  if (length(y) != expected_length) {
    PEcAn.logger::logger.error(
      "Standardized ensemble output has ", length(y),
      " values but expected ", expected_length,
      " for Sobol design size N = ", sobol_obj$N
    )
  }

  sobol_indices_result <- sensobol::sobol_indices(
    matrices = sobol_obj$matrices,
    Y = y,
    N = sobol_obj$N,
    params = sobol_obj$params,
    first = sobol_obj$first,
    total = sobol_obj$total,
    order = "first",
    boot = FALSE
  )

  sobol_results <- tibble::as_tibble(sobol_indices_result$results)
  if (!is.null(sobol_obj$factor_metadata)) {
    factor_metadata <- tibble::as_tibble(sobol_obj$factor_metadata) |>
      dplyr::rename(parameters = .data$factor)
    sobol_results <- dplyr::left_join(
      sobol_results,
      factor_metadata,
      by = "parameters"
    )
  }

  sobol_design <- sobol_obj
  save(
    sobol_design,
    file = file.path(
      outdir,
      sub("^ensemble\\.output\\.", "sobol.design.", basename(output_file))
    )
  )
  save(
    sobol_indices_result,
    sobol_results,
    file = file.path(
      outdir,
      sub("^ensemble\\.output\\.", "sobol.indices.", basename(output_file))
    )
  )

  return(sobol_results)
}
