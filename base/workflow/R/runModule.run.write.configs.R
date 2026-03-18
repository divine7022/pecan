#' Generate model-specific run configuration files for one or more PEcAn runs
#'
#' @param settings a PEcAn Settings or MultiSettings object
#' @param overwrite logical: Replace config files if they already exist?
#' @param input_design data.frame design matrix linking parameter draws and any
#'   sampled inputs across runs. Must include a `param` column selecting rows
#'   from `trait.samples`/`ensemble.samples`; may also include columns named for
#'   `settings$run$inputs` tags (e.g. `met`, `soil`) giving row indices into
#'   each input's `path` list. If supplied, `input_design` defines the runs to
#'   execute and the run count is `nrow(input_design)`, overriding
#'   `settings$ensemble$size`. If NULL, the design is generated from
#'   `settings$ensemble$size` using `generate_joint_ensemble_design()`.
#' @return A modified settings object, invisibly
#' @importFrom dplyr %>%
#' @export


runModule.run.write.configs <- function(settings,
                                        overwrite = TRUE,
                                        input_design = NULL) {
  if (PEcAn.settings::is.MultiSettings(settings)) {
    if (overwrite && file.exists(file.path(settings$rundir, "runs.txt"))) {
      PEcAn.logger::logger.warn("Existing runs.txt file will be removed.")
      unlink(file.path(settings$rundir, "runs.txt"))
    }
    if (is.null(input_design)) {
      ensemble_size <- settings$ensemble$size
      design_result <- PEcAn.uncertainty::generate_joint_ensemble_design(
        settings = settings[1],
        ensemble_size = ensemble_size
      )
      input_design <- design_result$X
    } else if (!is.null(settings$ensemble$size) &&
               settings$ensemble$size != nrow(input_design)) {
      PEcAn.logger::logger.warn(
        "Ignoring settings$ensemble$size = ", settings$ensemble$size,
        "; using nrow(input_design) = ", nrow(input_design), "."
      )
    }
    settings$ensemble$size <- nrow(input_design)

    return(PEcAn.settings::papply(settings,
                                  runModule.run.write.configs,
                                  overwrite = FALSE,
                                  input_design = input_design))
  } else if (PEcAn.settings::is.Settings(settings)) {
    # double check making sure we have method for parameter sampling
    if (is.null(settings$ensemble$samplingspace$parameters$method)) {
      settings$ensemble$samplingspace$parameters$method <- "uniform"
    }
    if (is.null(input_design)) {
      ensemble_size <- settings$ensemble$size
      design_result <- PEcAn.uncertainty::generate_joint_ensemble_design(
        settings = settings,
        ensemble_size = ensemble_size
      )
      input_design <- design_result$X
    } else if (!is.null(settings$ensemble$size) &&
               settings$ensemble$size != nrow(input_design)) {
      PEcAn.logger::logger.warn(
        "Ignoring settings$ensemble$size = ", settings$ensemble$size,
        "; using nrow(input_design) = ", nrow(input_design), "."
      )
    }
    ensemble_size <- nrow(input_design)
    settings$ensemble$size <- ensemble_size

    # check to see if there are posterior.files tags under pft
    posterior.files <- settings$pfts %>%
      purrr::map_chr("posterior.files", .default = NA_character_)

    return(PEcAn.workflow::run.write.configs(
      settings = settings,
      ensemble.size = ensemble_size,
      write = isTRUE(settings$database$bety$write), # treat null as FALSE
      posterior.files = posterior.files,
      overwrite = overwrite,
      input_design = input_design
    ))
  } else {
    PEcAn.logger::logger.error(
      "runModule.run.write.configs only works with Settings or MultiSettings"
    )
  }
}
