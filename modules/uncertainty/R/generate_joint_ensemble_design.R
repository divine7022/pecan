.trait_sample_bank_size <- function(trait.samples) {
  if (is.null(trait.samples) || length(trait.samples) == 0) {
    return(0L)
  }

  bank_sizes <- unlist(
    lapply(trait.samples, function(pft_traits) {
      if (is.null(pft_traits) || length(pft_traits) == 0) {
        return(integer(0))
      }

      vapply(
        pft_traits,
        function(trait_values) {
          if (is.null(trait_values) || length(trait_values) == 0) {
            return(NA_integer_)
          }

          as.integer(length(trait_values))
        },
        integer(1)
      )
    }),
    use.names = FALSE
  )

  bank_sizes <- bank_sizes[!is.na(bank_sizes) & bank_sizes > 0L]
  if (length(bank_sizes) == 0) {
    return(0L)
  }

  as.integer(min(bank_sizes))
}

.sobol_parameter_bank_size <- function(samples_file) {
  if (!file.exists(samples_file)) {
    return(0L)
  }

  samples <- new.env(parent = emptyenv())
  load(samples_file, envir = samples)

  if (is.null(samples$trait.samples) || length(samples$trait.samples) == 0) {
    return(0L)
  }

  .trait_sample_bank_size(samples$trait.samples)
}

.map_sobol_to_indices <- function(x, size) {
  indices <- floor(stats::qunif(x, min = 1, max = size + 1))
  as.integer(pmin(indices, size))
}

#' Generate joint ensemble design for parameter sampling
#' Creates a joint ensemble design that maintains parameter correlations across
#' all sites in a multi-site run. This function generates sample indices that
#' are shared across sites to ensure consistent parameter sampling.
#'
#' @param settings A PEcAn settings object containing ensemble configuration.
#' @param ensemble_size Integer specifying the number of ensemble members.
#'   When \code{sobol = TRUE}, this is the Sobol base sample size \code{N}, not
#'   the expanded number of model runs.
#' @param sobol Logical, generate a variance-based Sobol design using
#'   \code{sensobol}.
#'
#' @return A list with component \code{X}, a data frame design matrix
#'   describing PEcAn parameter and sampled-input indices. If \code{sobol = TRUE},
#'   the list also includes the metadata needed by
#'   \code{\link{compute_sobol_indices}}.
#'
#' @export

generate_joint_ensemble_design <- function(settings,
                                           ensemble_size,
                                           sobol = FALSE) {
  ens.sample.method <- settings$ensemble$samplingspace$parameters$method
  design_list <- list()
  sampled_inputs <- list()
  posterior.files <- settings$pfts %>%
    purrr::map_chr("posterior.files", .default = NA_character_)
  samp <- settings$ensemble$samplingspace
  parents <- lapply(samp, "[[", "parent")
  order <- names(samp)[
    lapply(parents, function(tr) which(names(samp) %in% tr)) %>%
      unlist()
  ]
  samp.ordered <- samp[c(order, names(samp)[!(names(samp) %in% order)])]

  if (sobol) {
    sobol_factors <- c(
      "param",
      names(samp.ordered)[
        names(samp.ordered) != "parameters" &
          vapply(
            samp.ordered,
            function(x) is.null(x$parent),
            logical(1)
          )
      ]
    )

    total_runs <- as.integer(ensemble_size) * (length(sobol_factors) + 2L)
    samples_file <- file.path(settings$outdir, "samples.Rdata")
    if (.sobol_parameter_bank_size(samples_file) < total_runs) {
      PEcAn.uncertainty::get.parameter.samples(
        settings = settings,
        ensemble.size = total_runs,
        posterior.files = posterior.files,
        ens.sample.method = ens.sample.method
      )
    }

    sobol_design <- sensobol::sobol_matrices(
      matrices = c("A", "B", "AB"),
      N = as.integer(ensemble_size),
      params = sobol_factors,
      order = "first",
      type = "QRN"
    )
    sobol_design <- as.data.frame(sobol_design)

    sobol_indices <- list()
    sobol_indices[["param"]] <- .map_sobol_to_indices(
      sobol_design[["param"]],
      total_runs
    )
    sampled_inputs[["parameters"]] <- list(ids = sobol_indices[["param"]])

    for (input_tag in setdiff(sobol_factors, "param")) {
      input_paths <- settings$run$inputs[[tolower(input_tag)]]$path
      if (is.null(input_paths) || length(input_paths) == 0) {
        PEcAn.logger::logger.error("Input ", sQuote(input_tag), " has no paths specified")
      }
      sobol_indices[[input_tag]] <- .map_sobol_to_indices(
        sobol_design[[input_tag]],
        length(input_paths)
      )
    }

    for (i in seq_along(samp.ordered)) {
      input_tag <- names(samp.ordered)[i]

      if (identical(input_tag, "parameters")) {
        next
      }

      parent_name <- samp.ordered[[i]]$parent
      if (!is.null(parent_name)) {
        input_result <- PEcAn.uncertainty::input.ens.gen(
          settings = settings,
          ensemble_size = total_runs,
          input = input_tag,
          method = samp.ordered[[i]]$method,
          parent_ids = sampled_inputs[[parent_name]]
        )
        sampled_inputs[[input_tag]] <- input_result
        design_list[[input_tag]] <- input_result$ids
      } else if (input_tag %in% names(sobol_indices)) {
        sampled_inputs[[input_tag]] <- list(ids = sobol_indices[[input_tag]])
        design_list[[input_tag]] <- sobol_indices[[input_tag]]
      }
    }

    design_list[["param"]] <- sobol_indices[["param"]]
    design_matrix <- data.frame(design_list)

    factor_metadata <- data.frame(
      factor = sobol_factors,
      source_type = sobol_factors,
      source_tag = ifelse(sobol_factors == "param", NA_character_, sobol_factors),
      stringsAsFactors = FALSE
    )

    return(list(
      X = design_matrix,
      N = as.integer(ensemble_size),
      params = sobol_factors,
      backend = "sensobol",
      matrices = c("A", "B", "AB"),
      first = "saltelli",
      total = "jansen",
      factor_metadata = factor_metadata
    ))
  }

  sampled_inputs[["parameters"]] <- list(ids = seq_len(ensemble_size))
  for (i in seq_along(samp.ordered)) {
    input_tag <- names(samp.ordered)[i]
    if (identical(input_tag, "parameters")) {
      next
    }

    parent_name <- samp.ordered[[i]]$parent
    parent_ids <- if (!is.null(parent_name)) {
      sampled_inputs[[parent_name]]
    } else {
      NULL
    }

    input_result <- PEcAn.uncertainty::input.ens.gen(
      settings = settings,
      ensemble_size = ensemble_size,
      input = input_tag,
      method = samp.ordered[[i]]$method,
      parent_ids = parent_ids
    )

    sampled_inputs[[input_tag]] <- input_result
    design_list[[input_tag]] <- input_result$ids
  }

  if (!file.exists(file.path(settings$outdir, "samples.Rdata"))) {
    PEcAn.uncertainty::get.parameter.samples(
      settings,
      ensemble.size = ensemble_size,
      posterior.files,
      ens.sample.method
    )
  }

  design_list[["param"]] <- sampled_inputs[["parameters"]]$ids
  design_matrix <- data.frame(design_list)
  return(list(X = design_matrix))
}
