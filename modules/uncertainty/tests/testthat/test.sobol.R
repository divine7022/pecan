make_sobol_settings <- function(outdir) {
  samples_file <- file.path(outdir, "samples.Rdata")
  trait.samples <- list(
    temperate = list(
      SLA = seq_len(40)
    )
  )
  save(trait.samples, file = samples_file)

  PEcAn.settings::Settings(
    outdir = outdir,
    pfts = list(list(name = "temperate", posterior.files = "post.distns.Rdata")),
    run = list(inputs = list(
      met = list(path = c("met1", "met2", "met3")),
      poolinitcond = list(path = c("ic1", "ic2")),
      events = list(path = c("evt1", "evt2", "evt3"))
    )),
    ensemble = list(
      samplingspace = list(
        parameters = list(method = "uniform"),
        met = list(method = "sampling"),
        poolinitcond = list(method = "looping"),
        events = list(method = "sampling", parent = "met")
      )
    )
  )
}

test_that("Sobol design expands to N * (k + 2) rows with metadata", {
  withr::with_tempdir({
    settings <- make_sobol_settings(getwd())

    result <- generate_joint_ensemble_design(
      settings = settings,
      ensemble_size = 4,
      sobol = TRUE
    )

    expect_equal(nrow(result$X), 20)
    expect_equal(result$N, 4)
    expect_identical(result$params, c("param", "met", "poolinitcond"))
    expect_identical(result$backend, "sensobol")
    expect_identical(result$matrices, c("A", "B", "AB"))
    expect_identical(result$first, "saltelli")
    expect_identical(result$total, "jansen")
    expect_true(all(c("param", "met", "poolinitcond", "events") %in% names(result$X)))
    expect_equal(result$X$events, result$X$met)
    expect_true(all(result$X$param >= 1))
    expect_true(all(result$X$param <= 20))
    expect_true(all(c("factor", "source_type", "source_tag") %in% names(result$factor_metadata)))
    expect_identical(result$factor_metadata$source_type, c("param", "met", "poolinitcond"))
  })
})

test_that("Non-Sobol design generation remains row-for-row", {
  withr::with_tempdir({
    settings <- make_sobol_settings(getwd())

    result <- generate_joint_ensemble_design(
      settings = settings,
      ensemble_size = 5,
      sobol = FALSE
    )

    expect_named(result, "X")
    expect_equal(nrow(result$X), 5)
    expect_false("backend" %in% names(result))
  })
})

test_that("compute_sobol_indices matches direct sensobol results", {
  withr::with_tempdir({
    sobol_obj <- list(
      N = 4L,
      params = c("param", "met"),
      backend = "sensobol",
      matrices = c("A", "B", "AB"),
      first = "saltelli",
      total = "jansen",
      factor_metadata = data.frame(
        factor = c("param", "met"),
        source_type = c("param", "met"),
        source_tag = c(NA_character_, "met"),
        stringsAsFactors = FALSE
      )
    )

    y <- seq_len(16)
    ensemble.output <- as.list(y)
    save(ensemble.output, file = file.path(
      getwd(),
      "ensemble.output.NOENSEMBLEID.GPP.NA.NA.Rdata"
    ))

    result <- compute_sobol_indices(
      outdir = getwd(),
      sobol_obj = sobol_obj,
      var = "GPP"
    )

    expected <- tibble::as_tibble(
      sensobol::sobol_indices(
        matrices = c("A", "B", "AB"),
        Y = y,
        N = 4L,
        params = c("param", "met"),
        first = "saltelli",
        total = "jansen",
        order = "first",
        boot = FALSE
      )$results
    )

    expect_equal(
      result[, c("parameters", "sensitivity", "original")],
      expected[, c("parameters", "sensitivity", "original")]
    )
    expect_true(file.exists(file.path(
      getwd(),
      "sobol.design.NOENSEMBLEID.GPP.NA.NA.Rdata"
    )))
    expect_true(file.exists(file.path(
      getwd(),
      "sobol.indices.NOENSEMBLEID.GPP.NA.NA.Rdata"
    )))
  })
})
