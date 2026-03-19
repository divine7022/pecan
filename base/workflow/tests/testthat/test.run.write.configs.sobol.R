test_that("runModule.run.write.configs uses input_design row count", {
  settings <- PEcAn.settings::Settings(
    ensemble = list(
      size = 3,
      samplingspace = list(parameters = list(method = "uniform"))
    ),
    database = list(bety = list(write = FALSE)),
    pfts = list(list(posterior.files = "post.distns.Rdata"))
  )
  input_design <- data.frame(param = seq_len(5))
  captured <- new.env(parent = emptyenv())

  mockery::stub(
    runModule.run.write.configs,
    "PEcAn.workflow::run.write.configs",
    function(settings, ensemble.size, input_design, write, posterior.files, overwrite) {
      captured$ensemble.size <- ensemble.size
      captured$input_design <- input_design
      list(
        ensemble = list(ensemble.id = 123),
        pfts = settings$pfts
      )
    }
  )

  result <- runModule.run.write.configs(
    settings = settings,
    input_design = input_design
  )

  expect_equal(captured$ensemble.size, 5)
  expect_identical(captured$input_design, input_design)
  expect_equal(result$ensemble$ensemble.id, 123)
})
