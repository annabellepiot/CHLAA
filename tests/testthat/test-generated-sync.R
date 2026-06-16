test_that("bundled odin/dust artefacts are installed", {
  odin_dir <- system.file("odin", package = "chlaa")
  dust_dir <- system.file("dust", package = "chlaa")

  expect_true(nzchar(odin_dir))
  expect_true(nzchar(dust_dir))

  odin_files <- list.files(odin_dir, pattern = "\\.R$", full.names = TRUE)
  dust_files <- list.files(dust_dir, pattern = "\\.cpp$", full.names = TRUE)

  expect_gt(length(odin_files), 0)
  expect_gt(length(dust_files), 0)
})
