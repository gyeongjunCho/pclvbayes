test_that("package loads", {
  expect_true("pglvbayes" %in% loadedNamespaces())
})

test_that("key exports exist", {
  expect_true(all(c(
    "fit_glv_pairwise",
    "cor_meta_resid",
    "compute_edge_weights_loo",
    "summarize_bayes_pglv"
  ) %in% getNamespaceExports("pglvbayes")))
})

test_that("options are sane (no CRAN-only heavy work)", {
  # 무거운 샘플링/외부 툴 필요 테스트는 여기서 하지 않음
  skip_on_cran()
  skip_if_not_installed("progressr")
  expect_true(is.character(getNamespaceName(asNamespace("pglvbayes"))))
})
