#' @keywords internal
#' 타깃(to) 하나에 대해 후보 엣지들의 공통 관측치로 정렬 후 LOO 가중치 계산
.edge_weights_by_target <- function(loo_items, method = c("stacking","pseudobma","pseudobma_bb")) {
  method <- match.arg(method)
  # 1) 공통 관측치 교집합
  keys_list <- lapply(loo_items, `[[`, "obs_key")
  common    <- Reduce(intersect, keys_list)
  if (!length(common)) return(NULL)

  # 2) 각 모델의 log_lik을 공통 관측치 순서로 서브셋
  re_loo <- lapply(loo_items, function(it){
    if (is.null(it$log_lik) || is.null(it$obs_key)) return(NULL)
    ord <- match(common, it$obs_key)
    if (anyNA(ord)) return(NULL)
    ll_sub <- it$log_lik[, ord, drop = FALSE]
    r_eff  <- try(loo::relative_eff(exp(ll_sub)), silent = TRUE)
    if (inherits(r_eff, "try-error")) r_eff <- NULL
    loo::loo(ll_sub, r_eff = r_eff, save_psis = TRUE)
  })
  re_loo <- Filter(Negate(is.null), re_loo)
  if (length(re_loo) < 2L) return(NULL)

  # 3) 가중치
  w <- loo::loo_model_weights(re_loo, method = method)
  data.frame(
    from = vapply(loo_items, `[[`, "", "from")[seq_along(re_loo)],
    to   = vapply(loo_items, `[[`, "", "to")[seq_along(re_loo)],
    method = method,
    weight = as.numeric(w),
    stringsAsFactors = FALSE
  )
}
