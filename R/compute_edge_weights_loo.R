#' @export
compute_edge_weights_loo <- function(res,
                                     method = c("stacking","pseudobma","pseudobma_bb"),
                                     sign_with = c("psp","mean_sign"),
                                     require_ok_diag = TRUE,
                                     ebfmi_min = 0.3) {
  method    <- match.arg(method)
  sign_with <- match.arg(sign_with)

  ls <- attr(res, "loo_store")
  if (is.null(ls)) stop("No 'loo_store' found in result. Run with save_loo=TRUE.")

  # 진단 필터를 메타 테이블로 만들어서 적용
  meta <- do.call(rbind, lapply(names(ls), function(k){
    it <- ls[[k]]
    dir <- it$dir; to <- it$to; from <- it$from
    if (dir == "j->i") {
      row <- subset(res, i == to & j == from)
      data.frame(key = k, to = to, from = from,
                 ok_diag = isTRUE(row$ok_diag_ij[1]),
                 ebfmi   = row$ebfmi_min_ij[1],
                 psp     = row$psp_ij[1])
    } else {
      row <- subset(res, i == from & j == to)
      data.frame(key = k, to = to, from = from,
                 ok_diag = isTRUE(row$ok_diag_ji[1]),
                 ebfmi   = row$ebfmi_min_ji[1],
                 psp     = row$psp_ji[1])
    }
  }))
  # 필터
  keep <- rep(TRUE, nrow(meta))
  if (require_ok_diag) keep <- keep & meta$ok_diag
  if (is.finite(ebfmi_min)) keep <- keep & (!is.finite(meta$ebfmi) | meta$ebfmi >= ebfmi_min)
  meta <- meta[keep, , drop = FALSE]

  # 타깃별로 가중치 계산
  by_to <- split(meta, meta$to)
  out <- lapply(by_to, function(m){
    items <- lapply(m$key, function(k) ls[[k]])
    .edge_weights_by_target(items, method = method)
  })
  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (is.null(out) || !nrow(out)) return(out)

  # 부호 반영 (기본: PSP)
  if (sign_with == "psp") {
    out$sign_score <- vapply(seq_len(nrow(out)), function(i){
      if (out$to[i] == res$i & out$from[i] == res$j) {
        NA_real_
      }
      # 안전하게 조회
      r1 <- subset(res, i == out$to[i] & j == out$from[i])
      r2 <- subset(res, i == out$from[i] & j == out$to[i])
      if (nrow(r1)) 2 * r1$psp_ij[1] - 1 else if (nrow(r2)) 2 * r2$psp_ji[1] - 1 else NA_real_
    }, numeric(1))
  } else {
    # 대안: 후방 평균 부호 사용
    out$sign_score <- vapply(seq_len(nrow(out)), function(i){
      r1 <- subset(res, i == out$to[i] & j == out$from[i])
      if (nrow(r1)) sign(r1$a_ij_mean[1])
      else {
        r2 <- subset(res, i == out$from[i] & j == out$to[i])
        if (nrow(r2)) sign(r2$a_ji_mean[1]) else NA_real_
      }
    }, numeric(1))
  }
  out$signed_weight <- out$weight * out$sign_score
  out
}
