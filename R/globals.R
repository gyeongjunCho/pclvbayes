#' @keywords internal
NULL

utils::globalVariables(c(
  "Sample","subject","time","xi","xj","xi_next","time_next",
  "n_pairs_ij","a_ij_mean","a_ij_sd","a_ij_q2.5","a_ij_q97.5","p_ij","q_ij",
  "rhat_ij","essb_ij","esst_ij","div_ij","tdhit_ij","ok_diag_ij","keep_ij_final",
  "n_pairs_ji","a_ji_mean","a_ji_sd","a_ji_q2.5","a_ji_q97.5","p_ji","q_ji",
  "rhat_ji","essb_ji","esst_ji","div_ji","tdhit_ji","ok_diag_ji","keep_ji_final",
  "xi_lag","xj_lag","alr_i_lag","time_lag","dt","dt_min","dt_adj",
  "xi_raw","xj_raw","rest_raw",
  ".data","xj_next","rest_now","rest_next","alr_i_next","alr_i_now","y",
  "i","j","k",
  "q_bayes","p_sign2","r_pooled",
  "taxa_vec","progress","kfold_K","kfold_R"
))
