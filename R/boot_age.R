#' Parametric bootstrap of the male age model
#'
#' Simulates data sets from a zero-inflated beta-binomial model including
#' MaleAge, refits the model to each, and returns percentile 95% CIs for the
#' MaleAge slope and for the mean hatch rate of each observed line x age cell.
#' Refits without a positive-definite Hessian are discarded; estimates are not
#' truncated.
#'
#' @param m_age glmmTMB model including MaleAge (hatchrate ~ -1 + Male + MaleAge)
#' @param dat data used to fit m_age
#' @param n_boot number of converged refits to keep (default 1000)
#' @return list with
#'   pred: hatch rate by line x age (n, est, lo, hi);
#'   slope: MaleAge slope (n, est, lo, hi, n_ok)
#' @export
boot_age <- function(m_age, dat, n_boot = 500) {
  nd <- dat |>
    dplyr::count(Male, CrossType, MaleAge) |>
    dplyr::mutate(sum_eggcount_corr = 1)

  sims <- stats::simulate(m_age, nsim = n_boot + 100, seed = 1)  # +100 spare for failed refits
  slope <- c()
  pred  <- NULL

  for (i in seq_along(sims)) {
    d <- dat
    y <- sims[[i]]
    d$hatchrate <- if (is.matrix(y)) y[, 1] / rowSums(y) else y
    fit <- try(stats::update(m_age, data = d), silent = TRUE)
    if (inherits(fit, "try-error") || !isTRUE(fit$sdr$pdHess)) next
    slope <- c(slope, glmmTMB::fixef(fit)$cond[["MaleAge"]])
    pred  <- rbind(pred, stats::predict(fit, newdata = nd, type = "response"))
    if (length(slope) == n_boot) break
  }

  if (length(slope) < n_boot) warning("Fewer than n_boot refits converged")

  nd$est <- stats::predict(m_age, newdata = nd, type = "response")
  nd$lo  <- apply(pred, 2, stats::quantile, 0.025)
  nd$hi  <- apply(pred, 2, stats::quantile, 0.975)

  list(pred = nd,
       slope = data.frame(n    = nrow(dat),
                          est  = glmmTMB::fixef(m_age)$cond[["MaleAge"]],
                          lo   = stats::quantile(slope, 0.025, names = FALSE),
                          hi   = stats::quantile(slope, 0.975, names = FALSE),
                          n_ok = length(slope)))
}
