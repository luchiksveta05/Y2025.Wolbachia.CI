#' Model selection rows: model without vs with male age
#'
#' Returns df, logLik, AICc, delta AICc, Akaike weight and the likelihood-ratio
#' test (with-age vs without-age model) for one analysis.
#'
#' @param m0 glmmTMB model without MaleAge
#' @param m1 glmmTMB model with MaleAge
#' @param analysis label of the analysis (first column of the table)
#' @return data.frame with two rows (No age, Age)
#' @export
model_sel_age <- function(m0, m1, analysis) {
  aicc <- c(MuMIn::AICc(m0), MuMIn::AICc(m1))
  data.frame(Analysis = analysis,
             Model    = c("No age", "Age"),
             df       = c(attr(stats::logLik(m0), "df"), attr(stats::logLik(m1), "df")),
             logLik   = c(as.numeric(stats::logLik(m0)), as.numeric(stats::logLik(m1))),
             AICc     = aicc,
             dAICc    = aicc - min(aicc),
             Weight   = exp(-(aicc - min(aicc)) / 2) / sum(exp(-(aicc - min(aicc)) / 2))
             )
}
