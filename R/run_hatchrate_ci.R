#' Fit a zero-inflated binomial model
#' @description Test for differences in hatch rate among the four types of crosses using a zero-inflated binomial model
#' @param data database
#' @param uninfstrain uninfected strain
#' @param infstrain infected strain
#' @param nsim number of replicates for bootstrap)
#'
#' @return A vector or a dataframe with the required output (AICc table, confintp1, confinP2, etc.)
#'
#' @export
#'
#' @examples
#' #extract_results(id=1, outputper_infstrain=results, parameter="confintp2")
## Define function
run_hatchrate_ci <- function(data,
                             uninfstrain,
                             infstrain,
                             nsim,
                             model = c("m1", "m3"),
                             seed = NULL) {
  model <- match.arg(model)
 
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # 1) Subset
  dat <- data |>
    dplyr::filter(
      Female %in% c(uninfstrain, infstrain),
      Male %in% c(uninfstrain, infstrain))
    
  
  # 2) Fit model
  if (model == "m1") {
    m_fit <- glmmTMB::glmmTMB(
      hatchrate ~ -1 + CrossType + (1 | CoupleID),
      data = dat,
      ziformula = ~ 1,
      family = binomial,
      weight = sum_eggcount_corr
    )
  } else {
    m_fit <- glmmTMB::glmmTMB(
      hatchrate ~ -1 + CrossType + (1 | CoupleID),
      data = dat,
      ziformula = ~ -1 + CrossType,
      family = binomial,
      weight = sum_eggcount_corr
    )
  }
  
  # 3) Point estimates on probability scale (p2)
  p2 <- 1 / (1 + exp(-glmmTMB::fixef(m_fit)$cond))
  
  # 4) Bootstrap refit
  b1 <- lme4::bootMer(
    m_fit,
    FUN = extract_model,
    nsim = nsim,
    .progress = "txt"
  )
  
  # 5) CI table -> keep only "cond.*"
  confint <- bCI_tab(b = b1)
  
  out <- confint[grepl("^cond\\.", rownames(confint)), , drop = FALSE] |>
    as.data.frame() |>
    tibble::rownames_to_column("type") |>
    dplyr::mutate(continent = "europe", uninfstrain = uninfstrain, infstrain = infstrain) |>
    dplyr::select(continent, uninfstrain, infstrain, type, Estimate, `2.5 %`, `97.5 %`) |>
    dplyr::rename(X2.5. = `2.5 %`, X97.5. = `97.5 %`)
  
  attr(out, "model") <- m_fit
  attr(out, "p2") <- p2
  out
}

