#' Rescue assay: zero-inflated beta-binomial models for one male line
#'
#' @description Same as fit_ZI_Beta_binomial, for a design where the male line is fixed
#' and the female line varies (uninfected Am females = test cross,
#' infected females = control crosses). Returns the same list as fit_ZI_Beta_binomial,
#' so extract_results() can be used unchanged.
#' Requires extract_model() and bCI_tab().
#'
#' @param crosses vector: continent, uninfected strain, infected male line, e.g. c("europe", "Am", "Bal")
#' @param data data already filtered to the rescue blocks and crosses
fit_rescue <- function(crosses = c("europe", "Am", "Bal"), data, nbootsim, dec_places = 3,
                       saturated_model_estimates = NULL, forms = NULL,
                       zis = list(ziCT = ~ -1 + CrossType,
                                  zi1  = ~ 1),
                       disps = list(dispCT = ~ -1 + CrossType,
                                    disp1  = ~ 1)) {

  ## Extract the name of the uninfected strain and of the infected male line
  strains <- crosses[2:3]
  test    <- paste(strains, collapse = "x")

  print(paste("Testing for rescue in crosses with", strains[2], "males"))

  ## Subset the data: all crosses with males of the infected line
  dat <- data[data$Male == strains[2], ]
  block <- unique(dat$Batch)

  print(paste("Use data from the following blocks:", paste(sort(block), collapse = " ")))

  ## Compute proportion of eggs that hatch into larvae
  dat$prop <- dat$sum_larvacount / dat$sum_eggcount_corr

  ## Test cross (Am females) vs control crosses (infected females)
  dat$CI_test <- factor(ifelse(dat$CrossType == test, "Test_cross", "Control_cross"))

  ## Order of the crosses: test cross first
  dat$CrossType <- factor(dat$CrossType,
                          levels = c(test, setdiff(sort(unique(dat$CrossType)), test)))

  datCrossType <- dat |>
    dplyr::group_by(CrossType, Batch, Year) |>
    dplyr::summarise(Sum_eggcount_corr = sum(sum_eggcount_corr),
                     Sum_larvacount = sum(sum_larvacount)) |>
    dplyr::mutate(hatchrate = Sum_larvacount / Sum_eggcount_corr)

  print(datCrossType)

  if (is.null(forms)) {
    ## List all models
    forms <- list(
      m0 = prop ~ 1,                 ## Null model
      m1 = prop ~ CI_test,           ## Test vs control crosses
      m2 = prop ~ -1 + CrossType)    ## Saturated model
  }

  ## Fit one model, storing literal formulas and data in the call (needed by bootMer)
  fit_one <- function(f, zi, disp) {
    m <- glmmTMB::glmmTMB(f, ziformula = zi, dispformula = disp,
                          family = glmmTMB::betabinomial(),
                          weights = sum_eggcount_corr, data = dat)
    m$call$formula     <- f
    m$call$ziformula   <- zi
    m$call$dispformula <- disp
    m$call$data        <- dat
    m
  }

  ## All combinations
  grid <- expand.grid(cond = names(forms), zi = names(zis), disp = names(disps),
                      stringsAsFactors = FALSE)
  ## Fit all combinations
  mods <- lapply(seq_len(nrow(grid)), function(i)
    fit_one(forms[[grid$cond[i]]], zis[[grid$zi[i]]], disps[[grid$disp[i]]]))
  names(mods) <- paste(grid$cond, grid$zi, grid$disp, sep = "_")

  ## Compare models based on AICc
  modelaicc <- data.frame(MuMIn::model.sel(mods))
  modelaicc <- modelaicc |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, dec_places)))

  ## Check convergence
  conv <- sapply(mods, \(x) x$fit$convergence)
  pd   <- sapply(mods, \(x) x$sdr$pdHess)

  labels_cond <- c(m0 = "null", m1 = "CI (test vs control)", m2 = "saturated (cross type)")
  labels_zi   <- c(ziCT = "zi ~ cross type", zi1 = "zi constant")
  labels_disp <- c(dispCT = "disp ~ cross type", disp1 = "disp constant")
  parts <- do.call(rbind, strsplit(rownames(modelaicc), "_"))

  modelaicc <- data.frame(model = rownames(modelaicc),
                          cond  = labels_cond[parts[, 1]],
                          zi    = labels_zi[parts[, 2]],
                          disp  = labels_disp[parts[, 3]],
                          modelaicc,
                          convergence = ifelse(conv[rownames(modelaicc)] == 0, "converged", "not_converged"),
                          Hessian     = ifelse(pd[rownames(modelaicc)], "positive_definite", "non_positive_definite"),
                          row.names = rownames(modelaicc))

  ## Best model: the user's choice, or the best-ranked model that converged cleanly
  bestmodelname <- if (!is.null(saturated_model_estimates)) {
    saturated_model_estimates
  } else {
    modelaicc$model[modelaicc$Hessian == "positive_definite"][1]
  }
  if (!bestmodelname %in% names(mods)) stop("Model '", bestmodelname, "' not found among: ",
                                            paste(names(mods), collapse = ", "))
  bestmodel <- mods[[bestmodelname]]
  bestmodel$call$data <- dat
  bestmodel$call$formula <- formula(bestmodel)

  ## Bootstrap: same estimates as extract_model, CI = 1 - p2(test) / p2(infected female), negatives set to 0
  boot_fun <- function(fit) {
    v  <- extract_model(fit, estimate_response = TRUE, estimate_overdisp = TRUE, newdata = dat)
    p2 <- v[paste0("cond.", levels(dat$CrossType))]
    CI_estimate <- 1 - p2[1] / p2[-1]
    CI_estimate[CI_estimate < 0] <- 0
    names(CI_estimate) <- paste0("TestCrossvs.", levels(dat$CrossType)[-1])
    c(v, CI_estimate)
  }

  ## Refit model nsim times
  b1 <- lme4::bootMer(bestmodel, FUN = boot_fun, nsim = nbootsim, .progress = "txt")

  ## Check whether bootstrap failed
  nfail <- sum(!complete.cases(b1$t))
  if (nfail == nrow(b1$t)) {
    stop("All bootstrap refits failed: ", paste(unique(attr(b1, "boot.fail.msgs")), collapse = "; "))
  }
  if (nfail > 0) warning(nfail, " of ", nrow(b1$t), " bootstrap refits failed and were dropped")

  ## Extract 95%-CIs
  confint <- as.data.frame(bCI_tab(b = b1))

  ## Add cross (any number of crosses: k crosses, k - 1 CI estimates)
  k   <- nlevels(dat$CrossType)
  nCI <- k - 1
  nrep <- (nrow(confint) - nCI) / k
  cross_lab <- c("Test cross W- x W+", paste0("Control cross ", seq_len(nCI), " W+ x W+"))
  ci_lab    <- c("Test_cross", rep("Control_cross", nCI))
  confint$Cross   <- c(rep(cross_lab, times = nrep), rep(NA, nCI))
  confint$CI_test <- c(rep(ci_lab,    times = nrep), rep(NA, nCI))

  ## Add sample size
  n <- tapply(dat$CrossType, dat$CrossType, length)
  confint$n <- c(rep(n, times = nrep), rep(NA, nCI))

  ## Add cross type label
  confint <- confint |>
    dplyr::relocate(n) |>
    dplyr::mutate(CrossTypeLabel = paste0(CrossType, " (", n, ")"))

  print(confint)
  cols <- c("Estimate", "2.5 %", "97.5 %")
  confint <- confint |>
    dplyr::mutate(dplyr::across(dplyr::all_of(cols), \(x) round(as.numeric(x), dec_places)))

  confintp1   <- confint[grepl("zi", rownames(confint)), ]
  confintp2   <- confint[grepl("cond", rownames(confint)), ]
  confintdisp <- confint[grepl("phi|rho|sd_beta_bin", rownames(confint)), ]
  confintCI   <- confint[grepl("TestCrossvs", rownames(confint)), ]

  confintCI <- confintCI |>
    dplyr::select(-c(CI_test, CrossTypeLabel, Cross, n))

  return(list(continent = as.vector(crosses)[1], uninf_strain = as.vector(strains)[1],
              inf_strain = as.vector(strains)[2], Observed_hatch_rate = datCrossType,
              modelaicc = modelaicc, confintp1 = confintp1, confintp2 = confintp2,
              confintdisp = confintdisp, confintCI = confintCI,
              block = paste(sort(block), collapse = "_")))
}
