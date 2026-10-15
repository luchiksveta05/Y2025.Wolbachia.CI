#' Fit a zero-inflated beta-binomial model
#' @description Test for differences in hatch rate among the four types of crosses using a zero-inflated binomial model
#' @param crosses vector with the continent and the names of the uninfected and infected strains (in that order) used for the four types of crosss
#' @param data dataframe with hatch rate data
#' @param female name of the column in data which corresponds to the female strain in the cross (default: "Female")
#' @param male name of the column in data which corresponds to the male strain in the cross (default: "Male")
#' @param nbootsim number of simulations for parametric bootstrapping (default: nbootsim=500)
#' @param dec_places number of decimal places for the output table (default: dec_places=3)
#' @param saturated_model_estimates whether to return estimates and 95%-confidence interval from one saturated model ("m2_zi1_disp1", "m2_zi1_dispCT", "m2_ziCT_disp1" or "m2_ziCT_dispCT") and not from the best model with lowest AICc  (default: saturated_model_estimates=NULL)
#' @param include_all_blocks whether all blocks should be included for the analysis (default: include_all_blocks=FALSE)
#' @return A list with the different outpus (AICc table, parameters estimate and confidence intervals for the best model, etc.)
#'
#' @export
#'
#' @examples


# Define the function
fit_ZI_Beta_binomial <- function(crosses=c("Europe", "Am", "Bal"), data=merged_data, female="Female", male="Male", nbootsim=500, dec_places=3, saturated_model_estimates=NULL, include_all_blocks=FALSE,
                                 zis = list(ziCT = ~ -1 + CrossType, ## zero-inflation differs among cross types
                                             zi1  = ~ 1),
                                 disps = list(
                                   dispCT = ~ -1 + CrossType,      ## dispersion differs among cross types
                                   disp1  = ~ 1)) {

  ## Extract the name of the uninfected and infected strains
  strains <- crosses[2:3]

  print(paste("Testing for differences in hatch rates among crosses between", strains[1], "and", strains[2] ))

  if(include_all_blocks){
    block <- unique(c(block, block2, block3, block4))
  }else{
    ## Find block for the test cross
    block <- unique(data$Batch[data[, female]==strains[1]&data[, male]==strains[2]])
    block2 <- unique(data$Batch[data[, female]==strains[2]&data[, male]==strains[1]])
    block3 <- unique(data$Batch[data[, female]==strains[1]&data[, male]==strains[1]])
    block4 <- unique(data$Batch[data[, female]==strains[2]&data[, male]==strains[2]])

    ## Test if blocks of control crosses are the same as the test cross
    ## If not, add one block to that of the test cross
    if(all(!block2%in%block)){
      block <- c(block, block2[1])
    }
    if(all(!block3%in%block)){
      block <- c(block, block3[1])
    }
    if(all(!block4%in%block)){
      block <- c(block, block4[1])
    }
    if(grepl("F1|F2", strains[2])){
      block <- unique(c(block, block2, block3, block4))
    }
  }

  print(paste("Use data from the following blocks:", paste(sort(block), collapse = " " )))

  ## Subset the data with a single type of cross
  dat <- data[tibble::deframe(data[, female])%in%strains & tibble::deframe(data[, male])%in%strains & data$Batch%in%block, ]
  ## Compute proportion of eggs that hatch into larvae
  dat$prop <- dat$sum_larvacount / dat$sum_eggcount_corr

  if(female!="Female"| male!="Male"){
    dat$CrossType <- dat$CrossType2
  }

  ## Reorder crosstypes
  dat$CrossType <- factor(dat$CrossType)
  lev <- levels(dat$CrossType)

  ## Name of the crosses
  uninf_inf <- paste(strains, collapse="x")
  inf_uninf <- paste(rev(strains), collapse="x")
  uninf_uninf <- paste(rep(strains[1], 2), collapse="x")
  inf_inf <- paste(rep(strains[2], 2), collapse="x")

  ## Take into account special character in names
  uninf_uninf_corr <- paste0(uninf_uninf, "$")
  inf_inf_corr <- gsub("\\^", "\\\\^", inf_inf)
  uninf_inf_corr <- gsub("\\^", "\\\\^", uninf_inf)
  inf_uninf_corr <- paste0(gsub("\\^", "\\\\^", inf_uninf), "$")

  ## Order of the crosses in the right way
  dat$CrossType <- factor(dat$CrossType, levels=lev[c(grep(uninf_uninf_corr, lev), grep(uninf_inf_corr, lev), grep(inf_uninf_corr, lev), grep(inf_inf_corr, lev))])

  datCrossType <- dat |>
    dplyr::group_by(CrossType, Batch, Year) |>
    dplyr::summarise(Sum_eggcount_corr = sum(sum_eggcount_corr),
                     Sum_larvacount = sum(sum_larvacount))|>
    dplyr::mutate(hatchrate = Sum_larvacount / Sum_eggcount_corr)

  print(datCrossType)

  ## List all models
  forms <- list(
    m0   = prop ~ 1, ## Null model
    m1   = prop ~ CI_test, ## CI model
    m2   = prop ~ -1 + CrossType) ## Saturated model

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

  ## Fit all models
  #mods <- lapply(forms, function(f)
  #  glmmTMB::glmmTMB(update(f, prop ~ .),
  #                   ziformula = ~-1 +CrossType, dispformula = ~ -1+ CrossType, family = glmmTMB::betabinomial(), weights = sum_eggcount_corr, data = dat))

  ## Compare models based on AICc
  modelaicc <- data.frame(MuMIn::model.sel(mods))
  modelaicc <- modelaicc |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, dec_places)))

  ## Check convergence
  conv <- sapply(mods, \(x) x$fit$convergence)
  pd   <- sapply(mods, \(x) x$sdr$pdHess)
  #modelaicc$convergence <- ifelse(conv[rownames(modelaicc)] == 0, "converged", "not_converged")
  #modelaicc$Hessian     <- ifelse(pd[rownames(modelaicc)], "positive_definite", "non_positive_definite")

  ## Add model variables
  #labels <- c(m0 = "null (intercept)", m1 = "CI (test vs control crosses)", m2 = "saturated (cross type)")
  #modelaicc <- data.frame(model = rownames(modelaicc), model_var = labels[rownames(modelaicc)], modelaicc)

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
                          Hessian     = ifelse(pd[rownames(modelaicc)], "positive_definite", "non_positive_definite"))

  ## Get name of best model
  bestmodelname <- ifelse(is.null(saturated_model_estimates), rownames(modelaicc)[1], saturated_model_estimates)
  bestmodel <- mods[[bestmodelname]]

  # This step helps bootMer find the data
  #assign("dat", dat, envir = .GlobalEnv)
  bestmodel$call$data <- dat
  bestmodel$call$formula <- formula(bestmodel)

  ## Refit model nsim times
  b1 <- lme4::bootMer(bestmodel, FUN=function(fit) extract_model(fit, estimateCI = TRUE, estimate_response=TRUE, estimate_overdisp=TRUE, newdata=dat), nsim=nbootsim, .progress="txt")

  ## Check wether bootstrap failed
  nfail <- sum(!complete.cases(b1$t))
  if (nfail == nrow(b1$t)) {
    stop("All bootstrap refits failed: ", paste(unique(attr(b1, "boot.fail.msgs")), collapse = "; "))
  }
  if (nfail > 0) warning(nfail, " of ", nrow(b1$t), " bootstrap refits failed and were dropped")

  ## Extract 95%-CIs
  confint <- as.data.frame(bCI_tab(b=b1))

  ## Add cross
  confint$Cross <- c(rep(c("Control cross 1 W- x W-", "Test cross W- x W+", "Control cross 2 W+ x W-", "Control cross 3 W+ x W+"), times=(nrow(confint)-3)/4), rep(NA, 3))
  confint$CI_test <- c(rep(c("Control_cross", "Test_cross", "Control_cross", "Control_cross"), times=(nrow(confint)-3)/4), rep(NA, 3))

  ## Add sample size
  n <- tapply(dat$CrossType, dat$CrossType, length)
  confint$n <- c(rep(n, times=(nrow(confint)-3)/4), rep(NA, 3))

  ## Add cross type label
  confint <- confint|>
    dplyr::relocate(n)|>
    dplyr::mutate(CrossTypeLabel = paste0(CrossType, " (", n, ")"))

  print(confint)
  cols <- c("Estimate", "2.5 %", "97.5 %")
  confint <- confint|>
    dplyr::mutate(dplyr::across(dplyr::all_of(cols), \(x) round(as.numeric(x), dec_places)))

  confintp1 <- confint[grepl("zi", rownames(confint)),]
  confintp2 <- confint[grepl("cond", rownames(confint)),]
  confintdisp<- confint[grepl("phi|rho|sd_beta_bin", rownames(confint)),]
  confintCI <- confint[grepl("TestCrossvs", rownames(confint)),]

  confintCI <- confintCI |>
    dplyr::select(-c(CI_test, CrossTypeLabel, Cross, n))

  return(list(continent=as.vector(crosses)[1], uninf_strain=as.vector(strains)[1], inf_strain=as.vector(strains)[2], Observed_hatch_rate=datCrossType, modelaicc=modelaicc, confintp1=confintp1, confintp2=confintp2, confintdisp=confintdisp, confintCI=confintCI, block=paste(sort(block), collapse = "_" )))

}

# Extract estimates from a zero-inflated binomial model
extract_model_logit_scale <- function(.){

  ## Extract p1
  logitp1 <- unlist(glmmTMB::fixef(.)$zi)
  names(logitp1) <- paste0("zi.", names(logitp1))

  ## Extract p2
  logitp2 <- unlist(glmmTMB::fixef(.)$cond)
  names(logitp2) <- paste0("cond.", names(logitp2))

  vec <- c(logitp1, logitp2)

  return(vec)
}


#
#' Extract estimates from a zero-inflated (beta-)binomial model
#' @description Function to extract estimates from a zero-inflated (beta-)binomial model when using lme4::bootMer for bootstrapping
#' @param estimateCI estimate CI from the four cross types
#' @param estimate_overdisp estimate overdispersion in the beta-binomial distribution
#' @param estimate_response make prediction on the response scale (1-p1)*p2
#' @param newdata name of the original dataframe to pass it to the function
#' @return A dataframe with the estimates of the different parameters
#'
#' @export
#'
#' @examples


extract_model <- function(., estimateCI=FALSE, estimate_overdisp=FALSE, estimate_response=FALSE, newdata){


  ## Extract p1
  p1 <- tapply(predict(., type="zprob"), newdata$CrossType, mean)
  names(p1) <- paste0("zi.", names(p1))

  ## Extract p2
  p2 <- tapply(predict(., type="conditional"), newdata$CrossType, mean)
  names(p2) <- paste0("cond.", names(p2))


  vec <- c(p1, p2)

  if(estimate_response) {
    ## Compute predicted hatch rate
    prop <- tapply(predict(., type="response"), newdata$CrossType, mean)
    names(prop) <- paste0("res.", names(prop))

    vec <- c(vec, prop)

  }

  if(estimate_overdisp) {
    ## Extract phi (precision parameter)
    phi <- tapply(predict(., type="disp"), newdata$CrossType, mean)

    ## Compute rho (overdispertion parameter)
    rho <- 1 / (phi + 1)

    ## Update names
    names(phi) <- paste0("phi.", names(phi))
    names(rho) <- paste0("rho.", names(rho))
    ## Compute between-couple standard deviation in hatch probability
    sd_beta_bin <- round(sqrt(p2 * (1 - p2) * rho), 3)
    names(sd_beta_bin) <- paste0("sd_beta_bin.", gsub("cond.", "", names(sd_beta_bin)))

    vec <- c(vec, phi, rho, sd_beta_bin)

  }

  if(estimateCI){
    if(length(p2)>1){
      ## Reorder estimates
      CI_estimate <- 1-p2[2]/p2[c(1, 3, 4)]
      CI_estimate[CI_estimate<0] <- 0
    }else{
      CI_estimate <- rep(0, 3)
    }
    ## Export results
    names(CI_estimate) <- c("TestCrossvs.CtrlCross1", "TestCrossvs.CtrlCross2", "TestCrossvs.CtrlCross3")
    vec <- c(vec, CI_estimate)
  }

  return(vec)
}

#' @title 95%CI compute
#'
#' @description  Compute 95% confidence interval of a mixed model
#' @param boo01 bootMer object
#'
#' @return
#' @export
#'
#' @examples
#'boo01 <- lme4::bootMer(m2, mySumm, nsim = 2, re.form = NA)
#'bCI.tab(boo01)
#'
bCI_tab <- function(b,ind=length(b$t0), type="perc", conf=0.95) {
  btab0 <- t(sapply(as.list(seq(ind)),
                    function(i)
                      boot::boot.ci(b, index=i, conf=conf, type=type)$percent))
  btab <- btab0[,4:5]
  rownames(btab) <- names(b$t0)
  a <- (1 - conf)/2
  a <- c(a, 1 - a)
  pct <- format.perc(a, 3)
  colnames(btab) <- pct
  return(cbind(CrossType=gsub("zi.|cond.|res.|phi.|rho.|sd_beta_bin.","", names(b$t0)), Estimate=b$t0, btab))
}

format.perc <- function (probs, digits) {
  return(paste(format(100 * probs, trim = TRUE,
                      scientific = FALSE, digits = digits),  "%"))
}


#' Extract the likelihood from glmmTMB model
#' @export
#'
#assignInNamespace(
#  "logLik.glmmTMB",
#  function(object, ...) {
#    structure(
#      -object$fit$objective,
#      df   = as.integer(length(object$fit$par)),
#      nobs = as.integer(nobs(object)),
#      class = "logLik"
#    )
#  },
#  ns = "glmmTMB"
#)
