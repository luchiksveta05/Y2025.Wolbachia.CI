#' Fit a zero-inflated binomial model
#' @description Test for differences in hatch rate among the four types of crosses using a zero-inflated binomial model
#' @param crosses vector with the continent and the names of the uninfected and infected strains (in that order) used for the four types of crosss
#' @param data dataframe with hatch rate data
#' @param female name of the column in data which corresponds to the female strain in the cross (default: "Female")
#' @param male name of the column in data which corresponds to the male strain in the cross (default: "Male")
#' @param nbootsim number of simulations for parametric bootstrapping (default: nbootsim=500)
#'
#' @return A list with the different outpus (AICc table, parameters estimate and confidence intervals for the best model, etc.)
#'
#' @export
#'
#' @examples


# Define the function
fit_ZIBinomial <- function(crosses=c("Europe", "Am", "Bal"), data=merged_data, female="Female", male="Male", nbootsim=500) {

  ## Extract the name of the uninfected and infected strains
  strains <- crosses[2:3]

  print(paste("Testing for differences in hatch rates among crosses between", strains[1], "and", strains[2] ))

  ## Find batch for the test cross
  batch <- unique(data$Batch[data[, female]==strains[1]&data[, male]==strains[2]])
  batch2 <- unique(data$Batch[data[, female]==strains[2]&data[, male]==strains[1]])
  batch3 <- unique(data$Batch[data[, female]==strains[1]&data[, male]==strains[1]])
  batch4 <- unique(data$Batch[data[, female]==strains[2]&data[, male]==strains[2]])

  if(all(!batch2%in%batch)){
    batch <- c(batch, batch2[1])
  }
  if(all(!batch3%in%batch)){
    batch <- c(batch, batch3[1])
  }
  if(all(!batch4%in%batch)){
    batch <- c(batch, batch4[1])
  }
  if(grepl("F1|F2", strains[2])){
    batch <- unique(c(batch, batch2, batch3, batch4))
  }
  print(paste("Use data from the following batches:", paste(sort(batch), collapse = " " )))

  ## Subset the data with a single type of cross
  dat <- data[tibble::deframe(data[, female])%in%strains & tibble::deframe(data[, male])%in%strains & data$Batch%in%batch, ]
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

  ## Fit null model
  m0 <- glmmTMB::glmmTMB(prop ~ 1 + (1| CoupleID), data=dat, ziformula=~1 , family=binomial, weight=sum_eggcount_corr)

  ## cross type effect on hatch rate
  m1 <- glmmTMB::glmmTMB(prop ~ -1 + CrossType + (1| CoupleID), data=dat, ziformula=~ 1 , family=binomial, weight=sum_eggcount_corr)

  ## cross type effect on hatch rate and mating probability
  m2 <- glmmTMB::glmmTMB(prop ~ 1 + (1| CoupleID), data=dat, ziformula=~ -1 + CrossType , family=binomial, weight=sum_eggcount_corr)

  ## cross type effect on hatch rate and mating probability
  m3 <- glmmTMB::glmmTMB(prop ~ -1 + CrossType + (1| CoupleID), data=dat, ziformula=~ -1 + CrossType , family=binomial, weight=sum_eggcount_corr)

  ## Compare model based on AICc
  modelaicc <- data.frame(MuMIn::model.sel(m0, m1, m2, m3))

  ## Add model variables
  modelaicc <- data.frame(model_var = as.factor(rownames(modelaicc)), modelaicc)
  levels(modelaicc$model_var) <- c("zi~int|bin~int", "zi~int|bin~CrossType", "zi~CrossType|bin~int", "zi~CrossType|bin~CrossType")

  ## Get name of best model
  bestmodelname <- rownames(modelaicc)[1]
  bestmodel <- get(bestmodelname)

  # This step helps bootMer find the data
  assign("dat", dat, envir = .GlobalEnv)

  ## Refit model nsim times
  b1 <- lme4::bootMer(bestmodel, FUN=extract_model, nsim=nbootsim, .progress="txt")

  ## Extract all CIs
  confint <- bCI_tab(b=b1)
  confintp1 <- confint[grepl("zi", rownames(confint)),]
  confintp2 <- confint[grepl("cond", rownames(confint)),]
  confintCI <- confint[!grepl("zi|cond", rownames(confint)),]

  return(list(continent=as.vector(crosses)[1], uninf_strain=as.vector(strains)[1], inf_strain=as.vector(strains)[2], Observed_hatch_rate=datCrossType, modelaicc=modelaicc, confintp1=confintp1, confintp2=confintp2, confintCI=confintCI))

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


# Extract estimates from a zero-inflated binomial model
extract_model <- function(., estimateCI=FALSE, estimate_hatche_rate=FALSE){

  ## Extract p1
  p1 <- unlist(1/(1+exp(-glmmTMB::fixef(.)$zi)))
  names(p1) <- paste0("zi.", names(p1))

  ## Extract p2
  p2 <- unlist(1/(1+exp(-glmmTMB::fixef(.)$cond)))
  names(p2) <- paste0("cond.", names(p2))

if(estimate_hatche_rate) {
  vec <- c(p1, p2, hatch_rate=(1-p1)*p2)
  names(vec)[grepl("hatch_rate", names(vec))] <- gsub(".cond.", ".", names(vec)[grepl("hatch_rate", names(vec))], fixed = TRUE)
}else{
  vec <- c(p1, p2)
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
