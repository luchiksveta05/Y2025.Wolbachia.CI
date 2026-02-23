#' Fit a zero-inflated binomial model
#' @description Test for differences in hatch rate among the four types of crosses using a zero-inflated binomial model
#' @param id number of the infected strain from which the data should be extracted
#' @param outputper_infstrain list output from function  fit_ZIBinomial
#' @param parameter parameter of interest (should be one of "confintp1", "confintp2", "confintCI" or "modelaicc")
#'
#' @return A vector or a dataframe with the required output (AICc table, confintp1, confinP2, etc.)
#'
#' @export
#'
#' @examples
#' #extract_results(id=1, outputper_infstrain=results, parameter="confintp2")
## Define function
extract_results <- function(id, outputper_infstrain, parameter){
  
  ## Extract continent and name of strain infected with Wolbahcia
  continent <- outputper_infstrain[[id]]["continent"][[1]]
  uninf_strain <- outputper_infstrain[[id]]["uninf_strain"][[1]]
  inf_strain <- outputper_infstrain[[id]]["inf_strain"][[1]]
  
  ## Extract output for parameter of interest
  output <- outputper_infstrain[[id]][parameter]
  
  if(length(output[parameter][[1]])==3){
    output <- unlist(output)
    names(output) <- gsub("confintp1.", "", names(output))
    names(output) <- gsub("confintp2.", "", names(output))
    
  }else{
    output <- as.data.frame(output[parameter][[1]])
    type <- rownames(output)
  }
  
  
  if(parameter=="modelaicc"){
    output2 <- data.frame(continent=rep(continent, nrow(output)), uninfstrain=rep(uninf_strain, nrow(output)), infstrain=rep(inf_strain, nrow(output)), bestmodel=c("Yes", rep("No", nrow(output)-1)), model=type, output)
  } else{
    
    if(is.vector(output)){
      output2 <- data.frame(continent=continent, uninfstrain=uninf_strain, infstrain=inf_strain, type=ifelse(parameter=="confintp1", "zi.(Intercept)", "cond.(Intercept)"), t(output))
      
    }else{
      output2 <- data.frame(continent=rep(continent, nrow(output)),  uninfstrain=rep(uninf_strain, nrow(output)), infstrain=rep(inf_strain, nrow(output)), type=type, output)
    }
    colnames(output2) <- gsub("X2.5..", "2.5%", colnames(output2))
    colnames(output2) <- gsub("X97.5..", "97.5%", colnames(output2))
  }
  
  rownames(output2) <- NULL
  
  return(output2)
  
}