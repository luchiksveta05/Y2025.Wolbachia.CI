#' Define description of the cross
#' @description Define the description of the cross (test or control)
#' @param infected - subset of infected lines
#' @param uninfected - subset of uninfected lines
#'
#' @return A vector with infected or uninfected status according given information for the each line
#'
#' @export
#'
#' @examples
#' 

# Define the cross description function
determine_cross <- function(male_infect, female_infect) {
  if (male_infect == "Infected" & female_infect == "Uninfected") {
    return("Test cross W- x W+")
  } else if (male_infect == "Uninfected" & female_infect == "Infected") {
    return("Control cross 2 W+ x W-")
  } else if (male_infect == "Uninfected" & female_infect == "Uninfected") {
    return("Control cross 1 W- x W-")
  } else if (male_infect == "Infected" & female_infect == "Infected") {
    return("Control cross 3 W+ x W+")
  } else {
    return(NA) # or "NA" if you want it as a string
  }
}
