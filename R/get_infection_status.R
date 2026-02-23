#' Define the infection status function
#' @description Define the infection status for individuals using information for lines
#' @param individual column in the database with information for female or male line 
#' @param infected - subset of infected lines
#' @param uninfected - subset of uninfected lines
#'
#' @return A vector with infected or uninfected status according given information for the each line
#'
#' @export
#'
#' @examples
#' 

# Define the infection status function
get_infection_status <- function(individual, infected, uninfected) {
  dplyr::case_when(
    individual %in% infected ~ "Infected",
    individual %in% uninfected ~ "Uninfected",
    TRUE ~ NA_character_
  )
}

