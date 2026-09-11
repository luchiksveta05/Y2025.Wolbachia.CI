#' Making clones function
#' @description Duplicate rows for specific cross
#' @param data datafram from where the data will be replicated 
#' @param crosstype - specific crosstype which will be replicated
#' @param repeats - number of clones needed
#'
#' @return Datafram with repeated rows from the database
#'
#' @export
#'
#' @examples
#' 
## Define the function
duplicate_rows_by_crosstype <- function(data, crosstype, repeats) {
  filtered_data <- data %>%
    dplyr::filter(CrossType == crosstype)
  
  # Create an empty data frame to store the duplicated rows
  duplicated_data <- data.frame()
  
  # Loop to duplicate rows and modify the CrossType
  for (i in 1:repeats) {
    temp_data <- filtered_data
    temp_data$CrossType <- paste0(crosstype, "_clone", i)
    duplicated_data <- rbind(duplicated_data, temp_data)
  }
  
  return(duplicated_data)
}
