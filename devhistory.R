# Create project on Github

## Create compendium
rrtools::use_compendium("C:/Users/sergasvi/Documents/Y2025.Wolbachia.CI", open = FALSE)

## Add to .gitignore
usethis::use_git_ignore(".DS_Store")
usethis::use_build_ignore(".DS_Store")
usethis::use_git(message = ":see_no_evil: Ban .DS_Store files")

## Modify DESCRIPTION file
usethis::edit_file("DESCRIPTION")
usethis::use_git(message = ":bulb: Update documentation")

dir.create("action2")
dir.create("action2/data/raw_data")
dir.create("action2/reports")

## Create a R directory and a file for functions
usethis::use_r("simulate_dynamics")

## Update DESCRIPTION file
usethis::use_package("here")
usethis::use_package("stringr")
usethis::use_package("ggplot2", type = "Depends")

## Update NAMESPACE file
devtools::document()

## Load all required packages
devtools::load_all()

