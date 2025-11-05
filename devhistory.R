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

dir.create("data")
dir.create("data/raw_data")
dir.create("data/derived_data")
dir.create("reports")
dir.create("plots")
dir.create("tutorials")

## Create a R directory and a file for functions
usethis::use_r("utils-pipe")

## Update DESCRIPTION file
usethis::use_package("here")
usethis::use_package("stringr")
usethis::use_package("ggplot2", type = "Depends")
usethis::use_package("dplyr")
usethis::use_package("lme4")
usethis::use_package("stringr")
usethis::use_package("ggpubr")
usethis::use_package("glmmTMB")
usethis::use_package("MuMIn")
usethis::use_package("DHARMa")
usethis::use_package("emmeans")
usethis::use_package("multcomp")

## Update NAMESPACE file
devtools::document()

## Load all required packages
devtools::load_all()

