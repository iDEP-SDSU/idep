# Building a Prod-Ready, Robust Shiny Application.
# 
# README: each step of the dev files is optional, and you don't have to 
# fill every dev scripts before getting started. 
# 01_start.R should be filled at start. 
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
# 
# 
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Engineering

## Dependencies ----
## Add one line by package you want to add as dependency
usethis::use_package("dplyr")
usethis::use_package("tidyselect")
usethis::use_package("limma")
usethis::use_package("edgeR")
usethis::use_package("DESeq2")
usethis::use_package("SummarizedExperiment")
usethis::use_package("BiocGenerics")
usethis::use_package("ggplot2")
usethis::use_package("tidyr")
usethis::use_package("circlize")
usethis::use_package("ComplexHeatmap")
usethis::use_package("grid")
usethis::use_package("stats")
usethis::use_package("factoextra")
usethis::use_package("InteractiveComplexHeatmap")
usethis::use_package("GetoptLong")
usethis::use_package("reshape2")
usethis::use_package("dendextend")
usethis::use_package("visNetwork")
usethis::use_package("gage")
usethis::use_package("PGSEA")
usethis::use_package("ReactomePA")
usethis::use_package("fgsea")
usethis::use_package("pathview")
usethis::use_package("png")
usethis::use_package("KEGGREST")
usethis::use_package("graphics")
usethis::use_package("plotly")
usethis::use_package("PREDA")
usethis::use_package("PREDAsampledata")
usethis::use_package("biclust")
usethis::use_package("WGCNA")
usethis::use_package("flashClust")
usethis::use_package("dynamicTreeCut")
usethis::use_package("igraph")
usethis::use_package("DBI")
usethis::use_package("RSQLite")
usethis::use_package("htmltools")
usethis::use_package("shinyjs")
usethis::use_package("DT")
usethis::use_package("shinybusy")
usethis::use_package("stringr")
usethis::use_package("QUBIC")
usethis::use_package("runibic")
usethis::use_package("grDevices")
usethis::use_package("loose.rock")
usethis::use_package("PCAtools")
usethis::use_package("shinyBS")
usethis::use_package("tippy")
usethis::use_package("ottoPlots")

## Add modules ----
## Create a module infrastructure in R/
golem::add_module( name = "10_doc" ) # Name of the module
#golem::add_module( name = "name_of_module2" ) # Name of the module

## Add helper functions ----
## Creates fct_* and utils_*
golem::add_fct( "10_doc" ) 
#golem::add_utils( "helpers" )

## External resources
## Creates .js and .css files at inst/app/www
# golem::add_js_file( "script" )
# golem::add_js_handler( "handlers" )
# golem::add_css_file( "custom" )

## Add internal datasets ----
## If you have data in your package
#usethis::use_data_raw( name = "my_dataset", open = FALSE ) 

## Tests ----
## Add one line by test you want to create
usethis::use_test( "app" )

# Documentation

## Vignette ----
usethis::use_vignette("DEG")
devtools::build_vignettes("Load Data")

## Code Coverage----
## Set the code coverage service ("codecov" or "coveralls")
# usethis::use_coverage()

# # Create a summary readme for the testthat subdirectory
# covrpage::covrpage()

# ## CI ----
# ## Use this part of the script if you need to set up a CI
# ## service for your application
# ## 
# ## (You'll need GitHub there)
# usethis::use_github()

# # GitHub Actions
# usethis::use_github_action() 
# # Chose one of the three
# # See https://usethis.r-lib.org/reference/use_github_action.html
# usethis::use_github_action_check_release() 
# usethis::use_github_action_check_standard() 
# usethis::use_github_action_check_full() 
# # Add action for PR
# usethis::use_github_action_pr_commands()

# # Travis CI
# usethis::use_travis() 
# usethis::use_travis_badge() 

# # AppVeyor 
# usethis::use_appveyor() 
# usethis::use_appveyor_badge()

# # Circle CI
# usethis::use_circleci()
# usethis::use_circleci_badge()

# # Jenkins
# usethis::use_jenkins()

# # GitLab CI
# usethis::use_gitlab_ci()

# # You're now set! ----
# # go to dev/03_deploy.R
# rstudioapi::navigateToFile("dev/03_deploy.R")

