# Launch the interactive estimator-histories learnr tutorial from the repository
# root. Uses simulated data only; writes nothing to any production directory.
rmarkdown::run("inst/tutorials/estimator-histories/estimator-histories.Rmd",
               shiny_args = list(launch.browser = TRUE))
