library(plumber)
r <- plumb("plumberImp.R")  # Where 'myfile.R' is the location of the file shown above
r$run(port=8000)