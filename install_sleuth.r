source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
biocLite("biomaRt")
install.packages("devtools")
devtools::install_github("pachterlab/sleuth")


source("http://bioconductor.org/biocLite.R") # This allows you to install Bioconductor packages
biocLite("seqplots")
library(seqplots)
run()