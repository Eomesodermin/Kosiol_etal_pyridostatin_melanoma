# pipeline variables
start.time <- Sys.time()

# Set working directory to source file location
setwd(dirname(getActiveDocumentContext()$path))

# since moving script from local to github - I want to adjust work dir to be main github dir - therefore 
setwd("..")

working.dir <- getwd()

# For reproducibility
set.seed(42)

# Explicitly load key packages
library("Seurat")
library("SeuratDisk")
library("dplyr")
library("scCustomize")
library("usefulfunctions")
library("ggplot2")
library("usefulfunctions")

# Establish colour scheme
source("scripts/variables/Colour_scheme_variable.R", local = knitr::knit_global())


###############################
# create output directories
###############################

# Common directories 
if(!dir.exists("saves")){dir.create("saves", recursive = T)}

# scRNAseq directories
if(!dir.exists("results/")){dir.create("results/", recursive = T)}


