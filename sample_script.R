### sample script

source("maxAbsDiffTopFeatures.R")
source("getAllConsensusMIWz.R")
source("findMultipleConsPSAttractors.R")

load("dataList_exprs_PFC.rda")
load("map_exprs.rda")
seeds <- read.csv("sample_seeds.csv", as.is = TRUE, header = FALSE)
seeds <- seeds$V1

cons_exprs <- findMultipleConsPSAttractors(dataMatrixList = dataList, 
					   map_short = map_short, seeds = seeds)
