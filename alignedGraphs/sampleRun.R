setwd("~/OneDrive\ -\ purdue.edu/Research/sna/alignedGraphs/")
source("https://raw.githubusercontent.com/kurtmaia/sna/master/sna_utils.R")
source('netCompMM.R')


# Synthetic data
source("generateSyntheticData.R")

mm <- netCompMM(graphs, side, family = "binom", link = "logit", heterogenity = TRUE, threshold = NULL, matrix.plot = TRUE, R =15, H = 12, Tsamples = 1500, burnIn = 500, by.print = 15, nCores = 8)	

# Brazil twitter dataset
dtset <- "dados.RDS"
dtset_side <- "brazilTwitter_userSide.RDS"

graphs <- readRDS(dtset)$graphs
side_dt <- readRDS(dtset_side)
side_dt <- readRDS(dtset_side)
side_dt <- dplyr::filter(side_dt,user %in% names(graphs))
xtabs(~side,side_dt)
side <- side_dt$side
names(side) <- side_dt$user

mm <- netCompMM(graphs, side, family = "binom", link = "logit",heterogenity = FALSE, threshold = NULL, R =10, H = 15, Tsamples = 1500, burnIn = 500, by.print = 10, nCores = 4)
