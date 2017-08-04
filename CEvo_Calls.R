#!/usr/bin/env RScript
rm (list = ls())

args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("$PWD variable must be given!", call.=FALSE)
}

#source(paste(args[1], "/CEvo_Analysis.R", sep = ""))
#source(paste(args[1], "/CEvo_Sequence.R", sep = ""))
source(paste("D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Dynamics_test/CEvo_Analysis.R", sep = ""))
source(paste("D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Dynamics_test/CEvo_Sequence.R", sep = ""))

#Mother folder for a whole daily session (*/Output/Euler/20170711/)
pa <- "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/20170712/"
#Path for finding the infection chain file (*/20170710/1000bp_5SNP_04Malus_seqReal_ImSys005_HighCap2/)
p <- list.dirs(pa, recursive = F)[1]
#p <- paste(pa, "/1000bpRS_5SNP_02Malus_ImSys001_HighCap_LowEvoRate/", sep = "")
#p <- paste(pa, "/1000bp_5SNP_04Malus_seqReal_ImSys005_HighCap3/", sep = "")
#Epidemics folder for the analysis of the sequences (*/SOME_NAME/Epidemics_1/)
pas <- paste(p, "/Epidemics_1/", sep = "")
#mother directory of the sample trees (before, after, infection, random) directories after BEAST analysis
pts <- "C:/Users/iotn_/Desktop/Toy_Bayesian/20170711_HighCap_LowEvo/Epi_2"
ptss <- paste(pts, "/Random_100Strains/", sep = "")



for (d in dir(pa, recursive = F))
{
  p <- paste(pa, "/", d, "/", sep = "")
  #  par_file <- paste(p, "parameter_search.sh", sep = "")
  nr.ep <- length(list.dirs(p, recursive = F))
  
  
  #read parameters
  
  #get Infection rate constant line and parse it
  # irc <- readLines(par_file)[3]
  #irc <- strsplit(irc, split = '[\\(\\)]')[[1]][2]
  #irc_label <- strsplit(irc, split = " ")[[1]]
  #irc <- as.numeric(strsplit(irc, split = " ")[[1]])
  
  #get fitness and parse it
  #fitness <- readLines(par_file)[4]
  #fitness <- strsplit(fitness, split = '[\\(\\)]')[[1]][2]
  #fitness_label <- strsplit(fitness, split = " ")[[1]]
  #fitness <- as.numeric(strsplit(fitness, split = " ")[[1]])
  
  #End of parameter reading
  
  
  
  for (e in 1:nr.ep) plot.ViralLoad(path = p, epidemics = e, per_vir = F, onlyHC = F)
  #for (e in 1:nr.ep) plot.ViralLoad(path = p, epidemics = e, per_vir = T, onlyHC = F)
  for (e in 1:nr.ep) plot.InfTree(path = p, epidemics = e)
}
#This to generate the multiple alignment and times of sampling
SampleAndConvert(path_to_epi = pas, which = "random", howMuchTime = 30, HowManyStrains = 100)
SampleAndConvert(path_to_epi = pas, which = "infection", howMuchTime = 30, HowManyStrains = 100)
SampleAndConvert(path_to_epi = pas, which = "transmitted", howMuchTime = 30, HowManyStrains = 100)
SampleAndConvert(path_to_epi = pas, which = "tips", howMuchTime = 30, HowManyStrains = 100)

#This to do after BEAST analysis
CompareTrees(path_to_ref = p, path_to_samples = pts, epidemics = 2, dist_method = "treeVec")
