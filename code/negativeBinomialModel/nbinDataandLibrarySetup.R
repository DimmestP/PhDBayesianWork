# Code to set up the libraries and format/import the data for the stan model

knitr::opts_chunk$set(echo = TRUE)

setupLibraries <- function(){
  library(reshape2)
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  library(tidyverse)
  library(cowplot)
  library(ggmcmc)
  library(shinystan)
  library(bayesplot)
  library(Biostrings)
}

setupLibraries()
# The wild type data needs to be extracted from its various files and imported into a R variable
wildTypeFileNames <- list.files("~/diffExpModel/simpleTwoMutModel/data/experimental/wildType", full.names=TRUE)
mutantFileNames <- list.files("~/diffExpModel/simpleTwoMutModel/data/experimental/snf2Mutant", full.names=TRUE)
wildTypeData <- lapply(wildTypeFileNames,read.table)
mutantData <- lapply(mutantFileNames,read.table)

# Merge replicates into same dataframe
combinedWildTypeData <-wildTypeData[[1]]
combinedMutantData <- mutantData[[1]]
for(i in 2:48){
  combinedWildTypeData <- merge(combinedWildTypeData,wildTypeData[[i]],by="V1") # replicates are merged horizontally by introducing new columns
  combinedMutantData <- merge(combinedMutantData,mutantData[[i]],by="V1")
}



# Give combinedWildTypeData columns suitable names
columnNames <- vector(mode="character", length=49)
columnNames[1] <- "transcriptName" 
for(i in 1:48) columnNames[i+1] <- paste0("Replicate_",i)
colnames(combinedWildTypeData) <- columnNames
colnames(combinedMutantData) <- columnNames

# convert to tidyverse tibble
combinedWildTypeData <- as_tibble(combinedWildTypeData)
combinedMutantData <- as_tibble(combinedMutantData)
combinedWildTypeData[[1]] <- as.character(combinedWildTypeData[[1]])
combinedMutantData[[1]] <- as.character(combinedMutantData[[1]])

# Remove any transcripts not mapped to open reading frames
Scer_ORF <- readDNAStringSet("https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_coding_all.fasta.gz")
Scer_ORFName <- as.tibble(names(Scer_ORF)) %>%
  separate(value,c("transcriptName",NA),extra="drop",sep=" ")
combinedWildTypeData <- filter(combinedWildTypeData,transcriptName %in% Scer_ORFName$transcriptName)
combinedMutantData <- filter(combinedMutantData,transcriptName %in% Scer_ORFName$transcriptName)

# Remove any transcripts with consistantly low count data
library(magrittr) # pipe R package
library(tidyverse) # dataframe manipulation package
usableWildTranscripts <- 
  combinedWildTypeData %>%
  gather(key ="Replicate",value ="value",columnNames[2:49]) %>%
  group_by(transcriptName) %>%
  summarise(sufficientCounts = sum(value > 5),keep=(sufficientCounts==48)) %>%
  filter(keep) %$%
  transcriptName

usableMutantTranscripts <- 
  combinedMutantData %>%
  gather(key ="Replicate",value ="value",columnNames[2:49]) %>%
  group_by(transcriptName) %>%
  summarise(sufficientCounts = sum(value > 5),keep=(sufficientCounts==48)) %>%
  filter(keep) %$%
  transcriptName

cleanWildTypeData <- 
  combinedWildTypeData %>%
  ungroup() %>%
  filter(transcriptName %in% usableWildTranscripts) %>%
  filter(transcriptName %in% usableMutantTranscripts)

cleanMutantData <- 
  combinedMutantData %>%
  ungroup() %>%
  filter(transcriptName %in% usableWildTranscripts) %>%
  filter(transcriptName %in% usableMutantTranscripts)
