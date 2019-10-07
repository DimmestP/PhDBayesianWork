## Code to parallelise the downsampling code in downSampledDESeqcomparison.Rmd so that multiple bootstrapping can occur
## Parallelisation occurs thanks to the genius of the Future package


# Import all experimental data
source("~/diffExpModel/PhDBayesianWork/code/negativeBinomialModel/nbinDataandLibrarySetup.R")

# Run multi-iteration, multi-downsampled Rstan code
library(doFuture)
registerDoFuture()
library(rstan)
library(Rcpp)
library(tictoc) # Measure the run time of a program


setupLibraries()
WTData <- cleanWildTypeData %>%
  select(-transcriptName)
MTData <- cleanMutantData %>%
  select(-transcriptName)
chosenRNABatchSize <- 5755

# Bootstrap replicate downsampling to have multiple runs of  randomly choosen replicates
chosenNumberReplicates <- c(2,4,8,12,24,48)

# Empty to lists to hold vectors of the columns of downsampled replicates, dynamic because the number of replicates to downsample to changes
WTcolNamesRep <- list()
MTcolNamesRep <- list()

# Logic to randomly select 4 different combination of replicates according the given value to downsample to
for(i in 1:length(chosenNumberReplicates)){
  currentReplicates <- chosenNumberReplicates[i]
  if(currentReplicates != 48){
    for(j in 1:4){
      WTcolNamesRep[[(((i-1) * 4) + j)]] <- sample(1:48, currentReplicates)
      MTcolNamesRep[[(((i-1) * 4) + j)]] <- sample(1:48, currentReplicates)
    }
  }
  else {
    WTcolNamesRep[[(((i-1) * 4) + 1)]] <- 1:ncol(WTData)
    MTcolNamesRep[[(((i-1) * 4) + 1)]] <- 1:ncol(MTData)
    for(j in 2:4){
      WTcolNamesRep[[(((i-1) * 4) + j)]] <- WTcolNamesRep[[(((2-1) * 4) + 1)]]
      MTcolNamesRep[[(((i-1) * 4) + j)]] <- WTcolNamesRep[[(((2-1) * 4) + 1)]]
    }
  }
}

# Function to prepare all of the different combinations of replicate data s.t. rstan accepts it
prepareData <- function(chosenRNABatchSize,WTcolNames,MTcolNames){
  countsArray <- array(
    cbind(
      data.matrix(
        WTData[1:chosenRNABatchSize, WTcolNames]
      ),
      data.matrix(
        MTData[1:chosenRNABatchSize, MTcolNames]
      )
    ),
    dim = c(chosenRNABatchSize,length(WTcolNames),2)
  )# manipulation to merge count data into a single array for use in Rstan
  
  stanData <- list(NRNA = chosenRNABatchSize, NReplicates = length(WTcolNames), counts = countsArray) 
}

# Function to save output of successful Rstan run for future downstream analysis
saveStanFit <- function(i,j,nbinlogStanFit) {
  fileName <- paste0("~/diffExpModel/PhDBayesianWork/data/simulated/autoDownSample",length(WTcolNamesRep[[(((i-1) * 4) + j)]]),"Replicates",j,"iter.rds")
  currentMetaData <- read_tsv("~/diffExpModel/PhDBayesianWork/data/simulated/fourRepdownSamplednBinCompMetaData.tsv")
  metaData <- currentMetaData %>%  
    bind_rows(tibble(file=fileName,repNo=length(WTcolNamesRep),wildRep=WTcolNamesRep[[(((i-1) * 4) + j)]],mutRep=MTcolNamesRep[[(((i-1) * 4) + j)]],maxRhat = max(summary(nbinlogStanFit)$summary[,10])))
  write_tsv(metaData,paste0("~/diffExpModel/PhDBayesianWork/data/simulated/fourRepdownSamplednBinCompMetaData.tsv"))
  save(nbinlogStanFit,file=fileName)
}

# set up future plan s.t. the parallel processes know exactly how many cores it has access to 
plan(list(sequential,tweak(multiprocess, workers = 4L)))

tic()

foreach(i = 1:6) %dopar% {
  foreach(j = 1:4) %dopar% {
    rstan_options(auto_write = TRUE)
    
    # Once bootstrapped downsampling has been completed, run the simulation
    # Using more cores to reduce the timing!
    stanDataCurrent <- prepareData(chosenRNABatchSize,WTcolNamesRep[[(((i-1) * 4) + j)]],MTcolNamesRep[[(((i-1) * 4) + j)]])

    nbinlogStanFit %<-% stan("~/diffExpModel/PhDBayesianWork/code/negativeBinomialModel/nbinModel.stan",data = stanDataCurrent,iter = 2000, chains = 4,cores = 4)
    
    saveStanFit(i,j,nbinlogStanFit)
  }
}  

toc()

plan(sequential)


