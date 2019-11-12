runDESeq <- function(e) {
  cds <- newCountDataSet(exprs(e), pData(e)$condition)
  cds <- DESeq::estimateSizeFactors(cds)
  cds <- DESeq::estimateDispersions(cds)
  suppressWarnings({capture.output({fit1 <- DESeq::fitNbinomGLMs(cds, count ~ condition)})})
  suppressWarnings({capture.output({fit0 <- DESeq::fitNbinomGLMs(cds, count ~ 1)})})
  pvals <- DESeq::nbinomGLMTest(fit1, fit0)
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=fit1$conditionB))
}

runDESeq2 <- function(e, retDDS=FALSE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,quiet=TRUE)
  res <- results(dds)
  beta <- res$log2FoldChange
  pvals <- res$pvalue
  padj <- res$padj
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=beta))
}

runDESeq2Outliers <- function(e, retDDS=FALSE) {
  dds <- DESeqDataSetFromMatrix(exprs(e), DataFrame(pData(e)), ~ condition)
  ddsDefault <- DESeq(dds, quiet=TRUE)
  ddsNoRepl <- ddsDefault
  if (ncol(e) >= 14) {
    # insert original maximum Cook's distances
    # so the rows with replacement will be filtered
    # this avoid re-running with minReplicateForReplace=Inf
    mcols(ddsNoRepl)$maxCooks <- apply(assays(ddsNoRepl)[["cooks"]], 1, max)
  }
  resDefault <- results(ddsDefault)
  resNoFilt <- results(ddsDefault, cooksCutoff=FALSE)
  resNoRepl <- results(ddsNoRepl)
  resList <- list("DESeq2"=resDefault, "DESeq2-noFilt"=resNoFilt, "DESeq2-noRepl"=resNoRepl)
  resOut <- lapply(resList, function(res) {
    pvals <- res$pvalue
    padj <- res$padj
    pvals[is.na(pvals)] <- 1
    pvals[rowSums(exprs(e)) == 0] <- NA
    padj <- p.adjust(pvals,method="BH")
    padj[is.na(padj)] <- 1
    list(pvals=pvals, padj=padj)
  })
  return(resOut)
}

runEdgeR <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  dgel <- estimateGLMCommonDisp(dgel, design)
  dgel <- estimateGLMTrendedDisp(dgel, design)
  dgel <- estimateGLMTagwiseDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  predbeta10 <- predFC(exprs(e), design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
       predbeta=predbeta[,"pData(e)$conditionB"], predbeta10=predbeta10[,"pData(e)$conditionB"])
}

runEdgeRRobust <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # settings for robust from robinson_lab/edgeR_robust/robust_simulation.R
  dgel <- estimateGLMRobustDisp(dgel, design, maxit=6)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predbeta <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pvals <- edger.lrt$table$PValue
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
       predbeta=predbeta[,"pData(e)$conditionB"])
}


runDSS <- function(e) {
  X <- as.matrix(exprs(e))
  colnames(X) <- NULL
  designs <- as.character(pData(e)$condition)
  seqData <- newSeqCountSet(X, designs)
  seqData <- estNormFactors(seqData)
  seqData <- estDispersion(seqData)
  result <- waldTest(seqData, "B", "A")
  result <- result[match(rownames(seqData),rownames(result)),]
  pvals <- result$pval
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=( log2(exp(1)) * result$lfc ))
}

runVoom <- function(e) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pvals <- tt$P.Value 
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals, padj=padj, beta=tt$logFC)
}

runSAMseq <- function(e) {
  set.seed(1)
  x <- exprs(e)
  y <- pData(e)$condition
  capture.output({samfit <- SAMseq(x, y, resp.type = "Two class unpaired")})
  beta <- log2(samfit$samr.obj$foldchange)
  pvals <- samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pvals,method="BH")
  padj[is.na(padj)] <- 1
  list(pvals=pvals,padj=padj,beta=beta)
}

runEBSeq <- function(e) {
  sizes <- MedianNorm(exprs(e))
  out <- capture.output({
    suppressMessages({
      res <- EBTest(Data = exprs(e),
                    Conditions = pData(e)$condition,
                    sizeFactors = sizes,
                    maxround = 5)
    })
  })
  padj <- rep(1, nrow(exprs(e)))
  # we use 1 - PPDE for the FDR cutoff as this is recommended in the EBSeq vignette
  padj[match(rownames(res$PPMat), rownames(e))] <- res$PPMat[,"PPEE"]
  beta <- rep(0, nrow(exprs(e)))
  list(pvals=padj, padj=padj, beta=beta)
}

runGFOLD <- function(e) {
  gfold.samples.dir <- "/path/to/gfold/samples"
  gfold.out.dir <- "/path/to/gfold/output"
  m <- ncol(e)
  n <- nrow(e)
  sample.fls <- replicate(m, tempfile(pattern = "sample",
                                      tmpdir = gfold.samples.dir, fileext = ""))
  out.fl <- tempfile(pattern = "out",
                     tmpdir = gfold.out.dir, fileext = "")
  for (i in seq_len(m)) {
    out <- cbind(1:n, rep(0,n), exprs(e)[,i], rep(0,n), rep(0,n))
    write.table(out, file=sample.fls[i],
                quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  }
  samples1 <- paste(sample.fls[pData(e)$condition == "A"], collapse=",")
  samples2 <- paste(sample.fls[pData(e)$condition == "B"], collapse=",")
  system(paste0("gfold diff -s1 ", samples1, " -s2 ", samples2, " -o ", out.fl),
         ignore.stdout = TRUE, ignore.stderr = TRUE)
  res <- read.table(out.fl)
  
  # clean up
  for (fl in c(sample.fls, out.fl, paste0(out.fl,".ext"))) {
    system(paste0("rm -f ",fl))
  }
  colnames(res) <- c("GeneSymbol","GeneName","GFOLD(0.01)",
                     "E-FDR","log2fdc","1stRPKM","2ndRPKM")
  list(pvals = res$"E-FDR",
       padj = res$"E-FDR",
       beta = res$log2fdc)
}

# function to use my Bayesian analysis given the simulated datasets
runBayes <- function(e, retDDS=FALSE) {
  numRNA = nrow(exprs(e))
  numReplicates = ncol(exprs(e)) / 2
  countsArray = array(exprs(e),dim=c(numRNA,numReplicates,2))
  
  stanDataCurrent <- list(NRNA = numRNA, NReplicates = numReplicates, counts =countsArray)
  
  nbinlogStanFit <- stan("~/diffExpModel/PhDBayesianWork/code/negativeBinomialModel/nbinModel.stan",data = stanDataCurrent,iter = 500, chains = 8,cores = 8)
  
  negBinPosteriorSample <- rstan::extract(nbinlogStanFit)
  
  negBinSampleMuWild <- negBinPosteriorSample$muWild     
  negBinSampleMuMut <- negBinPosteriorSample$muMut
  
  beta <- apply(negBinSampleMuMut,2,median)  - apply(negBinSampleMuWild,2,median) 
  diffMeanCounts <- 2^(negBinSampleMuMut) - 2^(negBinSampleMuWild)
  pvals <- apply(diffMeanCounts, 2, function(x) 1 - abs((sum(x < 0) * 2/ 2000) - 1))
  
  padj <- p.adjust(pvals,method = "fdr")
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=beta))
}

# function to extract pvalues from previous bayesian results
openBayes <- function(e, retDDS=FALSE) {
  
  numReplicates = ncol(exprs(e)) / 2
  k <- as.integer(annotation(e)) - 1
  filenameStanData <- paste0("~/diffExpModel/PhDBayesianWork/data/simulated/autoDownSample",numReplicates,"Replicates",k%%4 + 1,"iter.rds")
  
  load(filenameStanData)
  
  negBinPosteriorSample <- rstan::extract(nbinlogStanFit)
  
  negBinSampleMuWild <- negBinPosteriorSample$muWild     
  negBinSampleMuMut <- negBinPosteriorSample$muMut
  
  beta <- apply(negBinSampleMuMut,2,median)  - apply(negBinSampleMuWild,2,median) 
  diffMeanCounts <- 2^(negBinSampleMuMut) - 2^(negBinSampleMuWild)
  pvals <- apply(diffMeanCounts, 2, function(x) 1 - abs((sum(x < 0) * 2/ 2000) - 1))
  
  padj <- p.adjust(pvals,method = "fdr")
  pvals[is.na(pvals)] <- 1
  pvals[rowSums(exprs(e)) == 0] <- NA
  padj[is.na(padj)] <- 1
  return(list(pvals=pvals, padj=padj, beta=beta))
}