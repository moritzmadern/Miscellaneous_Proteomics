## DE testing function using limma ###
DE_limma <- function(m,
                     groupname1 = groupname1,
                     groupname2 = groupname2,
                     groups,
                     trend_limma=FALSE,
                     batch=batch){
  
  # convert groups to a factor
  factor_groups <- factor(groups)
  
  # convert batch to a factor
  factor_batch <- factor(batch)
  
  # define contrast vector x
  x <- paste0(groupname2, "-", groupname1)
  
  ## test limma when no batch is specified
  if(is.null(batch)){
    design <- model.matrix(~0 +factor_groups)
    colnames(design) <- sub(colnames(design), pattern="factor_groups", replacement = "")
    fit <- lmFit(m, design)
    contrast_matrix <- makeContrasts(contrasts = x, levels=colnames(design))
    fit_contrasts <- contrasts.fit(fit, contrast_matrix)
    fit_contrasts <- eBayes(fit_contrasts)
    LIMMAresults <- topTable(fit_contrasts, number = Inf, sort.by = "none")
  }
  
  # test limma when batch is specified
  if(!is.null(batch)){
    design <- model.matrix(~0 +factor_groups)
    colnames(design) <- sub(colnames(design), pattern="factor_groups", replacement = "")
    dupcor <- duplicateCorrelation(m,design,block=factor_batch)
    fit <- lmFit(m, design, block= factor_batch, correlation = dupcor$consensus)
    contrast_matrix <- makeContrasts(contrasts = x, levels=colnames(design))
    fit_contrasts <- contrasts.fit(fit, contrast_matrix)
    fit_contrasts <- eBayes(fit_contrasts, trend = trend_limma, robust = TRUE)
    LIMMAresults <- topTable(fit_contrasts, number = Inf, sort.by = "none")
  }
  
  # exctract limma results (fc, p.val and adj.p.val) for plotting
  fc <- LIMMAresults$logFC
  fc[is.na(fc)] <- 0
  p_val <- LIMMAresults$P.Value
  p_val[is.na(p_val)] <- 1
  
  # plot volcano plot (p-value vs fc)
  par(mfrow=c(1,1))
  par(mar=c(4,4,4,5))
  par(mgp=c(2.5,0.8,0))
  par(xpd=TRUE)
  plot(x= fc, y= -log(p_val,base=10), pch=19, cex=1, xaxt="n", yaxt="n", xlab= paste0(groupname2, " / ", groupname1, " \nfold change [log2]" ), ylab = "- log10 (p-value)",main= "Volcano Plot \n p-val", cex.lab=0.7, cex.main= 0.8, font.lab=2, xlim=c(-max(abs(fc)),max(abs(fc))), bty="n", col=rgb(red=100, green=100, blue=100, alpha=50, maxColorValue = 255))
  axis(side=1, cex.axis=0.7)
  axis(side=2, cex.axis=0.7, las=2, mgp=c(2.5, 0.8, 0))
  
  
  # change colnames of results table to indicate specific comparison
  colnames(LIMMAresults) <- paste0(colnames(LIMMAresults), "__", x)
  
  # return LIMMAresults
  return(LIMMAresults)
}







## Build example dataframe
df = data.frame(WT_Th1.1 = rnorm(50),
                WT_Th1.2 = rnorm(50),
                WT_Th1.3 = rnorm(50),
                KO_Th1.1 = rnorm(50),
                KO_Th1.2 = rnorm(50),
                KO_Th1.3 = rnorm(50),
                WT_Th2.1 = rnorm(50),
                WT_Th2.2 = rnorm(50),
                WT_Th2.3 = rnorm(50),
                KO_Th2.1 = rnorm(50),
                KO_Th2.2 = rnorm(50),
                KO_Th2.3 = rnorm(50))

## Define groups
groups = c("WT_Th1", "WT_Th1", "WT_Th1",
           "KO_Th1", "KO_Th1", "KO_Th1",
           "WT_Th2", "WT_Th2", "WT_Th2",
           "KO_Th2", "KO_Th2", "KO_Th2")
           
## Convert quantitative dataframe to matrix
m_from_df<- as.matrix(df)


## Artificially create one highly significant datapoint
m_from_df[1,1] <- 10.1
m_from_df[1,2] <- 10.0
m_from_df[1,3] <- 10.2


## Apply function to m_from_df to test WT_Th1 vs KO_Th1. Set trend = FALSE, and no batch (not necessary anymore after batch-correction)
library(limma)
limma_res <- DE_limma(m = m_from_df,
                       groupname1 = "WT_Th1", groupname2 = "KO_Th1",
                       groups = groups,
                       trend_limma = FALSE,
                       batch = NULL)


## Add results to your input dataframe of interest via cbind()
df_merged <- cbind(df, limma_res)





