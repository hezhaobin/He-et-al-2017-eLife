## R functions for DESeq -- RNAseq analysis

# define stat_sum_df
stat_sum_df <- function(fun, geom="crossbar", ...) {
  stat_summary(fun.data=fun, colour="red", geom=geom, width=0.2, ...)
}

# define myPlotMD
myPlotMD <- function(fit, sig, i, ylim = c(-6,12)){
  # annotate fit2 object
  genes <- fit$genes$GeneID
  fit$genes$Status <- ifelse( genes %in% scToCg[Xu.genes], "Xu.28", 
                              ifelse( genes %in% scToCg[all.80], "all80", "other") )
  # create plots
  up.down <- table(sig[,i])[c(1,3)]
  plotMD( fit, coef = i, main = i, sub = paste( c("down","up"), up.down, sep=":", collapse = " "), ylim = ylim )
}

# another version of plotMD
myPlotMD.hl <- function(fit, sig, coef, ...){
  up.down <- table( sig[,coef] )[c(1,3)]
  plotMD( fit, coef = coef, ..., status = sig[,coef], 
          hl.col=c("green","red"), hl.cex=0.8, sub = paste( c("down","up"), up.down, sep=":", collapse = " " ) )
}
# define my own plotCounts using ggplot2
myPlotCounts <- function(gene){
  d <- plotCounts( dds1, gene = gene, intgroup = c("Pho4","Pho2"), returnData = T )
  # drop Scer2
  d1 <- subset( d, Pho2 != "Scer2" )
  d1 <- droplevels( d1 )
  p <- ggplot( d1, aes(x=Pho2, y=count) )
  p <- p + geom_point( size=2 ) + ggtitle(label = use.genenames[gene]) + 
    stat_summary( fun.y="mean",geom="point",col="red",size=3 ) + 
    facet_wrap(~Pho4, scales = "free_x", nrow = 1)
  print( p )
  return( p )
}

# for printing in percent format
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x/sum(x), format = format, digits = digits, ...), "%")
}

# define myPlotPHO4()
myPlotPHO4 <- function(fit, sig, pho4, unannotated = F, ...){
  mf <- ifelse( unannotated, 2, 1)
  par(mfrow=c(1, mf), mar=c(4,4,4,2))
  pho4.alone <- paste(pho4,".alone",sep="")
  pho4.w.sc2 <- paste(pho4,".w.cg2",sep="")
  c1 <- abs(sig[,pho4.alone]) 
  c2 <- abs(sig[,pho4.w.sc2]) 
  dat <- as.data.frame(coef(fit))
  x <- which(names(dat) == pho4.alone); y <- which(names(dat) == pho4.w.sc2)
  names(dat)[c(x,y)] <- c("alone","w.sc2")
  xlab = "PHO4 alone"; ylab = "PHO4 w. ScerPho2"
  if( unannotated ){
    with(dat, 
         plot( x = alone, y = w.sc2, pch = 21, cex = 1.2,
               bg = rgb(0,0,0,0.1), col = rgb(0,0,0,0.1), 
               xlab = xlab, ylab = ylab, ... ))
    #abline(0,1,lty=2)
    abline(h=0,v=0,lty=3,col=rgb(0,0,0,0.5))
    lmFit <- lm( w.sc2 ~ alone, dat )
    slope <- coef(lmFit)[2]; r.sqr <- summary(lmFit)$adj.r.squared
    title( main = sprintf("slope = %.2f; adjusted r^2 = %.2f", slope, r.sqr) )
  }
  with(dat, 
       plot( x = alone, y = w.sc2, pch = 16, cex = 1.2,
             col = rgb(c1,0,c2,(c1+c2)/3+0.1), main = pho4, 
             xlab = xlab, ylab = ylab, ... ))
  legend("bottomright",pch=16,col=rgb(c(1,0,1,0),0,c(0,1,1,0)),bty="n",
         legend=c("sig. w/o pho2","sig. w Pho2","sig. in both","not sigificant"))
  #abline(0,1,lty=2)
  abline(h=0,v=0,lty=3,col=rgb(0,0,0,0.5))
}


# definee my own vennDiagram
myVennDiagram <- function(list){
  vennDiagram(list, include = c("up","down"), counts.col = c("red","blue"))
}

# Function for calculating DE genes given the data, design 
PairwiseDE <- function( countTable, design, cond.A, cond.B ){
  index.A <- which( design==cond.A )
  index.B <- which( design==cond.B )
  subTable <- countTable[, c(index.A,index.B)]
  conditions <- factor( c(rep(cond.A,length(index.A)), rep(cond.B, length(index.B))) )
  cds <- newCountDataSet( subTable, conditions )
  cds <- estimateSizeFactors( cds )
  cds <- estimateDispersions( cds )
  res = nbinomTest( cds, cond.A, cond.B )
  res1 <- res[ !is.na(res$padj), ]
  return( res )
}

# Function to report significant genes at a given fdr and log2FoldChange threshold
sigGene <- function( res, fdr=0.05, log2FoldChange=1, direction="up", baseMean=0 ) {
  par <- paste("sigGene --fdr=", fdr, " --foldchange=", 2^log2FoldChange, " --direction=", direction, " --baseMean=",baseMean, sep="", collapse="")
  if( direction != "up" & direction != "down" ){
    list = rep(FALSE, nrow(res))
  }
  else if( direction == "up"){
    list = res$padj < fdr & res$log2FoldChange > log2FoldChange & res$baseMean > baseMean
  }
  else if( direction == "down"){
    list = res$padj < fdr & res$log2FoldChange < -log2FoldChange & res$baseMean > baseMean
  }
  list = replace(list, is.na(list), FALSE)
  attr(list, "par") <- par
  return( list )
}

# Use DESeq2 to call significantly DE genes
PairwiseDE2 <- function( countTable, design, cond.A, cond.B ){
  index.A <- which( design==cond.A )
  index.B <- which( design==cond.B )
  subTable <- countTable[, c(index.A,index.B)]
  conditions <- factor( c(rep(cond.A,length(index.A)), rep(cond.B, length(index.B))), levels=c(cond.A,cond.B) )
  colData <- DataFrame( condition = conditions )
  dds <- DESeqDataSetFromMatrix( countData = subTable, colData = colData, design = ~condition)
  dds <- DESeq( dds )
  return( dds )
}

# Compare two lists
compList <- function(listA, listB, total=0, print=FALSE){
  res <- list(common = intersect(listA, listB), A.only = setdiff(listA, listB), B.only = setdiff(listB, listA))
  len <- sapply(res, length)
  #cat(paste(c("In both"," A only"," B only"), len, sep=":", collapse="\n"))
  print(len)
  if( total ) { # if total # of genes is given, calculate the enrichment factor and p-value from hypergeometric dist.
    expect <- length(listB) / total * length(listA)
    enrich <- round( len["common"] / expect, 1 )
    p.hypergeom <- phyper(q=len["common"],m=length(listB),n=total-length(listB),k=length(listA),lower.tail=F)
    cat("\n")
    cat( sprintf( "Enrichment factor: %.1f\np< %.2g", enrich, p.hypergeom ) )
  }
  return(res)
}
