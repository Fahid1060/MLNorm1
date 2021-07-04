


funCleaning<-function(datamatrix)
{
  library(M3Drop)
  library(Seurat)
  
ds<-as.matrix(datamatrix)
ds<-M3DropCleanData(datamatrix,is.counts = FALSE)
ds<-ds$data

ds<-CreateSeuratObject(counts  = ds)
ds <- NormalizeData(ds)
ds<-FindVariableFeatures(ds) 
ds<-ds@assays$RNA@var.features
ds<-as.matrix(ds)
ds<-as.data.frame(ds)

common <- intersect(ds$V1, rownames(datamatrix))
common<-datamatrix[common,]
common<-as.matrix(common)

return (common)
}


funLimma<-function(datamatrix){
  library(limma)
  n<-ncol(datamatrix)/2
  #Create limma design matrix (first 5 columns are tumor, last 5 columns are normal)
  designmatrix <- c(rep(1,n),rep(2,n))
  designmatrix <- model.matrix(~ 0+factor(designmatrix))
  colnames(designmatrix) <- c("group1", "group2")
  rownames(designmatrix) <- colnames(datamatrix)
  
  contrast.matrix <- makeContrasts("group1-group2", levels=designmatrix)
  
  fit2 <- lmFit(datamatrix,designmatrix)
  fit2 <- contrasts.fit(fit2, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  p.Values<-length(which(fit2$p.value<0.05))
  cat("p.Values :",p.Values)
  
  p<-as.numeric(fit2$p.value)
  adjusted<-as.data.frame(p.adjust(p,method="fdr",n=length(p)))
  adjusted<-as.data.frame(adjusted)
  
  q.Values<-length(which(adjusted<0.05))
  cat("\nq.Values :",q.Values)
  
  DE<-p.Values-q.Values
  cat("\nDifference :",DE,"\n")
}




