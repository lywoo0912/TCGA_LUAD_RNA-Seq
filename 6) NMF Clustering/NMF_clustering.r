final_survival_info <- read.csv("final_survival_info.csv", row.names = 1)
Data <- read.csv("count_matrix3.csv", row.names = 1)
deg_sig3 <- read.csv("deg_sig3.csv", row.names = 1)
colnames(Data) <- chartr(".", "-", colnames(Data))
time <- final_survival_info[colnames(Data), "time"]
status <- final_survival_info[colnames(Data), "event"]

## log(x+1)로 정규화
norm_Data <- log2(Data + 1)

## mad기준 feature selection
mads <- apply(norm_Data, 1, mad)
hist(mads, breaks=100, main="MAD Distribution")
mad_threshold <- 1.5
idx <- which(mads >= mad_threshold)
filtered_norm_Data <- norm_Data[idx, ] #2958
write.csv(filtered_norm_Data, "filtered_norm_Data.csv")

# ecotyper ver data
library(tibble)
filtered_norm_Data <- as.data.frame(filtered_norm_Data)
colnames(filtered_norm_Data) <- chartr("-", ".", colnames(filtered_norm_Data))
filtered_norm_Data2 <- rownames_to_column(filtered_norm_Data, "Gene")
write.table(filtered_norm_Data2, "filtered_norm_Data2.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## Clustering(Consensus Matrix & CDF)
library(ConsensusClusterPlus)

ExecuteCC <- function(clusterNum,
    d, maxK = 10, clusterAlg = "hc",
    distance = "pearson", title = "Consensus Clustering Result",
    reps=500, pItem=0.8, pFeature=1, plot="png",
    innerLinkage="average", finalLinkage="average", writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL, verbose=FALSE,corUse="everything", seed=42){
  if(is.list(d))
  {
    temp=NULL
    for(i in 1: length(d))
    {
      temp=rbind(temp,d[[i]])
    }
    temp=t(scale(t(temp)))
  }
  else
   temp=d
  originalResult=ConsensusClusterPlus(
      temp, maxK=maxK,clusterAlg=clusterAlg,
      distance=distance,title=title,
      reps=reps, pItem=pItem, pFeature=pFeature,plot=plot,
      innerLinkage=innerLinkage, finalLinkage=finalLinkage,
      writeTable=writeTable,weightsItem=weightsItem,weightsFeature=weightsFeature,
      verbose=verbose,corUse=corUse)
  group=originalResult[[clusterNum]][["consensusClass"]]
  distanceMatrix=originalResult[[clusterNum]][["consensusMatrix"]]
  attr(distanceMatrix,'class')="Similarity"
  #icl=calcICL(result,title =fileName,plot="png" )
  result=list(group=group,distanceMatrix=distanceMatrix,originalResult=originalResult)
  result
}

find_k <- ExecuteCC(clusterNum = 3, d = filtered_norm_Data, seed=42)
find_k2 <- ExecuteCC(clusterNum = 4, d = filtered_norm_Data, seed=42)

## silhouette valiation
silhouette_SimilarityMatrix<-function(group, similarity_matrix, seed=42) {
  similarity_matrix=as.matrix(similarity_matrix)
  similarity_matrix<-(similarity_matrix+t(similarity_matrix))/2
  diag(similarity_matrix)=0
  normalize <- function(X) X / rowSums(X)
  similarity_matrix<-normalize(similarity_matrix)
  
  n <- length(group)
  if(!all(group == round(group))) stop("'group' must only have integer codes")
  cluster_id <- sort(unique(group <- as.integer(group)))
  k <- length(cluster_id)
  if(k <= 1 || k >= n)
    return(NA)
  doRecode <- (any(cluster_id < 1) || any(cluster_id > k))
  if(doRecode)
    group <- as.integer(fgroup <- factor(group))
  cluster_id <- sort(unique(group))
  
  wds <- matrix(NA, n,3, dimnames =list(names(group), c("cluster","neighbor","sil_width")))  
  for(j in 1:k)
  { 
    index <- (group == cluster_id[j])
    Nj <- sum(index)
    wds[index, "cluster"] <- cluster_id[j]
    dindex <- rbind(apply(similarity_matrix[!index, index, drop = FALSE], 2,
                          function(r) tapply(r, group[!index], mean)))
    maxC <- apply(dindex, 2, which.max)
    wds[index,"neighbor"] <- cluster_id[-j][maxC]
    s.i <- if(Nj > 1) {
      a.i <- colSums(similarity_matrix[index, index])/(Nj - 1)
      b.i <- dindex[cbind(maxC, seq(along = maxC))]
      ifelse(a.i != b.i, (a.i - b.i) / pmax(b.i, a.i), 0)
    } else 0
    wds[index,"sil_width"] <- s.i
  }
  attr(wds, "Ordered") <- FALSE
  class(wds) <- "silhouette"
  wds
}

sil_k3 <- silhouette_SimilarityMatrix(find_k$group, find_k$distanceMatrix, seed=42)
plot(sil_k3)
sil_k4 <- silhouette_SimilarityMatrix(find_k2$group, find_k2$distanceMatrix, seed=42)
plot(sil_k4)

## Clustering(NMF)
library(NMF)

final_nmf_result <- nmf(filtered_norm_Data, rank=3, method="brunet", nrun=30)
best_result <- final_nmf_result
W_matrix <- basis(best_result) 
H_matrix <- coef(best_result) 

D <- apply(W_matrix, 2, max)
W_norm <- sweep(W_matrix, 2, D, FUN="/")
H_norm <- sweep(H_matrix, 1, D, FUN="*")
colnames(W_norm) <- c('c1', 'c2', 'c3')
rownames(H_norm) <- c('c1', 'c2', 'c3')

matrix_new <- as.data.frame(apply(H_norm, 2, which.max), colnames(H_norm))
names(matrix_new) <- c("cluster")
W_norm_T <- t(W_norm)
integ_final_survival_info <- read.csv("integ_final_survival_info.csv", row.names=1)
matrix_new$time <- integ_final_survival_info[rownames(matrix_new), ]$time
matrix_new$event <- integ_final_survival_info[rownames(matrix_new), ]$event
write.csv(matrix_new, "matrix_new.csv")

# consensusmap(
#   best_result,
#   annCol = matrix$cluster,
#   tracks = c("consensus"),
#   main = "Consensus Matrix(k=3)"
# )

H_norm_T <- t(H_norm)
nmf_cluster <- apply(H_norm_T, 1, which.max)
dist_samples <- dist(H_norm_T, method="euclidean")
sil_nmf <- silhouette(nmf_cluster, dist_samples)
plot(sil_nmf)
mean(sil_nmf[, "sil_width"])

## Survival Analysis(Kaplan Meier)
library(survival)
library(survminer)
library(dplyr)
ce_matrix_new <- matrix_new[rownames(eco_assign), ] # sample 334개

Surv(
  ce_matrix_new$time,
  ce_matrix_new$event == 1
)

survival_fit <- survfit(
  Surv(
  ce_matrix_new$time,
  ce_matrix_new$event == 1
) ~ ce_matrix_new$cluster,
data = ce_matrix_new
)

## Log rank test
survdiff(
  Surv(
  ce_matrix_new$time,
  ce_matrix_new$event == 1
) ~ ce_matrix_new$cluster,
data = ce_matrix_new
)

ggsurvplot(survival_fit, break.time.by = 365, risk.table = T, fun = 'pct', pval = T, pval.coord = c(0.2, 0.2))

# c1/2 차이 검증 
c1_c2_matrix <- subset(matrix_new, cluster == 1 | cluster == 2)

Surv(
  c1_c2_matrix$time,
  c1_c2_matrix$event == 1
)

survival_fit <- survfit(
  Surv(
  c1_c2_matrix$time,
  c1_c2_matrix$event == 1
) ~ c1_c2_matrix$cluster,
data = c1_c2_matrix
)

## Log rank test
survdiff(
  Surv(
  c1_c2_matrix$time,
  c1_c2_matrix$event == 1
) ~ c1_c2_matrix$cluster,
data = c1_c2_matrix
)
ggsurvplot(survival_fit, break.time.by = 365, risk.table = T, fun = 'pct', pval = T, pval.coord = c(0.2, 0.2))

## Pathway
library(clusterProfiler)
library(org.Hs.eg.db)

ef <- extractFeatures(
    best_result,
    method = c('kim', 'max'),
    format = c('list', 'combine', 'subset'),
    nodups = TRUE
)

c1_genes <- rownames(best_result)[ef[[1]]]
c2_genes <- rownames(best_result)[ef[[2]]]
c3_genes <- rownames(best_result)[ef[[3]]]

C_go <- enrichGO(
        gene = c3_genes,
        universe = rownames(best_result),
        OrgDb = org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )

dotplot(C_go, showCategory = 20, title = paste0("C3 NMF GO (Biological Process)"))

set.seed(42)
C_kegg <- enrichKEGG(
    gene = deg_sig3$Entrez[deg_sig3$symbol %in% c3_genes],
    universe = as.character(deg_sig3$Entrez[deg_sig3$symbol %in% rownames(filtered_norm_Data)]),
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)
dotplot(C_kegg, showCategory = 20, title = paste0("C3 NMF KEGG"))
