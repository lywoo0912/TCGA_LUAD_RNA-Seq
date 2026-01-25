library(DESeq2)
library(SummarizedExperiment)
library(sva)
library(ggplot2)
library(harmony)
library(uwot)
library(plyr)

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix2,
    colData = sample_info2,
    design = ~ sample_type
) #nrow: 60660 ncol: 599

# 저발현 유전자 필터링(counts(dds) >= 10 인 것만)
filtered_gene <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[filtered_gene, ]

dds <- DESeq(dds)
vst_dds <- vst(dds, blind = FALSE)
vst_matrix <- assay(vst_dds) #35909

# batch 효과 적용
batch <- sample_info2$preservation_method
mod <- model.matrix(~sample_type, data = sample_info2)

vst_combat <- ComBat(
    dat = vst_matrix,
    batch = batch,
    mod = mod,
    par.prior = TRUE,
    prior.plots = FALSE
)

write.csv(vst_combat, "vst_combat.csv")
variances <- apply(vst_combat, 1, var)
names(variances[variances == 0 | is.na(variances)])
filtered_vst <- vst_combat[variances > 0, ]
write.csv(filtered_vst, "filtered_vst.csv")

filtered_vst <- read.csv("filtered_vst.csv", row.names = 1)

# PCA 실행
res_pca <- prcomp(
    t(filtered_vst),
    center = TRUE,
    scale. = TRUE
)

# PCA 결과해석
screeplot(res_pca, type = "lines")

pca_df <- data.frame(
    PC1 = res_pca$x[, 1],
    PC2 = res_pca$x[, 2],
    sample_id = rownames(res_pca$x),
    sample_type = sample_info2$sample_type
) #599

# plot 작성  
ggplot(pca_df, aes(x=PC1, y=PC2, color= batch)) +
geom_point(size=2, alpha=0.7) +
labs(title = "PCA: PC1 vs PC2") +
theme_classic() +
labs(title = "PCA: PC1 vs PC2",
    x = paste0("PC1 (", round(summary(res_pca)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(res_pca)$importance[2, 2] * 100, 1), "%)")) +
    theme(plot.title = element_text(hjust = 0.5))


# umap

set.seed(42)

df_scaled <- scale(t(filtered_vst))

res_umap <- umap(
    df_scaled,
    n_neighbors = 15,
    min_dist = 0.1,
    metric = "cosine"
)

umap_df <- data.frame(
    V1 = res_umap[, 1],
    V2 = res_umap[, 2],
    sample_id = rownames(res_umap),
    sample_type = sample_info2$sample_type
)

ggplot(umap_df, aes(x = V1, y = V2, color = sample_type)) +
geom_point(size = 2, alpha = 0.7) +
labs(title = "Umap: Tumor vs Normal", x = "V1", y = "V2") +
theme_classic() +
theme(legend.position = "right")

# PCA 30차원 -> Umap
res_pca_30 <- res_pca$x[, 1:30]
res_umap_30 <- umap(res_pca_30)

umap_30_df <- data.frame(
    V1 = res_umap_30[, 1],
    V2 = res_umap_30[, 2],
    sample_id = rownames(res_umap_30),
    sample_type = sample_info2$sample_type,
    gender = sample_info2$gender,
    status = sample_info2$vital_status
)

# plot 작성
ggplot(umap_30_df, aes(x = V1, y = V2, color = batch)) +
geom_point(size = 2, alpha = 0.7) +
labs(title = "Umap_30: Tumor vs Normal", x = "V1", y = "V2") +
theme_classic() +
theme(legend.position = "right")

#======================================================
# DEG 리스트 이용한 PCA
deg_sig3 <- read.csv("deg_sig3.csv", row.names = 1)
deg_sig_names <- rownames(deg_sig3)
deg_sig_vst <- filtered_vst[deg_sig_names, ]
colnames(deg_sig_vst) <- chartr(".", "-", colnames(deg_sig_vst))
colnames(filtered_vst) <- chartr(".", "-", colnames(filtered_vst))

sample_info2 <- readRDS("sample_info2.rds")
tumor_samples <- colnames(filtered_vst)[sample_info2$sample_type == "Primary Tumor"]
deg_sig_vst_tumor <- deg_sig_vst[, tumor_samples]

deg_pca <- prcomp(
    t(deg_sig_vst_tumor),
    center = TRUE,
    scale. = TRUE
)


deg_pca_df <- data.frame(
    PC1 = deg_pca$x[, 1],
    PC2 = deg_pca$x[, 2],
    sample_id = rownames(deg_pca$x),
    sample_type = subset(sample_info2, sample_type == "Primary Tumor")$sample_type
) #540

ggplot(deg_pca_df, aes(x=PC1, y=PC2, color=sample_type)) +
geom_point(size=2, alpha=0.7) +
theme_classic() +
labs(title = "DEG_PCA: PC1 vs PC2",
    x = paste0("PC1 (", round(summary(deg_pca)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(deg_pca)$importance[2, 2] * 100, 1), "%)")) +
theme(plot.title = element_text(hjust = 0.5))

outlier_samples <- rownames(deg_pca$x)[deg_pca$x[, 2] < -60]
outlier_samples #13


# deg_pca 30차원 umap
deg_pca_30 <- deg_pca$x[, 1:30]
deg_umap_30 <- umap(deg_pca_30)

deg_umap_30_df <- data.frame(
    V1 = deg_umap_30[, 1],
    V2 = deg_umap_30[, 2],
    sample_id = rownames(deg_umap_30),
    sample_type = subset(sample_info2, sample_type == "Primary Tumor")$sample_type ,
    tobacco_status = subset(sample_info2, sample_type == "Primary Tumor")$tobacco_smoking_status 
)

deg_umap_30_df$outlier <- ifelse(
    rownames(deg_umap_30_df) %in% outlier_samples,
    "outlier",
    "Big cluster Tumor"
)

ggplot(deg_umap_30_df, aes(x = V1, y = V2, color = outlier)) +
geom_point(size = 2, alpha = 0.7) +
labs(title = "DEG Umap_30: Tumor", x = "V1", y = "V2") +
theme_classic() +
theme(legend.position = "right")

# UMAP 그래프 clustering
library(factoextra)

deg_pca_30 <- deg_pca$x[, 1:30]

# 최적 n_neighbors, min_dist 탐색
par(mfrow=c(2,2))
for(n in c(15, 50)) {
  for(d in c(0.01, 0.1)) {
    set.seed(42)
    umap_res <- umap(deg_pca_30, n_neighbors=n, min_dist=d)
    plot(umap_res, main=paste("n:", n, "dist:", d))
  }
}

set.seed(42)
deg_umap_30 <- umap(deg_pca_30, n_neighbors=15, min_dist=0.01)
umap_coords <- deg_umap_30[, 1:2] #umap 좌표값

# Determing and visualizing the optimal number of clusters
fviz_nbclust(umap_coords, kmeans, method="wss", k.max=10) #k=5
fviz_nbclust(umap_coords, kmeans, method="silhouette", k.max=10) #k=8
fviz_nbclust(umap_coords, kmeans, method="gap_stat", k.max=10) #k=1

set.seed(42)
kmeans_res <- kmeans(umap_coords, centers=8)
cluster_labels <- kmeans_res$cluster

umap_df <- as.data.frame(umap_coords)

colnames(umap_df) <- c("V1", "V2")
umap_df$smoking_status <- sample_info2[rownames(umap_df), ]$tobacco_smoking_status 

fviz_cluster(kmeans_res,
             data = umap_df,
             geom = "point",
             pointsize = 2,
             ellipse = TRUE,
             ellipse.type = "convex",
             main = "K-means Clustering on UMAP(silhouette, k=8)")

umap_df$cluster <- kmeans_res$cluster
write.csv(umap_df, "umap.df.csv")

find_hull <- function(df) {
    df[chull(df$V1, df$V2), ]
}

hulls <- ddply(umap_df, "cluster", find_hull)

ggplot(umap_df, aes(x = V1, y = V2, color = smoking_status)) +
    geom_point(size = 2.5, alpha = 0.7) +
    geom_polygon(data = hulls, aes(fill = as.factor(cluster)), 
                 alpha = 0.15, color = "black", size = 0.7, linetype = 1) +
    labs(title = "K-means Clustering on UMAP (Smoking Status)",
         x = "V1", y = "V2",
         color = "Smoking Status",
         fill = "Cluster") +
    theme_classic() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
