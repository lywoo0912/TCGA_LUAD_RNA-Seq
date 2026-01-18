library(DESeq2)
library(SummarizedExperiment)
library(sva)
library(ggplot2)
library(harmony)

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

variances <- apply(vst_combat, 1, var)
names(variances[variances == 0 | is.na(variances)])
filtered_vst <- vst_combat[variances > 0, ]

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
library(uwot)
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
deg_sig_names <- rownames(deg_sig3)
deg_sig_vst <- filtered_vst[deg_sig_names, ]

deg_pca <- prcomp(
    t(deg_sig_vst),
    center = TRUE,
    scale. = TRUE
)



deg_pca_df <- data.frame(
    PC1 = deg_pca$x[, 1],
    PC2 = deg_pca$x[, 2],
    sample_id = rownames(deg_pca$x),
    sample_type = sample_info2$sample_type
)

ggplot(deg_pca_df, aes(x=PC1, y=PC2, color= sample_type)) +
geom_point(size=2, alpha=0.7) +
theme_classic() +
labs(title = "DEG_PCA: PC1 vs PC2",
    x = paste0("PC1 (", round(summary(deg_pca)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(deg_pca)$importance[2, 2] * 100, 1), "%)")) +
    theme(plot.title = element_text(hjust = 0.5))

outlier_samples <- rownames(deg_pca$x)[deg_pca$x[, 2] > 100]
outlier_samples 

# DEG 리스트 이용한 UMAP
set.seed(42)
deg_sig_umap <- umap(
    scale(t(deg_sig_vst)),
    n_neighbors = 15,
    min_dist = 0.1,
    metric = "cosine"
)

deg_sig_umap_df <- data.frame(
    V1 = deg_sig_umap[, 1],
    V2 = deg_sig_umap[, 2],
    sample_id = rownames(deg_sig_umap),
    sample_type = sample_info2$sample_type
)

ggplot(deg_sig_umap_df, aes(x = V1, y = V2, color = sample_type)) +
geom_point(size = 2, alpha = 0.7) +
labs(title = "DEG Umap: Tumor vs Normal", x = "V1", y = "V2") +
theme_classic() +
theme(legend.position = "right")

# deg_pca 30차원 umap
deg_pca_30 <- deg_pca$x[, 1:30]
deg_umap_30 <- umap(deg_pca_30)

deg_umap_30_df <- data.frame(
    V1 = deg_umap_30[, 1],
    V2 = deg_umap_30[, 2],
    sample_id = rownames(deg_umap_30),
    sample_type = sample_info2$sample_type,
    tobacco_status = sample_info2$tobacco_smoking_status
)

ggplot(deg_umap_30_df, aes(x = V1, y = V2, color = tobacco_status)) +
geom_point(size = 2, alpha = 0.7) +
labs(title = "DEG Umap_30: Tumor vs Normal", x = "V1", y = "V2") +
theme_classic() +
theme(legend.position = "right")
