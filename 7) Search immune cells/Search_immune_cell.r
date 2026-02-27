install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")
library(immunedeconv)

## raw counts -> TPM normalized matrix
Data <- read.csv("count_matrix3.csv", row.names=1)
colnames(Data) <- chartr('.', '-', colnames(Data))

filtered_norm_Data <- read.csv("filtered_norm_Data.csv", row.names=1)
colnames(filtered_norm_Data) <- chartr('.', '-', colnames(filtered_norm_Data))

eco_filtered_norm_Data <- filtered_norm_Data[, rownames(eco_assign)]
filtered_Data <- Data[rownames(eco_filtered_norm_Data), colnames(eco_filtered_norm_Data)]

# gene length 가져오기
library(biomaRt)
library(dplyr)

mygenes <- rownames(filtered_Data)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

gene_info <- getBM(
  attributes = c("hgnc_symbol", "transcript_length"),
  filters = "hgnc_symbol",
  values = mygenes,
  mart = mart
)

gene_length <- gene_info %>% 
              group_by(hgnc_symbol) %>%
              summarise(length=median(transcript_length)) %>%
              filter(hgnc_symbol %in% mygenes)

gene_length <- gene_length[match(mygenes, gene_length$hgnc_symbol), "length", drop=FALSE]

## TPM normalization
library(DGEobj.utils)

filtered_Data_mat <- as.matrix(filtered_Data)
filtered_Data_mat <- apply(filtered_Data_mat, 2, as.numeric)
rownames(filtered_Data_mat) <- mygenes

tpm_conv <- convertCounts(
  countsMatrix = filtered_Data_mat,
  unit = "tpm",
  geneLength = gene_length$length,
  log = FALSE,
  normalize = "none"
)

# cell_type 6개
n <- ncol(tpm_conv)
res_timer <- deconvolute(tpm_conv, method = "timer", indications = rep("LUAD", n))

# quanTIseq: 11개
res_qt <- deconvolute(tpm_conv, method = "quantiseq", tumor = TRUE)

# EPIC: 면역 + 섬유아세포/내피 + uncharacterized; 8개
res_epic <- deconvolute(tpm_conv, method = "epic")

# xCell: input => numeric이어야 함
## "Num. of genes: 1151"
## "ERROR: not enough genes"
"""
library(xCell)
num_eco_filtered_norm_Data <- eco_filtered_norm_Data
num_eco_filtered_norm_Data <- apply(num_eco_filtered_norm_Data, 2, function(x) as.numeric(as.character(x)))
rownames(num_eco_filtered_norm_Data) <- rownames(eco_filtered_norm_Data)
num_eco_filtered_norm_Data <- as.matrix(num_eco_filtered_norm_Data)
res_xc   <- deconvolute(tpm_conv, method = "xcell")
res_xc   <- deconvolute(num_eco_filtered_norm_Data, method = "xcell")
res_xc <- xCellAnalysis(num_eco_filtered_norm_Data, rnaseq=FALSE)
"""

# MCP-counter: 9개
res_mcp  <- deconvolute(tpm_conv, method = "mcp_counter")

# Cibersort
library(CIBERSORT)

res_cib <- cibersort(
  sig_matrix = "/Users/lywoo/Desktop/winter/Renewal/New_norm_Ver/LM22.txt",
  mixture_file = tpm_conv,
  perm = 100,
  QN = FALSE
)

immune_result <- list(
  timer = res_timer,
  quantiseq = res_qt,
  epic = res_epic,
  mcp_counter = res_mcp,
  cibersort = res_cib
)

for (i in names(immune_result)) {
  write.csv(
    immune_result[[i]],
    paste0("res_", i, ".csv")
  )
}

#============================================================

## mycluster별 immune cells 분포 확인
library(dplyr)
library(ggplot2)
library(reshape2)
library(rstatix)
library(ggpubr)

eco_assign <- read.csv("eco_assign.csv", row.names=1)
eco_assign$sample <- rownames(eco_assign)
res_epic <- read.csv("/Users/lywoo/Desktop/winter/Renewal/New_norm_Ver/res_epic.csv", row.names=1)
colnames(res_epic) <- chartr('.', '-', colnames(res_epic))

#method 결과 res 변형 & [sample x (cell types, mycluster)]
res_epic <- as.data.frame(res_epic)
rownames(res_epic) <- res_epic$cell_type
res_epic$cell_type <- NULL
immune_data <- as.data.frame(t(res_epic))
immune_data$sample <- rownames(immune_data)

merge_df <- merge(immune_data, eco_assign[, c("mycluster", "sample")], by="sample")
rownames(merge_df) <- merge_df$sample
merge_df$sample <- NULL
merge_df <- merge_df[rownames(immune_data), ]

cell_types <- setdiff(colnames(merge_df), "mycluster")
df_long <- melt(
  merge_df,
  measure.vars = cell_types,
  variable.name = "Cell_type",
  value.name = "Proportion"
)
df_long$mycluster <- as.factor(df_long$mycluster)
write.csv(df_long, "df_long.csv")

#==============================================================

## Kruskal Wallis(그룹 간 차이 확인, 유의 cell type찾기)
kruskal_test <- df_long %>%
                group_by(Cell_type) %>%
                kruskal_test(Proportion ~ mycluster) %>%
                adjust_pvalue(method="BH") %>%
                add_significance("p.adj")

sig_celltypes <- kruskal_test %>%
                  filter(p.adj < 0.05) %>%
                  pull(Cell_type)

df_sig <- df_long %>% filter(Cell_type %in% sig_celltypes)

## Dunn test(사후검정, 다중비교, 어떤 그룹 간 차이가 있는지 확인)
library(dunn.test)

dunn_test <- df_sig %>%
  group_by(Cell_type) %>%
  dunn_test(Proportion ~ mycluster, p.adjust.method="BH") %>%
  add_xy_position(x="mycluster", scales="free_y")

#===============================================================

## 시각화 (Boxplot)

ggplot(df_sig, aes(x = mycluster, y = Proportion, fill = mycluster)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
    stat_pvalue_manual(
      dunn_test,
      label = "p.adj.signif",
      tip.length = 0.02,
      hide.ns = TRUE,
      size = 3.5
    ) +
    facet_wrap(~ Cell_type, scales = "free_y", ncol = 3) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    labs(
      title = "Significantly Different Immune Cells Across Clusters",
      subtitle = paste("Kruskal-Wallis p.adj < 0.05 (", length(sig_celltypes), "cell types)"),
      x = "Cluster", 
      y = "Cell Proportion",
      fill = "Cluster"
    ) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "lightblue"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )

#====================================================================
