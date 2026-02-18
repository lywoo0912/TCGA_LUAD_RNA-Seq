TCGAbiolinks:::getProjectSummary(("TCGA-LUAD"))
library(TCGAbiolinks)

query_TCGA = GDCquery(
    project = "TCGA-LUAD",
    data.category = 'Transcriptome Profiling',
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)

lihc_res = getResults(query_TCGA)
colnames(lihc_res)
head(lihc_res$sample_type)
summary(factor(lihc_res$sample_type))

# search data in GDC => STAR를 통해 align & count 됨
query_TCGA = GDCquery(
    project = "TCGA-LUAD",
    data.category = 'Transcriptome Profiling',
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
# download from GDC repository
GDCdownload(query_TCGA)

# make R object from the download data
data <- GDCprepare(query_TCGA)

# extract gene expression matrix
library(SummarizedExperiment)
eset <- assay(data) #행렬 추출

# save the matrix as csv
write.csv(eset, file="GE.csv")

sample_info <- as.data.frame(colData(data)) #샘플 메타데이터정보 추출
colnames(sample_info)

nrow(sample_info) # 599(barcode)
ncol(gedata) # 599(barcode)
ncol(sample_info) # 90

gedata <- read.csv('/Users/lywoo/Desktop/winter/GE.csv', row.names = 1)
count_matrix <- as.matrix(gedata) #행렬로 변환
nrow(count_matrix)

filtered_type <- sample_info$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
colnames(count_matrix) <- chartr(".", "-", colnames(count_matrix))

sample_info2 <- sample_info[filtered_type, ]
count_matrix2 <- count_matrix[, sample_info2$barcode] 

stopifnot(all(colnames(count_matrix2) == sample_info2$barcode))

# reference(기준): Solid Tissue Normal -> 0으로 설정
# contrast(대상): Primary Tumor -> 1로 설정
sample_info2$sample_type <- relevel(factor(sample_info2$sample_type), ref = "Solid Tissue Normal")

library(DESeq2)
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix2,
    colData = sample_info2,
    design = ~ sample_type
)

# 저발현 유전자 필터링(counts(dds) >= 10 인 것만)
filtered_gene <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[filtered_gene, ]

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name = "sample_type_Primary.Tumor_vs_Solid.Tissue.Normal")

str(res)

# 결과 테이블 정리/저장
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, "LUAD_DESeq2_results.csv")

# log2FC shrinkage for visualization and ranking
# log2FoldChange = log2(Tumor_평균 / Normal_평균)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="sample_type_Primary.Tumor_vs_Solid.Tissue.Normal", type="apeglm")
write.csv(as.data.frame(resLFC), "LUAD_DESeq2_results_log2FC.csv")
nrow(resLFC) #nrow: 35909

# DEG 리스트 만들기(cut 생성)
deg <- as.data.frame(resLFC)
deg$gene_id <- rownames(deg)
deg_sig <- subset(deg, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(deg_sig, "LUAD_DEG_Cut.csv", row.names=FALSE)
View(read.csv('LUAD_DEG_Cut.csv')) #nrow: 13764

# Visualization(volcano plot)
# 양수(+) -> Tumor > Normal
# 음수(-) -> Tumor < Normal
library(EnhancedVolcano)
draw_plot <- EnhancedVolcano(
    resLFC,
    lab = rownames(resLFC),
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1,
    title = "TCGA-LUAD: Tumor vs Normal",
    subtitle = 'DESeq2',
    caption = 'padj < 0.05, |log2FC| > 1'
)
draw_plot

# png -> pdf
library(ggplot2)
ggsave(
    "LUAD Tumor vs Normal.pdf",
    plot = draw_plot,
    device = "pdf")





