### Epic으로 CE를 설명
### EPIC 결과(세포의 양)를 이용해 CE(세포의 성격/상태)가 실제로 어떤 면역학적 환경을 만들어냈는지 설명하는 것이 핵심입니다.

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(corrplot)
library(RColorBrewer)

## 데이터 준비
res_epic <- read.csv("/Users/lywoo/Desktop/winter/Renewal/New_norm_Ver/res_epic.csv", row.names=1)
rownames(res_epic) <- res_epic$cell_type
res_epic$cell_type <- NULL
t_res_epic <- t(res_epic)
rownames(t_res_epic) <- chartr('.', '-', rownames(t_res_epic))

eco_assign <- read.csv("eco_assign.csv", row.names=1)

merge_eco_epic <- cbind(t_res_epic, eco_assign$Carcinoma.Ecotype)
colnames(merge_eco_epic)[9] <- "Ecotype"
merge_eco_epic <- as.data.frame(merge_eco_epic)


## Heatmap
immune_cells <- c(
    "B cell",
    "Cancer associated fibroblast",
    "T cell CD4+",
    "T cell CD8+",
    "Endothelial cell",
    "Macrophage",
    "NK cell",
    "uncharacterized cell"
)
merge_eco_epic <- merge_eco_epic %>%
                  mutate(across(all_of(immune_cells), as.numeric))

heatmap_df <- merge_eco_epic %>%
              group_by(Ecotype) %>%
              summarise(across(all_of(immune_cells), mean, na.rm=TRUE)) %>%
              column_to_rownames("Ecotype")
heatmap_df <- heatmap_df[c("CE1", "CE2", "CE3", "CE5", "CE6", "CE7", "CE8", "CE9", "CE10"), ]

scaled_heatmap_df <- scale(heatmap_df)

pheatmap(t(scaled_heatmap_df), 
         cluster_cols = FALSE,      # CE 순서 유지
         cluster_rows = TRUE,       # 세포끼리 클러스터링
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = "white",
         fontsize = 12,
         angle_col = 45,
         main = "Immune Signature by Carcinoma Ecotype (Z-score)",
         display_numbers = FALSE,
         legend = TRUE)

#=======================================================================

## Cluster2의 Significant CAF 유전자 찾기(caf fraction을 연속변수로 design)
c2s <- subset(eco_assign, mycluster=="2") # 104
colnames(res_epic) <- chartr(".", "-", colnames(res_epic))
res_epic_c2s <- res_epic[, rownames(c2s)]
write.csv(res_epic_c2s, "res_epic_c2s.csv")
caf_frac <- as.numeric(res_epic_c2s["Cancer associated fibroblast", ])
names(caf_frac) <- colnames(res_epic_c2s)


library(DESeq2)
library(EnhancedVolcano)
library(biomaRt)

gedata <- read.csv('/Users/lywoo/Desktop/winter/GE.csv', row.names = 1)
count_matrix <- as.matrix(gedata)
colnames(count_matrix) <- chartr(".", "-", colnames(count_matrix))

countdata <- count_matrix[, names(caf_frac)]
coldata <- data.frame(row.names = colnames(countdata))
coldata$caf <- caf_frac
coldata$caf_c <- scale(coldata$caf, center = TRUE, scale = FALSE)  

dds <- DESeqDataSetFromMatrix(
    countData = countdata,
    colData = coldata,
    design = ~ caf_c
)
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="caf_c", type="apeglm")

#resLFC ensembl id -> symbol
"""
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- gsub("\\..*", "", row.names(resLFC))

annotations <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = ensembl_ids,
  mart = ensembl
)

symbol_name <- annotations$hgnc_symbol
names(symbol_name) <- annotations$ensembl_gene_id

resLFC$symbol <- ifelse(
    gsub("\\..*", "", rownames(resLFC)) %in% names(symbol_name),
    symbol_name[gsub("\\..*", "", rownames(resLFC))],
    rownames(resLFC)
)

#symbol로 변환 안된 것 제외
symbolLFC <- resLFC[!grepl("^ENSG0", resLFC$symbol, ignore.case = TRUE), ]
nrow(resLFC) #60660
nrow(symbolLFC) #59402

#symbol == ""인 것 제외
count(symbolLFC, symbol == "")
symbolLFC <- symbolLFC |> filter(!symbol %in% "")
nrow(symbolLFC) 

#중복 symbol 제거
symbolLFC <- symbolLFC[!duplicated(symbolLFC$symbol), ]
nrow(symbolLFC) 
"""

draw_volcano <- EnhancedVolcano(
    resLFC,
    lab=resLFC$symbol,
    x="log2FoldChange",
    y="padj",
    pCutoff=0.05,
    FCcutoff=2,
    title="Significant genes of CAF",
    caption="padj < 0.05, |Log2fc| > 2"
)

draw_volcano


sig_deg <- subset(as.data.frame(resLFC), !is.na(padj) & padj <= 0.05 & abs(log2FoldChange) > 2) # 2788
sig_deg <- sig_deg[!(sig_deg$symbol == ""), ] # 2199

sig_high_caf_genes <- subset(sig_deg, log2FoldChange > 2) # 1055
sig_high_caf_genes <- sig_high_caf_genes[!(sig_high_caf_genes$symbol == ""), ] # 955

sig_low_caf_genes <- subset(sig_deg, log2FoldChange < 2) # 1733
sig_low_caf_genes <- sig_low_caf_genes[!(sig_low_caf_genes$symbol == ""), ] # 1244

## CAF genes pathway
library(msigdbr)
library(fgsea)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

# high&low caf hallmark
symbolLFC <- resLFC[!(resLFC$symbol == ""), ]
symbolLFC <- symbolLFC[!duplicated(symbolLFC$symbol), ]
symbolLFC <- symbolLFC[!grepl("^ENSG0", symbolLFC$symbol, ignore.case = TRUE), ]

gene_list <- symbolLFC$log2FoldChange
names(gene_list) <- symbolLFC$symbol
gene_list <- sort(gene_list, decreasing=TRUE)

res_fgsea <- fgsea(
    pathways = hallmark_list,
    stats = gene_list,
    minSize = 15,
    maxSize = 500,
    scoreType = "std",
    nPermSimple = 100000
)

high_caf_path <- subset(res_fgsea, NES > 0)
low_caf_path <- subset(res_fgsea, NES < 0)
head(high_caf_path[order(padj, -NES)], 15)
head(low_caf_path[order(padj, NES)], 15)

## high_caf_path의 공통 leadingEdge 탐색
high_top15 <- head(high_caf_path[order(padj, -NES)], 15)
leading_list <- high_top15$leadingEdge
all_genes <- unlist(leading_list)
gene_freq <- table(all_genes)
common_genes <- names(gene_freq[gene_freq >= 5]) #"CCND2"(5), "NRP1"(5), "SERPINE1"(6)
high_top15$pathway[
    sapply(high_top15$leadingEdge, function(x) "SERPINE1" %in% x)
]

## sig_high_caf_genes와의 overlap gene -> SERPINE1
common_genes %in% sig_high_caf_genes$symbol #SERPINE1

## SERPINE1의 예후와의 관계(cluster2 vs cluster3)
library(survival)
library(survminer)
library(dplyr)

c23s <- subset(eco_assign, mycluster == "2" | mycluster == "3")
res_epic_c23 <- res_epic[, rownames(c23s)]
c23_caf_frac <- as.numeric(res_epic_c23["Cancer associated fibroblast", ])
names(c23_caf_frac) <- colnames(res_epic_c23)

countdata2 <- count_matrix[, names(c23_caf_frac)]
coldata2 <- data.frame(row.names = colnames(countdata2))
coldata2$caf <- c23_caf_frac
coldata2$caf_c <- scale(coldata2$caf, center = TRUE, scale = FALSE) 

dds <- DESeqDataSetFromMatrix(
    countData=countdata2,
    colData=coldata2,
    design=~caf_c
)

dds <- DESeq(dds)
vst_norm <- vst(dds, blind=FALSE)
serpine1_norm <- assay(vst_norm)["ENSG00000106366.9", ]

matrix_new <- read.csv("matrix_new.csv", row.names=1)
c23_matrix <- matrix_new[rownames(subset(eco_assign, mycluster=="2" | mycluster=="3")), ]
c23_matrix$SERPINE1 <- serpine1_norm[rownames(c23_matrix)]

res.cut <- surv_cutpoint(
    c23_matrix,
    time="time",
    event="event",
    variables="SERPINE1"
) #cutpoint = 13.11695, statistic = 3.303263
surv_cat <- surv_categorize(res.cut)
c23_matrix$SERPINE1_level <- surv_cat$SERPINE1


## KM curve
survival_fit <- survfit(
    Surv(
        surv_cat$time,
        surv_cat$event == 1
    ) ~ surv_cat$SERPINE1,
    data=surv_cat
)

survdiff(
    Surv(
        surv_cat$time,
        surv_cat$event == 1
    ) ~ surv_cat$SERPINE1,
    data=surv_cat
)
ggsurvplot(survival_fit, break.time.by = 365, risk.table = T, fun = 'pct', pval = T, pval.coord = c(0.2, 0.2))

## Cox regression(SERPINE1 & stage와의 관계)
serpine_cox <- coxph(
    formula=Surv(time, event==1) ~ stage + SERPINE1_level,
    data=c23_matrix
)

summary(serpine_cox)
