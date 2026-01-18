library(EnhancedVolcano)
library(dplyr)

res_LFC <- read.csv("LUAD_DESeq2_results_log2FC.csv", row.names = 1)
deg_sig <- subset(as.data.frame(res_LFC), !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1) 

# Up-regulated top20 (padj기준)
up_top20 <- deg_sig |> filter(log2FoldChange > 1) |> arrange(padj) |> slice(1:20)

# Down-regulated top20 (padj기준)
down_top20 <- deg_sig |> filter(log2FoldChange < -1) |> arrange(padj) |> slice(1:20)
genes_top40 <- rbind(up_top20, down_top20)
genes_top40_id <- rownames(genes_top40)


# ID -> Symbol
library(biomaRt)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- gsub("\\..*", "", genes_top40_id)

annotations <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = ensembl_ids,
  mart = ensembl
)

# symbol mapping
symbol_name <- annotations$hgnc_symbol
names(symbol_name) <- annotations$ensembl_gene_id

# genes_top40에 symbol 추가
genes_top40$symbol <- ifelse(
    gsub("\\..*", "", rownames(genes_top40)) %in% names(symbol_name),
    symbol_name[gsub("\\..*", "", rownames(genes_top40))],
    rownames(genes_top40)
)

View(genes_top40$symbol)

EnhancedVolcano(
    genes_top40,
    lab = genes_top40$symbol,
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.05,
    FCcutoff = 1,
    title = "TCGA-LUAD: Tumor vs Normal",
    subtitle = "Top20 Up&Down"
)

#전체에 top20 plot======================================================================

# res_LFC에 symbol 컬럼 추가
library(EnhancedVolcano)
library(ggplot2)

res_LFC$symbol <- ifelse(
    gsub("\\..*", "", rownames(res_LFC)) %in% names(symbol_name),
    symbol_name[gsub("\\..*", "", rownames(res_LFC))],
    rownames(res_LFC)
)

keyvals <- ifelse(
  res_LFC$log2FoldChange < -1 & res_LFC$padj < 0.05, 'blue',
  ifelse(res_LFC$log2FoldChange > 1 & res_LFC$padj < 0.05, 'red',
    'grey'))

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'red'] <- 'Up-regulated'
names(keyvals)[keyvals == 'blue'] <- 'Down-regulated'
names(keyvals)[keyvals == 'grey'] <- 'NS'

EnhancedVolcano(
    res_LFC,
    lab = res_LFC$symbol,
    x = "log2FoldChange",
    y = "padj",
    selectLab = genes_top40$symbol,
    pCutoff = 0.05,
    FCcutoff = 1,
    ylim = c(0, 200),
    colCustom = keyvals,
    colAlpha = 0.6,
    drawConnectors = TRUE,
    widthConnectors = 0.3,
    colConnectors = "grey30",
    max.overlaps = Inf,
    boxedLabels = FALSE,
    title = "TCGA-LUAD: Tumor vs Normal",
    subtitle = "Top20"
)


