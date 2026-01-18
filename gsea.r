library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)

# res_LFC(Ensembl ID -> symbol)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- gsub("\\..*", "", row.names(res_LFC))

annotations <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = ensembl_ids,
  mart = ensembl
)

symbol_name <- annotations$hgnc_symbol
names(symbol_name) <- annotations$ensembl_gene_id

res_LFC$symbol <- ifelse(
    gsub("\\..*", "", rownames(res_LFC)) %in% names(symbol_name),
    symbol_name[gsub("\\..*", "", rownames(res_LFC))],
    rownames(res_LFC)
)

#symbol로 변환 안된 것 제외
symbol_LFC <- res_LFC[!grepl("^ENSG0", res_LFC$symbol, ignore.case = TRUE), ]
nrow(res_LFC) #35909
nrow(symbol_LFC) #35232

#symbol == ""인 것 제외
count(symbol_LFC, symbol == "")
symbol_LFC <- symbol_LFC |> filter(!symbol %in% "")
nrow(symbol_LFC) #27466

#중복 symbol 제거
symbol_LFC <- symbol_LFC[!duplicated(symbol_LFC$symbol), ]
nrow(symbol_LFC) #27462
#========================================================

# gene_list준비
gene_list <- symbol_LFC$log2FoldChange
names(gene_list) <- symbol_LFC$symbol
gene_list <- sort(gene_list, decreasing = TRUE)


# msigdb에서 hallmark gene set data 준비
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)

# fgsea 
res_fgsea <- fgsea(
    pathways = hallmark_list,
    stats = gene_list,
    minSize = 15,
    maxSize = 500,
    scoreType = "std",
    nPermSimple = 100000
)
head(res_fgsea[order(pval), ], 10)

topPathwaysUp <- res_fgsea[ES > 0][head(order(NES), n=10), pathway]
topPathwaysDown <- res_fgsea[ES < 0][head(order(NES), n=10), pathway]
topPathways <- c(topPathwaysDown, topPathwaysUp)


# 1) Top pathway만 추출
top_res <- res_fgsea[pathway %in% topPathways]


# 음수를 아래, 양수를 위로 재배치
top_res <- rbind(
  top_res[NES < 0],  # 절댓값 작은 음수가 맨 아래
  top_res[NES > 0]   # 절댓값 작은 양수가 0 위
)

# 4) Pathway 이름 정리 (HALLMARK_ 제거)
library(stringr)
top_res[, pathway_clean := str_remove(pathway, "HALLMARK_")]
top_res[, pathway_clean := str_replace_all(pathway_clean, "_", " ")]
top_res[, pathway_clean := factor(pathway_clean, levels = pathway_clean)]

# 5) Lollipop plot 그리기
ggplot(top_res, aes(x = NES, y = pathway_clean)) +
  # 선 (0 → NES)
  geom_segment(aes(x = 0, xend = NES, 
                   y = pathway_clean, yend = pathway_clean),
               color = "red", linewidth = 1) +
  # 점 (NES 위치)
  geom_point(color = "red", size = 3) +
  # 기준선 (x = 0)
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  # 축 라벨
  labs(
    x = "Normalized Enrichment Score",
    y = NULL,
    title = "GSEA - Hallmark Pathways"
  ) +
  # 테마
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )


# enrichment score plot of top3 abs(NES)

library(fgsea)
library(ggplot2)

top3 <- res_fgsea |> filter(padj < 0.05) |> arrange(desc(abs(NES))) |> head(3)

get_top5_label <- function(pathway_name, res_fgsea, gene_list, n = 5) {
  leading_genes <- res_fgsea[pathway == pathway_name]$leadingEdge[[1]]
  top_genes <- head(sort(gene_list[leading_genes]), n)   # 가장 음수/작은 쪽 n개
  label <- paste(names(top_genes), collapse = ", ")
  paste("Leading genes:", label)
}

for (i in 1:3) {
  pw <- top3$pathway[i]
  
  sub_txt <- get_top5_label(pw, res_fgsea, gene_list, n = 5)
  p <- plotEnrichment(hallmark_list[[pw]], gene_list) +
    labs(title = gsub("HALLMARK_", "", pw), subtitle = sub_txt) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 11, colour = "blue"),
      axis.title = element_text(size = 11)
    )
  
  print(p)
}
