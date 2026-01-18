#GO, KEGG -> ORA(Over Representation Analysis)
#deg_sig2 = padj < 0.05, abs(lfc) > 1 사용
#Entrez_id 컬럼 붙이기

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

deg_sig2 <- subset(symbol_LFC, padj < 0.05 & abs(log2FoldChange) > 1)
hs <- org.Hs.eg.db
my_symbols <- deg_sig2$symbol
entrez_map <- select(
    hs,
    keys = my_symbols,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "SYMBOL"
) |> distinct(SYMBOL, .keep_all = TRUE)

sum(is.na(entrez_map$ENTREZID)) #Symbol -> Entrez 변환 후 125개 NA

deg_sig2$Entrez <- entrez_map$ENTREZID[match(deg_sig2$symbol, entrez_map$SYMBOL)]
deg_sig2 #9690

up_deg_sig2 <- subset(deg_sig2, log2FoldChange > 0) #7447
down_deg_sig2 <- subset(deg_sig2, log2FoldChange < 0) #2243
universe_genes <- symbol_LFC$symbol

# Up-regulated GO 실행
up_go <- enrichGO(
    gene = up_deg_sig2$symbol,
    universe = universe_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)

dotplot(up_go, showCategory = 20, title = "Up-regulated GO(Cellular Component)")

# Down-regulated GO 실행
down_go <- enrichGO(
    gene = down_deg_sig2$symbol,
    universe = universe_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)

dotplot(down_go, showCategory = 20, title = "Down-regulated GO(Cellular Component)")

# KEGG====================================================

hs2 <- org.Hs.eg.db
my_symbols2 <- symbol_LFC$symbol
entrez_map2 <- select(
    hs2,
    keys = my_symbols2,
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "SYMBOL"
) |> distinct(SYMBOL, .keep_all = TRUE)

sum(is.na(entrez_map2$ENTREZID)) #res_LFC$symbol -> Entrez 변환 후 NA 261개

symbol_LFC$Entrez <- entrez_map2$ENTREZID[match(symbol_LFC$symbol, entrez_map2$SYMBOL)]

# NA 제외 symbol_LFC(27462 - 261)
fil_symbol_LFC <- subset(symbol_LFC, !is.na(Entrez)) #27201

deg_sig3 <- subset(deg_sig2, !is.na(Entrez)) #9565
up_deg_sig3 <- subset(deg_sig3, log2FoldChange > 0) #7342
down_deg_sig3 <- subset(deg_sig3, log2FoldChange < 0) #2223
universe_entrez <- fil_symbol_LFC$Entrez

# Up-regulated KEGG 실행
up_kegg <- enrichKEGG(
    gene = up_deg_sig3$Entrez,
    universe = universe_entrez,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)

dotplot(up_kegg, showCategory = 20, title = "Up-regulated KEGG")

# Down-regulated KEGG 실행
down_kegg <- enrichKEGG(
    gene = down_deg_sig3$Entrez,
    universe = universe_entrez,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
)

dotplot(down_kegg, showCategory = 20, title = "Down-regulated KEGG")

#pathway map작성
library(pathview)

up_kegg$ID[up_kegg$Description == "Cytokine-cytokine receptor interaction"]
pathway_id <- "hsa04060"

fc <- deg_sig3$log2FoldChange
names(fc) <- deg_sig3$Entrez

pathview(
    gene.data = fc,
    pathway.id = pathway_id,
    species = "hsa",
    limit = list(gene = 3),
    out.suffix = "cytokine-cytokine receptor interaction"
)

