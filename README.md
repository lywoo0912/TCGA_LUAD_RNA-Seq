# TCGA_LUAD_RNA-Seq


## test.r + test_dev.r(DEG Analysis)
1. GDCquery(TCGA)에서 workflow.type = "STAR-Counts", sample.type = (Primary Tumor, Solid Tissue Normal)로 다운로드
- STAR - Count: RNA-seq reads를 GRCh38에 STAR로 정렬한 뒤 gene 단위로 집계한 raw read count
2. reference(기준): Solid Tissue Normal -> 0, contrast(대상): Primary Tumor -> 1로 설정하여 DESeq2 진행 
3. 저발현 유전자 필터링(각 유전자가 최소 10개 샘플에서 count가 10 이상인 것만 남김)
4. Log2FC shrinkage for visualization & ranking
5. padj < 0.05 & |log2FoldChange| > 1으로 Volcano plot 작성
6. Up-regulated와 Down-regulted 각각 top20을 ensembl ID -> Symbol로 변환하고 labeling하여 최종 Volcano plot 작성
<br>
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/44c1696f-c7a2-43f4-94d4-d2210668e9b1" /><br>
- Tumor에서 발현된 유전자의 양이 Normal Tissue에서 발현된 유전자의 양보다 월등히 많음<br>
- Top20_Up-regulated 유전자: FAM83A, PYCR1, AFAP1-AS1, TEDC2 등<br>
- Top20_Down-regulated 유전자: OTUD1, EPAS1, STX11 등

## gsea.r(GSEA Analysis)
1. fgsea를 수행하기 위한 gene_list 준비; symbol컬럼 필터링(symbol로 변환 안된 것, ""인 것, 중복 symbol 제외)
2. Barplot of gene_list 작성
3. msigdb에서 hallmark gene set 으로 gene_list와 align
4. GSEA table plot 작성
<br>
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/2efecacb-d66c-437c-99ac-22895c81ab6d" /><br>
- gene_list: res_LFC - (ensembl ID → symbol 변환 안된 것 + symbol == “” + 중복 symbol)<br>
- x축: gene_list의 symbol name, y축: log2FoldChange 기준으로 유전자 발현량 차이를 보여줌<br>
- PSG1, RAX, F7, SIX1 등이 tumor에서 강하게 발현된 유전자들이다.<br>
- KLB, SPN 등이 normal에서 강하게 발현된 유전자들이다.<br>
- log2FC > 0: 같은 유전자가 tumor에서 normal일 때보다 발현율이 더 높다.<br>
- log2FC < 0: 같은 유전자가 normal에서 tumor일 때보다 발현율이 더 높다.<br>
<br>
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/aec411a3-9ed1-4d27-a2b3-df33497a0cae" /><br>
NES(Normalized Enrichment Score) 해석<br>
양수 NES: 해당 pathway의 유전자들이 ranked list의 왼쪽(상위, log2FC 양수)에 몰려 있음
- G2M_CHECKPOINT, E2F_TARGETS → 세포 증식/분열 관련
- PANCREAS_BETA_CELLS → 췌장 베타 세포 기능 관련
- SPERMATOGENESIS, GLYCOLYSIS → 대사 활성화
- MYC_TARGETS → 암관련 신호<br>

음수 NES: 해당 pathway의 유전자들이 ranked list의 오른쪽(하위, log2FC 음수)에 몰려 있음
- TGF_BETA_SIGNALING, APOPTOSIS → 정상세포 분화/죽음 경로 억제
- INTERFERON_ALPHA_RESPONSE, IL6_JAK_STAT3_SIGNALING → 면역 반응 억제<br>
pval, padj 낮을수록 유의하다.


