# Immunedeconv를 이용한 mycluster별 Immunce cell 분포 파악
### Immunedeconv란?<br>
: Ecotyper와 달리 면역세포 비율을 계산하는데 특화된 R 패키지, bulk RNA-seq에서 면역·기질 세포 타입별 fraction/score를 추정한다.<br>
---
1. Immunedeconv는 input으로 TPM normalized matrix를 요구함.<br>
2. method(timer, quanTIseq, epic, MCP-counter, cibersort)별로 탐색 수행<br>
3. 수행 후 null값이 가장 적은 method = "epic"방법 선택<br>
4. Epic 결과 8개의 cell types를 Kruskal-Wallis & Dunn test를 통해 유의한 cell types만 추출<br>
5. mycluster별 boxplot을 작성하여 cluster별 cell types 차이 확인<br>
<br>
<img width="800" height="600" alt="Image" src="https://github.com/user-attachments/assets/5da793af-74ab-4be7-ba0c-b3e1357f3f47" /><br>
<br>

### Epic으로 CE 해석(Heatmap)<br>
EPIC 결과(세포의 양)를 이용해 CE(세포의 성격/상태)가 실제로 어떤 면역학적 환경을 만들어냈는지 설명하는 것이 핵심입니다.<br>
<img width="800" height="600" alt="Image" src="https://github.com/user-attachments/assets/466eb200-3b4f-465f-8729-7ac46d062b1f" /><br>

### Cell types별 특징<br>
- Cancer associated fibroblast (cluster2 highest)<br>
: a cell type within the TME that promotes tumorigenic features by initiating the remodeling of the extracellular matrix or by secreting cytokines.<br>
- T cell CD4+ / T cell CD8+ (cluster2 lowest, cluster3 high)<br>
: key players in the immune response against both pathogenic infections and cancer.<br>
- Endothelial cell (cluster2 highest)<br>
: angiogenesis<br>
- Macrophage (cluster2 highest)<br>
: destroy germs, damaged cells and cancer cells; promote tissue repair and healing.<br>
- Nk(Natural killer) cell (cluster3 highest)<br>
: destroy infected and diseased cells, like cancer cells.<br>



