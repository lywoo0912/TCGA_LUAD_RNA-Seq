# Finding significant gene in CAF that affect poor prognosis of cluster2 
## Goal: cluster2의 예후가 가장 불량하였고 그 중 CAF 비율이 가장 컸으므로, CAF에서 cluster2의 예후를 악화시킨 중요 유전자를 탐색한다.<br>
---
1. Cluster2의 sample들에서 CAF fraction을 추출<br>
2. CAF fraction을 연속변수로 하여 DESeq 수행<br>
3. padj < 0.05 & |Log2FC| > 2로 filtering하여 volcano plot 작성, sig_high_caf_genes/sig_low_caf_genes로 구분<br>
4. High-regulated genes의 상위 15개 HALLMARK pathway 수행 후, 각 pathway에서 공통적으로 많이 들어있는 genes 추출("CCND2"(5), "NRP1"(5), "SERPINE1"(6))<br>
5. 이 genes와 sig_high_caf_genes와의 overlap을 확인("SERPINE1")<br>
6. 실제로 SERPINE1이 예후에 결정적인 영향을 끼치는 유전자인지 알아보기 위해 CAF 비율이 가장 높았던 cluster2와 CAF 비율이 가장 낮았던 cluster3 sample들에서 CAF fraction 추출 후 Survival analysis 진행<br>
7. Cox regression을 통해 SERPINE1과 stage와의 관계 파악
---
<img width="800" height="600" alt="Image" src="https://github.com/user-attachments/assets/7fc73a74-cb50-4b7b-885d-cc458d05c113" /><br>
<br>
<img width="800" height="600" alt="Image" src="https://github.com/user-attachments/assets/a06c04e5-2e37-4b22-9a28-4952431549fe" /><br>
- 실제로 Km curve에서 SERPINE1 high가 SERPINE1 low보다 예후가 급격히 나빠지는 경향을 보임.<br>
<br>

