# Cluster prediction using saliva exosome genes with Machine Learning
## Goal: 타액에서 나온 유전자를 통해 LUAD 예후 예측을 할 수 있다.
---
1. Feature selection된 genes와 saliva exosome genes의 공통 genes를 추출(65개 genes)<br>
2. 예후 양상이 유사했던 cluster1, cluster3을 합치고 cluster2와 DEG 수행<br>
3. padj <= 0.05, |Log2FC| > 2인 genes로 필터링 하여 그 중 65개 genes가 있는지 확인(High-regulated[cluster2]: 4개, Down-regulated[cluster1&3]: 5개)<br>
4. 9개의 feature genes로 예후 예측이 가능한 cluster 분류 machine learning 수행<br>
<br>



