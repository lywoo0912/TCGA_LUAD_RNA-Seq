# Log2(x+1) 정규화 & Consensus/NMF 방식 클러스터링
Goal: 기존 PCA/UMAP 방식 클러스터링과 다른 log2(x+1) normalization과 NMF 클러스터링 방식을 사용한다.<br>
---
1. Gene x Sample raw read counts를 log2(x+1)로 정규화한다.
2. 정규화한 데이터를 threshold = MAD(Mean Absolute Deviation) > 1.5로 필터링하여 feature gene selection을 한다.(2958 genes)
3. Consensus clustering 방식으로 consensus matrix, cdf curve, delta area를 이용하여 k후보들 중 최적 k 탐색
4. 최적 k를 NMF 방식으로 확정
5. 최적 k에 따른 Survival analysis 후 km curve 작성, pathway 탐색
