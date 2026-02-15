# Log2(x+1) 정규화 & Consensus/NMF 방식 클러스터링
Goal: 기존 PCA/UMAP 방식 클러스터링과 다른 log2(x+1) normalization과 NMF 클러스터링 방식을 사용한다.<br>
---
### Steps
1. Gene x Sample raw read counts를 log2(x+1)로 정규화한다.
2. 정규화한 데이터를 threshold = MAD(Mean Absolute Deviation) > 1.5로 필터링하여 feature gene selection을 한다.(2958 genes)
3. Consensus clustering 방식으로 consensus matrix, cdf curve, delta area를 이용하여 k후보들 중 최적 k 탐색
4. 최적 k를 NMF 방식으로 확정
5. 최적 k에 따른 Survival analysis 후 km curve 작성, pathway 탐색
<br>
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/98af21a1-cb83-41ae-862e-4631fe046adb" /><br>
- Cluster2은 초기부터 생존율이 급감하는 현상을 보이며 예후가 가장 좋지 않았음.<br>
- Cluster3의 생존율은 4~5년차까지는 가장 높았지만 그 이후로는 다소 감소하는 경향을 보임.<br>
- Cluster1의 생존율은 3년차까지는 Cluster2와 유사하게 생존율이 급감하는 현상을 보였으나 그 이후로는 안정적인 현상을 보임.<br>
<br>
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/e015cd30-4ef5-47fd-b1c1-2d697d3028c7" /><br>
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/a9cce195-5ae7-48de-af17-a6c378e61499" /><br>
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/deff99bd-0b1f-4a4a-ba33-4e173bca2ba1" /><br>
