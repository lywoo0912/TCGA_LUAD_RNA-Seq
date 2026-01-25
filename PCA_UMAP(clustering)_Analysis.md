## PCA Analysis
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/66df0c3d-fba5-4b3c-af44-9f98fb010808" />
<br>
- PC1(13.8%), PC2(7.2%)의 결과가 나옴<br>
- sample_type(tumor vs normal) 기준으로 분류했을 때 PC1: (-100 ~ -50) PC2: (- ~ 50) 사이에서 normal cluster가 형성됨<br>
- PC1: (-50 ~ 50), PC2: (- ~ 50) 구간에서 커다란 tumor cluster가 형성됨<br>
- PC2: (100 >) 구간에서 흩뿌려진 outlier tumor cluster가 형성됨을 확인<br>
- 두 tumor cluster를 구분하는 요인이 무엇인지 확인 필요<br>

## UMAP Analysis
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/0cd8fa8e-bc97-4c7f-9065-6cdcb137cf10" />
<br>
- PCA의 1~30차원으로만 UMAP을 수행<br>
- plot의 오른쪽 커다란 tumor cluster가 형성됨<br>
- plot의 왼쪽에 normal cluster가 형성됨<br>

---
DEG 리스트를 이용해서 deg_pca 30차원 UMAP의 clustering 수행<br>
1. umap 함수의 n_neighbors, min_dist 파라미터의 최적값 탐색(n_neighbors=15, min_dist=0.01)<br>
2. fviz_nbclust 함수의 kmeans, method=("wss", "silhouette", "gap_stat")을 이용해서 최적 k값 선정(k=8)<br>
3. UMAP clustering plot 작성


