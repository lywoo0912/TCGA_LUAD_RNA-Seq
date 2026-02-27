# Cluster prediction using saliva exosome genes with Machine Learning
## Goal: 타액에서 나온 유전자를 통해 LUAD 예후 예측을 할 수 있다.
---
1. Feature selection된 genes와 saliva exosome genes의 공통 genes를 추출(65개 genes)<br>
2. 예후 양상이 유사했던 cluster1, cluster3을 합치고 cluster2와 DEG 수행<br>
3. padj <= 0.05, |Log2FC| > 2인 genes로 필터링 하여 그 중 65개 genes가 있는지 확인(Up-regulated[cluster2]: 4개, Down-regulated[cluster1&3]: 5개)<br>
4. 9개의 feature genes로 예후 예측이 가능한 cluster 분류 machine learning 수행<br>
### DEG between Cluster1&3 / Cluster2
<img width="800" height="800" alt="Image" src="https://github.com/user-attachments/assets/2ca851b4-7170-451f-898b-f8bec298779c" /><br>

- Up-regulated genes(cluster2): "S100A9", "KRT16", "SBSN", "KRT6A"<br>
- Down-regulatd genes(cluster1&3): "SLPI", "AQP5", "ALB", "LRRC26", "MUC5AC"<br>
### Randomforest model classification report
<img width="500" height="500" alt="Image" src="https://github.com/user-attachments/assets/2178022e-55ab-4f93-af70-62134233f4bf" /><br>

- 모델의 최종 예측 정확도: 0.88<br>
- cluster1&3의 f1-score가 높은것이 cluster1&3 : cluster2 = 2 : 1 정도의 비율인 것을 감안해서 과적합 되었다고 판단해도 cluster2의 f1-score도 충분히 높은 값으로 나오는 것으로 보아 전체적인 모델의 성능은 준수한 것으로 확인됨.<br>
### AUC Curve & ROC Score
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/e5cef059-c5c3-4430-8c94-a3ccd12a9288" /><br>
### Barplot of feature importances
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/d38e2905-c59a-4eeb-b303-c6496fb4f0d0" /><br>

- Randomforest모델이 9개 genes의 가중치를 계산한 feature importance<br>
- SBSN: 약 0.28, KRT16: 약 0.166으로 전체의 약 44%차지함.<br>
- 이 두 genes는 Cluster2에서 Up-regulated된 유전자였음
---
## Cluster1&3, Cluster2 비율 1:1 & 임상데이터 추가 Ver
