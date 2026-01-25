# Survival Analysis
## Kaplan-Meier Analysis & Log rank test
<img width="600" height="600" alt="Image" src="https://github.com/user-attachments/assets/8610c188-ee6a-4674-9ce6-c9bf16c94adc" /><br>
- cluster 2, cluster 7의 기울기가 완만함 -> 예후 더 좋음<br>
- 나머지 cluster는 기울기가 급격히 감소함 -> 예후 안 좋음<br>
## Cox Regression & Forest plot
<figure>
    <img src='https://github.com/user-attachments/assets/dc4c9683-753f-49d3-9f7c-e35dad58047b' alt='missing' width="600" height="600"/>
    <figcaption>[🎯Only cluster]</figcaption>
</figure>
<br>
<br>

<figure>
    <img src='https://github.com/user-attachments/assets/6ca849d8-1477-4256-ac62-1d7e3e46fa16' alt='missing' width="600" height="600"/>
    <figcaption>[🎯Cluster with age]</figcaption>
</figure>
<br>
- cluster 2, 3, 4, 5, 7이 reference인 cluster 1보다 HR < 1이고, 특히 2와 7은 HR ≈ 0.25~0.45 수준이라 "명확한 low-risk 그룹"<br>
- p-value도 2, 3, 4, 5, 7은 유의함(특히 2, 7은 매우 작음)<br>
- age는 HR ≈ 1이고 p-value > 0.1 이므로 생존과 유의한 관계가 없음<br>
<br>
  
<figure>
    <img src='https://github.com/user-attachments/assets/7b502e3f-1dfc-44a5-9842-21c72ae408da' alt='missing' width="600" height="600"/>
    <figcaption>[🎯Cluster with smoking_status]</figcaption>
</figure>
<br>
- smoking_status의 p-value값들이 대체로 높다 -> 생존에 유의한 영향을 끼치지 않음<br>
- 오히려 Lifelong Non Smoker의 HR이 Current Smoker보다 높음<br>
- UMAP에서 확인했듯이 only cluster의 HR값과 cluster with smoking_status의 HR값의 각 cluster별 HR값 차이가 거의 없으므로 smoking_status와 cluster는 관계가 없음

