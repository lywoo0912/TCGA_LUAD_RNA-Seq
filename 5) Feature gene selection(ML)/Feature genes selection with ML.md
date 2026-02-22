# Feature genes selection using Machine Learning
## Goal: Cluster를 형성했던 sample들의 각 cluster별 주요 영향 genes를 파악한다.
1. Genes set에서 DESeq결과 padj < 0.05 & abs(log2FoldChange) > 1 인 genes들 중 padj 기준 오름차순 정렬 후 상위 top20 genes selection
2. Sample의 각 genes들에 대한 read counts matrix를 VST로 정규화한 데이터로 sample x (gene+cluster) 형태로 ML 데이터 준비(458 rows x 21 cols)
3. RandomForestClassifier를 분류모델로 하여 GridSearch로 세부 파라미터 값 결정(Optuna도 사용해 보았음; score 거의 차이없음) -> prediction accuracy score: 0.66
4. 모델의 Feature Importance을 출력하여 내림차순으로 정렬 후 상위 top10 genes selection
5. Top10 selected genes와 cluster의 z-score heatmap 작성하여 상관관계 파악
6. Top10 selected genes와 cluster의 Kaplan Meier curve & Cox를 통해 생존예후와의 관계 파악
---
- Top10 selected genes: ['EPAS1', 'B3GNT3', 'TOP2A', 'PECAM1', 'TEDC2', 'EMP2', 'S1PR1', 'FAM83A', 'AFAP1-AS1', 'PYCR1']
<br>
<img width="800" height="700" alt="Image" src="https://github.com/user-attachments/assets/ec4e2bf8-5571-4ae6-a8a4-296664676f55" /><br>
<img width="566" height="493" alt="Image" src="https://github.com/user-attachments/assets/dba59c1b-2484-4dc7-86da-0d3009dd6c56" /><br>

- Cluster 1(가장 나쁜 예후):<br>
  - 'FAM83A'(1.6) -> EGFR/RAS/MAPK pathway 활성화, 항암제 저항성과 연관<br>
  - 'TEDC2'(1.34) -> 중심체 구성, 세포분열 시 미세소관 조직 조절(과발현하면 분열 속도 증가, 중심체 기능 이상 -> 염색체 불안정성 증가)<br>
  - 'B3GNT3'(1.14) -> 혈관 내 침윤(intravasation) 증가<br>
  - 'AFAP1-AS1'(-1.88) -> 전이에 의존하지 않는 혈관 접근성이 높은 고증식 종양(proliferation-driven cancer; "전이형 암"이 아닌 "국소에서 폭주하는 증식형 암"일 가능성)<br>


- Cluster 7(가장 좋은 예후):<br>
  - 'PECAM1'(1.92) -> 세포 접착 유지, 혈관 장벽 안정성 유지<br>
  - 'S1PR1'(1.84) -> 혈관 안정화, Immune cell migration/Angiogenesis regulation pathway<br>
  - 'TOP2A'(-1.06) -> low proliferation
  - 'PYCR1'(-1.42) -> 침윤성 및 증식 증가 유전자의 억제

---

Limitation: 약 9000개의 genes 중 top20개만 모델 훈련에 사용했다는 점에서 모델결과를 해석/단정짓기는 어려움.




