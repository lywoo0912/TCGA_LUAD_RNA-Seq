# Feature genes selection using Machine Learning
## Goal: Cluster를 형성했던 sample들의 각 cluster별 주요 영향 genes를 파악한다.
1. Genes set에서 DESeq결과 padj < 0.05 & abs(log2FoldChange) > 1 인 genes들 중 padj 기준 내림차순 정렬 후 상위 top20 genes selection
2. Sample의 각 genes들에 대한 read counts matrix를 VST로 정규화한 데이터로 sample x (gene+cluster) 형태로 ML 데이터 준비(458 rows x 21 cols)
3. RandomForestClassifier를 분류모델로 하여 GridSearch로 세부 파라미터 값 결정(Optuna도 사용해 보았음; score 거의 차이없음) -> prediction accuracy score: 0.66
4. 모델의 Feature Importance을 출력하여 내림차순으로 정렬 후 상위 top10 genes selection
5. Top10 selected genes와 cluster의 z-score heatmap 작성하여 상관관계 파악
6. Top10 selected genes와 cluster의 Kaplan Meier curve & Cox를 통해 생존예후와의 관계 파악
---
- Top10 selected genes: ['EPAS1', 'B3GNT3', 'TOP2A', 'PECAM1', 'TEDC2', 'EMP2', 'S1PR1', 'FAM83A', 'AFAP1-AS1', 'PYCR1']
<br>

<img width="566" height="493" alt="Image" src="https://github.com/user-attachments/assets/dba59c1b-2484-4dc7-86da-0d3009dd6c56" /><br>
- 




