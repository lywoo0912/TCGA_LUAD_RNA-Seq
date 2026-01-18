## GO-ORA
### Up-regulated DEG
- x축: GeneRatio -> (해당 GO term에 속한 DEG 수) / (전체 DEG 수)
- 값이 클수록 "이 GO term에 속한 up-regulated 유전자 비율이 높다"
- y축: GO term 이름
- 그래프 해석:
  1. 세포 분열/증식 관련(organelle fission, nuclear division, meiotic cell cycle process)<br>
    : LUAD 종양 세포가 빠르게 증식하고 있으며, 세포주기와 염색체 분리/핵산 구조 조립 과정이 활발히 up-regulated 되어 있음. 암세포의 비정상적인 세포분열을 나타냄
  2. 감각 자극 감지(detection of stimulus involved in sensory perception, sensory perception of chemical stimulus 등)<br>
    : 후각/화학 감각 수용체 유전자들이 종양에서 비정상적으로 발현되는 현상을 반영할 가능성이 있음
  3. 염색체/염색질 구조(intermediate filament organization, protein localization to chromosome, centromeric region)<br>
    : 세포 분열 관련 단백질의 염색체 배치/구조 조절이 활발, 역시 빠른 증식과 관련

### Down-regulated DEG
- down-GO ORA -> "Normal에서 높던 기능이 Tumor에서 억제되었다"
- 그래프 해석:
  1. 정상 조직의 상처 치유와 항상성(wound healing, chemotaxis, regulation of body fluid levels)<br>
    : 정상 폐 조직은 손상 시 재생, 복구 기능과 정상적인 세포 이동, 체액 조절을 함
  2. 정상 조직의 혈관 조절(regulation of angiogenesis, endothelial cell proliferation, vascular process 등)<br>
    : 정상 폐 조직은 정상적인 혈관 발달과 항상성 유지를 위한 조절 프로그램을 작동시킴
  3. 정상 폐 상피의 섬모 기능(cilium movement, cilium-dependent cell motility, axoneme assembly 등)<br>
    : 정상 기관지 상피 세포는 섬모를 가지고 활발히 움직여서 점액과 이물질을 제거함
  4. 정상 조직의 세포외기질 조직화(extracellular matrix organization, external encapsulating structure organization)<br>
    : 정상 폐 조직은 구조화된 ECM을 유지, 조직
     
