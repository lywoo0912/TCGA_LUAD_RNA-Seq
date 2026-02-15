### Goal: my cluster와 ecotype cluster(CEs) 비교
library(dplyr)
library(tidyr)
library(ggplot2)

matrix_new <- read.csv("matrix_new.csv", row.names=1) #sample 개수: 458
eco_assign <- read.delim("/Users/lywoo/Desktop/ecotyper_output/Carcinoma_Ecotypes/Ecotype_Assignment.txt", header=TRUE, sep='\t', dec='.', row.names=1) #sample 개수: 334
rownames(eco_assign) <- chartr(".", "-", rownames(eco_assign))
eco_assign$mycluster <- matrix_new[rownames(eco_assign), ]$cluster
sample_info4 <- readRDS("sample_info2.rds")
matrix_new$stage <- sample_info4[rownames(matrix_new), ]$ajcc_pathologic_stage
write.csv(matrix_new, "matrix_new.csv")
eco_assign$stage <- sample_info4[rownames(eco_assign), ]$ajcc_pathologic_stage
write.csv(eco_assign, "eco_assign.csv")

## CEs & mycluster 간 비율
eco_prop <- eco_assign |>
  dplyr::count(mycluster, Carcinoma.Ecotype, name = "n") |>
  dplyr::group_by(mycluster) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::ungroup() 

eco_prop <- as.data.frame(eco_prop)

eco_prop_sorted <- eco_prop %>%
  dplyr::group_by(mycluster) %>%
  dplyr::arrange(dplyr::desc(prop), .by_group = TRUE) %>%
  dplyr::ungroup()

heat_df <- eco_prop %>% 
            filter(Carcinoma.Ecotype %in% paste0("CE", 1:10)) %>%
            dplyr::select(Carcoma=Carcinoma.Ecotype, mycluster, prop) %>%
            pivot_wider(
                names_from = mycluster,
                values_from = prop,
                values_fill = list(prop=0)
            ) %>%
            tibble::column_to_rownames("Carcoma") %>%
            as.data.frame()

CE_list <- c("CE1", "CE2", "CE3", "CE5", "CE6", "CE7", "CE8", "CE9", "CE10")
heat_df <- heat_df[CE_list, ]

## heatmap
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales) # 퍼센트 변환용

# 1. 데이터 준비 (Long format 변환)
plot_df <- eco_prop %>%
  # (1) CE1~CE10만 남기기
  filter(Carcinoma.Ecotype %in% CE_list) %>%
  mutate(Carcinoma.Ecotype = factor(Carcinoma.Ecotype, levels = CE_list)) %>%
  # (2) mycluster 문자열로 변환 (필요시)
  mutate(mycluster = as.character(mycluster)) %>%
  mutate(
      y_label = case_when(
      mycluster == "1" ~ "1 (N=63)",
      mycluster == "2" ~ "2 (N=104)",
      mycluster == "3" ~ "3 (N=167)",
      TRUE ~ mycluster
    )
  ) %>%
  # (5) 빠진 조합(N=0, Prop=0) 채우기
  complete(y_label, Carcinoma.Ecotype, fill = list(n = 0, prop = 0)) %>%
  # (6) 텍스트 라벨 만들기 (n값과 % 함께 표시)
  mutate(text_label = ifelse(n == 0, 0, paste0(n, "\n(", percent(prop, accuracy = 0.1), ")")))

# 2. 히트맵 그리기
ggplot(plot_df, aes(x = Carcinoma.Ecotype, y = y_label)) +
  # 타일 그리기 (색상은 prop 기준)
  geom_tile(aes(fill = prop), color = "white", lwd = 0.5) +
  # 텍스트 넣기 (n과 prop)
  geom_text(aes(label = text_label), size = 3.5, color = "black") +
  # 색상 스케일 설정 (흰색 -> 빨간색)
  scale_fill_gradient(low = "white", high = "#FF6B6B", labels = percent) +
  # 축 및 범례 설정
  labs(x = "Carcinoma Ecotype", 
       y = "My Cluster (Sum per cluster)",
       title = "Proportion of mycluster & CEs",
       fill = "Proportion") +
  # 테마 깔끔하게
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),       # 격자 제거
    axis.text.x = element_text(angle = 45, hjust = 1) # X축 글자 기울기
  )


