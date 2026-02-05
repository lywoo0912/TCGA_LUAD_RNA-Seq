library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)
sample_info2 <- readRDS("sample_info2.rds")

## 데이터 준비
only_tumor_sample_info2 <- sample_info2[rownames(umap_df), ]
survival_info <- only_tumor_sample_info2[, c("vital_status", "days_to_last_follow_up", "days_to_death")]

survival_info$time <- case_when(
    (survival_info$vital_status == "Alive") ~ survival_info$days_to_last_follow_up,
    (survival_info$vital_status == "Dead") ~ survival_info$days_to_death
) #540

nrow(subset(survival_info, is.na(time))) #78
nrow(subset(survival_info, time == 0)) #4
survival_info <- subset(survival_info, !is.na(time) & time > 0) #458
survival_info$event <- ifelse(
    survival_info$vital_status == "Alive",
    0,
    1
)

final_survival_info <- survival_info[, c("time", "event")]
final_survival_info$cluster <- umap_df[rownames(final_survival_info), ]$cluster
write.csv(final_survival_info, "final_survival_info.csv")

# cluster integrate ver.
integ_cluster <- case_when(
    final_survival_info$cluster == 2 ~ "2/4/6",
    final_survival_info$cluster == 4 ~ "2/4/6",
    final_survival_info$cluster == 6 ~ "2/4/6",
    .default = final_survival_info$cluster
)
integ_final_survival_info <- cbind(final_survival_info, integ_cluster)
write.csv(integ_final_survival_info, "integ_final_survival_info.csv")

## Kaplan Meier
Surv(final_survival_info$time, integ_final_survival_info$event == 1)
survival_fit <- survfit(
    Surv(integ_final_survival_info$time, integ_final_survival_info$event == 1) ~ integ_final_survival_info$integ_cluster, data=integ_final_survival_info)

## Log rank test
survdiff(Surv(integ_final_survival_info$time, integ_final_survival_info$event == 1) ~ integ_final_survival_info$integ_cluster, data=integ_final_survival_info) 

## plot 그리기
ggsurvplot(survival_fit, break.time.by = 365, risk.table = T, fun = 'pct', pval = T, pval.coord = c(0.2, 0.2))

## Cox
age <- sample_info2[rownames(final_survival_info), ]$age_at_index
smoking_status <- sample_info2[rownames(final_survival_info), ]$tobacco_smoking_status


final_plus_age <- cbind(final_survival_info, age)
final_plus_age <- subset(final_plus_age, !is.na(age))
final_plus_smoking <- cbind(final_survival_info, smoking_status)
final_plus_smoking <- subset(final_plus_smoking, !is.na(smoking_status))
final_plus_smoking <- subset(
    final_plus_smoking,
    !smoking_status %in% c(
        "Current Reformed Smoker, Duration Not Specified",
        "Not Reported",
        "Unknown"
    )
)
write.csv(final_plus_age, "final_plus_age.csv") #450
write.csv(final_plus_smoking, "final_plus_smoking.csv") #442

final_plus_age <-read.csv("final_plus_age.csv", row.names=1)
final_plus_age$cluster <- as.factor(final_plus_age$cluster)
# final_plus_age$smoking_status <- as.factor(final_plus_age$smoking_status)

survival_cox <- coxph(
    formula = Surv(time, event == 1) ~ cluster + age,
    data = final_plus_age
)

summary(survival_cox)

ggforest(survival_cox, data=final_plus_age)
