library(optiSolve)
library(WeightIt)
library(ggplot2)
library(dplyr)
library(tidyr)

Denom_Linear= 0.3126122
Numer_Linear= 0.3127464
LATE_Linear=1.000429

Denom_nonLinear= 0.31246
Numer_nonLinear= 0.3124363
LATE_nonLinear=0.9999239

results <- read.csv("Result_merge.csv",header = TRUE)

mismatches <- results %>%
  select(kernel, DGP.type, N, SEED, Inference, num, denom, late) %>%
  pivot_wider(
    names_from = Inference, 
    values_from = c(num, denom, late)
  ) %>%
  filter(
    abs(num_SS - num_Boot) > 1e-6 | 
      abs(denom_SS - denom_Boot) > 1e-6 | 
      abs(late_SS - late_Boot) > 1e-6
  )
results=anti_join(results, mismatches, by = c("kernel", "DGP.type", "N", "SEED"))

results <- results %>% mutate(across(where(is.numeric), ~ round(.x, 7)))
results <- unique(results)

freq <- table(results$SEED)
results <- results[results$SEED %in% names(freq[freq == 100]), ]

results=results%>%
  group_by(N,DGP.type,Inference,kernel) %>%
  slice_head(n = 500) %>%
  ungroup()

results <- results %>%
  mutate(
    Lower.num=sqrt(m_num/N)*(Lower.num-num)+num,
    Upper.num=sqrt(m_num/N)*(Upper.num-num)+num,
    Lower.denom=sqrt(m_denom/N)*(Lower.denom-denom)+denom,
    Upper.denom=sqrt(m_denom/N)*(Upper.denom-denom)+denom,
    Lower.late=sqrt(m_late/N)*(Lower.late-late)+late,
    Upper.late=sqrt(m_late/N)*(Upper.late-late)+late
  )


num1=results %>%
  filter(DGP.type==1)%>%
  group_by(kernel, N,Inference) %>%
  mutate(
    cover = (Numer_Linear >= Lower.num & Numer_Linear <= Upper.num),
    length_vec = Upper.num - Lower.num
  ) %>%
  summarise(
    Bias = 100*mean(num - Numer_Linear),
    SE=100*sd(num),
    MSE=Bias^2+SE^2,
    Coverage = mean(cover),
    length = mean(length_vec),
    .groups = "drop"
  )

denom1=results %>%
  filter(DGP.type==1)%>%
  group_by(kernel, N,Inference) %>%
  summarise(
    Bias = 100*mean(denom - Denom_Linear),
    SE=100*sd(denom),
    MSE=Bias^2+SE^2,
    Coverage = mean(Denom_Linear >= Lower.denom & Denom_Linear <= Upper.denom),
    length = mean(Upper.denom - Lower.denom),
    .groups = "drop"
  )

late1=results %>%
  filter(DGP.type==1)%>%
  group_by(kernel, N,Inference) %>%
  summarise(
    Bias = 100*mean(late - LATE_Linear),
    SE=100*sd(late),
    MSE=Bias^2+SE^2,
    Coverage = mean(LATE_Linear >= Lower.late & LATE_Linear <= Upper.late),
    length = mean(Upper.late - Lower.late),
    .groups = "drop"
  )


num2=results %>%
  filter(DGP.type==2)%>%
  group_by(kernel, N,Inference) %>%
  summarise(
    Bias = 100*mean(num - Numer_nonLinear),
    SE=100*sd(num),
    MSE=Bias^2+SE^2,
    Coverage = mean(Numer_nonLinear >= Lower.num & Numer_nonLinear <= Upper.num),
    length = mean(Upper.num - Lower.num),
    .groups = "drop"
  )

denom2=results %>%
  filter(DGP.type==2)%>%
  group_by(kernel, N,Inference) %>%
  summarise(
    Bias =100* mean(denom - Denom_nonLinear),
    SE=100*sd(denom),
    MSE=Bias^2+SE^2,
    Coverage = mean(Denom_nonLinear >= Lower.denom & Denom_nonLinear <= Upper.denom),
    length = mean(Upper.denom - Lower.denom),
    .groups = "drop"
  )

late2=results %>%
  filter(DGP.type==2)%>%
  group_by(kernel, N,Inference) %>%
  summarise(
    Bias = 100*mean(late - LATE_nonLinear),
    SE=100*sd(late),
    MSE=Bias^2+SE^2,
    Coverage = mean(LATE_nonLinear >= Lower.late & LATE_nonLinear <= Upper.late),
    length = mean(Upper.late - Lower.late),
    .groups = "drop"
  )