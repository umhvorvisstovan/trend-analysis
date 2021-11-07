#testing tacont() function
library(tidyverse)

df <- read_csv("data/demo_dataset.csv", col_types = cols())
source("taContFunctions.R", encoding = "utf-8")

ylab <- df$contaminant %>% unique()
results <- taCont(df$sampling_year, y = df$measval_num_calc, censored = df$measval_cen_calc)

results[[1]] + ylab(ylab)
