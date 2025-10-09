# Librairies ####
library(tidyverse)

# Parameters ####
setwd(dir = "D:/Louis/Transcriptome/Regeneration_analysis/")

# Functions ####
load_matrix <- function(filename){
  return(read_csv(file = filename, col_names = TRUE,))
}

# Main ####

matrix_0dpa = load_matrix("Table/0dpa/0dpa_TMM_matrix.csv")

matrix_1dpa = load_matrix("Table/1dpa/1dpa_TMM_matrix.csv")

matrix_2dpa = load_matrix("Table/2dpa/2dpa_TMM_matrix.csv")

matrix_3dpa = load_matrix("Table/3dpa/3dpa_TMM_matrix.csv")

matrix_5dpa = load_matrix("Table/5dpa/5dpa_TMM_matrix.csv")

matrix_NA = load_matrix("Table/NA/NA_TMM_matrix.csv")

prop_0_1 = c()
prop_1_5 = c()
prop_5_10 = c()
prop_10_inf = c()
conditions = c("0 dpa", "1 dpa", "2 dpa", "3 dpa", "5 dpa", "NA")

list_matrix = list(matrix_0dpa, matrix_1dpa, matrix_2dpa, matrix_3dpa, matrix_5dpa, matrix_NA)

for (i in 1:6) {
  i = list_matrix[[i]]
  median_matrix = i[,ncol(i)]
  total = nrow(median_matrix)
  prop_0_1 = c(prop_0_1, length(which(median_matrix < 1))/total*100)
  prop_1_5 = c(prop_1_5, length(which(median_matrix >= 1 & median_matrix < 5))/total*100)
  prop_5_10 = c(prop_5_10, length(which(median_matrix >= 5 & median_matrix < 10))/total*100)
  prop_10_inf = c(prop_10_inf, length(which(median_matrix >= 10))/total*100)
}

df = tibble(Conditions = conditions, `> 0 TMM` = prop_0_1, `> 1 TMM` = prop_1_5,
            `> 5 TMM` = prop_5_10, `> 10 TMM` = prop_10_inf)
df = df %>% pivot_longer(c(`> 0 TMM`, `> 1 TMM`, `> 5 TMM`, `> 10 TMM`),
                         names_to = 'Limits', values_to = "Proportion")

df$Limits = factor(df$Limits, levels = c("> 0 TMM", "> 1 TMM", "> 5 TMM", "> 10 TMM"))

ggplot(data = df, mapping = aes(x = Conditions, y = Proportion, fill = Limits)) +
  geom_bar(position = "stack", stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("#696969", "#088ac7", "#8708c7", "#c70b08")) + 
  theme_light()+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(size = 0.4)
  ) +
  xlab("") +
  ylab("Proportion of genes (%)") +
  labs(fill = "Threshold")
  scale_y_continuous(limits = c(0,100),labels = scales::label_percent(scale = 1))

#