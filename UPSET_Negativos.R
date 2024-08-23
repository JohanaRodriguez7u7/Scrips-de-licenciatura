# Instala y carga las librerías necesarias
if (!require("UpSetR")) install.packages("UpSetR")
library(UpSetR)
library(ggplot2)

# Establece el directorio de trabajo y carga los archivos CSV
setwd("C:/Users/Johana Banana/Documents/GraficasCorelacionadas/Negativos")
PGLS_TMM <- read.csv("PGLSTMMSue_et_al.csv")
PGLS_TMM_Fantom <- read.csv("PGLSTMM_Fantom.csv")
PGLS_Cuartiles <- read.csv("PGLSsuetalCuartiles.csv")
PGLS_Cuartiles_Fantom <- read.csv("PGLSCuartiles_Fantom.csv")
LM_TMM <- read.csv("LMTMMsueetal.csv")
LM_TMM_Fantom <- read.csv("LMTMM_Fantom.csv")
LM_Cuartiles <- read.csv("LMCuartilessueetal.csv")
LM_Cuartiles_Fantom <- read.csv("LMCuartiles_Fantom.csv")

# Unir todos los genes en un solo vector
all_genes <- unique(c(
  LM_TMM$P_Value.positivo, LM_TMM$RT.Value.Negativo,
  LM_TMM_Fantom$P_Value.positivo, LM_TMM_Fantom$RT.Value.Negativo,
  PGLS_TMM$P_Value.positivo, PGLS_TMM$RT.Value.Negativo,
  PGLS_TMM_Fantom$P_Value.positivo, PGLS_TMM_Fantom$RT.Value.Negativo,
  LM_Cuartiles$P_Value.positivo, LM_Cuartiles$RT.Value.Negativo,
  LM_Cuartiles_Fantom$P_Value.positivo, LM_Cuartiles_Fantom$RT.Value.Negativo,
  PGLS_Cuartiles$P_Value.positivo, PGLS_Cuartiles$RT.Value.Negativo,
  PGLS_Cuartiles_Fantom$P_Value.positivo, PGLS_Cuartiles_Fantom$RT.Value.Negativo
))

# Crear una matriz binaria
upset_data <- data.frame(
  Gene = all_genes,
  LM_TMM_P_Value = ifelse(all_genes %in% LM_TMM$P_Value.positivo, 1, 0),
  LM_TMM_RT_Value_Negativo = ifelse(all_genes %in% LM_TMM$RT.Value.Negativo, 1, 0),
  LM_TMM_Fantom_P_Value = ifelse(all_genes %in% LM_TMM_Fantom$P_Value.positivo, 1, 0),
  LM_TMM_Fantom_RT_Value_Negativo = ifelse(all_genes %in% LM_TMM_Fantom$RT.Value.Negativo, 1, 0),
  PGLS_TMM_P_Value = ifelse(all_genes %in% PGLS_TMM$P_Value.positivo, 1, 0),
  PGLS_TMM_RT_Value_Negativo = ifelse(all_genes %in% PGLS_TMM$RT.Value.Negativo, 1, 0),
  PGLS_TMM_Fantom_P_Value = ifelse(all_genes %in% PGLS_TMM_Fantom$P_Value.positivo, 1, 0),
  PGLS_TMM_Fantom_RT_Value_Negativo = ifelse(all_genes %in% PGLS_TMM_Fantom$RT.Value.Negativo, 1, 0),
  LM_Cuartiles_P_Value = ifelse(all_genes %in% LM_Cuartiles$P_Value.positivo, 1, 0),
  LM_Cuartiles_RT_Value_Negativo = ifelse(all_genes %in% LM_Cuartiles$RT.Value.Negativo, 1, 0),
  LM_Cuartiles_Fantom_P_Value = ifelse(all_genes %in% LM_Cuartiles_Fantom$P_Value.positivo, 1, 0),
  LM_Cuartiles_Fantom_RT_Value_Negativo = ifelse(all_genes %in% LM_Cuartiles_Fantom$RT.Value.Negativo, 1, 0),
  PGLS_Cuartiles_P_Value = ifelse(all_genes %in% PGLS_Cuartiles$P_Value.positivo, 1, 0),
  PGLS_Cuartiles_RT_Value_Negativo = ifelse(all_genes %in% PGLS_Cuartiles$RT.Value.Negativo, 1, 0),
  PGLS_Cuartiles_Fantom_P_Value = ifelse(all_genes %in% PGLS_Cuartiles_Fantom$P_Value.positivo, 1, 0),
  PGLS_Cuartiles_Fantom_RT_Value_Negativo = ifelse(all_genes %in% PGLS_Cuartiles_Fantom$RT.Value.Negativo, 1, 0)
)

# Crear el gráfico de UpSet
upset_plot <- upset(
  upset_data,
  sets = c(
    "LM_TMM_P_Value", "LM_TMM_RT_Value_Negativo",
    "LM_TMM_Fantom_P_Value", "LM_TMM_Fantom_RT_Value_Negativo",
    "PGLS_TMM_P_Value", "PGLS_TMM_RT_Value_Negativo",
    "PGLS_TMM_Fantom_P_Value", "PGLS_TMM_Fantom_RT_Value_Negativo",
    "LM_Cuartiles_P_Value", "LM_Cuartiles_RT_Value_Negativo",
    "LM_Cuartiles_Fantom_P_Value", "LM_Cuartiles_Fantom_RT_Value_Negativo",
    "PGLS_Cuartiles_P_Value", "PGLS_Cuartiles_RT_Value_Negativo",
    "PGLS_Cuartiles_Fantom_P_Value", "PGLS_Cuartiles_Fantom_RT_Value_Negativo"
  ),
  sets.bar.color = "#CD5B45",
  main.bar.color = "#8B3E2F",
  matrix.color = "black",
  mainbar.y.label = "Número de Genes Compartidos",
  sets.x.label = "Número de Genes en Cada Análisis",
  order.by = "freq",
  empty.intersections = "on"
)

# Mostrar el gráfico utilizando ggplot2
print(upset_plot)
