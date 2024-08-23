# Instala y carga el paquete limma
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")

library(limma)
library(Biobase)

# Establece el directorio de trabajo
setwd("C:/Users/Johana Banana/Documents/Atahualpa/Listas")

# Lee los datos de expresión
expresion <- read.csv("promedioadultoSU.csv", row.names = 1)

# Convierte el data frame a una matriz
expresion_rounded <- round(as.matrix(expresion))

# Crea un objeto ExpressionSet
eset <- ExpressionSet(assayData = expresion_matrix)

# Normaliza los datos con voom() y normalizeQuantiles()
normalized_eset <- voom(eset, normalize.method = "quantile")

# Ahora, los datos están normalizados en 'normalized_eset'.
# Puedes acceder a los datos normalizados usando:

normalized_exprs <- normalized_eset$E

# Y guardar los datos normalizados si lo deseas
write.csv(normalized_exprs, "datos_normalizados_limma.csv")

# Además, puedes continuar con tu análisis de expresión diferencial utilizando limma u otros métodos disponibles en el paquete.
