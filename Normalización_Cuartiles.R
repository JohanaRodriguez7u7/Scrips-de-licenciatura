# Ajusta la ruta y el nombre del archivo según tu caso
setwd("C:/Users/Johana Banana/Documents/fantom/Modifiaciones/Normalizacion")

data <- read.csv("promediofilasfantomconsumatoriamayor_a_diez.csv", row.names = 1)

# Cargar la librería limma si aún no está cargada
if (!requireNamespace("limma", quietly = TRUE)) {
  install.packages("limma")
}
library(limma)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

# Diseño de experimentos (asumo que 'tejidos' es una columna en tu conjunto de datos)


# Antes de la normalización
par(mfrow = c(1, 2))
boxplot(data, main = "Antes de la Normalización", col = c("red", "blue"))
abline(h = 0, col = "black", lty = 2)

# Normalizar los datos usando el método de cuantiles
normalized_data <- normalizeBetweenArrays(
  object = data,
    # Ajusta según tu caso
  method = "quantile",
  
)

# Después de la normalización
boxplot(normalized_data, main = "Después de la Normalización", col = c("red", "blue"))
abline(h = 0, col = "black", lty = 2)
write.csv(normalized_data,"Briernormalizado.csv")
