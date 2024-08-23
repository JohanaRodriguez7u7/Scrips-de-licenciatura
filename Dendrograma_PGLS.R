
# Instalar y cargar el paquete ape (si no está instalado)
if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}

library(ape)

# Leer tus datos (ajusta la ruta y formato según tu caso)
#setwd("C:/Users/Johana Banana/Documents/Atahualpa/DENDROGRAMA")
setwd("C:/Users/Johana Banana/Documents/fantom/Modifiaciones/Dendrograma")

data <- read.csv("Normalizacion_Cuartiles_Fantom.csv", row.names = 1)


# Transponer los datos si es necesario
data_transpuesta <- t(data)

# Calcular la matriz de distancias
dist_matrix <- dist(data_transpuesta)

# Aplicar clustering jerárquico
hclust_result <- hclust(dist_matrix)

# Cortar el dendrograma para obtener grupos (ajusta la altura según tus necesidades)
cut_tree_result <- cutree(hclust_result, h = 100)

# Ordenar las columnas según el resultado del clustering
data_ordenada <- data_transpuesta[, order(cut_tree_result)]

# Graficar el dendrograma y el mapa de calor
pdf("C:/Users/Johana Banana/Documents/fantom/Modifiaciones/Dendrograma.pdf")
plot(hclust_result, main = "Dendrograma de Clustering Jerárquico", xlab = "Tejidos", cex = 0.3)

# Convertir el resultado del clustering jerárquico a un objeto "phylo"
phy_tree <- as.phylo(hclust_result)

# Guardar el dendrograma en un archivo Newick (.nwk) (ajusta la ruta y el nombre del archivo según tu caso)
write.tree(phy_tree, file = "C:/Users/Johana Banana/Documents/fantom/Modifiaciones/Dendrograma.nwk")

dev.off()
