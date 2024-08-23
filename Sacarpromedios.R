# Carga de datos (reemplaza 'tu_archivo.csv' con el nombre de tu archivo)
setwd("C:/Users/Johana Banana/Documents/Brier/Modificados")

# Leer el archivo CSV con identificadores de m
datos <- read.csv("promedio.csv", check.names = FALSE)

# Limpiar los nombres de las columnas
colnames(datos) <- make.names(colnames(datos))

# Función para promediar columnas con el mismo nombre
promediar_columnas <- function(data) {
  columnas <- names(data)
  columnas_unicas <- unique(columnas)
  
  for (col in columnas_unicas) {
    columnas_con_nombre <- which(columnas == col)
    
    if (length(columnas_con_nombre) > 1) {
      promedio <- rowMeans(data[, columnas_con_nombre])
      data <- cbind(data, promedio)
      colnames(data)[ncol(data)] <- paste("Promedio_", col, sep = "")
    }
  }
  
  return(data)
}

# Aplicar la función al dataframe
datos_con_promedios <- promediar_columnas(datos)
write.csv(datos_con_promedios, "datos_con_promedios3.csv", row.names = FALSE)

