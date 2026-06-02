
################################################################################
# Simulación de las trayectorias de un modelo de WF de tres alelos
# Autor: Elías González Nieto
# Afil : Facultad de Ciencias - UNAM
# Tesis: Distribuciones Tipo Fase para un Modelo de Coalescencia de Tres Alelos
################################################################################

# Librerías
library(PhaseTypeR)
library(ggplot2)
library(tidyr)
library(dplyr)
################################################################################

# Función para generar una trayectoria de WF con3 alelos
WF_tres_alelos <- function(N, p0, generations) {
  # p0: vector c(pm_1, pm_2, pm_3)
  # devuelve una matriz de generations x 3 con la trayectoria
  
  # La matriz donde guardaremos la trayectoria
  P <- matrix(0, nrow = generations, ncol = 3)
  P[1, ] <- p0 # la inicial es lo que le dimos al inicio (el acomodo inicial)
  extinction_time <- NA
  fixation_time <- NA
  # Por cada generación hacemos la muestra multinomial y actualizamos
  for (t in 2:generations) {
    # muestra multinomial de la generación actual
    counts <- rmultinom(1, size = N, prob = P[t-1, ])
    # revisa si hay extinción
    if (is.na(extinction_time) && any(counts == 0)) {
      extinction_time <- t
    }
    
    # revisa si hay fijación
    if (is.na(fixation_time) && any(counts == N)) {
      fixation_time <- t
    }
    
    # Convertimos a frecuencias
    P[t, ] <- counts / N
    
  }
  
  colnames(P) <-  c(expression(m[1]), expression(m[2]), expression(m[3]))
  WF <- N*P
  print(extinction_time)
  print(fixation_time)
  print(fixation_time - extinction_time)
  return(P)
}

# Ejemplo de uso
set.seed(123)
N <- 100      # tamaño poblacional
p0 <- c(1/3, 1/3, 1/3)
gens <- 50
result <- WF_tres_alelos(N, p0, gens)
print(result)

# Plots de las frecuencias de cada alelo
matplot(result, type = "l", lwd = 2,
        col = c("red", "blue", "darkgreen"),
        ylab = "Frecuencia", xlab = "Generación")

legend("topright",
       legend = c(expression(m[1]), expression(m[2]), expression(m[3])),
       col = c("red", "blue", "darkgreen"),
       lwd = 2)

graficar_trayectoria <- function(tray){
  matplot(tray, type = "l", lwd = 2,
          col = c("red", "blue", "darkgreen"),
          ylab = "Frecuencia", xlab = "Generación")
  
  legend("topright",
         legend = c(expression(m[1]), expression(m[2]), expression(m[3])),
         col = c("red", "blue", "darkgreen"),
         lwd = 2)
  
}
# La misma función para graficar trayectoria pero con ggplot
graficar_trayectoria <- function(tray) {
  
  df <- as.data.frame(tray)
  df$Generacion <- 1:nrow(df)
  colnames(df)[1:3] <- c("m1", "m2", "m3")
  
  df_long <- pivot_longer(df, cols = c("m1","m2","m3"),
                          names_to = "Alelo",
                          values_to = "Frecuencia")
  
  ggplot(df_long, aes(x = Generacion, y = Frecuencia, color = Alelo)) +
    geom_line(size = 0.8) +
    scale_color_manual(values = c("red", "blue", "darkgreen"),
                       labels = c(expression(m[1]), expression(m[2]), expression(m[3]))) +
    labs(x = "Generación", y = "Frecuencia",
         title = "Trayectoria Wright–Fisher de 3 alelos") +
    theme_minimal(base_size = 14)
}

# Ejemplo (Figura 3.1 en la tesis)
set.seed(81)
tesis <- WF_tres_alelos(5000, p0, 6000)
graficar_trayectoria(tesis)
