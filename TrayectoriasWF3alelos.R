library(PhaseTypeR)
library(ggplot2)
library(tidyr)
library(dplyr)

subintensity_matrix <- matrix(c(-1.5, 0, 0,
                                1.5, -1, 0,
                                0, 1, -0.5), ncol = 3)
initial_probabilities <- c(0.9, 0.1, 0)
ph <- PH(subintensity_matrix, initial_probabilities)
ph

matrix2 <- matrix(c(-1, 0, 0.5,
                    0.7, -1, 0,
                    0.2, 0, -0.5), ncol = 3)
pi2 <- c(0.3, 0.7, 0)
ph2 <- PH(matrix2, pi2)
summary(ph2)
summary(ph)

recompensa <- matrix(c(1, 1, 1,
                       1, 1, 1,
                       1, 1, 1), ncol = 3)

media1 = mean(ph)
media2 = mean(ph2)

rFullPH(ph)
rFullPH(ph2)

matrix3 <- matrix(c(-1, 0.6, 0.1,0.2,
                    0, -1, 0.9, 0,
                    0.2, 0.7, -1, 0,
                    0.9, 0, 0, -1), ncol = 4, byrow = TRUE)
ph3 <- PH(matrix3)
rFullPH(ph3)
summary(ph3)
media3 <- mean(ph3)


contador <- 0
for(i in 1:N){tray1 <- rFullPH(ph3)
tau <- sum(unlist(tray1$time))
contador <- contador + tau}
contador / N
# ok esto funciona bien para aproximar la media (es monte carlo)
N <- 1000
contador2 <- 0
media2 <- mean(ph2)
media3 <- mean(ph3)
for(i in 1:N){
  tray2 <- rFullPH(ph2)
  tray3 <- rFullPH(ph3)
  tau2 <- sum(unlist(tray2$time))
  tau3 <- sum(unlist(tray3$time))
  contador2 <- contador2 + (tau2-media2)*(tau3-media3)
}
contador2/N

# La idea sería tener la misma trayectoria del proceso en 3 dimensiones

# LA VARIABLE ALEATORIA I

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
    # check extinction
    if (is.na(extinction_time) && any(counts == 0)) {
      extinction_time <- t
    }
    
    # check fixation
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
#set.seed(123)
N <- 100      # tamaño poblacional
p0 <- c(1/3, 1/3, 1/3)
gens <- 50

result <- WF_tres_alelos(N, p0, gens)
print(result)


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
# La misma función para graficar tray pero con ggplot
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

set.seed(81)
tesis_example <- WF_tres_alelos(5000, p0, 6000)
graficar_trayectoria(tesis_example)

# EJEMPLO ÚTIL PARA TESIS #
res <- WF_tres_alelos_actualizado(5000, c(0.3,0.3,0.4), 12000)

attr(res, "extinction_time")
attr(res, "fixation_time")
graficar_trayectoria(res$P)
# Todo esto está en el último capítulo de la tesis, sirve bastante

# Todo lo anterior está excelente, nos va a servir mucho
N = 3
gens = 100
p0 <- c(1/3, 1/3, 1/3) # esto es como 1/3(1,1,1)
pruebaN3 <- WF_tres_alelos(N, p0, gens)

encontrarI <- function(trayectorias){
  for(i in 1:nrow(trayectorias)){
    fila <- trayectorias[i, ]
    
    if(any(fila == 0)){
      Ii <- which(fila == 0)
      taotilde <- i
      return(data.frame(
        I = Ii,
        ttau = taotilde
      ))  # la PRIMERA población que se vuelve 0
    }
  }
  
  return(data.frame(
    I = NA,
    ttau = NA
  ))  # ningún alelo se extinguió
}
gens = 100
tray1 <- WF_tres_alelos(N, p0, gens)
tray1
encontrarI(tray1)
a <- encontrarI(tray1)$I
b <- encontrarI(WF_tres_alelos(N, p0, gens))$I
mean(c(a,b))
graficar_trayectoria(tray1)
# Lo logramos!!!! R está bien confuso después de tanto usar python

# Simulación de I
trays <- list()
n_tray <- 1000
for(i in 1:n_tray){
  trays[[i]] <- WF_tres_alelos(N, p0, gens)
}

Is <- numeric(n_tray)

for (i in 1:n_tray) {
  Is[i] <- encontrarI(trays[[i]])$I
}

# Estimamos la función de masa de I

contador1 <- 0
contador2 <- 0
contador3 <- 0
for (i in Is){
  if(i==1){contador1 = contador1 + 1}
  if(i==2){contador2 = contador2 + 1}
  if(i==3){contador3 = contador3 + 1}
}
proba1 <- contador1/n_tray
proba2 <- contador2/n_tray
proba3 <- contador3/n_tray

masa <- c(proba1, proba2, proba3)


informacion <- data.frame(
  alelo = id,
  prob_extincion = masa
)

barplot(informacion$prob_extincion,
        names.arg = informacion$alelo,
        col = c("red", "blue", "green"),
        ylab = "Probabilidad de extinción",
        main = "Alelos que se extinguen primero")


#hist(unlist(Is), 
#     breaks = 3,
#     col = "lightblue",
#     main = "Distribución del alelo que se extingue primero")

# Ahora todo esto será una función que devuelve la masa!

masa_wf3 <- function(N, p0, gens, n_tray){
  trays <- vector("list", n_tray)
  
  for(i in 1:n_tray){
    trays[[i]] <- WF_tres_alelos(N, p0, gens)
  }
  
  Is <- numeric(n_tray)
  Taus <- numeric(n_tray)
  
  for (i in 1:n_tray) {
    res <- encontrarI(trays[[i]])
    Is[i] <- res$I
    Taus[i] <- res$ttau
  }
  
  # Quitar NA
  Is <- Is[!is.na(Is)]
  Taus <- Taus[!is.na(Taus)]
  
  # Tablas de frecuencia como densidades
  masaI <- prop.table(table(Is))
  masaTau <- prop.table(table(Taus))
  
  # Convertir a data frame
  df_I <- as.data.frame(masaI)
  colnames(df_I) <- c("valor", "densidad")
  
  df_tau <- as.data.frame(masaTau)
  colnames(df_tau) <- c("valor", "densidad")
  df_I$valor <- as.numeric(as.character(df_I$valor))
  df_tau$valor <- as.numeric(as.character(df_tau$valor))
  
  return(list(
    masa_I = df_I,
    masa_tau = df_tau
  ))
}

graficar_masa <- function(f, titulo){
  ggplot(f, aes(x = valor, y = densidad)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    scale_y_continuous(limits = c(0, max(f$densidad))) +
    labs(
      x = "Valor",
      y = "Probabilidad",
      title = titulo) +
    theme_minimal()
}

# Probemos las funciones
gens = 400
taut1 <- masa_wf3(100, p0, gens, n_tray)$masa_tau
taut1
library(glue)

graficar_masa(
  taut1,
  glue("tautilde con p0 = 1/3(1,1,1), {gens} gens y {n_tray} trayectorias")
)
lambda_hat <- 1 / mean(taut1$densidad)
1 / lambda_hat
ks.test(taut1, "pexp", rate = lambda_hat)

library(nortest)
# Prueba KS
ks.test(taut1, "pexp", rate = lambda_hat)
install.packages("kSamples")
library(kSamples)

ad.test(taut1, null = "exponential")

# Ahora voy a hacer una prueba de bondad de ajuste para ver si
# la I se distribuye igual que p_0

bondadI <- function(N, p0, proba_bondad, gens, n_tray){
    trays <- list()
    for(i in 1:n_tray){
      trays[[i]] <- WF_tres_alelos(N, p0, gens)
    }
    
    Is <- numeric(n_tray)
    
    for (i in 1:n_tray) {
      Is[i] <- encontrarI(trays[[i]])$I
    }
    Is <- Is[!is.na(Is)]   # quitar los NA
    
    tabla <- table(factor(Is, levels = 1:3))
    p <- chisq.test(tabla, p = proba_bondad)$p.value
    
    return(p)
  }

p1 <- c(0.9, 0.05, 0.05)
p0 <- c(1/3, 1/3, 1/3)
genera <- 10000
N <- 250
n_tray <- 500 # el tamaño de la muestra de I
m <- masa_wf3(N, p1, genera, n_tray)$masa_I # esto necesita de muchas trayectorias
graficar_masaI(m)
bondadI(N, p1, p0, genera, n_tray)  
bondadI(N, p1, p1, genera, n_tray)
p2 <- c(0.1, 0.45, 0.45)
bondadI(N, p1, p2, genera, n_tray)
p3 <- c(0.05, 0.475, 0.475)
bondadI(N, p1, p3, genera, n_tray)
p4 <- c(0.01, 0.495, 0.495) # Esta es la mejor aproximación!!! mucho mejor
bondadI(N, p1, p4, genera, n_tray)
p5 <- c(0.001, 0.4995, 0.4995)
bondadI(N, p1, p5, genera, n_tray) 


# Simulaciones para el caso N=3
N3 <- 3
gen3 <- 300
# Ojo porque esto no puede ser de la misma trayectoria
m3 <- masa_wf3(3, p0, gen3, n_tray)$masa_I# esto necesita de muchas trayectorias
graficar_masaI(m3)
tray3 <- WF_tres_alelos(N, p0, gen3)
graficar_trayectoria(tray3) # acá solo se necesita una trayectoria



