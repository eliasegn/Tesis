################################################################################
# Simulaciones de WF de tres alelos para N=7
# Autor: Elías González Nieto
# Afil : Facultad de Ciencias - UNAM
# Tesis: Distribuciones Tipo Fase para un Modelo de Coalescencia de Tres Alelos
################################################################################

# Librerías
library(igraph)
install.packages("ape")
library(ape)
library(PhaseTypeR)
library(ggplot2)
library(plotly)
library(akima)
library(scales)
set.seed(123) #(Esta semilla se cambia al aplicar Monte Carlo)

#######################################################################
# Construcción de WF de 3 alelos como cadena de Markov

matriz_transicion <- function(N) {
  # N es el tamaño de la poblacion 
  malla <- expand.grid(rep(list(0:N), 3))
  estados_totales <- as.matrix(malla[rowSums(malla) == N, , drop = FALSE])
  cardinal_S <- nrow(estados_totales)
  
  estados_absorbentes <- which(apply(estados_totales, 1, function(x) any(x == N)))
  estados_transientes  <- setdiff(seq_len(cardinal_S), estados_absorbentes)
  
  nuevo_orden <- c(estados_transientes, estados_absorbentes)
  orden_estados <- estados_totales[nuevo_orden, , drop = FALSE]
  nombres <- apply(orden_estados, 1, paste, collapse = ",")
  
  P <- matrix(0, nrow = cardinal_S, ncol = cardinal_S,
              dimnames = list(nombres, nombres))
  
  for (i in seq_len(cardinal_S)) {
    p_vec <- as.numeric(orden_estados[i, ]) / N
    probs <- apply(orden_estados, 1, function(y) dmultinom(as.integer(y), prob = p_vec))
    probs[probs < 0] <- 0
    probs <- probs / sum(probs)
    P[i, ] <- probs
  }
  
  list(
    P = P,
    estados = orden_estados,
    estados_transientes = seq_len(length(estados_transientes)),
    estados_absorbentes = (length(estados_transientes) + 1):cardinal_S
  )
}

# Descomposicion para tau (tiempo de fijacion total)
# la descomposición de tipo fase
separacion_matrices <- function(f) {
  P <- f$P
  num_estados <- nrow(P)
  num_trans <- length(f$estados_transientes)
  
  T <- P[1:num_trans, 1:num_trans]
  t <- matrix(rowSums(P[1:num_trans, (num_trans + 1):num_estados]), nrow = num_trans)
  list(T = T, t = t)
}

# Descomposicion para tau_tilde (tiempo de la primera extincion)
separacion_tau_tilde <- function(f) {
  P <- f$P
  estados <- f$estados
  
  absorbentes <- which(apply(estados, 1, function(x) any(x == 0)))
  transientes <- setdiff(seq_len(nrow(P)), absorbentes)
  
  P2 <- P[c(transientes, absorbentes), c(transientes, absorbentes)]
  k <- length(transientes)
  
  list(T = P2[1:k, 1:k],
       t = matrix(rowSums(P2[1:k, (k + 1):nrow(P2)]), nrow = k),
       estados_trans = estados[transientes, ])
}

# Descomposicion para tau_K (tiempo de competencai de dos alelos)
separacion_tau_k <- function(f, N) {
  P <- f$P
  estados <- f$estados
  
  dos_alelos <- which(apply(estados, 1, function(x) sum(x == 0) == 1))
  fijacion   <- which(apply(estados, 1, function(x) any(x == N)))
  transientes <- setdiff(dos_alelos, fijacion)
  
  P2 <- P[c(transientes, fijacion), c(transientes, fijacion)]
  k <- length(transientes)
  
  list(T = P2[1:k, 1:k],
       t = matrix(rowSums(P2[1:k, (k + 1):nrow(P2)]), nrow = k),
       estados_trans = estados[transientes, ])
}

# Indicadora del conjunto interior E, sobre los estados transientes de T (tau)
indicadora_E <- function(estados_trans, N) {
  as.integer(apply(estados_trans, 1, function(v) all(v >= 1 & v <= N - 1)))
}

#######################################################################
# Superifices de esperanza y varianza condicional

# Calcula E_i[Y] o Var_i[Y] para cada estado transiente i
momentos_condicionales <- function(T, R = NULL, tipo = c("media", "varianza")) {
  # R es por si tiene recompensa, para que el código pueda ser replicado
  tipo <- match.arg(tipo)
  k <- nrow(T)
  out <- numeric(k)
  for (i in 1:k) {
    pi_i <- rep(0, k); pi_i[i] <- 1
    dph <- DPH(T, pi_i)
    obj <- if (is.null(R)) dph else reward_phase_type(dph, R)
    out[i] <- if (tipo == "media") mean(obj) else var(obj)
  }
  out
}

# Curvas de nivel sobre el simplejo
graficar_curvas_nivel <- function(estados_trans, valores, N, titulo,
                                  etiquetas_borde = FALSE) {
  interp_out <- akima::interp(
    x = estados_trans[, 1], y = estados_trans[, 2], z = valores,
    xo = seq(0, N, length = 300), yo = seq(0, N, length = 300), linear = TRUE
  )
  image(interp_out$x, interp_out$y, interp_out$z,
        col = hcl.colors(60, "Plasma"),
        xlab = "Alelo 1", ylab = "Alelo 2", axes = FALSE, main = titulo)
  polygon(x = c(0, N, 0), y = c(0, 0, N), border = "#8B1C62", lwd = 4)
  axis(1); axis(2)
  contour(interp_out$x, interp_out$y, interp_out$z, add = TRUE, nlevels = 15)
  abline(h = 0:N, v = 0:N, col = "black")
  if (etiquetas_borde) {
    text(N / 2, -0.25, "Extincion de Alelo 2")
    text(-0.25, N / 2, "Extincion de Alelo 1", srt = 90)
    text(N / 2, N - N / 2 + 0.25, "Extincion de Alelo 3")
  }
  invisible(interp_out)
}

# Superficie 3D sobre el simplejo
graficar_superficie_3d <- function(estados_trans, valores, N, titulo_z) {
  interp_out <- akima::interp(x = estados_trans[, 1], y = estados_trans[, 2],
                              z = valores, linear = TRUE)
  mask <- outer(interp_out$x, interp_out$y, function(a, b) a + b <= N)
  z_base <- matrix(NA, length(interp_out$x), length(interp_out$y))
  z_base[mask] <- 0
  
  plot_ly() %>%
    add_surface(x = interp_out$x, y = interp_out$y, z = z_base,
                showscale = FALSE, opacity = 0.4) %>%
    add_surface(x = interp_out$x, y = interp_out$y, z = interp_out$z) %>%
    layout(scene = list(xaxis = list(title = "Alelo 1"),
                        yaxis = list(title = "Alelo 2"),
                        zaxis = list(title = titulo_z)))
}

#######################################################################
# Construcción de cadena y figuras para la tesis

# Construcción de la cadena
N <- 7
f <- matriz_transicion(N)

sep_tau       <- separacion_matrices(f)
sep_tau_tilde <- separacion_tau_tilde(f)
sep_tau_k     <- separacion_tau_k(f, N)

estados_trans_tau       <- f$estados[f$estados_transientes, ]
estados_trans_tau_tilde <- sep_tau_tilde$estados_trans
estados_trans_tau_k     <- sep_tau_k$estados_trans

# Esperanzas
Etau       <- momentos_condicionales(sep_tau$T,       tipo = "media")
Etau_tilde <- momentos_condicionales(sep_tau_tilde$T,  tipo = "media")
Etau_k     <- momentos_condicionales(sep_tau_k$T,      tipo = "media")

graficar_superficie_3d(estados_trans_tau, Etau, N, "E_x[tau]")
graficar_curvas_nivel(estados_trans_tau, Etau, N, "E[tau]", etiquetas_borde = TRUE)

graficar_curvas_nivel(estados_trans_tau_tilde, Etau_tilde, N, "E[tau_tilde]")

graficar_curvas_nivel(estados_trans_tau_k, Etau_k, N, "E[tau_K]")
graficar_superficie_3d(estados_trans_tau_k, Etau_k, N, "E_x[tau_K]")

# Varianzas
Vartau_tilde <- momentos_condicionales(sep_tau_tilde$T, tipo = "varianza")
Vartau_k     <- momentos_condicionales(sep_tau_k$T,     tipo = "varianza")
Vartau       <- momentos_condicionales(sep_tau$T,       tipo = "varianza")

par(mfrow = c(2, 2))
graficar_curvas_nivel(estados_trans_tau_tilde, Vartau_tilde, N, "Var(tau_tilde)")
graficar_curvas_nivel(estados_trans_tau_k,     Vartau_k,     N, "Var(tau_K)")
graficar_curvas_nivel(estados_trans_tau,       Vartau,       N, "Var(tau)")
par(mfrow = c(1, 1))

#######################################################################
# Covarianzas teórica y estimación Monte Carlo

# Estimación Monte Carlo
estimar_covarianza_mc <- function(f, N, num_simulaciones = 5000, estado_inicial = NULL) {
  P <- f$P
  estados <- f$estados
  
  idx_interior    <- which(apply(estados, 1, function(x) all(x > 0)))
  idx_frontera    <- which(apply(estados, 1, function(x) sum(x == 0) == 1))
  idx_absorbentes <- which(apply(estados, 1, function(x) sum(x == 0) == 2))
  
  if (is.null(estado_inicial)) {
    distancias <- apply(estados[idx_interior, ], 1, function(x) sum((x - N / 3)^2))
    estado_inicial <- idx_interior[which.min(distancias)]
  }
  
  n_estados <- nrow(P)
  sim_tau_tilde <- numeric(num_simulaciones)
  sim_tau_k     <- numeric(num_simulaciones)
  
  for (i in 1:num_simulaciones) {
    estado_actual <- estado_inicial
    t_tilde <- 0; t_k <- 0
    
    while (!(estado_actual %in% idx_absorbentes)) {
      if (estado_actual %in% idx_interior) {
        t_tilde <- t_tilde + 1
      } else if (estado_actual %in% idx_frontera) {
        t_k <- t_k + 1
      }
      estado_actual <- sample(1:n_estados, size = 1, prob = P[estado_actual, ])
    }
    sim_tau_tilde[i] <- t_tilde
    sim_tau_k[i]     <- t_k
  }
  
  list(
    estado_inicial = estados[estado_inicial, ],
    covarianza = cov(sim_tau_tilde, sim_tau_k),
    correlacion = cor(sim_tau_tilde, sim_tau_k),
    datos = data.frame(Tau_Tilde = sim_tau_tilde, Tau_K = sim_tau_k)
  )
}

recomp_tilde <- indicadora_E(estados_trans_tau, N)
recomp_K     <- 1L - recomp_tilde
R_mat        <- cbind(recomp_tilde, recomp_K)

idx_interior_tau <- which(recomp_tilde == 1)
distancias_i0 <- apply(estados_trans_tau[idx_interior_tau, ], 1, function(x) sum((x - N / 3)^2))
i0  <- idx_interior_tau[which.min(distancias_i0)]
pi0 <- rep(0, nrow(sep_tau$T)); pi0[i0] <- 1

mdph <- MDPH(sep_tau$T, init_probs = pi0, reward_mat = R_mat)

cov_exacta_mdph <- var(mdph)[1, 2]
muestras_mdph   <- rMDPH(10000, mdph)

resultado_teorico <- cov_exacta_mdph
set.seed(22)
tamanos_simulacion <- c(100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 750000, 1000000)
covarianzas_empiricas <- sapply(tamanos_simulacion, function(n) {
  estimar_covarianza_mc(f, N, num_simulaciones = n)$covarianza
})
 
datos_plot <- data.frame(Iteraciones = tamanos_simulacion, Covarianza = covarianzas_empiricas)
datos_plot
grafico_convergencia <- ggplot(datos_plot, aes(x = Iteraciones, y = Covarianza)) +
  geom_hline(yintercept = resultado_teorico, color = "firebrick",
             linetype = "dashed", linewidth = 1) +
  geom_line(color = "dodgerblue3", linewidth = 1) +
  geom_point(color = "dodgerblue3", size = 3) +
  scale_x_log10(labels = label_comma()) +
  labs(title = "Convergencia de Monte Carlo al Valor Teorico Exacto",
       subtitle = sprintf("Covarianza Teorica Exacta: %.4f", resultado_teorico),
       x = "Iteraciones (Escala Logaritmica)",
       y = "Covarianza") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(color = "firebrick", face = "italic", hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 15), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), face = "bold"),
        panel.grid.minor = element_blank())

print(grafico_convergencia)

# Resultados de Covarianza
cat("Cov. exacta:", cov_exacta_mdph, "\n")
cat("Cov. Monte Carlo:", estimar_covarianza_mc(f, N, 500000, estado_inicial = i0)$covarianza, "\n")

#######################################################################
# Interavalos de Confianza

resultados_mc_masivo <- estimar_covarianza_mc(f, N, num_simulaciones = 500000, estado_inicial = i0)
datos_sim <- resultados_mc_masivo$datos
datos_sim$Tau_Total <- datos_sim$Tau_Tilde + datos_sim$Tau_K
n_sim <- nrow(datos_sim)

obtener_ic <- function(variable, conf_level = 0.95) {
  media <- mean(variable)
  error_est <- sd(variable) / sqrt(n_sim)
  z <- qnorm(1 - (1 - conf_level) / 2)
  c(Media = media, L_Inf = media - z * error_est, L_Sup = media + z * error_est)
}

ic_tilde <- obtener_ic(datos_sim$Tau_Tilde)
ic_k     <- obtener_ic(datos_sim$Tau_K)
ic_total <- obtener_ic(datos_sim$Tau_Total)

cat("Tabla\n")
cat(sprintf("tilde{tau}: Media = %.3f, I.C. 95%% = [%.3f, %.3f]\n", ic_tilde[1], ic_tilde[2], ic_tilde[3]))
cat(sprintf("tau_K:      Media = %.3f, I.C. 95%% = [%.3f, %.3f]\n", ic_k[1], ic_k[2], ic_k[3]))
cat(sprintf("tau:        Media = %.3f, I.C. 95%% = [%.3f, %.3f]\n", ic_total[1], ic_total[2], ic_total[3]))