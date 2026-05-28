################################################################################
# Grafos de modelos de coalescencia
# Autor: Elías González Nieto
# Afil : Facultad de Ciencias - UNAM
# Tesis: Distribuciones Tipo Fase para un Modelo de Coalescencia de Tres Alelos
################################################################################

# Librerías
library(igraph)
install.packages("ape")
library(ape)

#######################################################################
# Simulación del Modelo de Wright-Fisher

# Función para simular un proceso de Wright-Fisher bialélico
WF_grafo_completo <- function(N, p0 = 0.5) {
  # Parámetros iniciales
  gen <- rbinom(N, 1, p0)
  history <- list(gen)
  edges <- c()
  
  # Simulación hasta la fijación
  while (!(all(gen == 0) || all(gen == 1))) {
    next_gen <- numeric(N)
    t_actual <- length(history) - 1
    t_siguiente <- length(history)
    
    for (i in 1:N) {
      # Simulación de la siguiente generación
      parent <- sample(1:N, 1)
      next_gen[i] <- gen[parent]
      
      # Nombres para los nodos
      from <- paste0("g", t_actual, "_", parent)
      to   <- paste0("g", t_siguiente, "_", i)
      edges <- c(edges, from, to)
    }
    gen <- next_gen
    history[[length(history) + 1]] <- gen
    
    # Límite por si tarda mucho en fijar
    if(length(history) > 50) break
  }
  
  # Crear el grafo con la librería igraph
  g <- graph_from_edgelist(matrix(edges, ncol = 2, byrow = TRUE), directed = TRUE)
  
  # Consideraciones de nodos
  all_node_names <- unlist(lapply(0:(length(history)-1), function(t) paste0("g", t, "_", 1:N)))
  missing_nodes <- setdiff(all_node_names, V(g)$name)
  if(length(missing_nodes) > 0) g <- add_vertices(g, nv = length(missing_nodes), name = missing_nodes)
  
  # Sincronización de atributos para igraph
  node_data <- data.frame(name = V(g)$name)
  node_data$gen <- as.numeric(gsub("g(\\d+)_\\d+", "\\1", node_data$name))
  node_data$ind <- as.numeric(gsub("g\\d+_(\\d+)", "\\1", node_data$name))
  
  # Coordenadas 
  coords <- matrix(0, nrow = vcount(g), ncol = 2)
  coords[, 1] <- node_data$ind
  coords[, 2] <- -node_data$gen
  
  # Colores y edición
  vertex_colors <- sapply(1:vcount(g), function(i) {
    t <- node_data$gen[i] + 1 # +1 porque R empieza en 1
    id <- node_data$ind[i]
    if(history[[t]][id] == 1) "darkolivegreen" else "cadetblue"
  })
  
  # Dibujo
  plot(g,
       layout = coords,
       vertex.color = vertex_colors,
       vertex.size = 15,
       vertex.label = NA,
       edge.arrow.size = 0.2,
       edge.color = "gray10",
       edge.curved = 0,
       main = paste("Modelo de Wright-Fisher"))
  
  return(history)
}

set.seed(111) # Para reproducibilidad
historial <- WF_grafo_completo(6) # Simulación (figura 1.2 de la tesis)

# Función para hacer un grafo de la finámica con k alelos
set.seed(123)
WF_grafo_k_alelos <- function(N, k = 3) {
  
  # Definimos una paleta de colores preestablecida para k alelos
  colores_alelos <- rainbow(k)
  
  # Usamos como distribución inicial la uniforme
  gen <- sample(1:k, N, replace = TRUE)
  history <- list(gen)
  edges <- c()
  
  # Simulamos hasta fijación
  while (length(unique(gen)) > 1) {
    next_gen <- numeric(N)
    t_prev <- length(history) - 1
    t_curr <- length(history)
    
    for (i in 1:N) {
      # Simulación de WF
      parent <- sample(1:N, 1)
      next_gen[i] <- gen[parent]
      
      # Creación del grafo
      from <- paste0("g", t_prev, "_", parent)
      to   <- paste0("g", t_curr, "_", i)
      edges <- c(edges, from, to)
    }
    
    gen <- next_gen
    history[[length(history) + 1]] <- gen
    
    # Límite para la convergencia
    if(length(history) > 150) break
  }
  
  # Creación y Edición del grafo
  g <- graph_from_edgelist(matrix(edges, ncol = 2, byrow = TRUE), directed = TRUE)
  all_nodes <- unlist(lapply(0:(length(history)-1), function(t) paste0("g", t, "_", 1:N)))
  missing <- setdiff(all_nodes, V(g)$name)
  if(length(missing) > 0) g <- add_vertices(g, nv = length(missing), name = missing)
  nombres <- V(g)$name
  gen_num <- as.numeric(gsub("g(\\d+)_.*", "\\1", nombres))
  ind_num <- as.numeric(gsub(".*_(\\d+)", "\\1", nombres))
  coords <- matrix(0, nrow = vcount(g), ncol = 2)
  coords[, 1] <- ind_num
  coords[, 2] <- -gen_num
  vertex_colors <- sapply(1:vcount(g), function(i) {
    t <- gen_num[i] + 1
    id <- ind_num[i]
    alelo_val <- history[[t]][id]
    colores_alelos[alelo_val]
  })
  
  # Display del grafo
  plot(g,
       layout = coords,
       vertex.color = vertex_colors,
       vertex.size = 10,
       vertex.label = NA,
       edge.arrow.size = 0.2,
       edge.color = "gray10",
       edge.curved = 0, # Líneas rectas
       main = paste("Wright-Fisher con k =", k, "alelos y N =", N))
  
  return(history)
}

# Ejemplo
historial <- WF_grafo_k_alelos(N = 7, k = 4)

#######################################################################
# Modelo de Moran

Moran_grafo_continuo <- function(N, p0 = 0.5) {
  
  # Inicializamos con la probabilidad dada
  gen <- rbinom(N, 1, p0)
  history <- list(gen)
  times <- c(0) # Registro del tiempo por ser a tiempo continuo
  edges <- c()
  
  # Simulación hasta fijación
  while (!(all(gen == 0) || all(gen == 1))) {
    
    # Simulación con tasas exponenciales
    dt <- rexp(1, rate = N)
    times <- c(times, tail(times, 1) + dt)
    t_actual_idx <- length(history)
    next_gen <- gen
    
    # Muerte aleatoria
    muere <- sample(1:N, 1)
    # Nacimiento por elección
    nace_de <- sample(1:N, 1)
    
    # El individuo que nace reemplaza al que muere
    next_gen[muere] <- gen[nace_de]
    
    # Creación y Edición del grafo
    for (i in 1:N) {
      from_node <- paste0("s", t_actual_idx - 1, "_", i)
      to_node   <- paste0("s", t_actual_idx, "_", i)
      
      # Si i es el que murió, su ancestro es el padre (nace_de)
      if (i == muere) {
        from_node <- paste0("s", t_actual_idx - 1, "_", nace_de)
      }
      
      edges <- c(edges, from_node, to_node)
    }
    
    gen <- next_gen
    history[[length(history) + 1]] <- gen
    
    # Límite para convergencia
    if(length(history) > 200) break 
  }
  
  # Display del grafo
  g <- graph_from_edgelist(matrix(edges, ncol = 2, byrow = TRUE), directed = TRUE)
  
  # Leyendas y ediciones
  nombres <- V(g)$name
  step_num <- as.numeric(gsub("s(\\d+)_.*", "\\1", nombres))
  ind_num  <- as.numeric(gsub(".*_(\\d+)", "\\1", nombres))
  
  coords <- matrix(0, nrow = vcount(g), ncol = 2)
  coords[, 1] <- ind_num
  # Usamos el tiempo acumulado para el eje Y en lugar de pasos discretos
  coords[, 2] <- -sapply(step_num, function(s) times[s + 1])
  
  vertex_colors <- sapply(1:vcount(g), function(i) {
    s <- step_num[i] + 1
    id <- ind_num[i]
    if(history[[s]][id] == 1) "darkolivegreen" else "cadetblue"
  })
  
  # 5. Plot
  plot(g,
       layout = coords,
       vertex.color = vertex_colors,
       vertex.size = 8,
       vertex.label = NA,
       edge.arrow.size = 0.15,
       edge.color = "gray10",
       edge.curved = 0,
       main = paste("Modelo de Moran"))
  
  return(list(history=history, times=times))
}


# Ejemplo
set.seed(61) # Semilla que funciona para tener una buena representación
resultado <- Moran_grafo_continuo(N = 5) # N pequeño para apreciar los saltos de tiempo (figura 1.3 en la tesis)

#######################################################################
# Árboles Coaelescentes

par(mfrow = c(1, 1))

# Función para generar un árbol coalescente donde la muestra es de población n 
arbol_coalescente <- function(n){
  # Generar el árbol usando la librería ape
  tree <- rcoal(n)
  
  # Etiquetas para los nodos
  tree$tip.label <- 1:n
  
  # Gráfico del árbol
  plot(tree, 
       main = paste("Árbol Coalescente con n =", n),
       edge.width = 2, 
       direction = "downwards",
       show.tip.label = FALSE) # Ocultamos el texto original "t..."
  
  # Edición de nodos y eje del tiempo
  tiplabels(pch = 21, bg = "black", col = "white", cex = 1)
  axisPhylo(side = 2)
}
# Árboles coalescentes con poblaciones de 50, 100, 500 y 1000 (figura 1.5 en la tesis)
set.seed(42)
arbol_coalescente(15)
poblaciones <- c( 50, 100,500, 1000)
for(p in poblaciones){arbol_coalescente(p)}
arbol_coalescente(500)

