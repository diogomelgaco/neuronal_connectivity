######################################################################################################################## o
####  Adaptação do Estimador de conectividade neuronal pelo método dos pares de neurônios (Figueiredo 2025)
######################################################################################################################## o
####  Autor: Diogo Maia de Figueiredo
####  Orientadora: Dra. Aline Duarte
####  Programa de Pós Graduação do IME-USP

######################################################################################################################## o
######################################################################################################################## o
#####  1: Contagens todas as janelas  ----

SA_star <- function(X, i, aux = NULL) {
  max_aux <- ncol(X) - 1
  if (is.null(aux)) aux <- max_aux
  ks <- seq(1, min(aux, max_aux))
  sum(X[i, ks] == 1)
}

SB_star <- function(X, i, aux = NULL) {
  max_aux <- ncol(X) - 1
  if (is.null(aux)) aux <- max_aux
  ks <- seq(1, min(aux, max_aux))
  sum(X[i, ks] == 1 & X[i, ks + 1] == 1)
}

SC_star <- function(X, i, j, aux = NULL) {
  max_aux <- ncol(X) - 2
  if (is.null(aux)) aux <- max_aux
  ks <- seq(1, min(aux, max_aux))
  sum(X[i, ks] == 1 & X[j, ks + 1] == 1)
}

SD_star <- function(X, i, j, aux = NULL) {
  max_aux <- ncol(X) - 2
  if (is.null(aux)) aux <- max_aux
  ks <- seq(1, min(aux, max_aux))
  sum(X[i, ks] == 1 & X[j, ks + 1] == 1 & X[i, ks + 2] == 1)
}

######################################################################################################################## o
######################################################################################################################## o
#####  2: Estimador Pares Extendido  ----

estimador_redeneural_pares_vetorizado_star_fast <- function(X, u, Delta = 10.33,csi1=0.005, csi2=0.10,  verbose = FALSE) {
  n <- ncol(X)
  m <- nrow(X)
  
  # Inicializa matrizes de saída
  theta <- matrix(0, m, m)
  diag(theta) <- 0
  R_list <- numeric(m)
  G_matrix <- matrix(0, m, m)
  
  # Parâmetros da sigmoide
  p_min <- phi(-1e6, u)
  p_max <- phi(1e6, u)
  delta <- phi(1, u) - phi(0, u)
  a <- p_min / p_max
  b <- delta / p_max
  
  # Cotas amostrais
  n_aux <- floor(n / 3)
  n_lim <- 2*ceiling(p_min * Delta * n_aux)
  n_m <- ceiling((19 / 20) * (p_min^2 * Delta^2) * (1 - ((b / 10) * sqrt(p_min * Delta))) * n_aux)
  
  if (verbose) {
    cat("Parâmetros:\n")
    cat("p_min =", p_min, "p_max =", p_max, "delta =", delta, "\n")
    cat("a =", a, "b =", b, "\n")
    cat("n_aux =", n_aux, "n_lim =", n_lim, "n_m =", n_m, "\n")
  }
  
  # Pré-criação de janelas para SA/SB/SC/SD
  ks_SB <- 1:(n - 1)
  ks_SD <- 1:(n - 2)
  
  for (i in 1:m) {
    # Pré-calcular SA e SB para todas as janelas possíveis
    sa_vec <- cumsum(X[i, ks_SB] == 1)
    sb_vec <- cumsum((X[i, ks_SB] == 1) & (X[i, ks_SB + 1] == 1))
    
    # Encontrar o menor K tal que SA ≥ n_m
    K <- which(sa_vec >= n_m)[1]
    if (is.na(K)) K <- length(sa_vec)
    
    # Calcular R
    if (n_lim < K) {
      SA_val <- sa_vec[n_lim]
      SB_val <- sb_vec[n_lim]
      R <- SB_val / SA_val
    } else {
      SA_val <- sa_vec[K]
      SB_val <- sb_vec[K]
      R <- SB_val / n_m
    }
    
    R_list[i] <- R
    
    if (verbose) {
      cat("i =", i, "| K =", K, "| SA =", SA_val, "| SB =", SB_val, "| R =", round(R, 4), "\n")
    }
    
    for (j in setdiff(1:m, i)) {
      # Pré-calcular SC e SD vetorizadamente
      sc_vec <- cumsum((X[i, ks_SD] == 1) & (X[j, ks_SD + 1] == 1))
      sd_vec <- cumsum((X[i, ks_SD] == 1) & (X[j, ks_SD + 1] == 1) & (X[i, ks_SD + 2] == 1))
      
      # Encontrar o menor H tal que SC ≥ n_m
      H <- which(sc_vec >= n_m)[1]
      if (is.na(H)) H <- length(sc_vec)
      
      SC_val <- sc_vec[H]
      SD_val <- sd_vec[H]
      
      # Calcular G
      G <- if (SC_val == 0) 0 else SD_val / SC_val
      
      G_matrix[j, i] <- G
      theta[j, i] <- G - R
      
      if (verbose) {
        cat("   j =", j, "| H =", H, "| SC =", SC_val, "| SD =", SD_val, "| G =", round(G, 4), "| theta =", round(theta[j, i], 4), "\n")
      }
    }
  }
  
  theta <- round(theta, 4)
  
  V_est <- gerar_W_est(theta,csi1, csi2)
  
  return(list(R_list = R_list, G_matrix = G_matrix, theta = theta,V_est=V_est))
}



######################################################################################################################## o
######################################################################################################################## o
#####  3: Gráfico comparativo ----

plot_comparacao_V_V_est_excit_inib_star <- function(V_list, V_est_list, n_list) {
  library(ggplot2)
  library(gridExtra)
  library(grid)
  
  # 🛠️ Garantir que todos os inputs sejam listas
  if (!is.list(V_list)) V_list <- list(V_list)
  if (!is.list(V_est_list)) V_est_list <- list(V_est_list)
  if (!is.list(n_list)) n_list <- list(n_list)
  
  n_list <- lapply(n_list, as.numeric)
  plot_list <- list()
  
  # 🎨 Paleta de cores diferenciando excitatórios e inibitórios
  colors <- c(
    "Conexão excitatória existente (acerto)" = "steelblue",
    "Conexão inibitória existente (acerto)"  = "mediumpurple3",
    "Sem conexão (acerto)"                   = "white",
    "Conexão não existente (falso positivo)"= "yellow",
    "Conexão não identificada (falso negativo)"= "sandybrown",
    "Conexão existente invertida"                = "lightgrey"
  )
  
  titulos <- c(
    expression("Método por pares" ~ hat(V)[n]^{"''"}),
    expression("Método por pares estendido" ~ hat(V)[n]^{"''*"})
  )
  
  for (i in seq_along(V_list)) {
    n <- n_list[[i]]
    V <- V_list[[i]]
    V_est <- V_est_list[[i]]
    
    heatmap_data <- data.frame(
      Presynaptic = rep(1:nrow(V), each = ncol(V)),
      Postsynaptic = rep(1:ncol(V), times = nrow(V)),
      V = as.vector(t(V)),
      V_est = as.vector(t(V_est))
    )
    
    heatmap_data$Condition <- with(heatmap_data, ifelse(
      is.na(V_est), "Inconclusivo",
      ifelse(V_est == 1 & V == 1, "Conexão excitatória existente (acerto)",
             ifelse(V_est == -1 & V == -1, "Conexão inibitória existente (acerto)",
                    ifelse(V_est == 0 & V == 0, "Sem conexão (acerto)",
                           ifelse(V_est == 0 & V != 0, "Conexão não identificada (falso negativo)",
                                  ifelse(V_est != 0 & V == 0, "Conexão não existente (falso positivo)",
                                         ifelse(V_est != 0 & V != 0 & V_est != V, "Conexão existente invertida", "Outro")))))))
    )
    
    subtitle <- titulos[[i]]
    
    # 🔥 Criar heatmap individual
    p <- ggplot(heatmap_data, aes(x = Postsynaptic, y = Presynaptic, fill = Condition)) +
      geom_tile(color = "black") +
      scale_fill_manual(values = colors) +
      labs(
        x = "Neurônios pós-sinápticos", 
        y = "Neurônios pré-sinápticos", 
        title = subtitle
      ) +
      scale_x_continuous(breaks = 1:ncol(V), position = "top") +
      scale_y_reverse(breaks = 1:nrow(V)) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "none"
      )
    
    plot_list[[i]] <- p
  }
  
  # 🎨 Preparar legenda única
  legenda_data <- data.frame(
    Condition = factor(names(colors), levels = names(colors))
  )
  
  legenda_plot <- ggplot(legenda_data, aes(x = Condition, fill = Condition)) +
    geom_bar() +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(nrow = 2)) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key = element_rect(colour = "black", size = 0.5),
      legend.key.size = unit(0.7, "cm"),
      legend.spacing.x = unit(0, "cm"),
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0)
    )
  
  # 📦 Função auxiliar para extrair legenda
  g_legend <- function(a.gplot) {
    tmp <- ggplotGrob(a.gplot)
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    tmp$grobs[[leg]]
  }
  
  legenda_unica <- g_legend(legenda_plot)
  
  # 🔢 Ajuste automático no número de colunas
  n_graficos <- length(plot_list)
  ncol_plot <- ifelse(n_graficos == 1, 1, ifelse(n_graficos == 2, 2, 3))
  
  # 🖼️ Combinar gráficos + legenda
  grid.arrange(
    arrangeGrob(grobs = plot_list, ncol = ncol_plot),
    legenda_unica,
    heights = c(10, 1.5)
  )
}

######################################################################################################################## o
######################################################################################################################## o
#####  4: Aplicacao nos exemplos ----

#### Exemplo 3 ----


est_pares3_3star <-estimador_redeneural_pares_vetorizado_star_fast(mod3_3$X, mod3_3$u)

# Gráfico comparando V com V_est
V_list3_comp <- list(mod3_3$W,mod3_3$W)
V_list_est_par1 <- list(est_pares3_3$V_est,est_pares3_3star$V_est)
n_list3_comp <- list(ncol(mod3_3$X),ncol(mod3_3$X))

plot_comparacao_V_V_est_excit_inib_star(V_list3_comp,V_list_est_par1,n_list3_comp)


