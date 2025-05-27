######################################################################################################################## o
####  Estimador de conectividade neuronal pelo método dos pares de neurônios (Santis 2022)
######################################################################################################################## o
####  Autor: Diogo Maia de Figueiredo
####  Orientadora: Dra. Aline Duarte
####  Programa de Pós Graduação do IME-USP

######################################################################################################################## o
######################################################################################################################## o
#####  1: Funções para gerar eventos necessarios para o estimador  ----

######  1.1 Estimativa de u ----
media_u_sobre_mil_W <- function(m, n_amostras = 10000) {
  us <- numeric(n_amostras)
  
  for (i in 1:n_amostras) {
    W <- gerador_matriz_W(m)
    u <- max(1, mean(colSums(abs(W))))
    us[i] <- u
  }
  
  return(mean(us))
}

u_real=media_u_sobre_mil_W(5)

######  1.2 Contagens A,B,C,D ----

#SA conta quantos disparos do neurônio i aconteceram em janelas fixas a cada 2 instantes: X_t(i)=1.
SA <- function(X, i, aux = floor(ncol(X) / 2)) {                      
  max_aux <- floor(ncol(X) / 2)
  aux <- if (aux > max_aux) max_aux else aux
  ks <- seq(1, aux)
  sum(X[i, 2 * ks - 1] == 1)
}

#SB conta quantos disparos em sequência do neurônio i aconteceram em janelas fixas a cada 2 instantes: X_t(i)=1, X_(t+1)=1.
SB <- function(X, i, aux = floor(ncol(X) / 2)) {
  max_aux <- floor(ncol(X) / 2)
  aux <- if (aux > max_aux) max_aux else aux
  ks <- seq(1, aux)
  sum(X[i, 2 * ks - 1] == 1 & X[i, 2 * ks] == 1)
}

#SC conta quantos disparos do neuronio i, em sequencia do j aconteceram em janelas fixas a cada 2 instantes: X_t(i)=1, X_(t+1)(j)=1.
SC <- function(X, i, j, aux = floor(ncol(X) / 3)) {
  max_aux <- floor(ncol(X) / 3)
  aux <- if (aux > max_aux) max_aux else aux
  ks <- seq(1, aux)
  sum(X[i, 3 * ks - 2] == 1 & X[j, 3 * ks - 1] == 1)
}


#SD conta quantos disparos do neuronio i, em sequencia do j aconteceram em janelas fixas a cada 2 instantes: X_t(i)=1, X_(t+1)(j)=1.
SD <- function(X, i, j, aux = floor(ncol(X) / 3)) {
  max_aux <- floor(ncol(X) / 3)
  aux <- if (aux > max_aux) max_aux else aux
  ks <- seq(1, aux)
  sum(X[i, 3 * ks - 2] == 1 & X[j, 3 * ks - 1] == 1 & X[i, 3 * ks] == 1)
}

######################################################################################################################## o
######################################################################################################################## o
#####  2: Construção do Estimador pelo método por pares  ----  


gerar_W_est <- function(matriztheta, par1, par2) {
  W_est <- matrix(0, nrow = nrow(matriztheta), ncol = ncol(matriztheta))
  W_est[matriztheta >= par2] <- 1
  W_est[matriztheta <= -par1] <- -1
  return(W_est)
}

estimador_redeneural_pares<- function(X, u, Delta = 10.33,csi1=0.005, csi2=0.10, verbose = FALSE) {
  n <- ncol(X)
  m <- nrow(X)
  
  d <- n-1        #maxima cardinalidade    
  
  # Inicializa matrizes de saída
  theta <- matrix(0, m, m)
  diag(theta) <- 0
  R_list <- numeric(m)
  G_matrix <- matrix(0, m, m)
  
  # Calcula parâmetros da sigmoide para definição de cotas
  p_min <- phi(-1e6, u)
  p_max <- phi(1e6, u)
  delta <- phi(1, u) - phi(0, u)
  a <- p_min / p_max
  b <- delta / p_max
  
  # Cotas amostrais
  n_aux <- floor(n / 3)
  n_lim <- ceiling(p_min * Delta * n_aux)
  n_m <- ceiling((19 / 20) * (p_min^2 * Delta^2) *
                   (1 - ((b / 10) * sqrt(p_min * Delta))) * n_aux)
  
  if (verbose) {
    cat("Parâmetros:\n")
    cat("p_min =", p_min, "p_max =", p_max, "delta =", delta, "\n")
    cat("a =", a, "b =", b, "\n")
    cat("n_aux =", n_aux, "n_lim =", n_lim, "n_m =", n_m, "\n")
  }
  
  ### 📌 Pré-processamento vetorizado das janelas de SC e SD
  aux_max_SC <- floor(n / 3)
  ks <- 1:aux_max_SC
  inds_t   <- 3 * ks - 2
  inds_t1  <- 3 * ks - 1
  inds_t2  <- 3 * ks
  
  X_t   <- X[, inds_t, drop = FALSE]    # X[i, t]
  X_t1  <- X[, inds_t1, drop = FALSE]   # X[j, t+1]
  X_t2  <- X[, inds_t2, drop = FALSE]   # X[i, t+2]
  
  ### 🚀 Loop por neurônio i (pré-sináptico)
  for (i in 1:m) {
    ### Encontrar K tal que SA ≥ n_m
    max_aux_K <- floor(n / 2)
    K <- NA
    for (aux in 3:max_aux_K) {
      SA_val <- sum(X[i, 2 * (1:aux) - 1])
      if (SA_val >= n_m) {
        K <- aux
        break
      }
    }
    if (is.na(K)) K <- max_aux_K
    
    ### Calcular R[i]
    if (n_lim < K) {
      ks_lim <- 1:n_lim
      inds1 <- 2 * ks_lim - 1
      inds2 <- 2 * ks_lim
      x1 <- X[i, inds1]
      x2 <- X[i, inds2]
      SA_val <- sum(x1)
      SB_val <- sum(x1 & x2)
      R <- SB_val / SA_val
    } else {
      ks_K <- 1:K
      inds1 <- 2 * ks_K - 1
      inds2 <- 2 * ks_K
      x1 <- X[i, inds1]
      x2 <- X[i, inds2]
      SB_val <- sum(x1 & x2)
      R <- SB_val / n_m
    }
    
    R_list[i] <- R
    if (verbose) {
      cat("i =", i, "| K =", K, "| SA =", SA_val, "| SB =", SB_val, "| R =", round(R, 4), "\n")
    }
    
    ### 🔁 Loop por neurônio j (pós-sináptico)
    for (j in setdiff(1:m, i)) {
      x1 <- X_t[i, 1:aux_max_SC]
      x2 <- X_t1[j, 1:aux_max_SC]
      x3 <- X_t2[i, 1:aux_max_SC]
      
      SC_vec <- x1 & x2
      SD_vec <- x1 & x2 & x3
      
      # Encontrar H tal que SC acumulado ≥ n_m
      H <- which(cumsum(SC_vec) >= n_m)[1]
      if (is.na(H)) H <- aux_max_SC
      
      SC <- sum(SC_vec[1:H])
      SD <- sum(SD_vec[1:H])
      
      G <- ifelse(SC == 0, 0, SD / ifelse(H <= aux_max_SC, SC, n_m))
      
      G_matrix[j, i] <- G
      theta[j, i] <- G - R
      
      if (verbose){ 
        cat("   j =", j, "| H =", H, "| SC =", SC, "| SD =", SD, "| G =", round(G, 4), "| theta =", round(theta[j,i], 4), "\n")
    }
  }
  
  theta <- round(theta, 4)
 
  }
  
  if (verbose){ 
  print(theta)
  }
  
  V_est <- gerar_W_est(theta,csi1, csi2)
    
  return(list(
    R_list = R_list,
    G_matrix = G_matrix,
    theta = theta,
    V_est = V_est
  ))
  
} 


W_est <- gerar_W_est(theta = est$theta, csi1 = 0.005, csi2 = 0.11)

######################################################################################################################## o
######################################################################################################################## o
#####  3: Matriz de comparacao V com V_est ----
#diferente pq agora nas cores tem excitatorio e iniborio

plot_comparacao_V_V_est_excit_inib <- function(V_list, V_est_list, n_list) {
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
  
  # 🔁 Loop para gerar um gráfico por matriz
  for (i in seq_along(V_list)) {
    n <- n_list[[i]]
    V <- V_list[[i]]
    V_est <- V_est_list[[i]]
    
    # 📦 Montar dataframe para plot
    heatmap_data <- data.frame(
      Presynaptic = rep(1:nrow(V), each = ncol(V)),
      Postsynaptic = rep(1:ncol(V), times = nrow(V)),
      V = as.vector(t(V)),
      V_est = as.vector(t(V_est))
    )
    
    # 🏷️ Classificação das condições
    heatmap_data$Condition <- with(heatmap_data, ifelse(
      is.na(V_est), "Inconclusivo",
      ifelse(V_est == 1 & V == 1, "Conexão excitatória existente (acerto)",
             ifelse(V_est == -1 & V == -1, "Conexão inibitória existente (acerto)",
                    ifelse(V_est == 0 & V == 0, "Sem conexão (acerto)",
                           ifelse(V_est == 0 & V != 0, "Conexão não identificada (falso negativo)",
                                ifelse(V_est != 0 & V == 0, "Conexão não existente (falso positivo)",
                                  ifelse(V_est != 0 & V != 0 & V_est != V, "Conexão existente invertida", "Outro")))))))
    )
    
    # 📌 Título automático baseado em n
    subtitle <- bquote("n" == 10^.(log10(n)))
    
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

#### Exemplo 1 ----

# Algoritmos de estimação para n=10^d, d=4,5,6
est_pares1_1<-estimador_redeneural_pares(mod1_1$X,mod1_1$u)
est_pares1_2<-estimador_redeneural_pares(mod1_2$X,mod1_2$u)
est_pares1_3<-estimador_redeneural_pares(mod1_3$X,mod1_3$u)

# Gráfico comparando V com V_est
V_list1 <- list(mod1_1$W,mod1_2$W,mod1_3$W)
V_list_est_par1 <- list(est_pares1_1$V_est,est_pares1_2$V_est,est_pares1_3$V_est)
n_list1 <- list(ncol(mod1_1$X),ncol(mod1_2$X),ncol(mod1_3$X))

plot_comparacao_V_V_est_excit_inib(V_list1,V_list_est_par1,n_list1)

#### Exemplo 2 ----

# Algoritmos de estimação para n=10^d, d=4,5,6
est_pares2_1<-estimador_redeneural_pares(mod2_1$X,mod2_1$u)
est_pares2_2<-estimador_redeneural_pares(mod2_2$X,mod2_2$u)
est_pares2_3<-estimador_redeneural_pares(mod2_3$X,mod2_3$u)

# Gráfico comparando V com V_est
V_list2 <- list(mod2_1$W,mod2_2$W,mod2_3$W)
V_list_est_par2 <- list(est_pares2_1$V_est,est_pares2_2$V_est,est_pares2_3$V_est)
n_list2 <- list(ncol(mod2_1$X),ncol(mod2_2$X),ncol(mod2_3$X))

plot_comparacao_V_V_est_excit_inib(V_list2,V_list_est_par2,n_list2)

#### Exemplo 3 ----

# Algoritmos de estimação para n=10^d, d=4,5,6
est_pares3_1<-estimador_redeneural_pares(mod3_1$X,mod3_1$u)
est_pares3_2<-estimador_redeneural_pares(mod3_2$X,mod3_2$u)
est_pares3_3<-estimador_redeneural_pares(mod3_3$X,mod3_3$u)

# Gráfico comparando V com V_est
V_list3 <- list(mod3_1$W,mod3_2$W,mod3_3$W)
V_list_est_par3 <- list(est_pares3_1$V_est,est_pares3_2$V_est,est_pares3_3$V_est)
n_list3 <- list(ncol(mod3_1$X),ncol(mod3_2$X),ncol(mod3_3$X))

plot_comparacao_V_V_est_excit_inib(V_list3,V_list_est_par3,n_list3)

######################################################################################################################## o
######################################################################################################################## o
#######  5: Grafico percentual de acerto ----

# Criar dataframe
dados <- data.frame(
  n = factor(rep(c("10.000", "100.000", "1.000.000"), each = 3),
             levels = c("10.000", "100.000", "1.000.000")),
  Exemplo = factor(rep(c("Exemplo 1", "Exemplo 2", "Exemplo 3"), times = 3)),
  Acertos_total = c(80, 71, 61, 90, 100, 84, 95, 100, 93),
  Acerto_conexao = c(100, 85, 55, 89, 100, 66, 100, 100, 82)
)

# Gráfico
ggplot(dados, aes(x = n)) +
  geom_bar(
    aes(y = Acerto_conexao, fill = Exemplo),
    stat = "identity",
    position = position_dodge(width = 0.7),
    width = 0.6
  ) +
  geom_point(
    aes(y = Acertos_total, color = Exemplo),
    position = position_dodge(width = 0.7),
    size = 4,
    shape = 21,
    fill = "white",
    alpha = 0.8,        # Transparência dos pontos
    stroke = 1.3,
    show.legend = FALSE
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0, 110), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Tamanho da amostra (n)",
    y = "Acertos de conexão pelo estimador (%)",
    fill = NULL,
    color = "Exemplo"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),  # força a linha horizontal em y=0
    axis.ticks = element_line(color = "black")
  )