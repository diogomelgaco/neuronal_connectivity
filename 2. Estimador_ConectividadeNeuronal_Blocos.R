######################################################################################################################## o
####  Estimador de Conectividade Neuronal - Método por Blocos (Duarte 2019)
######################################################################################################################## o
####  Autor: Diogo Maia de Figueiredo
####  Orientadora: Dra. Aline Duarte
####  Programa de Pós Graduação do IME-USP

######################################################################################################################## o
######################################################################################################################## o
#######  1: Construindo o Estimador ----

estimador_redeneural_bloco <- function(X, xi = 0.001, epsilon = 0.05) {
  # X é uma matriz m x n composta por 1 e 0, onde 1 indica que o neurônio i disparou no instante t.
  
  # Definir número de linhas (neurônios) e colunas (instantes)
  n <- ncol(X)       
  m <- nrow(X)
  
  # Listas para armazenar submatrizes tau e probabilidades associadas para cada neurônio
  tau_list <- list()  
  tau_p_list <- list()  
  lambda <- matrix(NA, nrow = m, ncol = m)
  diag(lambda) <- 0
  lambda_total <- matrix(NA, nrow = m, ncol = m)
  diag(lambda_total) <- 0
  
  for (i in 1:m) {    
    tau_full <- list()  
    tau_p_full <- numeric()  
    l <- 1  
    
    repeat {
      # Criar lista de submatrizes w para cada instante considerando uma janela l
      w_list <- lapply(2:(n - l), function(w) {
        sub_matrix <- X[, w:(w + l - 1), drop = FALSE]
        return(sub_matrix)
      })
      
      # Criar lista de disparos anteriores do neurônio i
      i_list_ant <- lapply(1:(n - l - 1), function(s) {
        sub_matrix <- X[i, s:(s + l), drop = FALSE]
        return(sub_matrix)
      })
      
      # Criar lista de disparos posteriores do neurônio i
      i_list_pos <- if (l + 2 <= n) X[i, (l + 2):n, drop = FALSE] else NULL
      
      # Função para filtrar as listas mantendo somente i_list_ant com padrão 10^l
      filtrar_listas <- function(w_list, i_list_ant, i_list_pos) {
        w_list_filtrado <- list()
        i_list_ant_filtrado <- list()
        i_list_pos_filtrado <- list()
        
        for (k in 1:length(i_list_ant)) {
          i_ant <- i_list_ant[[k]]
          if (i_ant[1] == 1 && all(i_ant[-1] == 0)) {
            w_list_filtrado[[length(w_list_filtrado) + 1]] <- w_list[[k]]
            i_list_ant_filtrado[[length(i_list_ant_filtrado) + 1]] <- i_list_ant[[k]]
            if (!is.null(i_list_pos)) {
              i_list_pos_filtrado[[length(i_list_pos_filtrado) + 1]] <- i_list_pos[[k]]
            }
          }
        }
        return(list(w_list = w_list_filtrado, i_list_ant = i_list_ant_filtrado, i_list_pos = i_list_pos_filtrado))
      }
      
      # Aplicar a função de filtragem
      listas_filtradas <- filtrar_listas(w_list, i_list_ant, i_list_pos)
      i_list_ant <- listas_filtradas$i_list_ant
      i_list_pos <- listas_filtradas$i_list_pos
      w_list     <- listas_filtradas$w_list
      
      # Construir conjunto tau e contar frequência de ocorrências
      w_list_str <- sapply(w_list, function(mat) paste(as.vector(mat), collapse = "_"))
      aux_tau <- unique(w_list_str)
      aux_tau_n <- sapply(aux_tau, function(q) sum(w_list_str == q))
      
      # Conferir quantidade de w
      sum_disparo_i <- sum(X[i, 1:(n - l - 1)])
      sum_w <- sum(aux_tau_n)
      
      # Determinar quantidade mínima para tau
      corte <- n^(0.5 + xi)
      
      # Filtrar tau com base na quantidade mínima
      tau <- aux_tau[aux_tau_n > corte]
      
      # Sair do loop se tau estiver vazio
      if (length(tau) <= 1) break
      
      # Definir probabilidade de disparo após w ∈ tau
      tau_n <- aux_tau_n[aux_tau_n > corte]
      tau_i <- sapply(tau, function(q) {
        count <- 0
        for (K in 1:length(i_list_pos)) {
          if (i_list_pos[K] == 1 && w_list_str[K] == q) count <- count + 1
        }
        return(count)
      })
      
      tau_p <- tau_i / tau_n  # Probabilidades associadas a cada w
      #print(tau)
     # print(tau_n)
      #print(tau_i)
      #print(tau_p)
      
      # Converter tau em matrizes
      tau_matrices <- lapply(tau, function(t_str) {
        mat_values <- as.numeric(unlist(strsplit(t_str, "_")))
        return(matrix(mat_values, nrow = m, byrow = FALSE))
      })
      
      # Armazenar valores em listas completas
      tau_full <- c(tau_full, tau_matrices)
      tau_p_full <- c(tau_p_full, as.numeric(tau_p))
      
      # Incrementar o tam1anho do bloco l
      l <- l + 1
    }
    
    tau_list[[i]] <- tau_full
    tau_p_list[[i]] <- tau_p_full
    
    # Calcular lambda para todos os pares j -> i
    for (j in setdiff(1:m, i)) {
      aux_lambda <- numeric(0)
      
      if (length(tau_list[[i]]) > 1) {
        # Comparar submatrizes de tau entre si
        for (w in 1:(length(tau_list[[i]]) - 1)) {
          tau_w <- tau_list[[i]][[w]][-j, ]
          
          list_p_w_v <- sapply((w + 1):length(tau_list[[i]]), function(v) {
            tau_v <- tau_list[[i]][[v]][-j, ]
            if (identical(tau_w, tau_v)) return(tau_p_list[[i]][v]) else return(NA)
          })
          
          aux_lambda_tau <- abs(list_p_w_v - tau_p_list[[i]][w])
          aux_lambda <- c(aux_lambda, max(aux_lambda_tau, na.rm = TRUE))
        }
        
        # Atualizar lambda para j -> i
        lambda[j, i] <- if (length(aux_lambda) > 0) max(aux_lambda, na.rm = TRUE) else NA
      }
      
      lambda_total[j, i] <- lambda[j, i]
    }
  }
  
  lambda_total[is.infinite(lambda_total)] <- NA
  
  # Construir matriz de conectividade estimada V_est
  V_est <- matrix(NA, nrow = m, ncol = m)
  V_est[!is.na(lambda_total) & lambda_total > epsilon] <- 1
  V_est[!is.na(lambda_total) & lambda_total <= epsilon] <- 0
  
  # Retornar resultados
  return(list(
    X = X,
    xi = xi,
    epsilon = epsilon,
    tau_list = tau_list,
    tau_p_list = tau_p_list,
    lambda = lambda,
    lambda_total = lambda_total,
    V_est = V_est
  ))
}

#  . -----
######################################################################################################################## o
######################################################################################################################## o
#######  2: Matriz de Comparação (Real x Est) ----

plot_comparacao_V_V_est <- function(V_list, V_est_list, n_list) {
  
  if (!is.list(V_list)) V_list <- list(V_list)
  if (!is.list(V_est_list)) V_est_list <- list(V_est_list)
  if (!is.list(n_list)) n_list <- list(n_list)
  
  n_list <- lapply(n_list, as.numeric)
  plot_list <- list()
  
  # Paleta de cores com legenda atualizada
  colors <- c(
    "Inconclusivo" = "grey", 
    "Conexão existente (acerto)" = "mediumseagreen", 
    "Sem conexão (acerto)" = "white",
    "Conexão não existente (falso positivo)" = "yellow",
    "Conexão não identificada (falso negativo)" = "sandybrown"
  )
  
  for (i in 1:length(V_list)) {
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
      ifelse(V_est == 1 & V == 1, "Conexão existente (acerto)",
             ifelse(V_est == 0 & V == 0, "Sem conexão (acerto)",
                    ifelse(V_est == 1 & V == 0, "Conexão não existente (falso positivo)", 
                           "Conexão não identificada (falso negativo)")))))
    
    subtitle <- bquote("n" == 10^.(log10(n)) )
    
    p <- ggplot(heatmap_data, aes(x = Postsynaptic, y = Presynaptic, fill = Condition)) +
      geom_tile(color = "black") +
      scale_fill_manual(values = colors) +
      labs(
        x = "Neurônios pós-sinápticos", 
        y = "Neurônios pré-sinápticos", 
        title = subtitle
      ) +
      scale_x_continuous(breaks = 1:nrow(V), position = "top") +
      scale_y_reverse(breaks = 1:ncol(V)) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "none"
      )
    
    plot_list[[i]] <- p
  }
  
  # Criar legenda única baseada nas categorias
  legenda_data <- data.frame(
    Presynaptic = 1:5,
    Postsynaptic = 1:5,
    Condition = factor(names(colors), levels = names(colors))
  )
  
  legenda_plot <- ggplot(legenda_data, aes(x = Postsynaptic, y = Presynaptic, fill = Condition)) +
    geom_tile(color = "black") +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(nrow = 2))+
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.7, "cm")
    )
  
  g_legend <- function(a.gplot) {
    tmp <- ggplotGrob(a.gplot)
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    tmp$grobs[[leg]]
  }
  
  legenda_unica <- g_legend(legenda_plot)
  
  # 🔁 Ajuste dinâmico no número de colunas
  n_graficos <- length(plot_list)
  if (n_graficos == 1) {
    ncol_plot <- 1
  } else if (n_graficos == 2) {
    ncol_plot <- 2
  } else if (n_graficos <= 3) {
    ncol_plot <- 3
  } else {
    ncol_plot <- 3
  }
  
  # Combinar gráficos + legenda
  grid.arrange(
    arrangeGrob(grobs = plot_list, ncol = ncol_plot),
    legenda_unica,
    #top = textGrob("Comparação entre matriz de conectividade neuronal real e o estimador de blocos", 
                  # gp = gpar(fontsize = 16, fontface = "bold")),
    heights = c(10, 1.5)
  )
}



#  . -----
######################################################################################################################## o
######################################################################################################################## o
#######  3: Aplicação nos Exemplos ----

#### Exemplo 1 ----
# Algoritmos de estimação para n=10^d, d=4,5,6
est_bloco1_1<-estimador_redeneural_bloco(mod1_1$X)
est_bloco1_2<-estimador_redeneural_bloco(mod1_2$X)
est_bloco1_3<-estimador_redeneural_bloco(mod1_3$X)

# Gráfico comparando V com V_est
V_list1 <- list(mod1_1$V,mod1_2$V,mod1_3$V)
V_list_est_bloco1 <- list(est_bloco1_1$V_est,est_bloco1_2$V_est,est_bloco1_3$V_est)
n_list1 <- list(ncol(mod1_1$X),ncol(mod1_2$X),ncol(mod1_3$X))

plot_comparacao_V_V_est(V_list1,V_list_est_bloco1,n_list1)
est_bloco1_2$lambda_total

#### Exemplo 2 ----
# Algoritmos de estimação para n=10^d, d=4,5,6
est_bloco2_1<-estimador_redeneural_bloco(mod2_1$X)
est_bloco2_2<-estimador_redeneural_bloco(mod2_2$X)
est_bloco2_3<-estimador_redeneural_bloco(mod2_3$X)

# Gráfico comparando V com V_est
V_list2 <- list(mod2_1$V,mod2_2$V,mod2_3$V)
V_list_est_bloco2 <- list(est_bloco2_1$V_est,est_bloco2_2$V_est,est_bloco2_3$V_est)
n_list2 <- list(ncol(mod2_1$X),ncol(mod2_2$X),ncol(mod2_3$X))

plot_comparacao_V_V_est(V_list2,V_list_est_bloco2,n_list2)

#### Exemplo 3 ----
# Algoritmos de estimação para n=10^d, d=4,5,6
est_bloco3_1<-estimador_redeneural_bloco(mod3_1$X)
est_bloco3_2<-estimador_redeneural_bloco(mod3_2$X)
est_bloco3_3<-estimador_redeneural_bloco(mod3_3$X)

# Gráfico comparando V com V_est
V_list3 <- list(mod3_1$V,mod3_2$V,mod3_3$V)
V_list_est_bloco3 <- list(est_bloco3_1$V_est,est_bloco3_2$V_est,est_bloco3_3$V_est)
n_list3 <- list(ncol(mod3_1$X),ncol(mod3_2$X),ncol(mod3_3$X))

plot_comparacao_V_V_est(V_list3,V_list_est_bloco3,n_list3)

######################################################################################################################## o
######################################################################################################################## o
#######  4: Grafico percentual de acerto ----

# Criar dataframe
dados <- data.frame(
  n = factor(rep(c("10.000", "100.000", "1.000.000"), each = 3),
             levels = c("10.000", "100.000", "1.000.000")),
  Exemplo = factor(rep(c("Exemplo 1", "Exemplo 2", "Exemplo 3"), times = 3)),
  Acertos_total = c(30, 29, 2, 70, 74, 64, 95, 100, 96),
  Acerto_conexao = c(22, 8, 3, 56, 77, 48, 89, 100, 88)
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