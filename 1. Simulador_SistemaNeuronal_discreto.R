######################################################################################################################## o
####  Simulador para gerar uma sequência de disparos aleatórios com uma rede de m neurônios
######################################################################################################################## o
####  Autor: Diogo Maia de Figueiredo
####  Orientadora: Dra. Aline Duarte
####  Programa de Pós Graduação do IME-USP

######################################################################################################################## o
######################################################################################################################## o
#######  0: Instalação de Pacotes ----

pacotes_necessarios <- c(
  "dplyr",
  "ggplot2",
  "grid",
  "gridExtra",
  "igraph",
  "reshape2",
  "scales",
  "tidyr"
)

# Função completa para instalar, carregar e informar no console
instalar_e_carregar <- function(pacotes) {
  for (pacote in pacotes) {
    if (!require(pacote, character.only = TRUE)) {
      cat(sprintf("📦 Instalando pacote: %s...\n", pacote))
      install.packages(pacote, dependencies = TRUE)
      cat(sprintf("✅ Pacote %s instalado.\n", pacote))
    } else {
      cat(sprintf("✅ Pacote %s já instalado.\n", pacote))
    }
    library(pacote, character.only = TRUE)
  }
  cat("\n🎯 Todos os pacotes foram carregados com sucesso!\n")
}

# Executar
instalar_e_carregar(pacotes_necessarios)

#  . -----
######################################################################################################################## o
######################################################################################################################## o
#######  1: Funções para definição de parâmetros ----

######  1.1 Gerar Matrizes de peso sináptico ----
# Cada neurônio tem prob de 0.2 de ser excitatório e 0.05 de ser inibitório

####  Matriz de conectividade neuronal

gerador_matriz_W <- function(m, max_tentativas = 1000) {
  for (tentativa in 1:max_tentativas) {
    W <- matrix(0, m, m)
    for (i in 1:m) {
      for (j in 1:m) {
        if (i != j) {
          p <- runif(1)
          if (p < 0.2) {
            W[i, j] <- 1
          } else if (p < 0.25) {
            W[i, j] <- -1
          } else {
            W[i, j] <- 0
          }
        }
      }
    }
    linhas_com_zeros <- apply(W, 1, function(row) all(row == 0))
    colunas_com_zeros <- apply(W, 2, function(col) all(col == 0))
    if (!any(linhas_com_zeros) && !any(colunas_com_zeros)) {
      return(W)
    }
  }
  stop("Não foi possível gerar matriz W válida após muitas tentativas")
}


####  Grafo de peso sináptico

plot_neural <- function(W) {
  # Criar grafo com pesos
  graph <- graph_from_adjacency_matrix(W != 0, mode = "directed", weighted = NULL, diag = FALSE)
  E(graph)$weight <- W[as_edgelist(graph)]
  n <- vcount(graph)
  
  angles <- seq(from = pi/2, to = pi/2 - 2*pi + 2*pi/n, length.out = n)
  layout <- cbind(cos(angles), sin(angles))
  
  edge_list <- as_edgelist(graph)
  edge_curvatures <- sapply(1:nrow(edge_list), function(i) {
    from_node <- edge_list[i, 1]
    to_node <- edge_list[i, 2]
    is_bidirectional <- any(edge_list[,1] == to_node & edge_list[,2] == from_node)
    if (is_bidirectional) 0.3 else 0
  })
  
  edge_lty <- ifelse(E(graph)$weight > 0, 1, 2)
  edge_labels <- ifelse(abs(E(graph)$weight) == 1, "", as.character(E(graph)$weight))
  
  coords <- layout
  edges <- as_edgelist(graph, names = FALSE)
  edge_centers <- t(apply(edges, 1, function(e) {
    (coords[e[1], ] + coords[e[2], ]) / 2
  }))
  
  offset <- 0.08
  edge_label_x <- edge_centers[, 1] + offset * cos(atan2(edge_centers[, 2], edge_centers[, 1]))
  edge_label_y <- edge_centers[, 2] + offset * sin(atan2(edge_centers[, 2], edge_centers[, 1]))
  
  op <- par(no.readonly = TRUE)  # salvar parâmetros
  layout(matrix(1:2, ncol = 1), heights = c(1.5, 8))  # define layout antes
  
  # Parte 1: legenda no topo
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", 
         legend = c("Conexão excitatória", "", "Conexão inibitória"),  # item vazio entre os dois
         lty = c(1, NA, 2), 
         lwd = c(2, NA, 2),
         col = c("black", NA, "black"),
         pt.cex = 1.5,
         horiz = TRUE,
         x.intersp = 1.2,
         y.intersp = 1.5,
         bty = "n",           # remove contorno
         seg.len = 2,
         cex = 1.2,
         xpd = TRUE)
  
  # Parte 2: grafo
  par(mar = c(1, 1, 1, 1))
  plot(graph, 
       layout = layout,
       vertex.size = 30,
       vertex.color = "white",
       vertex.label.color = "black",
       vertex.label.cex = 1.6,
       edge.color = "black",
       edge.width = 1,
       edge.arrow.size = 0.6,
       edge.curved = edge_curvatures,
       edge.lty = edge_lty,
       edge.label = edge_labels,
       edge.label.cex = 1.1,
       edge.label.color = "black",
       edge.label.family = "sans",
       edge.label.x = edge_label_x,
       edge.label.y = edge_label_y
  )
  par(op)  # restaurar layout original
}

####  Matriz de peso sináptico

plot_matrix <- function(Matrix) {
  n <- nrow(Matrix)
  
  # Transformar a matriz em formato longo
  matrix_data <- data.frame(
    Row = rep(1:n, each = n),
    Col = rep(1:n, times = n),
    Value = as.vector(t(Matrix))  # Transpor para alinhar com os eixos
  )
  
  # Converter para fator com labels descritivos
  matrix_data$Value[matrix_data$Value == 1] <- "Conexão excitatória"
  matrix_data$Value[matrix_data$Value == -1] <- "Conexão inibitória"
  matrix_data$Value[matrix_data$Value == 0] <- "Sem conexão"
  
  matrix_data$Value <- factor(matrix_data$Value,
                              levels = c("Conexão excitatória",
                                         "Conexão inibitória",
                                         "Sem conexão"))
  
  # Definir cores
  colors <- c(
    "Conexão excitatória" = "lightblue",
    "Conexão inibitória" = "lightpink",
    "Sem conexão" = "white"
  )
  
  # Plotar
  ggplot(matrix_data, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "black") +
    scale_fill_manual(
      values = colors,
      name = NULL,
      guide = guide_legend(title.position = "top", nrow = 1)
    ) +
    scale_x_continuous(breaks = 1:n, position = "top") +
    scale_y_reverse(breaks = 1:n) +
    labs(x = "Neurônios pós-sinápticos", 
         y = "Neurônios pré-sinápticos") +
         #title = "Matriz de conectividade neuronal") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      legend.position = "top",
      legend.text = element_text(size = 10),
      legend.key.size   = unit(0.9, "cm"), 
      legend.box.spacing = unit(0.3, "cm"),
      legend.key.width = unit(0.8, "cm"),
      axis.text.x = element_text(size = 12, margin = margin(t = 2)),
      axis.text.y = element_text(size = 12, margin = margin(r = 2)),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 10))
    )
}

######################################################################################################################## o
######  1.2 Função de Disparo ----
# Probabilidade do neuronio i disparar no instante t com base em seu potencial de membrana em t-1

#####  Funcao sigmoide ajustada

sigmoid <- function(x, b, c, p_min, p_max) {
  # Função sigmoidal ajustada com limites
  p_min + (p_max - p_min) / (1 + exp(-b * (x - c)))
}

# Funcao que calcula a probabilidade baseada na funcao sigmoide ajustada
phi<- function(x, u) {
  
  # Condicoes
  p_min <- 0.02  # f(U -> -inf) = 0.02
  p0 <- 0.05     # f(U < 0) = 0.05
  pu <- 0.7      # f(U < u) = 0.8
  p_max <- 1     # Limite superior da sigmoide
  
  # Calculo dos parametros b e c
  b <- log((pu - p_min) / (p_max - pu) * (p_max - p0) / (p0 - p_min)) / u
  c <- -log((p0 - p_min) / (p_max - p0)) / b
  
  # Retorna a probabilidade ajustada
  return(sigmoid(x, b, c, p_min, p_max))
}

#####  Gráfico da Função de disparo

plot_phi <- function(u_values) {
  # Definir o intervalo para o eixo x e criar os valores x
  x_values <- seq(-6, 21, by = 0.01)
  colors <- c("blue", "orange", "green")
  
  # Configurar o grafico base com eixos e linhas de referencia, sem contorno
  plot(x_values, phi(x_values, u_values[1]), type = "n", 
       xlab = "", ylab = "", xlim = c(-6, 20), ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n")
  
  # Adicionar linhas pretas para o plano cartesiano (eixos x = 0 e y = 0)
  abline(h = 0, col = "black", lwd = 1.5)
  abline(v = 0, col = "black", lwd = 1.5)
  
  # Adicionar linhas horizontais para 0.02 e 0.8 (todas em vermelho)
  abline(h = 0.02, col = "red", lty = 2)
  abline(h = 0.7, col = "red", lty = 2)
  
  # Loop para cada valor de u, com cores diferentes para cada linha
  for (i in seq_along(u_values)) {
    u <- u_values[i]
    y_values <- phi(x_values, u)
    lines(x_values, y_values, col = colors[i], lwd = 2)
    # Adiciona a linha vertical para x = u
    abline(v = u, col = colors[i], lty = 2)
  }
  
  # Adicionando a legenda no canto superior esquerdo com contorno preto
  legend(x = par("usr")[1], y = par("usr")[4], legend = paste("u =", u_values), col = colors, lwd = 2, 
         inset = 0.05, box.lty = 1, box.col = "white", bg = "white", title = "Valores de u",cex = 1.2)
  
  # Configurar os eixos para exibir os numeros diretamente nas linhas do plano cartesiano
  axis(1, at = seq(-6, 20, by = 2), pos = 0, lwd = 0, lwd.ticks = 1, labels = seq(-6, 20, by = 2))  # Eixo x pares
  axis(2, at = seq(0.1, 1, by = 0.1), pos = 0, lwd = 0, lwd.ticks = 1, labels = seq(0.1, 1, by = 0.1))  # Eixo y sem o 0
  
  # Adicionar os rotulos dos eixos
  mtext(expression(U[t](i)), side = 1, line = 2.5)  # Rotulo para o eixo x
  mtext(expression(varphi* "'" *(U[t](i)*","*u)), side = 2, line = 2.5)  # Rotulo para o eixo y
}

#  . -----
######################################################################################################################## o
######################################################################################################################## o

#######  2: Simulação da Rede ----

#####  2.1 Gerador de Disparos ----

simulador_redeneural_discreto<-function(m,n,W=matrix(0, m, m)){  # m=neuronios na rede, n= instantes gerados
  
  # Ler a matrix de conectividade
  if (is.matrix(W) && all(W == 0)){
    W <- gerador_matriz_W(m)
  }
  V <- ifelse(W == 0, 0, 1)
  
  # Ler canal de vazamento
  g <- rep(0.8, m)
  
  # Definir função de disparo
  u<-max(1,mean(colSums(abs(W))))
  
  # Estado inicial U_0  
  gerador_U0 <- function(m, u) {
    u <- round(u)  # Garantir que u é inteiro
    U_0 <- sample(0:u, m, replace = TRUE)
    while (all(U_0 == 0)) {
      U_0 <- sample(0:u, m, replace = TRUE)
    }
    return(U_0)
  }
  
  #Vetor U_0
  U_0 <- gerador_U0(m, u)
  
  # Calculo de U e X
  # U[i,t] -> potencial de disparo do neuronio i no instante s
  # phi_t[i,t] -> prob. de i disparar no proximo instante
  #X[i,t] -> 1 se i disparou  no instante t, 0, cc.
  
  U <- matrix(0, nrow = m, ncol = n)
  phi_t <- matrix(0, nrow = m, ncol = n)
  X <- matrix(0, nrow = m, ncol = n)
  
  # Calculo para t = 1
  for (i in 1:m) {
    X[i, 1] <- rbinom(1, size = 1, prob = phi(U_0[i], u))
  }
  
  for(i in 1:m){
    if (X[i, 1] == 1) {
      U[i, 1] <- 0
    } else {
      U[i, 1] <- U_0[i] * g[i] + sum(X[, 1] * W[, i])
    }
    phi_t[i, 1] <- phi(U[i, 1], u)
  }
  
  # Calculo para t > 1
  for (t in 2:n) {
    for (i in 1:m) {
      X[i, t] <- rbinom(1, size = 1, prob = phi_t[i, t-1])
    }
    for (i in 1:m){
      if (X[i, t] == 1) {
        U[i, t] <- 0
      } else {
        U[i, t] <- U[i, t-1] * g[i] + sum(X[, t] * W[, i])
      }
      phi_t[i, t] <- phi(U[i, t], u)
    }
  }
  
  return(list(W = W, u = u, U_0 = U_0, U = U, X = X, g = g , phi_t = phi_t,n=n,m=m,V=V))
}


#####  2.2 L_sup ----

L_sup <- function(X, t = ncol(X)) {
  # Verifica se t e valido
  if (t > ncol(X)) {
    stop("Erro: 't' deve ser menor ou igual ao numero de colunas de X.")
  }
  
  # Inicializa o vetor de saida
  L_sup <- numeric(nrow(X))
  
  # Para cada linha de X
  for (aux in 1:nrow(X)) {
    # Encontra os indices onde X[a,t] == 1 para t < s
    idx <- which(X[aux, 1:t] == 1)
    
    # Se houver pelo menos um 1, pegue o maior valor de t
    if (length(idx) > 0) {
      L_sup[aux] <- max(idx)
    } else {
      # Se todos forem 0, preenche com -Inf
      L_sup[aux] <- -Inf
    }
  }
  return(L_sup)
}

#  . -----
######################################################################################################################## o
######################################################################################################################## o

#######  3: Representações Gráficas ----

#### 3.1 Gráficos Disparos + phi ----

plot_redeneural_discreto <- function(X, phi_t) {
  m <- nrow(X)
  n <- ncol(X)
  aux <- matrix(0, nrow = m, ncol = n)
  
  # Preencher a matriz auxiliar de acordo com a condicao
  for (i in 1:m) {
    for (t in 1:n) {
      if (X[i, t] == 1) {
        aux[i, t] <- X[i, t]
      } else {
        aux[i, t] <- phi_t[i, t]
      }
    }
  }
  
  # Converter a matriz aux para um data frame para facilitar a plotagem
  df_aux <- expand.grid(Neuronio = 1:m, Tempo = 1:n)
  df_aux$Valor <- as.vector(aux)
  
  # Criar o grafico com ggplot2
  ggplot(df_aux, aes(x = Tempo, y = Neuronio)) +
    # Desenhar os tracos verticais onde aux[x, y] = 1
    geom_segment(data = subset(df_aux, Valor == 1), aes(x = Tempo, xend = Tempo, y = Neuronio - 0.2, yend = Neuronio + 0.2), 
                 color = "black", size = 0.4) +
    # Desenhar as bolinhas coloridas onde aux[x, y] != 1
    geom_point(data = subset(df_aux, Valor != 1), aes(color = Valor), size = 3) +
    # Definir a escala de cores de branco para azul
    scale_color_gradient(low = "white", high = "blue", name = expression(varphi* "'" *(U[t](i)*","*u))) +
    # Configurar o eixo x para mostrar apenas numeros inteiros, comecar exatamente em 0 e remover expansao
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, n), expand = c(0, 0)) +
    # Configurar o eixo y para mostrar apenas numeros inteiros e inverter o eixo y
    scale_y_reverse(breaks = seq(1, m, by = 1)) +
    # Adicionar as linhas dos eixos
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # Remove as grades
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, size = 0.2),  # Adiciona as linhas dos eixos
      aspect.ratio = 0.25  # Ajuste da proporcao do aspecto para um grafico mais retangular
    ) +
    labs(x = "Tempo (t)", y = "Neuronios (i)")
}


#### 3.2 Gráficos Disparos (só traço) ----
plot_redeneuraldisparo_discreto <- function(X) {
  m <- nrow(X)
  n <- ncol(X)
  
  # Converter X para um data frame, mantendo apenas os pontos onde X == 1
  df_X <- expand.grid(Neuronio = 1:m, Tempo = 1:n)
  df_X$Valor <- as.vector(X)
  df_X <- subset(df_X, Valor == 1)  # Filtrar para manter apenas os valores onde X == 1
  
  # Criar o gráfico com ggplot2
  ggplot(df_X, aes(x = Tempo, y = Neuronio)) +
    # Desenhar os traços verticais onde X[i, t] = 1
    geom_segment(aes(x = Tempo, xend = Tempo, y = Neuronio - 0.2, yend = Neuronio + 0.2), 
                 color = "black", size = 0.4) +
    # Configurar o eixo x para mostrar apenas números inteiros, começar exatamente em 0 e remover expansão
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0, n), expand = c(0, 0)) +
    # Configurar o eixo y para mostrar apenas números inteiros e inverter o eixo y
    scale_y_reverse(breaks = seq(1, m, by = 1)) +
    # Adicionar as linhas dos eixos
    theme_minimal() +
    theme(
      panel.grid = element_blank(),  # Remove as grades
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, size = 0.2),  # Adiciona as linhas dos eixos
      aspect.ratio = 0.4  # Ajuste da proporção do aspecto para um gráfico mais retangular
    ) +
    labs(x = "Tempo (t)", y = "Neurônios (i)")
}

#  . -----
######################################################################################################################## o
######################################################################################################################## o

#######  4: Geração de Imagens ----

# Grafo ilustrativo

W0 <- matrix(0, 5, 5)
W0[1,2] <- 1
W0[2,3] <- 1
W0[3,1] <- -2
W0[5,1] <-2
W0[4,3] <- 3
W0[3,5] <- -1
W0[2,4] <- 1
plot_neural(W0)

# Funções de disparo
plot_phi(c(3, 8, 12))


#### Exemplo 1 ----

set.seed(844105)
W1 <- gerador_matriz_W(5)
W1
plot_neural(W1)
plot_matrix(W1)

mod1_1<-simulador_redeneural_discreto(5,10000,W1)
mod1_2<-simulador_redeneural_discreto(5,100000,W1)
mod1_3<-simulador_redeneural_discreto(5,1000000,W1)

#### Exemplo 2 ----

set.seed(326142)
W2 <- gerador_matriz_W(7)
W2
plot_neural(W2)
plot_matrix(W2)

mod2_1<-simulador_redeneural_discreto(7,10000,W2)
mod2_2<-simulador_redeneural_discreto(7,100000,W2)
mod2_3<-simulador_redeneural_discreto(7,1000000,W2)

plot_redeneural_discreto(mod2_1$X[,4450:4550],mod2_1$phi_t[,4450:4550])
plot_redeneuraldisparo_discreto(mod2_1$X[,4450:4550])

#### Exemplo 3 ----
set.seed(492562)
W3 <- gerador_matriz_W(10)
W3
plot_neural(W3)
plot_matrix(W3)

mod3_1<-simulador_redeneural_discreto(10,10000,W3)
mod3_2<-simulador_redeneural_discreto(10,100000,W3)
mod3_3<-simulador_redeneural_discreto(10,1000000,W3)

plot_redeneuraldisparo_discreto(mod3_1$X[,4450:4850])