#### plot funcao phi(u)

plot_phi <- function(u_values) {
  # Definir o intervalo para o eixo x e criar os valores x
  x_values <- seq(-5, 15, by = 0.01)
  colors <- c("blue", "orange", "green")  # Cores para os valores de u = 1, 4 e 7
  
  # Configurar o grafico base com eixos e linhas de referencia, sem contorno
  plot(x_values, phi(x_values, u_values[1]), type = "n", 
       xlab = "", ylab = "", xlim = c(-5, 10), ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n")
  
  # Adicionar linhas pretas para o plano cartesiano (eixos x = 0 e y = 0)
  abline(h = 0, col = "black", lwd = 1.5)
  abline(v = 0, col = "black", lwd = 1.5)
  
  # Adicionar linhas horizontais para 0.02 e 0.8 (todas em vermelho)
  abline(h = 0.02, col = "red", lty = 2)
  abline(h = 0.8, col = "red", lty = 2)
  
  # Loop para cada valor de u, com cores diferentes para cada linha
  for (i in seq_along(u_values)) {
    u <- u_values[i]
    y_values <- phi(x_values, u)
    lines(x_values, y_values, col = colors[i], lwd = 2)
    # Adiciona a linha vertical para x = u
    abline(v = u, col = colors[i], lty = 2)
  }
  
  # Adicionando a legenda no canto superior esquerdo com contorno preto
  legend("topleft", legend = paste("u =", u_values), col = colors, lwd = 2, 
         inset = 0.05, box.lty = 1, box.col = "black", bg = "white", title = "Valores de u")
  
  # Configurar os eixos para exibir os numeros diretamente nas linhas do plano cartesiano
  axis(1, at = seq(-4, 10, by = 2), pos = 0, lwd = 0, lwd.ticks = 1, labels = seq(-4, 10, by = 2))  # Eixo x pares
  axis(2, at = seq(0.2, 1, by = 0.2), pos = 0, lwd = 0, lwd.ticks = 1, labels = seq(0.2, 1, by = 0.2))  # Eixo y sem o 0
  
  # Adicionar os rotulos dos eixos
  mtext(expression(U[i](t)), side = 1, line = 2.5)  # Rotulo para o eixo x
  mtext(expression(phi(U[i](t))), side = 2, line = 2.5)  # Rotulo para o eixo y
}

sigmoid <- function(x, b, c, p_min, p_max) {
  # Função sigmoidal ajustada com limites
  p_min + (p_max - p_min) / (1 + exp(-b * (x - c)))
}

phi <- function(x, u) {
  # Novas condições
  p_min <- 0.02  # f(U -> -inf) = 0.02
  p0 <- 0.05     # f(U < 0) = 0.05
  pu <- 0.8      # f(U < u) = 0.8
  p_max <- 1     # Limite superior da sigmoide
  
  # Cálculo dos parâmetros b e c
  b <- log((pu - p_min) / (p_max - pu) * (p_max - p0) / (p0 - p_min)) / u
  c <- -log((p0 - p_min) / (p_max - p0)) / b
  
  # Retorna a probabilidade ajustada
  return(sigmoid(x, b, c, p_min, p_max))
}

plot_phi(c(1, 4, 7))


#### plot do grafo de conectividade neuronal W 

library(igraph)

# Funcao para criar o grafo de peso sinaptico da rede neuronal com arestas paralelas e sentido horário
plot_neural <- function(W) {
  # Criar grafo a partir da matriz de adjacencia
  graph <- graph_from_adjacency_matrix(W, mode = "directed", weighted = NULL, diag = FALSE)
  
  # Definir layout circular com rotacao para garantir o no 1 no topo e sentido horario
  layout <- layout_in_circle(graph)
  layout <- layout[order(V(graph)), ]  # Ordenar para garantir o no 1 no topo
  layout <- layout %*% matrix(c(cos(pi/2), -sin(pi/2), sin(pi/2), cos(pi/2)), ncol = 2)  # Rotacao de 90 graus
  
  # Identificar arestas bidirecionais e aplicar curvatura personalizada
  edge_curvatures <- ifelse(apply(get.edgelist(graph), 1, function(x) any(E(graph)[from(x[2]) & to(x[1])])),
                            0.3, 0)  # Curvatura apenas para arestas bidirecionais
  
  # Plotar o grafo com os ajustes
  plot(graph, 
       layout = layout,
       vertex.size = 30,                  # Aumentar o tamanho dos nós
       vertex.color = "white",            # Cor dos nós
       vertex.label.color = "black",      # Cor dos rótulos dos nós
       vertex.label.cex = 1.2,            # Tamanho dos rótulos
       edge.color = "black",              # Cor das arestas
       edge.width = 1,                    # Espessura das arestas
       edge.arrow.size = 0.6,             # Tamanho das setas
       edge.curved = edge_curvatures      # Curvatura para bidirecionais
  )
  
  # Adicionar titulo ao grafico
  title("Grafo de peso sinaptico da rede neuronal", cex.main = 1.2, line = 1.5)
}

W <- exemplo1$W

# Chamando a funcao para exibir o grafo
plot_neural(W)


####Funcao para a Matriz de Conectividade Neuronal V

library(ggplot2)

plot_matrix <- function(Matrix) {
  n <- nrow(Matrix)
  
  # Transformar a matriz em formato longo para ggplot
  matrix_data <- data.frame(
    Row = rep(1:n, each = n),
    Col = rep(1:n, times = n),
    Value = as.vector(t(Matrix))  # Transpor a matriz para alinhar com os eixos do grafico
  )
  
  # Definir as cores para a conexao
  colors <- c("Connection" = "lightblue", "No connection" = "white")
  
  # Converter valores para strings para mapeamento de cor
  matrix_data$Value <- as.character(matrix_data$Value)
  matrix_data$Value[matrix_data$Value == "0"] <- "No connection"  # Substituir 0 por "No connection"
  matrix_data$Value[matrix_data$Value == "1"] <- "Connection"     # Substituir 1 por "Connection"
  
  # Plotar a matriz
  ggplot(matrix_data, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "black") +
    scale_fill_manual(
      values = colors, 
      name = NULL,  # Remover título da legenda
      labels = c("Conexao entre os neuronios", "Sem conexao entre os neuronios"),
      guide = guide_legend(title.position = "top", nrow = 1)
    ) +
    scale_x_continuous(breaks = 1:n, position = "top") +
    scale_y_reverse(breaks = 1:n) +  # Inverter para que 1 fique no topo
    labs(x = "Neuronios pos-sinapticos", y = "Neuronios pre-sinapticos", title = "Matriz de conectividade neuronal") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),    # Afastar titulo do eixo x
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),    # Afastar titulo do eixo y
      legend.position = "top",                                            # Legenda acima do grafico
      legend.text = element_text(size = 8),
      legend.box.spacing = unit(0.2, "cm"),                               # Espacamento entre legenda e grafico
      legend.key.width = unit(0.8, "cm"),                                 # Largura compacta da legenda
      axis.text.x = element_text(margin = margin(t = 2)),                 # Aproximar rotulos dos eixos
      axis.text.y = element_text(margin = margin(r = 2)),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 10))     # Centralizar titulo e ajustar margem
    )
}

plot_matrix(exemplo1$W)


##### 3 graficos juntos
library(ggplot2)
library(gridExtra)

plot_matrix_list <- function(matrix_list, titulo_list) {
  n <- nrow(matrix_list[[1]])  # Assume que todas as matrizes têm o mesmo tamanho
  
  # Função interna para criar um gráfico individual com base na matriz e no título
  create_plot <- function(Matrix, subtitle) {
    # Transformar a matriz em formato longo para ggplot
    matrix_data <- data.frame(
      Row = rep(1:n, each = n),
      Col = rep(1:n, times = n),
      Value = as.vector(t(Matrix))
    )
    
    # Definir as cores para cada condição
    colors <- c("Connection" = "lightblue", "No connection" = "white", "Inconclusivo" = "grey")
    
    # Converter valores para strings para mapeamento de cor
    matrix_data$Value <- as.character(matrix_data$Value)
    matrix_data$Value[matrix_data$Value == "0"] <- "No connection"
    matrix_data$Value[matrix_data$Value == "1"] <- "Connection"
    matrix_data$Value[is.na(matrix_data$Value)] <- "Inconclusivo"
    
    # Criar o gráfico
    p <- ggplot(matrix_data, aes(x = Col, y = Row, fill = Value)) +
      geom_tile(color = "black") +
      scale_fill_manual(
        values = colors,
        name = NULL,  # Título da legenda removido
        labels = c("Conexao entre os neuronios", "Sem conexao entre os neuronios", "Inconclusivo")
      ) +
      scale_x_continuous(breaks = 1:n, position = "top") +
      scale_y_reverse(breaks = 1:n) +
      labs(x = "Neuronios pos-sinapticos", y = "Neuronios pre-sinapticos", subtitle = subtitle) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 10, margin = margin(t = 5)),
        axis.title.y = element_text(size = 10, margin = margin(r = 5)),
        legend.position = "bottom",                                           # Legenda abaixo do gráfico
        legend.text = element_text(size = 8),
        legend.box.spacing = unit(0.1, "cm"),
        legend.key.width = unit(0.8, "cm"),
        axis.text.x = element_text(margin = margin(t = 2)),
        axis.text.y = element_text(margin = margin(r = 2)),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 10)
      )
    return(p)
  }
  
  # Gerar gráficos individuais e armazenar em uma lista
  plot_list <- mapply(create_plot, matrix_list, titulo_list, SIMPLIFY = FALSE)
  
  # Exibir os gráficos lado a lado com uma única legenda centralizada
  grid.arrange(grobs = plot_list, ncol = 3, top = textGrob("Aplicação do estimador de blocos conectividade neuronal real em dados reais", gp = gpar(fontsize = 14)))
}


matrix_list <- list(est_basereal$V_est, est_basereal1$V_est, est_basereal2$V_est)
titulo_list <- c("Base inteira", "1ª Metade da Base", "2ª Metade da Base")

plot_matrix_list(matrix_list, titulo_list)


library(ggplot2)

plot_matrix <- function(Matrix) {
  n <- nrow(Matrix)
  
  # Transformar a matriz em formato longo para ggplot
  matrix_data <- data.frame(
    Row = rep(1:n, each = n),
    Col = rep(1:n, times = n),
    Value = as.vector(t(Matrix))  # Transpor a matriz para alinhar com os eixos do gráfico
  )
  
  # Definir as cores para cada condição
  colors <- c("Connection" = "lightblue", "No connection" = "white", "Inconclusivo" = "grey")
  
  # Converter valores para strings para mapeamento de cor
  matrix_data$Value <- as.character(matrix_data$Value)
  matrix_data$Value[matrix_data$Value == "0"] <- "No connection"  # Substituir 0 por "No connection"
  matrix_data$Value[matrix_data$Value == "1"] <- "Connection"     # Substituir 1 por "Connection"
  matrix_data$Value[is.na(matrix_data$Value)] <- "Inconclusivo"   # Mapear NA como "Inconclusivo"
  
  # Plotar a matriz
  ggplot(matrix_data, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "black") +
    scale_fill_manual(
      values = colors, 
      name = NULL,  # Remover título da legenda
      labels = c("Conexao entre os neuronios", "Sem conexao entre os neuronios", "Inconclusivo"),
      guide = guide_legend(title.position = "top", nrow = 1)
    ) +
    scale_x_continuous(breaks = 1:n, position = "top") +
    scale_y_reverse(breaks = 1:n) +  # Inverter para que 1 fique no topo
    labs(x = "Neuronios pos-sinapticos", y = "Neuronios pre-sinapticos", title = "Matriz de conectividade neuronal") +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),    # Afastar título do eixo x
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),    # Afastar título do eixo y
      legend.position = "top",                                            # Legenda acima do gráfico
      legend.text = element_text(size = 8),
      legend.box.spacing = unit(0.2, "cm"),                               # Espaçamento entre legenda e gráfico
      legend.key.width = unit(0.8, "cm"),                                 # Largura compacta da legenda
      axis.text.x = element_text(margin = margin(t = 2)),                 # Aproximar rótulos dos eixos
      axis.text.y = element_text(margin = margin(r = 2)),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 10))     # Centralizar título e ajustar margem
    )
}

# Exemplo de uso com uma matriz de exemplo
# exemplo1$W <- matrix(c(1, 0, NA, 1, 0, 1, 1, NA, 0), 3, 3)
plot_matrix(est_basereal$V_est)

