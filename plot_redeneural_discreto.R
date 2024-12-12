# Plot de observacoes de disparo com a prob phi

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
    scale_color_gradient(low = "white", high = "blue", name = expression(phi(U[t](i)))) +
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

# Plot de observacoes de disparo somente com os tracos

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