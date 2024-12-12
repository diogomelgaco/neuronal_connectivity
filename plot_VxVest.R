# Carregar bibliotecas necessárias
library(ggplot2)
library(gridExtra)

# Função para gerar gráficos de comparação
plot_multiple_comparison <- function(V_list, V_est_list, epsilon_list, n_list) {
  
  # Garantir que todos os valores de n são numéricos
  n_list <- lapply(n_list, as.numeric)
  epsilon_list <- lapply(epsilon_list, as.numeric)
  
  # Criar uma lista para armazenar os gráficos
  plot_list <- list()
  
  for (i in 1:length(V_list)) {
    n <- n_list[[i]]
    epsilon <- epsilon_list[[i]]
    V <- V_list[[i]]
    V_est <- V_est_list[[i]]
    
    # Preparar os dados para o gráfico
    heatmap_data <- data.frame(
      Presynaptic = rep(1:nrow(V), each = ncol(V)),
      Postsynaptic = rep(1:ncol(V), times = nrow(V)),
      V = as.vector(t(V)),
      V_est = as.vector(t(V_est))
    )
    
    heatmap_data$Condition <- with(heatmap_data, ifelse(
      is.na(V_est), "Inconclusivo",
      ifelse(V_est == 1 & V == 1, "Conexão existente",
             ifelse(V_est == 0 & V == 0, "Sem conexão",
                    ifelse(V_est == 1 & V == 0, "Conexão falsa (falso positivo)", "Conexão não identificada (falso negativo)")))))
    
    colors <- c("Inconclusivo" = "grey", 
                "Conexão existente" = "blue", 
                "Sem conexão" = "white",
                "Conexão falsa (falso positivo)" = "yellow",
                "Conexão não identificada (falso negativo)" = "red")
    
    # Subtítulo específico para cada gráfico
    subtitle <- bquote("n" == 10^.(log10(n)) ~ "," ~ epsilon == .(epsilon))
    
    # Criar o gráfico
    p <- ggplot(heatmap_data, aes(x = Postsynaptic, y = Presynaptic, fill = Condition)) +
      geom_tile(color = "black", aes(fill = Condition)) +
      scale_fill_manual(values = colors, guide = guide_legend(title = NULL)) + # Remove o título da legenda
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
        legend.position = "bottom" # Coloca a legenda embaixo de cada gráfico
      )
    
    plot_list[[i]] <- p
  }
  
  # Exibir todos os gráficos juntos com o título único
  grid.arrange(
    arrangeGrob(grobs = plot_list, ncol = 3),
    top = textGrob("Comparação entre matriz de conectividade neuronal real e o estimador de blocos", 
                   gp = gpar(fontsize = 16, fontface = "bold"))
  )
}

#Grafico Exemplo 1
V_list1<-list(mod1$V,mod2$V,mod3$V)
V_list_est1<-list(est1$V_est,est2$V_est,est3$V_est)
n_list1<-list(ncol(mod1$X),ncol(mod2$X),ncol(mod3$X))
epsilon_list1<-list(est1$epsilon, est2$epsilon, est3$epsilon)

plot_multiple_comparison(V_list1,V_list_est1, epsilon_list1,n_list1)

#Grafico Exemplo 2
V_list2<-list(mod4$V,mod5$V,mod6$V)
V_list_est2<-list(est4$V_est,est5$V_est,est6$V_est)
n_list2<-list(ncol(mod4$X),ncol(mod5$X),ncol(mod6$X))
epsilon_list2<-list(est4$epsilon, est5$epsilon, est6$epsilon)

plot_multiple_comparison(V_list2,V_list_est2, epsilon_list2,n_list2)