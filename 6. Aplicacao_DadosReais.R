######################################################################################################################## o
####  Aplicação dos Estimadores de conectividade neuronal nos dados reais usado em Brochinni (2017)
######################################################################################################################## o
####  Autor: Diogo Maia de Figueiredo
####  Orientadora: Dra. Aline Duarte
####  Programa de Pós Graduação do IME-USP

######################################################################################################################## o
######################################################################################################################## o
#######  1: Extração e Discretização dos dados----

carregar_ou_extrair_LocusData <- function(caminho_RDS = "LocusData.RDS") {
  
  # Se o arquivo RDS já existir, apenas carregar e devolver
  if (file.exists(caminho_RDS)) {
    cat("✅ Arquivo RDS encontrado. Carregando dados...\n")
    dados <- readRDS(caminho_RDS)
    return(dados)
  }
  
  cat("📦 Arquivo RDS não encontrado. Extraindo dados dos arquivos .txt do GitHub...\n")
  
  # IDs dos neurônios que vamos buscar
  neuronios_ids <- c(1, 2, 3, 4, 7)
  
  # Base do link RAW correto
  base_url <- "https://raw.githubusercontent.com/diogomelgaco/neuronal_connectivity/main/LocusData/"
  
  # Função auxiliar para extrair números de um texto
  extract_numbers <- function(text) {
    as.numeric(unlist(regmatches(text, gregexpr("\\d+\\.?\\d*", text))))
  }
  
  # Inicializar lista dos neurônios
  neuronios <- list()
  
  # Loop para ler os 5 arquivos
  for (id in neuronios_ids) {
    file_url <- paste0(base_url, "locust20010217_spont_tetD_u", id, ".txt")
    
    file_lines <- tryCatch(
      readLines(url(file_url)),  # <- agora corretamente usando url(file_url)
      error = function(e) {
        stop(paste("Erro ao acessar o arquivo do neurônio", id, ":", e$message))
      }
    )
    
    numbers <- unlist(lapply(file_lines, extract_numbers))
    neuronios[[paste0("neuronio_", id)]] <- numbers
  }
  
  # Definir parâmetros de aquisição
  taxa_aquisicao <- 15000  # Hz
  janela_discretizacao <- 155  # Em unidades de amostras
  
  # Calcular número de colunas
  max_time <- max(unlist(neuronios))
  num_colunas <- ceiling(max_time / janela_discretizacao)
  cat("✅ Número de colunas para a matriz X_real:", num_colunas, "\n")
  
  # Criar a matriz x_real
  x_real <- matrix(0, nrow = length(neuronios_ids), ncol = num_colunas)
  
  for (i in 1:length(neuronios_ids)) {
    neuronio_nome <- paste0("neuronio_", neuronios_ids[i])
    spikes <- neuronios[[neuronio_nome]]
    
    for (spike in spikes) {
      coluna <- ceiling(spike / janela_discretizacao)
      if (coluna > 0 && coluna <= num_colunas) {
        x_real[i, coluna] <- 1
      }
    }
  }
  
  # Salvar os dados extraídos no RDS
  saveRDS(list(x_real = x_real,
               neuronios = neuronios,
               taxa_aquisicao = taxa_aquisicao,
               janela_discretizacao = janela_discretizacao), 
          file = caminho_RDS)
  
  cat("💾 Dados extraídos e salvos como", caminho_RDS, "\n")
  
  # Retornar os dados para o R
  return(list(
    x_real = x_real,
    neuronios = neuronios,
    taxa_aquisicao = taxa_aquisicao,
    janela_discretizacao = janela_discretizacao
  ))
}

#  . -----
######################################################################################################################## o
######################################################################################################################## o
#####  2: Descritivas do Banco de Dados ----

# Carregando os dados

dados_locus <- carregar_ou_extrair_LocusData()

# Acessar os componentes:
x_real <- dados_locus$x_real

ncol(x_real)
rowSums(x_real)

rowSums(x_real)/ncol(x_real)

plot_redeneuraldisparo_discreto(x_real[,1000:1500])

#  . -----
######################################################################################################################## o
######################################################################################################################## o
#####  3: Graficos de conexão ----


####  3.1: Sem diferenciar Excitatoria e Inibitório ----

plot_estimadores_v <- function(matrix_list, titulo_list) {
  n <- nrow(matrix_list[[1]])  # Assume que todas as matrizes têm o mesmo tamanho
  
  # Paleta de cores
  colors <- c("Connection" = "palegreen3", "No connection" = "white", "Inconclusivo" = "gray")
  
  # Função interna para criar um gráfico individual com base na matriz e no título
  create_plot <- function(Matrix, subtitle) {
    # Transformar a matriz em formato longo para ggplot
    matrix_data <- data.frame(
      Row = rep(1:n, each = n),
      Col = rep(1:n, times = n),
      Value = as.vector(t(Matrix))
    )
    
    # Converter valores para strings para mapeamento de cor
    matrix_data$Value <- as.character(matrix_data$Value)
    matrix_data$Value[matrix_data$Value == "0"] <- "No connection"
    matrix_data$Value[matrix_data$Value == "1"] <- "Connection"
    matrix_data$Value[is.na(matrix_data$Value)] <- "Inconclusivo"
    
    # Criar o gráfico SEM legenda
    p <- ggplot(matrix_data, aes(x = Col, y = Row, fill = Value)) +
      geom_tile(color = "black") +
      scale_fill_manual(
        values = colors
      ) +
      scale_x_continuous(breaks = 1:n, position = "top") +
      scale_y_reverse(breaks = 1:n) +
      labs(
        x = "Neurônios pós-sinápticos",
        y = "Neurônios pré-sinápticos",
        subtitle = subtitle
      ) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 10, margin = margin(t = 5)),
        axis.title.y = element_text(size = 10, margin = margin(r = 5)),
        legend.position = "none",  # REMOVE a legenda individual
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
  
  # --------- Criar a legenda única (baseada em todas as categorias) ----------
  
  # Criar um dataframe falso para garantir que todas as categorias apareçam
  legenda_data <- data.frame(
    Row = 1:3,
    Col = 1:3,
    Value = factor(c("Connection", "No connection", "Inconclusivo"), 
                   levels = c("Connection", "No connection", "Inconclusivo"))
  )
  
  legenda_plot <- ggplot(legenda_data, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "black") +   # Borda preta nas caixinhas
    scale_fill_manual(
      values = colors,
      labels = c("Conexão estimada", "Sem conexão estimada", "Inconclusivo")
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.7, "cm")
    )
  
  # Função para extrair a legenda
  g_legend <- function(a.gplot) {
    tmp <- ggplotGrob(a.gplot)
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    tmp$grobs[[leg]]
  }
  
  legenda_unica <- g_legend(legenda_plot)
  
  # --------- Ajuste dinâmico no número de colunas ------------
  
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
  
  # --------- Combinar gráficos + legenda ----------
  
  grid.arrange(
    arrangeGrob(grobs = plot_list, ncol = ncol_plot),
    legenda_unica,
    #top = textGrob("Aplicação do estimador de blocos à conectividade neuronal real em dados reais", 
                  # gp = gpar(fontsize = 14, fontface = "bold")),
    heights = c(10, 1.5)
  )
}

####  3.2: Diferenciando Excitatoria e Inibitório ----

plot_estimadores_v2 <- function(matrix_list, titulo_list = NULL) {
  n <- nrow(matrix_list[[1]])
  
  # Paleta de cores
  colors <- c(
    "Conexão excitatória estimada" = "skyblue",
    "Conexão inibitória estimada"  = "thistle",
    "Sem conexão estimada"         = "white",
    "Inconclusivo"                 = "gray"
  )
  
  # Função para criar gráfico
  create_plot <- function(Matrix, subtitle = "") {
    matrix_data <- data.frame(
      Row = rep(1:n, each = n),
      Col = rep(1:n, times = n),
      Value = as.vector(t(Matrix))
    )
    
    matrix_data$Condition <- with(matrix_data, ifelse(
      is.na(Value), "Inconclusivo",
      ifelse(Value == 1, "Conexão excitatória estimada",
             ifelse(Value == -1, "Conexão inibitória estimada", "Sem conexão estimada"))
    ))
    
    ggplot(matrix_data, aes(x = Col, y = Row, fill = Condition)) +
      geom_tile(color = "black") +
      scale_fill_manual(values = colors) +
      scale_x_continuous(breaks = 1:n, position = "top") +
      scale_y_reverse(breaks = 1:n) +
      labs(x = "Neurônios pós-sinápticos", y = "Neurônios pré-sinápticos", subtitle = subtitle) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 10, margin = margin(t = 5)),
        axis.title.y = element_text(size = 10, margin = margin(r = 5)),
        legend.position = "none",
        axis.text.x = element_text(margin = margin(t = 2)),
        axis.text.y = element_text(margin = margin(r = 2)),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.subtitle = element_text(hjust = 0.5, face = "bold", size = 10)
      )
  }
  
  if (is.null(titulo_list)) {
    titulo_list <- rep("", length(matrix_list))
  }
  
  plot_list <- mapply(create_plot, matrix_list, titulo_list, SIMPLIFY = FALSE)
  
  legenda_data <- data.frame(
    Row = 1:3, Col = 1:3,
    Condition = factor(c("Conexão excitatória estimada", 
                         "Conexão inibitória estimada", 
                         "Sem conexão estimada"),
                       levels = c("Conexão excitatória estimada", 
                                  "Conexão inibitória estimada", 
                                  "Sem conexão estimada"))
  )
  
  legenda_plot <- ggplot(legenda_data, aes(x = Col, y = Row, fill = Condition)) +
    geom_tile(color = "black") +
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.7, "cm")
    )
  
  g_legend <- function(a.gplot) {
    tmp <- ggplotGrob(a.gplot)
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    tmp$grobs[[leg]]
  }
  
  legenda_unica <- g_legend(legenda_plot)
  
  n_graficos <- length(plot_list)
  ncol_plot <- ifelse(n_graficos == 1, 1, ifelse(n_graficos == 2, 2, 3))
  
  gridExtra::grid.arrange(
    arrangeGrob(grobs = plot_list, ncol = ncol_plot),
    legenda_unica,
    heights = c(10, 1.5)
  )
}




#  . -----
######################################################################################################################## o
######################################################################################################################## o
#####  4: Aplicação dos Estimadores ----

#####  4.1: Aplicação do Estimador por Blocos  ----

### Base completa
est_bloco_real<-estimador_redeneural_bloco(x_real)

# Primeira metade da base
est_bloco_real_1metade<-estimador_redeneural_bloco(x_real[,1:137839])

# Segunda metade da base
est_bloco_real_2metade<-estimador_redeneural_bloco(x_real[,137840:275678])

### Gerar gráfico
matrix_list <- list(est_bloco_real$V_est, est_bloco_real_1metade$V_est, est_bloco_real_2metade$V_est)
titulo_list <- c("Base inteira", "1ª Metade da Base", "2ª Metade da Base")

plot_estimadores_v(matrix_list, titulo_list)


#####  4.2: Aplicação do Estimador por Pares  ----

u_real=media_u_sobre_mil_W(5)

### Base completa
est_pares_real <-estimador_redeneural_pares(x_real,u_real,,,,TRUE)

# Primeira metade da base
est_pares_real_1metade<-estimador_redeneural_pares(x_real[,1:137839],u_real)

# Segunda metade da base
est_pares_real_2metade<-estimador_redeneural_pares(x_real[,137840:275678],u_real)

### Gerar gráfico
matrix_list2 <- list(est_pares_real$V_est, est_pares_real_1metade$V_est, est_pares_real_2metade$V_est)
titulo_list2 <- c("Base inteira", "1ª Metade da Base", "2ª Metade da Base")

plot_estimadores_v2(matrix_list2, titulo_list2)


#####  4.2: Aplicação do Estimador por Pares Estendido  ----

### Base completa
est_pares2_real <-estimador_redeneural_pares_vetorizado_star_fast(x_real,u_real,,,,TRUE)

# Primeira metade da base
est_pares_real2_1metade<-estimador_redeneural_pares_vetorizado_star_fast(x_real[,1:137839],u_real)

# Segunda metade da base
est_pares_real2_2metade<-estimador_redeneural_pares_vetorizado_star_fast(x_real[,137840:275678],u_real,,,,TRUE)

### Gerar gráfico
matrix_list3 <- list(est_pares2_real$V_est, est_pares_real2_1metade$V_est, est_pares_real2_2metade$V_est)
titulo_list3 <- c("Base inteira", "1ª Metade da Base", "2ª Metade da Base")

plot_estimadores_v2(matrix_list3, titulo_list3)

#  . -----
######################################################################################################################## o
######################################################################################################################## o
#####  5: Conclusão ----

### Simulação Exemplo 4
set.seed(63436)
W4 <- gerador_matriz_W(15)
plot_matrix(W4)
mod4_3<-simulador_redeneural_discreto(15,1000000,W4)


## Método por Blocos
est_bloco4_3<-estimador_redeneural_bloco(mod4_3$X)

V_list4 <- list(mod4_3$V)
V_list_est_bloco4 <- list(est_bloco4_3$V_est)
n_list4 <- list(ncol(mod4_3$X))

plot_comparacao_V_V_est(V_list4,V_list_est_bloco4,n_list4)

## Método por Pares
est_pares4_3<-estimador_redeneural_pares(mod4_3$X,mod4_3$u)

V_list4_2 <- list(mod4_3$W)
V_list_est_par4 <- list(est_pares4_3$V_est)

plot_comparacao_V_V_est_excit_inib(V_list4_2,V_list_est_par4,n_list4)


### Simulação Exemplo 5
set.seed(33463)
W5 <- gerador_matriz_W(30)
plot_matrix(W5)
mod5_3<-simulador_redeneural_discreto(30,1000000,W5)

## Método por Blocos
est_bloco5_3<-estimador_redeneural_bloco(mod5_3$X)

V_list5 <- list(mod5_3$V)
V_list_est_bloco5 <- list(est_bloco5_3$V_est)
n_list5 <- list(ncol(mod5_3$X))

plot_comparacao_V_V_est(V_list5,V_list_est_bloco5,n_list5)

## Método por Pares
est_pares5_3<-estimador_redeneural_pares_vetorizado_star_fast(mod5_3$X,mod5_3$u,Delta = 10.33, verbose = TRUE)
V_est_pares5_3<-gerar_W_est(est_pares5_3$lambda,0.005,0.1)


V_list5_2 <- list(mod5_3$W)
V_list_est_par5 <- list(V_est_pares5_3)

plot_comparacao_V_V_est_excit_inib(V_list5_2,V_list_est_par5,n_list5)

## Método por Pares Corrigido
tune_limiar(est_pares5_3$lambda,mod5_3$W,csi1_seq,csi2_seq)

V_est_pares5_3_v2<-gerar_W_est(est_pares5_3$lambda,0.02,0.03)
V_list_est_par5_v2 <- list(V_est_pares5_3_v2)
plot_comparacao_V_V_est_excit_inib(V_list5_2,V_list_est_par5_v2,n_list5)
