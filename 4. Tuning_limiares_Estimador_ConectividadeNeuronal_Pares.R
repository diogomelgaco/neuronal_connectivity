######################################################################################################################## o
####  Tuning para Estimador de conectividade neuronal pelo método dos pares de neurônios em tempo discreto (Figueiredo 2025)
######################################################################################################################## o
####  Autor: Diogo Maia de Figueiredo
####  Orientadora: Dra. Aline Duarte
####  Programa de Pós Graduação do IME-USP

######################################################################################################################## o
######################################################################################################################## o
#####  1: Funções encontrar o melhor limiar (otimizado) ----

tune_limiar <- function(lambda, W_real, csi1_seq, csi2_seq) {
  melhor <- list(score = -Inf)
  
  for (csi1 in csi1_seq) {
    for (csi2 in csi2_seq) {
      W_pred <- matrix(0, nrow = nrow(lambda), ncol = ncol(lambda))
      W_pred[lambda >  csi2] <-  1
      W_pred[lambda < -csi1] <- -1
      
      acertos <- sum(W_pred == W_real)
      
      if (acertos > melhor$score) {
        melhor <- list(
          score = acertos,
          csi1 = csi1,
          csi2 = csi2,
          V_est = W_pred
        )
      }
    }
  }
  
  return(melhor)
}


######################################################################################################################## o
######################################################################################################################## o
#####  2: Etapa para encontrar o intervalo para percorrer o tunning ----

### 2.1: Função para encontrar o tunning  ----
rodar_experimentos <- function(n = 100, m = 10, t = 1000,
                               csi1_seq = seq(0.000, 0.2, by = 0.005),
                               csi2_seq = seq(0.000, 0.2, by = 0.005),
                               verbose = TRUE) {
  
  resultados <- data.frame(
    Amostra = integer(n),
    Semente = integer(n),
    Csi1 = numeric(n),
    Csi2 = numeric(n),
    Acerto = numeric(n),
    P_inib = numeric(n)
  )
  
  for (i in 1:n) {
    semente <- sample(1e5:1e6, 1)
    set.seed(semente)
    
    mod <- simulador_redeneural_discreto(m = m, n = t)
    est <- estimador_redeneural_pares(mod$X, mod$u)
    melhor <- tune_limiar(est$lambda, mod$W, csi1_seq, csi2_seq)
    
    # Máscara fora da diagonal
    off_diag <- lower.tri(mod$W, diag = FALSE) | upper.tri(mod$W, diag = FALSE)
    
    # Cálculo de acertos e total
    acertos <- sum(melhor$V_est[off_diag] == mod$W[off_diag])
    total <- m * (m - 1)
    
    # % de conexões inibitórias
    p_inib <- round(100 * sum(mod$W[off_diag] < 0) / sum(off_diag), 2)
    
    # Armazenar resultados
    resultados[i, ] <- c(i, semente, melhor$csi1, melhor$csi2, round(100 * acertos / total, 2), p_inib)
    
    # Exibir progresso
    if (verbose) {
      cat("✓ Amostra", i, "concluída:\n")
      print(resultados[i, , drop = FALSE])
      cat("--------\n")
    }
  }
  
  return(list(
    dados = resultados,
    parametros = list(n = n, m = m, t = t)
  ))
}


### 2.2: Função para rodar o tunning muitas vezes  ----
loop_simulacoes_divididas <- function(num_blocos = 10, n_por_bloco = 5, m = 10, t = 1000) {
  for (N in 1:num_blocos) {
    cat("🚀 Rodando bloco", N, "com", n_por_bloco, "simulações...\n")
    
    res <- rodar_experimentos(n = n_por_bloco, m = m, t = t, verbose = TRUE)
    
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- paste0("resultados_tunning_", N, "_", timestamp, ".RDS")
    
    saveRDS(res, file = filename)
    cat("💾 Bloco", N, "salvo em:", filename, "\n\n")
    
    rm(res)
    gc()
  }
  
  cat("✅ Todas as simulações foram concluídas e salvas com sucesso.\n")
}

### 2.3: Função para salvar o tunning  ----
juntar_resultados_rds <- function(padrao_arquivo = "^resultados_tunning_\\d+_.*\\.RDS$") {
  arquivos <- list.files(pattern = padrao_arquivo)
  
  if (length(arquivos) == 0) {
    stop("❌ Nenhum arquivo RDS encontrado com o padrão fornecido.")
  }
  
  lista_blocos <- list()
  
  for (arq in arquivos) {
    # Ler cada .RDS
    bloco <- readRDS(arq)
    
    # Pegar número do bloco a partir do nome do arquivo
    numero_bloco <- as.integer(sub("resultados_tunning_(\\d+)_.*", "\\1", arq))
    
    # Pegar somente o data.frame com os dados
    df <- bloco$dados
    
    # Adicionar coluna "Bloco"
    df$Bloco <- numero_bloco
    
    lista_blocos[[length(lista_blocos) + 1]] <- df
  }
  
  # Juntar tudo em uma única tabela
  resultado_final <- do.call(rbind, lista_blocos)
  
  # Reordenar colunas
  resultado_final <- resultado_final[, c("Bloco", "Amostra", "Semente", "Csi1", "Csi2", "Acerto", "P_inib")]
  
  # Ordenar por bloco e amostra
  resultado_final <- resultado_final[order(resultado_final$Bloco, resultado_final$Amostra), ]
  
  return(resultado_final)
}


### 2.4: Resultados   ----

#Rodando a função
loop_simulacoes_divididas(num_blocos = 10, n_por_bloco = 10, m = 5, t = 300000)
tabela_tunning <- juntar_resultados_rds()
write.csv(tabela_tunning, file = "base_tunning_melhor_limiar.csv", row.names = FALSE)

# Lendo o arquivo CSV
caminho_arquivo <- "C:/Users/maril/OneDrive/Área de Trabalho/Mestrado/Tunning/base_tunning_melhor_limiar.csv"
tabela_tunning <- read.csv(caminho_arquivo, encoding = "UTF-8")  # ou "Latin1" se tiver acentuação errada

     
######################################################################################################################## o
######################################################################################################################## o
#####  3: Etapa para encontrar o intervalo para percorrer o tunning ----

# 3.1 Função para rodar simulações em blocos com csi1 e csi2 fixos e salvar os resultados ----
simular_todos_os_limiares_dados_fixo <- function(num_blocos = 10, n_por_bloco = 10, m = 5, t = 300000) {
  
  # Definir as combinações de csi1 e csi2
  csi1_vals <- c(0.005,0.001)
  csi2_vals <- c(0.06,0.07, 0.08,0.09, 0.10, 0.11, 0.12)
  combinacoes <- expand.grid(csi1 = csi1_vals, csi2 = csi2_vals)
  
  for (bloco in 1:num_blocos) {
    cat("\n🚀 Iniciando bloco", bloco, "com", n_por_bloco, "simulações fixas\n")
    
    resultados_bloco <- data.frame()
    
    for (i in 1:n_por_bloco) {
      semente <- sample(1e5:1e6, 1)
      set.seed(semente)
      
      # Simular rede e estimar
      mod <- simulador_redeneural_discreto(m = m, n = t)
      W_real <- mod$W
      X <- mod$X
      u <- mod$u
      
      qtd_inib <- sum(W_real == -1)
      qtd_exc  <- sum(W_real == 1)
      qtd_conex <- qtd_inib + qtd_exc
      
      est <- estimador_redeneural_pares_vetorizado_star_fast(X, u)
      lambda <- est$theta
      
      off_diag <- lower.tri(W_real, diag = FALSE) | upper.tri(W_real, diag = FALSE)
      
      # Criar tabela temporária para esta amostra
      resultados_amostra <- data.frame()
      
      for (k in 1:nrow(combinacoes)) {
        csi1 <- combinacoes$csi1[k]
        csi2 <- combinacoes$csi2[k]
        
        W_pred <- matrix(0, nrow = m, ncol = m)
        W_pred[lambda >= csi2] <- 1
        W_pred[lambda <= -csi1] <- -1
        
        ac_inib <- if (qtd_inib > 0) {
          sum(W_real == -1 & W_pred == -1) / qtd_inib * 100
        } else NA
        
        ac_exc <- if (qtd_exc > 0) {
          sum(W_real == 1 & W_pred == 1) / qtd_exc * 100
        } else NA
        
        ac_conex <- if (qtd_conex > 0) {
          sum(W_pred[W_real != 0] == W_real[W_real != 0]) / qtd_conex * 100
        } else NA
        
        ac_total <- sum(W_pred[off_diag] == W_real[off_diag]) / sum(off_diag) * 100
        
        resultados_amostra <- rbind(resultados_amostra, data.frame(
          Bloco = bloco,
          Amostra = i,
          Semente = semente,
          Qtd_inib = qtd_inib,
          Qtd_exc = qtd_exc,
          Qtd_conex = qtd_conex,
          csi1 = csi1,
          csi2 = csi2,
          Acerto_inib = round(ac_inib, 2),
          Acerto_exc = round(ac_exc, 2),
          Acerto_conex = round(ac_conex, 2),
          Acerto_total = round(ac_total, 2)
        ))
      }
      
      # Adiciona os resultados da amostra ao bloco
      resultados_bloco <- rbind(resultados_bloco, resultados_amostra)
      
      # Print amigável da amostra
      cat(sprintf("\n📊 **Resultados %02d_%02d** | Semente: %d | Inib: %d | Exc: %d | Conex: %d\n",
                  bloco, i, semente, qtd_inib, qtd_exc, qtd_conex))
      cat(" Cs1   | Cs2   | AcInib | AcExc | AcConex | AcTotal\n")
      cat("-------|--------|--------|--------|----------|---------\n")
      
      for (j in 1:nrow(resultados_amostra)) {
        linha <- resultados_amostra[j, ]
        cat(sprintf(" %.3f | %.3f | %6s | %6s | %8s | %7.2f\n",
                    linha$csi1,
                    linha$csi2,
                    ifelse(is.na(linha$Acerto_inib), " NA", sprintf("%5.2f", linha$Acerto_inib)),
                    ifelse(is.na(linha$Acerto_exc),  " NA", sprintf("%5.2f", linha$Acerto_exc)),
                    ifelse(is.na(linha$Acerto_conex), " NA", sprintf("%7.2f", linha$Acerto_conex)),
                    linha$Acerto_total))
      }
      
      rm(mod, est)
      gc()
    }
    
    # Salvar o bloco
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    nome_arquivo <- paste0("resultados_dadosfixos_bloco", bloco, "_", timestamp, ".RDS")
    saveRDS(resultados_bloco, file = nome_arquivo)
    
    cat("\n💾 Bloco", bloco, "salvo em:", nome_arquivo, "\n")
  }
  
  cat("\n🏁 Todas as simulações foram finalizadas com sucesso!\n")
}


# 3.2 Função para armazenar e salvar os resultados ----
juntar_ultimos_blocos_dadosfixos <- function(qtd_blocos = 10, salvar_csv = TRUE, nome_base = "tabela_ultimos_blocos") {
  # Listar todos os arquivos que seguem o padrão
  arquivos <- list.files(pattern = "^resultados_dadosfixos_bloco\\d+_\\d{8}_\\d{6}\\.RDS$")
  
  if (length(arquivos) == 0) {
    stop("❌ Nenhum arquivo encontrado com padrão 'resultados_dadosfixos_bloco'.")
  }
  
  # Ordenar por timestamp decrescente e pegar os mais recentes
  arquivos_ordenados <- arquivos[order(sub("^.*_(\\d{8}_\\d{6})\\.RDS$", "\\1", arquivos), decreasing = TRUE)]
  arquivos_selecionados <- head(arquivos_ordenados, qtd_blocos)
  
  cat("📂 Lendo os", length(arquivos_selecionados), "últimos blocos salvos:\n")
  print(arquivos_selecionados)
  
  # Ler arquivos e juntar
  lista <- lapply(arquivos_selecionados, readRDS)
  tabela_junta <- do.call(rbind, lista)
  
  # Ordenar por Bloco e Amostra
  tabela_junta <- tabela_junta[order(tabela_junta$Bloco, tabela_junta$Amostra), ]
  
  # Salvar em .csv, se solicitado
  if (salvar_csv) {
    timestamp_atual <- format(Sys.time(), "%Y%m%d_%H%M%S")
    nome_csv <- paste0(nome_base, "_", timestamp_atual, ".csv")
    write.csv(tabela_junta, file = nome_csv, row.names = FALSE)
    cat("💾 Tabela salva como:", nome_csv, "\n")
  }
  
  return(tabela_junta)
}

### 3.3: Resultados   ----
simular_todos_os_limiares_dados_fixo(num_blocos = 20, n_por_bloco = 10, m = 5, t = 300000)
tabela <- juntar_ultimos_blocos_dadosfixos(qtd_blocos = 10)


# Preparar função auxiliar para criar os boxplots agrupados
criar_boxplot <- function(tabela, var_y, titulo = NULL, csi1_incluir = NULL, csi2_incluir = NULL) {
  tabela_filtrada <- tabela %>%
    { if (!is.null(csi1_incluir)) filter(., csi1 %in% csi1_incluir) else . } %>%
    { if (!is.null(csi2_incluir)) filter(., csi2 %in% csi2_incluir) else . } %>%
    mutate(
      csi1 = factor(csi1, levels = sort(unique(csi1))),
      csi2 = factor(csi2, levels = sort(unique(csi2))),
      grupo = paste0("xi[1] == ", csi1)
    )
  
  ggplot(tabela_filtrada, aes(x = csi2, y = .data[[var_y]], fill = csi2)) +
    geom_boxplot(outlier.size = 0.8, width = 0.6) +
    facet_wrap(~ grupo, nrow = 1, scales = "free_x", labeller = label_parsed) +
    scale_fill_brewer(palette = "Blues", direction = 1) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 10),
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      title = titulo,
      x = expression(xi[2]),
      y = "Percentual de Conexões identificadas (%)",
      fill = expression(xi[2])
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 0),
      panel.grid = element_blank(),
      panel.spacing = unit(1.2, "lines")  # espaço entre os grupos de xi[1]
    )
}

# Criar os 4 gráficos
plot1 <- criar_boxplot(tabela_tunning_2, "Acerto_conex", ,csi1_incluir = setdiff(unique(tabela_tunning_2$csi1), 0))
grid.arrange(plot1, ncol = 1)

# Garantir que as colunas são numéricas
tabela$csi1 <- as.numeric(tabela$csi1)
tabela$csi2 <- as.numeric(tabela$csi2)
tabela$taxa_fp <- as.numeric(tabela$Taxa_FP)
tabela$Acerto_total <- as.numeric(tabela$Acerto_total)

# Filtrar os valores desejados
tabela_filtrada <- tabela %>%
  filter(
    csi1 %in% c(0.005, 0.01),
    csi2 %in% c(0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12)
  )

# Calcular médias por grupo (sem multiplicar por 100, pois a base já está em %)
resumo <- tabela_filtrada %>%
  group_by(csi1, csi2) %>%
  summarise(
    media_acerto = mean(Acerto_total, na.rm = TRUE),
    media_fp = mean(taxa_fp, na.rm = TRUE),
    .groups = "drop"
  )

# Converter xi2 para fator ordenado
resumo$csi2 <- factor(resumo$csi2, levels = c(0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12))

# Paleta coerente para barras e linhas
cores <- c("0.005" = "#006400", "0.01" = "#66CDAA")

# Gráfico final
ggplot(resumo, aes(x = csi2)) +
  # Linha pontilhada de referência em y = 98%
  geom_hline(yintercept = 85, linetype = "dotted", color = "black", size = 0.4) +
  
  # Barras para taxa de falso positivo
  geom_col(aes(y = media_fp, fill = factor(csi1)),
           position = position_dodge(width = 0.6), width = 0.5) +
  
  # Linhas para acerto de conexão
  geom_line(aes(y = media_acerto, color = factor(csi1), group = csi1),
            size = 0.9) +
  
  # Pontos nas linhas
  geom_point(aes(y = media_acerto, color = factor(csi1), group = csi1),
             size = 2.2, position = position_dodge(width = 0.6)) +
  
  # Escalas manuais para manter cores idênticas
  scale_color_manual(
    values = cores,
    labels = c("0.005" = expression(xi[1] == 0.005), "0.01" = expression(xi[1] == 0.01))
  ) +
  scale_fill_manual(
    values = cores,
    guide = "none"
  ) +
  
  # Labels e tema
  scale_y_continuous(
    breaks = seq(0, 100, by = 10),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05)),
    name = "Porcentagem (%)"
  ) +
  labs(
    x = expression(xi[2]),
    color = expression(xi[1])
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_blank(),   # remove as linhas principais do fundo
    panel.grid.minor = element_blank(),   # remove as linhas menores do fundo
    axis.line = element_line(color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )
