# Estimador de conectividade neuronal pelo metodo de blocos
estimador_redeneural_bloco <- function(X, xi = 0.001, epsilon = 0.05) {
  # X é uma matriz m x n composta por 1 e 0, onde 1 indica que o neurônio i disparou no instante t.
  
  # Definir número de linhas (neurônios) e colunas (instantes)
  n <- ncol(X)       
  m <- nrow(X)
  
  # Listas para armazenar submatrizes tau e probabilidades associadas para cada neurônio
  tau_list <- list()  
  tau_p_list <- list()  
  delta <- matrix(NA, nrow = m, ncol = m)
  diag(delta) <- 0
  delta_total <- matrix(NA, nrow = m, ncol = m)
  diag(delta_total) <- 0
  
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
      print(tau_i)
      print(tau_n)
      
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
    
    # Calcular delta para todos os pares j -> i
    for (j in setdiff(1:m, i)) {
      aux_delta <- numeric(0)
      
      if (length(tau_list[[i]]) > 1) {
        # Comparar submatrizes de tau entre si
        for (w in 1:(length(tau_list[[i]]) - 1)) {
          tau_w <- tau_list[[i]][[w]][-j, ]
          
          list_p_w_v <- sapply((w + 1):length(tau_list[[i]]), function(v) {
            tau_v <- tau_list[[i]][[v]][-j, ]
            if (identical(tau_w, tau_v)) return(tau_p_list[[i]][v]) else return(NA)
          })
          
          aux_delta_tau <- abs(list_p_w_v - tau_p_list[[i]][w])
          aux_delta <- c(aux_delta, max(aux_delta_tau, na.rm = TRUE))
        }
        
        # Atualizar delta para j -> i
        delta[j, i] <- if (length(aux_delta) > 0) max(aux_delta, na.rm = TRUE) else NA
      }
      
      delta_total[j, i] <- delta[j, i]
    }
  }
  
  # Construir matriz de conectividade estimada V_est
  V_est <- matrix(NA, nrow = m, ncol = m)
  V_est[!is.na(delta_total) & delta_total > epsilon] <- 1
  V_est[!is.na(delta_total) & delta_total <= epsilon] <- 0
  
  # Retornar resultados
  return(list(
    X = X,
    xi = xi,
    epsilon = epsilon,
    tau_list = tau_list,
    tau_p_list = tau_p_list,
    delta = delta,
    delta_total = delta_total,
    V_est = V_est
  ))
}



## Gerando o Exemplo 1:
exemplo1 <- matrix(0, nrow = 5, ncol = 5)

# Definição das posições específicas com valor 1
posicoes <- list(
  c(1, 2), c(1, 3), c(2, 3), c(2, 4),
  c(3, 1), c(3, 5), c(4, 5), c(5, 4)
)

# Atribuindo o valor 1 nas posições especificadas
for (pos in posicoes) {
  exemplo1[pos[1], pos[2]] <- 1
}


# Rodando os estimadores
set.seed(2512)
mod1<-simulador_redeneural_discreto(5,10000,exemplo1)
mod2<-simulador_redeneural_discreto(5,100000,exemplo1)
mod3<-simulador_redeneural_discreto(5,1000000,exemplo1)

exemplo1$W
est1<-estimador_redeneural_bloco(mod1$X)
est1$V
est1$delta_total
est2<-estimador_redeneural_bloco(mod2$X)
est2$V
est2$delta_total
est3<-estimador_redeneural_bloco(mod3$X)
est3$V
est3$delta_total


mod4<-simulador_redeneural_discreto(15,10000,exemplo2$W)
mod5<-simulador_redeneural_discreto(15,100000,exemplo2$W)
mod6<-simulador_redeneural_discreto(15,1000000,exemplo2$W)

exemplo2$W
est4<-estimador_redeneural_bloco(mod4$X)
est4$V
est4$delta_total
est5<-estimador_redeneural_bloco(mod5$X)
est5$V
est5$delta_total
est6<-estimador_redeneural_bloco(mod6$X)
est6$V
est6$delta_total

#dados reais
est_basereal<-estimador_redeneural_bloco(x_real)
est_basereal$V_est
est_basereal$delta_total

# primeira metade
est_basereal1<-estimador_redeneural_bloco(x_real[,1:137839])
est_basereal1$V_est
est_basereal1$delta_total

# segunda metade
est_basereal2<-estimador_redeneural_bloco(x_real[,137840:275678])
est_basereal2$V_est
est_basereal2$delta_total
