#Instacao de Pacotes Necessarios ----

install.packages("ggplot2")
install.packages("igraph")
install.packages("reshape2")

require(ggplot2)
require(igraph)
require(reshape2)


simulador_redeneural_discreto<-function(m,n,W=matrix(0, m, m)){  # m=neuronios na rede, n= instantes gerados
  
  # Declaracao de Parametros
  
    # 1. Matrizes W= peso sinaptico e V= conectividade neuronal
      gerador_matriz_W <- function(m) {
        
        # Probabilidade de conexao entre um par de neuronios
        prob_conexao <- 0.2
  
        # Inicializar a matriz
        W <- matrix(0, m, m)
        repeat {
          for (i in 1:m) {
            for (j in 1:m) {
              if (i != j) { # Ignorar a diagonal principal
                if (runif(1) < prob_conexao) {
                  W[i, j] <- 1
                }
              }
            }
          }
    
          # Validacao para todo neuronio da rede tenha conexao
          linhas_com_zeros <- apply(W, 1, function(row) all(row == 0))
          colunas_com_zeros <- apply(W, 2, function(col) all(col == 0))
    
            # Se nenhuma linha ou coluna tiver so zeros, saia do loop
            if (!any(linhas_com_zeros) && !any(colunas_com_zeros)) {
              break
            }
        }
        
        return(W)
      }

      # Gerar a matriz W caso nao tenha sido um input
      if (is.matrix(W) && all(W == 0)){
      W <- gerador_matriz_W(m)
      }
      V <- ifelse(W == 0, 0, 1)
      
    # 2. Funcao de vazamento g
      g <- rep(0.8, m)
  
    # 3. Funcao Phi
      u<-floor(mean(rowSums(W)))
      
      # Funcao sigmoide ajustada
      sigmoid <- function(x, b, c) {
        1 / (1 + exp(-b * (x - c)))
      }
      
      # Funcao que calcula a probabilidade baseada na funcao sigmoide ajustada
      phi<- function(x, u) {
        
        # Condicoes
        p0 <- 0.02  # f(U < 0) = 0,02
        pu <- 0.8   # f(U < u) = 0,8
    
        # Calculo dos parametros b e c
        b <- log(pu / (1 - pu) / (p0 / (1 - p0))) / u
        c <- -log(p0 / (1 - p0)) / b
        
        # Retorna a probabilidade para o valor x dado
        return(sigmoid(x, b, c))
      }
  
    # 4. Estado inicial U_0  
      
      # Funcao para gerar o vetor inicial U_0
      gerador_U0 <- function(m, u) {
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


