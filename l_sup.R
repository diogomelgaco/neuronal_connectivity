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