##### Extracao e discretizacao do banco de dados real
# Definir o diretório onde os arquivos .txt estão localizados
dir_path <- "C:/Users/User/Desktop/Mestrado/LocusData"

# Verificar se o diretorio existe
if (!dir.exists(dir_path)) {
  stop("O diretório especificado não existe.")
}

# Listar todos os arquivos .txt no diretorio que seguem o padrao especifico
file_list <- list.files(path = dir_path, pattern = "locust20010217_spont_tetD_u[12347]\\.txt$", full.names = TRUE)

# Verificar se ha arquivos correspondentes
if (length(file_list) == 0) {
  stop("Nenhum arquivo correspondente foi encontrado no diretório.")
}

# Funcao para extrair numeros de uma string
extract_numbers <- function(text) {
  as.numeric(unlist(regmatches(text, gregexpr("\\d+\\.?\\d*", text))))
}

# Inicializar uma lista para armazenar os vetores
neuronios <- list()

# Ler os arquivos .txt e extrair os números
for (file_path in file_list) {
  # Obter o numero do neuronio com base no nome do arquivo
  neuronio_num <- gsub(".*_u(\\d+)\\.txt$", "\\1", basename(file_path))
  
  # Ler o arquivo linha por linha
  file_lines <- readLines(file_path)
  
  # Extrair todos os numeros de cada linha e concatenar em um vetor
  numbers <- unlist(lapply(file_lines, extract_numbers))
  
  # Tambem armazenar no vetor neuronios para facil acesso
  neuronios[[paste0("neuronio_", neuronio_num)]] <- numbers
}

# Definir a taxa de aquisicao (Hz), a mesma usada no codigo Python
taxa_aquisicao <- 15000  # Hz

# Definir a janela de discretizacao em amostras (conforme codigo Python)
janela_discretizacao <- 155  # Em unidades de amostras, nao em segundos

# Encontrar o tempo maximo para definir o numero de colunas da matriz X
max_time <- max(unlist(neuronios))

# Definir o numero de colunas baseado no tempo maximo e a janela de discretizacao
# Aqui usamos diretamente a janela em amostras
num_colunas <- ceiling(max_time / janela_discretizacao)
print(paste("Número de colunas:", num_colunas))

# Verificar se o número de colunas e valido
if (is.na(num_colunas) || num_colunas <= 0) {
  stop("Número de colunas inválido. Verifique os dados.")
}

# Criar a matriz X com 5 linhas (um para cada neurônio) e 'num_colunas' colunas
x_real <- matrix(0, nrow = 5, ncol = num_colunas)

# Preencher a matriz X
for (i in 1:5) {
  neuronio_nome <- paste0("neuronio_", c(1, 2, 3, 4, 7)[i])  # Acessar os vetores de neuronios
  spikes <- neuronios[[neuronio_nome]]  # Obter os disparos do neurônio i
  
  # Iterar sobre cada disparo do neuronio
  for (spike in spikes) {
    # Encontrar a coluna correspondente ao intervalo de tempo para o disparo
    coluna <- ceiling(spike / janela_discretizacao)
    
    # Verificar se a coluna é válida antes de marcar o valor na matriz
    if (coluna > 0 && coluna <= num_colunas) {
      # Marcar com 1 a posição correspondente na matriz X
      x_real[i, coluna] <- 1
    }
  }
}

## Estimador pelo método dos blocos
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
