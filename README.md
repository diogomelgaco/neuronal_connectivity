# neuronal_connectivity
#### Maio-2025

Os códigos contidos nessa parte são partes da dissertação *"Comparação de estimadores de conectividade neuronal a partir de sequências de disparos discretizados"* elaborada por **Diogo Maia de Figueiredo** no programa de pós-graduação da USP do departamento de Matemática e Estatística do IME, sob a orientação da prof. **Dra. Aline Duarte de Oliveira**.

Todos os códigos foram desenvolvidos utilizando o R 4.4.2 e devem ser consultados na seguinte ordem:

1. *Simulador_SistemaNeuronal_discreto.R*

   Código para gerar um sistema neuronal tendo como entrada a quantidade de neurônios e a quantidade de instantes em que se deseja observar o comportamento de disparos destes neurônios. A simulação segue o modelo descrito por Galves (2013) em "Infinite systems of interacting chains with memory of variable length - a stochastic model for biological neural nets". Nesse simulador é possível adaptar os parâmetros W, de conectividade neuronal, a função de disparo e a função de perda. Além do simulador, há função para calcular o último instante de disparo para cada neurônios, geração dos exemplos de rede e todas as representações gráficas e visuais apresentadas no Capítulo 2 na dissertação.
   Todos os pacotes necessários para a execução de todos os códigos abaixo estão nesse arquivo.

3. *Estimador_ConectividadeNeuronal_Blocos.R*
   
   Algoritmo para estimar a conectividade entre neurônios com base em uma sequência de disparos baseado no método proposto por Duarte (2019) em "Estimating the interaction graph of stochastic neural dynamics". Além de todos as representações gráficas apresentadas no capítulo 3 da dissertação.

4. *Estimador_ConectividadeNeuronal_Pares.R*
   
   Algoritmo para estimar a conectividade entre neurônios de forma adaptada ao descrito em Santis (2022) em "Estimating the interaction graph of stochastic neuronal dynamics by observing only pairs of neurons", considerando agora uma sequência de disparos discretizados. Além de todos as representações gráficas apresentadas no capítulo 4 da dissertação.

5. *Estimador_ConectividadeNeuronal_ParesExtendido.R*
   
    Algoritmo para estimar a conectividade entre neurônios considerando pares de neurônios, porém, com uma extensão das janelas de tempo definidas no estimador de Santis (2022).
   
6. *Aplicacao_DadosReais.R*
   
   Extração dos dados reais utilizados em  "Estimation of neuronal interaction graph from spike train data", por Brochini (2017), seguidos os mesmos passos de discretização e em seguida, a aplicação dos estimadores por Blocos, Pares e Pares Extendido, visando identificar conexões entre os neurônios com base em uma sequência de disparos real.
