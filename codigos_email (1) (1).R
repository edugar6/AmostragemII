#email1

# A PRAE irá extrair uma amostra de tamanho 1500 (na qual questionaremos esse número mais para frente) e te enviou o e-mail abaixo.

# Prezad@ xxxxxx,

# Estamos felizes em contar com sua consultoria para aprimorarmos nosso serviço oferecido. Temos algumas dúvidas sobre como essa pesquisa será realizada:

# Precisamos coletar dados de todos os seis restaurantes? Se não precisa, de quantos restaurantes serão necessários?

# Dos 1500 clientes que irão compor a amostra, quantos serão de cada restaurante?

tamanhos <- RU_censo2 %>%
  count(RU) %>%
  mutate(proporcao = n / sum(n),
         amostra_por_ru = round(proporcao * 1000))
amostra_AASc <- RU_censo2 %>%
  group_by(RU) %>%
  group_modify(~ {
    ru_atual <- .y$RU
    tamanho_amostra <- tamanhos$amostra_por_ru[tamanhos$RU == ru_atual]
    sample_n(.x, size = tamanho_amostra, replace = TRUE) # AASc: replace = TRUE
  }) %>%
  ungroup()
amostra_AASc

#email2

# Prezad@ xxxxxx,

# Nossa equipe gostou da ideia de coletar dados de todos os 6 restaurante.

# Precisamos então definir quantos clientes serão entrevistados em cada restaurante, pode nos auxiliar com essa questão?

# Em tempo, lembro que pretendemos coletar uma amostra total de 1500 alunos.

library(tidyverse)
library(ggplot2)

#população por RU: tamanho, média e variância
RU_censo2 <- RU_censo2 %>%
  mutate(
    nota_satisfacao = ifelse(nota_satisfacao == "cinquenta", "50", nota_satisfacao),
    nota_satisfacao = as.numeric(nota_satisfacao)
  )

et_pop <- RU_censo2 %>%
  group_by(RU) %>%
  summarise(
    N_h = n(),
    media_h = mean(nota_satisfacao, na.rm = TRUE),
    var_h = var(nota_satisfacao, na.rm = TRUE),
    .groups = "drop"
  )


#tamanho total da população
N <- sum(et_pop$N_h)
#número de estratos (RUs)
H <- n_distinct(RU_censo2$RU)
#tamanho total da amostra
n <- 1500

# --- Métodos de alocação ---
# AE uniforme
et_pop <- et_pop %>%
  mutate(n_uniforme = round(n / H))

# AE proporcional
et_pop <- et_pop %>%
  mutate(W_h = N_h / N,
         n_proporcional = round(W_h * n))

# AE de Neyman
denom <- sum(et_pop$N_h * sqrt(et_pop$var_h), na.rm = TRUE)
et_pop <- et_pop %>%
  mutate(
    n_neyman = round((N_h * sqrt(var_h) / denom) * n)
  )

#resultado final
alocacoes <- et_pop %>%
  select(RU, N_h, media_h, var_h, n_uniforme, n_proporcional, n_neyman)

alocacoes

# --- depois de calcular n_uniforme, n_proporcional e n_neyman ---

#tamanho total da amostra (já fixado antes)
n_total <- 1500 

#adicionando pesos relativos W_h
et_pop <- et_pop %>%
  mutate(
    W_h = N_h / N,
    #desvio-padrão do estrato
    sd_h = sqrt(var_h)
  )

#variância da AAS simples (referência)
var_aas <- sum(et_pop$W_h^2 * et_pop$var_h / (n_total / nrow(et_pop)))

# --- Funções para cada plano de AE ---

# --- Cálculo da variância da média global considerando os diferentes planos de amostragem ---

# 1) AAS simples: amostragem sem estratos, cada unidade tem a mesma chance.
#referência (baseline). 
var_aas <- sum(et_pop$W_h^2 * et_pop$var_h / (n_total / nrow(et_pop)))

# 2) AE uniforme: cada RU recebe o mesmo tamanho de amostra (n_h = n_total / H).
#    Vantagem: simples de implementar, mas não leva em conta os tamanhos reais dos RUs.
var_uniforme <- sum(et_pop$W_h^2 * et_pop$var_h / et_pop$n_uniforme)

# 3) AE proporcional: cada RU recebe uma amostra proporcional ao seu tamanho (N_h).
#    Vantagem: garante representatividade da população.
#    Resultado esperado: variância menor que AAS simples e uniforme, pois respeita N_h.
var_proporcional <- sum(et_pop$W_h^2 * et_pop$var_h / et_pop$n_proporcional)

# 4) AE de Neyman: cada RU recebe amostra proporcional a N_h * sigma_h.
#    Vantagem: minimiza a variância global (mais eficiente estatisticamente).
#    Resultado esperado: variância ainda menor que a proporcional, pois também considera a heterogeneidade dos RUs.
var_neyman <- sum(et_pop$W_h^2 * et_pop$var_h / et_pop$n_neyman)

# Comparação em tabela
variancias <- tibble(
  Metodo = c("AAS simples", "AE uniforme", "AE proporcional", "AE Neyman"),
  Variancia = c(var_aas, var_uniforme, var_proporcional, var_neyman)
)

variancias

# Eficiência relativa de cada método em relação ao AAS simples
eficiencia <- variancias %>%
  mutate(Eficiencia_relativa = variancias$Variancia[1] / Variancia)

eficiencia

# Eficiência relativa = Var(AAS simples) / Var(método)
# Valores > 1 indicam que o método é mais eficiente (menor variância) que o AAS simples.

#gráfico de análise exploratória
ggplot(RU_censo2, aes(x = RU, y = nota_satisfacao)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  theme_minimal() +
  labs(
    title = "Distribuição da Satisfação por Restaurante Universitário (RU)",
    x = "Restaurante Universitário (RU)",
    y = "Nota de Satisfação"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#email3

# Submetemos nosso estudo ao Comitê de Ética em Pesquisa da Universidade. O projeto voltou com algumas pendências, uma relacionada ao tamanho amostral:


# “É necessário justificar o tamanho amostral de 1500 clientes. Forneça uma justificativa ou faça o cálculo do tamanho amostral necessário com base em algum critério”

#variância populacional da nota
S2 <- var(RU_censo2$nota_satisfacao, na.rm = TRUE)

z <- qnorm(0.975)   #não caberia ao pesquisador definir isso, 
                    #mas vamos considerar 95% de confiança
                
E <- 2              #erro máximo tolerado (margem de erro de 2 pontos)

#tamanho necessário para AAS simples
n_aas <- (z^2 * S2) / (E^2)

#efeito do delineamento
DEFF_neyman <- var_neyman / var_aas

#ajustado pelo DEFF (plano de Neyman, com base na última atividade)
n_neyman_final <- n_aas * DEFF_neyman

#precisão real para n prévio de 1500
SE_neyman <- sqrt(var_neyman)          #erro padrão
E_real <- z * SE_neyman                #margem de erro (95%)

n_aas
DEFF_neyman
n_neyman_final
E_real

# Interpretação:
#n_aas ≈ 442 → seria suficiente para erro = ±2 pontos na média sob AAS simples.
#n_neyman_final ≈ 260 → sob AE Neyman, ainda menor é suficiente (mais eficiente).
#Mas n = 1500 dá margem de erro ≈ ±0.9 pontos → estimativa da satisfação muito precisa.
#Ou seja, se a média verdadeira for 73, o intervalo ficará entre 72,1 e 73,9, por exemplo.
