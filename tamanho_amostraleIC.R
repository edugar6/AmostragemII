
library("survey")
library("dplyr")
data(api)

dados <- apistrat %>%
  mutate(api99 = as.numeric(api99)) %>%
  filter(!is.na(api99))

# Estatísticas por estrato
et_pop <- dados %>%
  group_by(stype) %>%
  summarise(
    N_h = n(),
    mu_h = mean(api99),
    # variância populacional (corrigida para N, não para n-1)
    var_pop_h = var(api99) * (n() - 1) / n(),
    .groups = "drop"
  )

# Estatísticas globais com pesos (estratos)
N <- sum(et_pop$N_h)   # população total

et_pop <- et_pop %>%
  mutate(
    sd_h = sqrt(var_pop_h),
    W_h = N_h / N,          # peso relativo
    total_h = N_h * mu_h    # total do estrato
  )

print(et_pop)

# Estatísticas globais
N <- sum(et_pop$N_h)   # população total
mu_pop <- sum((et_pop$N_h / N) * et_pop$mu_h)  # média ponderada global

# decomposição da variância
sigma2_d <- sum((et_pop$N_h / N) * et_pop$var_pop_h)                   # dentro
sigma2_e <- sum((et_pop$N_h / N) * (et_pop$mu_h - mu_pop)^2)           # entre
sigma2_total <- sigma2_d + sigma2_e                                    # total

# tabela global
globais <- tibble(
  N_total = N,
  mu_pop = mu_pop,
  sigma2_dentro = sigma2_d,
  sigma2_entre = sigma2_e,
  sigma2_total = sigma2_total
)

print(globais)

# -------------------
# Tamanho amostral
# -------------------
n_total <- 500  # escolha do tamanho da amostra

# -------------------
# Alocações
# -------------------
# uniforme
et_pop <- et_pop %>%
  mutate(n_uniforme = round(n_total / H))

# proporcional
et_pop <- et_pop %>%
  mutate(n_proporcional = round(W_h * n_total))

# Neyman
denom <- sum(et_pop$N_h * et_pop$sd_h, na.rm = TRUE)
et_pop <- et_pop %>%
  mutate(n_neyman = round((N_h * sd_h / denom) * n_total))

# Funções de variância (dentro de cada estrato)

# -------------------
# Variâncias do estimador da média (Variância do estimador da média estratificada)
# Fórmula: soma dos pesos² * var(h)/n_h
# -------------------
var_mean_AASc <- function(n_h) {
  sum((et_pop$W_h^2) * (et_pop$var_pop_h / n_h), na.rm = TRUE)
}
var_mean_AASs <- function(n_h) {
  sum((et_pop$W_h^2) * ((1 - n_h/et_pop$N_h) * (et_pop$var_pop_h / n_h)), na.rm = TRUE)
}

# -------------------
# Variâncias do estimador do total estratificado
# -------------------
var_total_AASc <- function(n_h) {
  sum((et_pop$N_h^2) * (et_pop$var_pop_h / n_h), na.rm = TRUE)
}
var_total_AASs <- function(n_h) {
  sum((et_pop$N_h^2) * ((1 - n_h/et_pop$N_h) * (et_pop$var_pop_h / n_h)), na.rm = TRUE)
}

# -------------------
# Variâncias do estimador da proporção estraficado
# -------------------
# Criar variável indicadora
dados <- dados %>%
  mutate(ybin = ifelse(api99 > 700, 1, 0))

# Proporção por estrato
prop_h <- dados %>%
  group_by(stype) %>%
  summarise(p_h = mean(ybin), .groups = "drop") %>%
  pull(p_h)

# Funções de variância da proporção
var_prop_AASc <- function(n_h, p_h) {
  sum((et_pop$W_h^2) * (p_h * (1 - p_h) / n_h), na.rm = TRUE)
}
var_prop_AASs <- function(n_h, p_h) {
  sum((et_pop$W_h^2) * ((1 - n_h/et_pop$N_h) * (p_h * (1 - p_h) / n_h)), na.rm = TRUE)
}

# -------------------
# Caso de referência: AASc e AASs sem estratos
# -------------------
# ⚠️ Aqui simulamos "como se fosse AAS simples" aplicando n_h iguais
# em cada estrato só para ter um comparativo de eficiência.

#Então o resultado não é a variância do AASs de verdade, porque o AAS simples não estratifica.
#O que você está calculando aqui é a variância da média estratificada sob uma alocação artificial uniforme
#(mesmo n_h para todos os estratos)

n_h_ref <- rep(n_total / H, H)

# Média
var_AASc_mean <- var_mean_AASc(n_h_ref)
var_AASs_mean <- var_mean_AASs(n_h_ref)

# Total
var_AASc_total <- var_total_AASc(n_h_ref)
var_AASs_total <- var_total_AASs(n_h_ref)

# Proporção
var_AASc_prop <- var_prop_AASc(n_h_ref, prop_h)
var_AASs_prop <- var_prop_AASs(n_h_ref, prop_h)

# -------------------
# Variâncias por alocação estratificada
# -------------------
# Média
var_uniforme_mean     <- var_mean_AASc(et_pop$n_uniforme)
var_proporcional_mean <- var_mean_AASc(et_pop$n_proporcional)
var_neyman_mean       <- var_mean_AASc(et_pop$n_neyman)

# Total
var_uniforme_total     <- var_total_AASc(et_pop$n_uniforme)
var_proporcional_total <- var_total_AASc(et_pop$n_proporcional)
var_neyman_total       <- var_total_AASc(et_pop$n_neyman)

# Proporção
var_uniforme_prop     <- var_prop_AASc(et_pop$n_uniforme, prop_h)
var_proporcional_prop <- var_prop_AASc(et_pop$n_proporcional, prop_h)
var_neyman_prop       <- var_prop_AASc(et_pop$n_neyman, prop_h)

# -------------------
# --- ADIÇÕES: calcular Vc, Vpr, Vun, Vney e efeitos relativos
# -------------------

# Vc: variância AASc "pura" (sem estratos) = sigma^2 / n
Vc <- sigma2_total / n_total

# Vpr: AE proporcional (com alocação n_h = n * W_h) -> Vpr = sigma_dentro^2 / n
# onde sigma2_d é a variância dentro (decomposição já calculada)
Vpr <- sigma2_d / n_total

# Vun: AE uniforme (n_h = n / H)
k <- n_total / H
Vun <- sum((et_pop$W_h^2) * (et_pop$var_pop_h / k), na.rm = TRUE)
# (equivalente a var_mean_AASc(rep(k,H)))

# Vney: AE Neyman ótima (V = 1/n * (sum_h W_h * sd_h)^2 )
Vney <- (1 / n_total) * (sum(et_pop$W_h * et_pop$sd_h))^2

# Efeitos em relação a AASc (redução percentual) e razão de eficiência
efeito_Vpr_pct  <- (Vc - Vpr) / Vc * 100
efeito_Vun_pct  <- (Vc - Vun) / Vc * 100
efeito_Vney_pct <- (Vc - Vney) / Vc * 100

# Razão (Vc / Vplan) — se >1, plano é mais eficiente que AASc
razao_Vpr  <- Vc / Vpr
razao_Vun  <- Vc / Vun
razao_Vney <- Vc / Vney

# Qual é o melhor (menor variância)
variancias_list <- c(AASc = Vc, Vpr = Vpr, Vun = Vun, Vney = Vney)
melhor_plano <- names(which.min(variancias_list))

# -------------------
# Comparação em tabelas
# -------------------
variancias_mean <- tibble(
  Metodo = c("AASc", "AASs", "AE Uniforme", "AE Proporcional", "AE Neyman"),
  Variancia = c(var_AASc_mean, var_AASs_mean, var_uniforme_mean,
                var_proporcional_mean, var_neyman_mean)
) %>%
  mutate(Eficiencia_relativa = Variancia[2] / Variancia) # usando AASs como referência

variancias_total <- tibble(
  Metodo = c("AASc", "AASs", "AE Uniforme", "AE Proporcional", "AE Neyman"),
  Variancia = c(var_AASc_total, var_AASs_total, var_uniforme_total,
                var_proporcional_total, var_neyman_total)
) %>%
  mutate(Eficiencia_relativa = Variancia[2] / Variancia)

variancias_prop <- tibble(
  Metodo = c("AASc", "AASs", "AE Uniforme", "AE Proporcional", "AE Neyman"),
  Variancia = c(var_AASc_prop, var_AASs_prop, var_uniforme_prop,
                var_proporcional_prop, var_neyman_prop)
) %>%
  mutate(Eficiencia_relativa = Variancia[2] / Variancia) # usando AASs como referência

# -------------------
# Resultados
# -------------------
cat("\n--- Variâncias do Estimador da Média ---\n")
print(variancias_mean)

cat("\n--- Variâncias do Estimador do Total ---\n")
print(variancias_total)

cat("\n--- Variâncias do Estimador da Proporção ---\n")
print(variancias_prop)

cat("\n--- Valores pedidos (n = ", n_total, ") ---\n", sep = "")
cat(sprintf("Vc (AASc, sigma^2 / n)        = %g\n", Vc))
cat(sprintf("Vpr (AE proporcional)         = %g\n", Vpr))
cat(sprintf("Vun (AE uniforme)             = %g\n", Vun))
cat(sprintf("Vney (AE Neyman ótima)        = %g\n", Vney))

cat("\n--- Efeitos em relação a AASc (redução percentual) ---\n")
cat(sprintf("Proporcional: %.2f %%\n", efeito_Vpr_pct))
cat(sprintf("Uniforme:     %.2f %%\n", efeito_Vun_pct))
cat(sprintf("Neyman:       %.2f %%\n", efeito_Vney_pct))

cat("\n--- Razões de eficiência (Vc / Vplan) ---\n")
cat(sprintf("Vc/Vpr  = %.4f\n", razao_Vpr))
cat(sprintf("Vc/Vun  = %.4f\n", razao_Vun))
cat(sprintf("Vc/Vney = %.4f\n", razao_Vney))

cat("\nMelhor plano (menor variância):", melhor_plano, "\n\n")

print(et_pop)

# -----------------------------
# --- EXERCÍCIO FINAL (apenas o necessário) ---
# Usa: et_pop, mu_pop, sigma2_total, N, H
# -----------------------------

# -----------------------------
# Intervalo de Confiança com n calculado pelo AAS simples
# -----------------------------

# Parâmetros globais
conf <- 0.95
z <- qnorm((1 + conf)/2)
E_choice <- 20   # margem de erro desejada (ex: ±30 pontos)

# Variância populacional total (já calculada no script)
S2_pop <- sigma2_total

# n requerido pelo AAS simples (com FPC)
n_aas_req <- {
  n0 <- (z^2 * S2_pop) / (E_choice^2)
  ceiling((N * n0) / (N + n0 - 1))   # com FPC
}

# variância do estimador da média sob AAS simples
var_AAS <- (1 - n_aas_req/N) * S2_pop / n_aas_req

# intervalo de confiança
ci_AAS <- c(mu_pop - z*sqrt(var_AAS), mu_pop + z*sqrt(var_AAS))

cat("\n--- Intervalo de Confiança (AAS simples, n calculado) ---\n")
print(tibble::tibble(
  Plano = "AAS simples (com FPC)",
  n = n_aas_req,
  Var = var_AAS,
  SE = sqrt(var_AAS),
  Lower = ci_AAS[1],
  Upper = ci_AAS[2]
))

#partindo para calculo tamanho amostral
# garantir H
H <- nrow(et_pop)

# Parâmetros de interesse
conf <- 0.95
z <- qnorm((1 + conf) / 2)
E <- 30    # margem de erro desejada (ajuste aqui se quiser)

# Variância populacional total
S2_pop <- sigma2_total

# 1) n sob AAS simples (sem FPC; pode aplicar depois se quiser)
n_aas <- (z^2 * S2_pop) / (E^2)
n_aas <- ceiling(n_aas)

# Referência
n_ref <- n_aas
Var_ref <- S2_pop / n_ref

# --- Variâncias teóricas ---
Vpr <- sigma2_d / n_ref

k <- n_ref / H
# garantir que k >= 1 para evitar divisão por zero
k <- max(1, k)
Vun <- sum((et_pop$W_h^2) * (et_pop$var_pop_h / k), na.rm = TRUE)

Vney <- (1 / n_ref) * (sum(et_pop$W_h * et_pop$sd_h))^2

# --- DEFFs teóricos ---
DEFF_proporcional <- Vpr / Var_ref
DEFF_uniforme    <- Vun / Var_ref
DEFF_neyman      <- Vney / Var_ref

# --- Ajuste dos tamanhos ---
n_proporcional_adj <- ceiling(n_aas * DEFF_proporcional)
n_uniforme_adj     <- ceiling(n_aas * DEFF_uniforme)
n_neyman_adj       <- ceiling(n_aas * DEFF_neyman)

# --- DEFFs empíricos (com n_total definido antes no script) ---
Var_ref_at_n_total <- S2_pop / n_total
DEFF_prop_emp <- var_proporcional_mean / Var_ref_at_n_total
DEFF_unif_emp <- var_uniforme_mean / Var_ref_at_n_total
DEFF_neym_emp <- var_neyman_mean / Var_ref_at_n_total

n_proporcional_adj_emp <- ceiling(n_aas * DEFF_prop_emp)
n_uniforme_adj_emp     <- ceiling(n_aas * DEFF_unif_emp)
n_neyman_adj_emp       <- ceiling(n_aas * DEFF_neym_emp)

# --- Resultados resumidos ---
cat("\n--- Resumo: cálculo de n ajustado pelo DEFF ---\n")
cat(sprintf("Parâmetros: conf = %.2f, E = %g\n", conf, E))
cat(sprintf("S2_pop (var total) = %g   => n_aas (sem FPC) = %d\n\n", S2_pop, n_aas))

res_deff <- tibble::tibble(
  Metodo = c("AE proporcional (teo)", "AE uniforme (teo)", "AE Neyman (teo)"),
  Var_plan = c(Vpr, Vun, Vney),
  DEFF = c(DEFF_proporcional, DEFF_uniforme, DEFF_neyman),
  n_ajustado = c(n_proporcional_adj, n_uniforme_adj, n_neyman_adj)
)
print(res_deff)

cat("\n--- DEFFs empíricos (para n_total usado anteriormente) e n ajustados (emp) ---\n")
res_emp <- tibble::tibble(
  Metodo = c("AE proporcional (emp)", "AE uniforme (emp)", "AE Neyman (emp)"),
  Var_plan_emp = c(var_proporcional_mean, var_uniforme_mean, var_neyman_mean),
  DEFF_emp = c(DEFF_prop_emp, DEFF_unif_emp, DEFF_neym_emp),
  n_ajust_emp = c(n_proporcional_adj_emp, n_uniforme_adj_emp, n_neyman_adj_emp)
)
print(res_emp)

# -----------------------------
# Função geral pedida pela questão
# -----------------------------
n_amostral <- function(E, conf = 0.95, 
                       S2,          # variância total populacional
                       Ns,          # tamanhos dos estratos
                       variancias,  # variâncias por estrato
                       sds,         # desvios-padrão por estrato
                       use_fpc = TRUE) {
  
  Npop <- sum(Ns)
  W <- Ns / Npop
  z <- qnorm((1 + conf) / 2)
  
  # AAS simples
  n0 <- (z^2 * S2) / (E^2)
  if (use_fpc) {
    n_aas <- ceiling((Npop * n0) / (Npop + n0 - 1))
  } else {
    n_aas <- ceiling(n0)
  }
  Var_ref <- S2 / n_aas
  
  # Proporcional
  Vpr <- sum(W * variancias) / n_aas
  DEFF_proporcional <- Vpr / Var_ref
  n_proporcional <- ceiling(n_aas * DEFF_proporcional)
  
  # Uniforme
  H <- length(Ns)
  k <- max(1, n_aas / H)  # evitar k=0
  Vun <- sum(W^2 * (variancias / k))
  DEFF_uniforme <- Vun / Var_ref
  n_uniforme <- ceiling(n_aas * DEFF_uniforme)
  
  # Neyman
  Vney <- (1 / n_aas) * (sum(W * sds))^2
  DEFF_neyman <- Vney / Var_ref
  n_neyman <- ceiling(n_aas * DEFF_neyman)
  
  tibble::tibble(
    Metodo = c("AAS simples", "AE proporcional", "AE uniforme", "AE Neyman"),
    n = c(n_aas, n_proporcional, n_uniforme, n_neyman),
    DEFF = c(1, DEFF_proporcional, DEFF_uniforme, DEFF_neyman)
  )
}

# Exemplo de uso
res_final <- n_amostral(E = 10, conf = 0.95, 
                        S2 = sigma2_total,
                        Ns = et_pop$N_h,
                        variancias = et_pop$var_pop_h,
                        sds = et_pop$sd_h)

print(res_final)
