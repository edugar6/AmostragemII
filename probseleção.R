# -------------------
# Amostra de tamanho total
# -------------------
n_total <- 1500

# Alocações (com mínimo de 1 por estrato)
et_pop <- et_pop %>%
  mutate(
    n_uniforme     = pmax(1, round(n_total / H)),
    n_proporcional = pmax(1, round(W_h * n_total))
  )

# Neyman
denom <- sum(et_pop$N_h * et_pop$sd_h, na.rm = TRUE)
et_pop <- et_pop %>%
  mutate(
    n_neyman = pmax(1, round((N_h * sd_h / denom) * n_total))
  )

# -------------------
# Probabilidades de seleção
# -------------------
# AASc (sem estratos, população toda)
pi_AASc <- 1 - (1 - 1/N)^n_total

# AASs (sem estratos, população toda)
pi_AASs <- n_total / N

# Funções auxiliares para AE
prob_AE_AASc <- function(n_h, N_h) {
  1 - (1 - 1/N_h)^n_h
}
prob_AE_AASs <- function(n_h, N_h) {
  n_h / N_h
}

# Probabilidades por escola (sname = identificador)
dados_probs <- dados %>%
  left_join(et_pop %>% select(stype, N_h, n_uniforme, n_proporcional, n_neyman),
            by = "stype") %>%
  mutate(
    # AE proporcional
    pi_proporcional = prob_AE_AASs(n_proporcional, N_h),
    # AE uniforme
    pi_uniforme = prob_AE_AASs(n_uniforme, N_h),
    # AE Neyman
    pi_neyman = prob_AE_AASs(n_neyman, N_h),
    # AASc e AASs globais
    pi_AASc = pi_AASc,
    pi_AASs = pi_AASs
  )

# -------------------
# Resultados
# -------------------
cat("\n--- Probabilidade de seleção de cada escola (n = 1500) ---\n")
print(
  dados_probs %>%
    select(sname, stype, pi_AASc, pi_AASs, pi_proporcional, pi_uniforme, pi_neyman) %>%
    head(100)
)

# Estatísticas resumo por estrato
cat("\n--- Probabilidade média de seleção por estrato ---\n")
print(
  dados_probs %>%
    group_by(stype) %>%
    summarise(
      mean_pi_proporcional = mean(pi_proporcional),
      mean_pi_uniforme     = mean(pi_uniforme),
      mean_pi_neyman       = mean(pi_neyman),
      .groups = "drop"
    )
)
