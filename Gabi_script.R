# Modelagem ----

## Pacotes ----
pacman::p_load(tidyverse, readxl,broom,aod)

## Leitura das bases ----

base_modelagem <- read_excel("Amostra_g04_Davi_GabrielaLobo_NataliaOliveira.xlsx") %>% 
  rename_all(~ c("numero_identficacao", "resultado_radiografia", "estagio_tumor", "nivel_fostatase", "env_nodal")) 

full_size <- nrow(base_modelagem)
n_cols <- ncol(base_modelagem)

base_teste <- read_excel("Amostra_VALIDACAO.xlsx") %>% 
  rename_all(~ c("numero_identficacao", "resultado_radiografia", "estagio_tumor", "nivel_fostatase", "env_nodal"))

## Teste de razao de verossimilhança ----
# H0: modelo nulo é adequadro
# H1: modelo saturado é adequado

### modelo completo/saturado 
modelo_saturado <- glm(env_nodal~resultado_radiografia+estagio_tumor+nivel_fostatase, data = base_modelagem, family = "binomial")
summary(modelo_completo)

## modelo somente com o intercepto 
modelo_nulo <- glm(env_nodal~1, data=base_modelagem, family = "binomial")
summary(modelo_nulo)

# Extrair log-likelihood dos modelos
logLik_nulo <- logLik(modelo_nulo)
logLik_saturado <- logLik(modelo_saturado)

# Calcular a estatística do teste de Razão de Verossimilhança
G2 <- -2 * (logLik_nulo - logLik_saturado)
G2_value <- as.numeric(G2)

# Calcular os graus de liberdade (diferença no número de parâmetros)
df <- df.residual(modelo_nulo) - df.residual(modelo_saturado)

# Calcular o p-valor
p_value <- pchisq(G2_value, df = df, lower.tail = FALSE)
## rejeitou h0, logo pelo menos um dos coeficientes de regressão é não nulo

## Selecao de variaveis/modelo ----

model_matrix_metrics <- function(data, variable_matrix) {
  models_formula <- apply(variable_matrix, MARGIN = 2, function(names_col) {
    glue::glue("env_nodal ~ {paste(names_col, collapse = ' + ')}")
  })
  return(purrr::map(
    models_formula, ~ glm(data = data, formula = .x, family = "binomial")
  ))
}

## nomes das variaveis explicativas
covariables <- names(base_modelagem)[names(base_modelagem) != "env_nodal"]

## 15 modelos
all_models <- map(
  1:(n_cols - 1),
  ~ model_matrix_metrics(
    data = base_modelagem,
    combn(covariables, .x)
  )
) %>%
  reduce(c)

model_info_to_plot <- data.frame(matrix(ncol = 8, nrow = n_models)) %>%
  rename_all(~ c(
    "Parametros", "n_parametros",
    "G^2", "g.l.", "p-valor",
    "AIC","BIC", "-2logL", "ACC", "round_threshold"
  ))

best_round_threshold <- function(data, model, par) {
  yhat <- predict(model, data, type = "response") >= par
  return(sum(data$Poupanca == yhat) / nrow(data))
}

for (i in 1:n_models) {
  current_model <- all_models[[i]]
  model_info_to_plot[i, "Parametros"] <- current_model$formula
  model_info_to_plot[i, "n_parametros"] <- length(
    current_model$coefficients
  )
  g2 <- current_model$deviance - modelo_saturado$deviance
  g2_gl <- current_model$df.residual - modelo_saturado$df.residual
  model_info_to_plot[i, "G^2"] <- g2
  model_info_to_plot[i, "g.l."] <- g2_gl
  model_info_to_plot[i, "p-valor"] <- pchisq(
    g2, g2_gl,
    lower.tail = FALSE
  )
  model_info_to_plot[i, "AIC"] <- current_model$aic
  
  model_info_to_plot[i, "BIC"] <- BIC(current_model)
  
  model_info_to_plot[i, "-2logL"] <- logLik(current_model)
    
  threshold <- optim(
    par = 0, fn = best_round_threshold, data = base_teste, model = current_model,
    method = "Brent", lower = 0, upper = 1
  )$value
  model_info_to_plot[i, "round_threshold"] <- threshold
  ytrue <- base_teste$Poupanca
  yhat <- predict(current_model, base_teste, type = "response") >= threshold
  
  model_info_to_plot[i, "ACC"] <- sum(ytrue == yhat) / nrow(base_teste)
}

###  plot AIC ----
ggplot(
  data = model_info_to_plot,
  aes(x = `n_parametros`, y = `AIC`)
) +
  geom_point(size = 5) +
  xlab("Quantidade de parâmetros") +
  ylab("AIC") +
  theme_bw()

### Plot BIC ----
ggplot(
  data = model_info_to_plot,
  aes(x = `n_parametros`, y = `BIC`)
) +
  geom_point(size = 5) +
  xlab("Quantidade de parâmetros") +
  ylab("BIC") +
  theme_bw()

### Plot Log vero ----
ggplot(
  data = model_info_to_plot,
  aes(x = `n_parametros`, y = `-2logL`)
) +
  geom_point(size = 5) +
  xlab("Quantidade de parâmetros") +
  ylab("BIC") +
  theme_bw()

### Modelos acima de 0.05 ----
## acima de 0.05
model_info_to_plot %>%
  select(-c(n_parametros, round_threshold, ACC)) %>%
  filter(`p-valor` > 0.05)

## Como só ficaram dois modelos no final, vou fazer  o teste de razao de verossimlhança para ver a diferença
# H0: modelo sem fosfatase é adequado
# H1: modelo saturado é adequado

## modelo com resultado da radiografia, estagio do tumor
modelo_semfosfatase <- glm(env_nodal~resultado_radiografia+estagio_tumor, data=base_modelagem, family = "binomial")

# Extrair log-likelihood dos modelos
logLik_semfosfatase <- logLik(modelo_semfosfatase)
logLik_saturado <- logLik(modelo_saturado)

# Calcular a estatística do teste de Razão de Verossimilhança
G2 <- -2 * (logLik_semfosfatase - logLik_saturado)
G2_value <- as.numeric(G2)

# Calcular os graus de liberdade (diferença no número de parâmetros)
df <- df.residual(modelo_semfosfatase) - df.residual(modelo_saturado)

# Calcular o p-valor
p_value <- pchisq(G2_value, df = df, lower.tail = FALSE)
## Não rejeitou h0, logo não há indicios de que o modelo sem fosfatase

## Efeito do nível de fosfatase ácida na predição para envolvimento nodal ----
## teste de WALD
## HO: Não há efeito de fosfatase acida
## H1: Há efeito de fosfatase acida

wald_result <- wald.test(b = coef(modelo_saturado), Sigma = vcov(modelo_saturado), Terms = 4)
## não rejeitou ho, logo não há indicios de que há efeito do nivel de fosfatase no modelo
## faço um intervalo de confianca??

## Analise de residuos ----
## vou fazer a analise de residuos do modelo sem fosfatase!!! Mas podemos mudar depois