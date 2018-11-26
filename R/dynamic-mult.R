library(tidyr)
library(ggplot2)

parse_number_special <- function(string){
  parsed <- parse_number(string)

  if(is.na(parsed)){
    return(0)
  }

  parsed
}

regex_dyn_model <- function(var, type){
  if(type == "lag"){
    return(str_glue("^{var}|{type}\\({var}"))
  }

  str_glue("diff\\({var}")
}

dynamic_mult <- function(model, mults = 100){
  # INIT VARIABLES
  dep_var <- attr(model, "dep_var")
  ex_var <- attr(model, "indep_var")
  type_model <- attr(model, "type")

  # SIZE OF VARS
  n_ex <- length(ex_var)

  # Parámetros
  alfa <- vector("numeric", mults) #ENDÓGENA
  teta <- vector("numeric", mults) #MA TERMS
  beta <- matrix(0, nrow = mults, ncol = n_ex) #EXÓGENAS
  colnames(beta) <- ex_var

  # Multiplicadores
  mult_aislado <- matrix(0, nrow = mults, ncol = n_ex)
  mult_acumulado <- matrix(0, nrow = mults, ncol = n_ex)
  mult_shock <- vector("numeric", mults)

  # COEFFICIENTS
  coefs <- coef(model)
  names_vars <- names(coefs)

  # Llenado de las variables
  # ENDOGENOUS
  endog <- coefs[str_subset(names_vars, regex_dyn_model(dep_var, "lag"))]
  alfa[parse_number(names(endog))] <- endog
  #print(alfa)

  #MAS
  ma_vars <- coefs[str_subset(names_vars, "^ma\\d$")]
  teta[1:length(ma_vars)] <- ma_vars
  #print(teta)

  # EXOGENOUS
  adjustment <- 1
  if(type_model == "mce"){
    adjustment <- 0
  }
  for(var in ex_var){
      exog <- coefs[str_subset(names_vars, regex_dyn_model(var, "lag"))]
      beta[map_dbl(names(exog), parse_number_special) + adjustment, var] <- exog
  }
  #print(beta)

  # OTHER MODELS
  if(type_model %in% c("mdgd", "mce")){
    # INIT VARIABLES
    d_alfa <- vector("numeric", mults) # ENDOG DIFF
    d_beta <- matrix(0, nrow = mults, ncol = n_ex) # EXOG DIFF
    colnames(d_beta) <- ex_var

    # UPDATE VARS
    if(type_model == "mce"){
      alfa[1] <- alfa[1] + 1
    }
    # ENDOGENOUS
    d_endog <- coefs[str_subset(names_vars, regex_dyn_model(dep_var, "diff"))]
    d_alfa[parse_number(names(d_endog))] <- d_endog

    for(i in 1:mults){
      alfa[i] <- alfa[i] + d_alfa[i]
      alfa[i + 1] <- alfa[i + 1] - d_alfa[i]
    }

    # EXOGENOUS
    for(var in ex_var){
      d_exog <- coefs[str_subset(names_vars, regex_dyn_model(var, "diff"))]
      #print(d_exog)
      d_beta[map_dbl(names(d_exog), parse_number_special) + 1, var] <- d_exog
    }

    for(i in 1:(mults-1)){
      for(j in 1:n_ex){
        beta[i, j] <- beta[i, j] + d_beta[i, j]
        beta[i+1, j] <- beta[i+1, j] + d_beta[i, j]
      }
    }

  }

  # SHOCK
  mult_shock[1] <- 1

  # Shock multipliers
  for(j in 2:mults){
    mult_shock[j] <- teta[j-1]
    for(i in 1:(j-1)){
      mult_shock[j] <- mult_shock[j] + mult_shock[j-i]*alfa[i]
    }
  }

  #print(beta)
  # Mult exogenous
  for(k in 1:n_ex){
    mult_aislado[1, k] <- beta[1, k]
    mult_acumulado[1, k] <- beta[1, k]
    for(j in 2:mults){
      mult_aislado[j, k] <- beta[j, k]
      for(i in 1:(j-1)){
        mult_aislado[j, k] <- mult_aislado[j, k] + mult_aislado[j-i, k]*alfa[i]
      }
      mult_acumulado[j, k] <- mult_aislado[j, k] + mult_acumulado[j-1, k]
    }
  }

  colnames(mult_aislado) <- ex_var
  colnames(mult_acumulado) <- ex_var

  # DATA FRAME
  # AISLADO
  mult_aislado <- as_tibble(mult_aislado)
  mult_aislado <- mutate(mult_aislado, tipo = "isolated", mult = 1:mults)
  mult_aislado <- gather(mult_aislado, variable, valor, -tipo, -mult)

  # ACUMULADO
  mult_acumulado <- as_tibble(mult_acumulado)
  mult_acumulado <- mutate(mult_acumulado, tipo = "accumulated", mult = 1:mults)
  mult_acumulado <- gather(mult_acumulado, variable, valor, -tipo, -mult)

  # FINAL MULT
  mult <- bind_rows(mult_aislado, mult_acumulado)

  #SHOCK
  mult_shock <- as_tibble(list(valor = mult_shock))
  mult_shock <- mutate(mult_shock, tipo = "shock", mult = 1:mults)
  mult_shock <- map(ex_var, ~mutate(mult_shock, variable = .))
  mult_shock <- reduce(mult_shock, bind_rows)

  # FINAL
  final <- bind_rows(mult, mult_shock)
  #print(final)

  ggplot(final, aes(mult, valor)) +
    geom_line() +
    facet_wrap(variable ~ tipo, scales = "free_y", ncol = 3) +
    ggtitle("Dynamic multipliers") +
    xlab("") +
    ylab("") +
    theme_bw()
}
