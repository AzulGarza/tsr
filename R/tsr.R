# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

require(purrr)
require(lmtest)
require(broom)
require(dplyr)
require(tseries)
require(stringr)
require(data.table)
require(readr)

stabilize <-
  function(
    df,
    lag_adf = 12,
    max_diff = 4,
    p_value_adf_tol = 0.01,
    p_value_det_tol = 0.01
  ){
    # Dicky Fuller
    adf_list <- adf.test(na.omit(df), alternative = "s", k = lag_adf)
    p_value_adf <- adf_list$p.value

    # Tendencia determinística
    time <- 1:nrow(df)
    summary_time <- summary(lm(df ~ time))
    p_value_det <- summary_time$coefficients[2, 4]

    # Variables de utilidad
    ts <- df
    n_diffs <- 0
    while(p_value_adf > p_value_adf_tol | p_value_det < p_value_det_tol | n_diffs <= max_diff){
      # Si existe alguna tendencia, diferenciamos
      ts <- diff(ts, differences = 1)
      n_diffs <- n_diffs + 1

      # DF
      adf_list <- adf.test(na.omit(ts), alternative = "s", k = lag_adf)
      p_value_adf <- adf_list$p.value

      #Determinista
      time <- 1:nrow(ts)
      summary_time <- summary(lm(ts ~ time))
      p_value_det <- summary_time$coefficients[2, 4]
    }

    summary <-
      tibble(
        diferences = n_diffs,
        lags = lag_adf,
        p_val_adf = p_value_adf,
        p_val_det = p_value_det,
        mean = mean(ts, na.rm = T),
        sd = sd(ts, na.rm = T)
      )

    list_r <-
      list(
        original_df = df,
        diff_df = ts,
        summary = summary
      )

    return(list_r)
}



auto_arma <-
  function(
    df,
    max_ar = 12,
    max_ma = 12,
    p_val_tol = 0.05,
    intercept = T
  ){

    if(intercept){
      fixed <- c(ar1 = NA, intercept = NA)
    } else {
      fixed <- c(ar1 = NA, intercept = 0)
    }

    for(ar in 1:max_ar){
      model <-
        arima(
          df,
          c(ar, 0, 0),
          fixed = fixed,
          transform.pars = F
        )

      # Update model only if another iteration
      tidy_model <- tidy(coeftest(model))

      drop_model <- drop_no_sign(tidy_model)


      if(ar < max_ar){
        terms <- actual_fixed_terms(drop_model, size = list("ar" = ar, "ma" = 0))
        fixed <- update_fixed_terms(terms, add = list("ar" = 1, "ma" = 0))
      }
    }

    # Updating for ma
    fixed <- update_fixed_terms(terms, add = list("ar" = 1, "ma" = 1))


    for(ma in 1:max_ma){
      model <-
        arima(
          df,
          c(max_ar, 0, ma),
          fixed = fixed,
          transform.pars = F
        )

      # Update model only if another iteration
      tidy_model <- tidy(coeftest(model))

      drop_model <- drop_no_sign(tidy_model)

      terms <- actual_fixed_terms(drop_model, size = list("ar" = max_ar, "ma" = ma))

      if(ma < max_ma){
        fixed <- update_fixed_terms(terms, add = list("ar" = 0, "ma" = 1))
      } else {
        fixed <- update_fixed_terms(terms, add = list("ar" = 0, "ma" = 0))
      }

    }

    arima(df, c(max_ar, 0, max_ma), fixed = fixed, transform.pars = F)

  }

auto_mdg <-
  function(
    df,
    xreg_df,
    max_ar = 12,
    max_cont = 12,
    max_ma = 12,
    p_val_tol = 0.05,
    intercept = T
  ){
    # Same arguments of auto_arma,
    # adding xreg (independent variables) in a data frame
    # (or matrix)

    #CONTROLS PART
    # THIS PART INCLUDE THE LAGS of dependent variable

    #DEP PART
    dep_var <- names(df)
    dep_df <- lag_data(df, max_cont)

    # INDEPENDET PART
    controls <- c(dep_var, names(xreg_df))

    max_terms <- list("ar" = 0, "ma" = 0)

    my_xreg <- bind_cols(dep_df, lag_data(xreg_df, max_cont))

    iter_vars <- NULL
    for(control in controls){
      for(ind in 1:max_cont){
        if(control == dep_var && ind == 1){
          next
        }

        iter_vars <- c(iter_vars, paste0(control, ind))

        model <-
          arima_error(
            df,
            c(0, 0, 0),
            fixed = NULL,
            xreg = my_xreg[,iter_vars] # Seleccionar variables
          )
        #print(model)

        #print(model)
        if(model == "error"){
        } else {
          # Tidy model for p.values
          tidy_model <- tidy(coeftest(model))
          # Get terms with sufficient significance
          drop_model <- drop_no_sign(tidy_model)
          terms_survived <- drop_model[["term"]]
          terms_survived <- terms_survived[!terms_survived == "intercept"]

        }

        # Updating new variables

        iter_vars <- terms_survived
      }
    }

    #print("-----MA PART----")

    # Original fixed
    # Intecept + iter_vars
    fixed <- c(NA, rep(NA, length(iter_vars)))

    #fixed_no_errors <- NULL
    #last_fixed <- NULL
    # ARMA PART

    #MA PART!
    # add_list[["ma"]] <- 1
    # #terms <- actual_fixed_terms(drop_model, size = max_terms, terms = names(max_terms))
    # fixed <- update_fixed_terms(terms, terms = names(max_terms), add = add_list)
    # fixed_no_errors <- fixed
    # #print(fixed)
    #
    fixed  <- NULL
    iter_max_ma <- 0
    for(ma in 1:max_ma){
      #print(fixed)
      #print(iter_vars)
      model <-
        arima_error(
          df,
          c(0, 0, ma),
          fixed = fixed,
          xreg = my_xreg[,iter_vars]
        )
      #print(fixed)
      #print(model)

      if(model == "error"){
        # En caso de error el modelo es el último que no posee errores
        # Actualizar los terms
        #terms[["ma"]][[ma]] <- 0
        #max_terms[["ma"]] <- ma
      } else {
        #fixed_no_errors <- fixed
        # Tidy model
        tidy_model <- tidy(coeftest(model))
        if(any(is.na(tidy_model))){
          #Si el ma regresa nas, pasamos al siguiente ma
          next
        }
        # Obtener los términos con suficiente significancia
        drop_model <- drop_no_sign(tidy_model)
        new_terms <- drop_model[["term"]]
        new_terms <- new_terms[!new_terms == "intercept"]

        #print(new_terms)

        # Update model
        iter_vars <- intersect(iter_vars, new_terms)

        #print(iter_vars)


        # MAS PART
        ma_terms <- str_subset(new_terms, "^ma\\d$")
        #print("MA TERMS")
        #print(ma_terms)

        fixed_ma <- vector("integer", ma)
        fixed_ma[na.omit(parse_number(ma_terms, "\\d"))] <- NA
        #print(fixed_ma)



        #print("fixed")
        #print(fixed)

        if(ma < max_ma){
          fixed <-
            c(
              fixed_ma, #MAS
              NA, #Otro MA
              NA, #Intercept
              rep(NA, length(iter_vars)) #Independent
            )
          # Otro MA
        } else {
          fixed <-
            c(
              fixed_ma, #MAS
              NA, #Intercept
              rep(NA, length(iter_vars)) #Independent
            )
        }

      }

    }
    #print("last_fixed")
    #print(fixed)
    arima_error(
      df,
      c(0, 0, max_ma),
      fixed = fixed,
      xreg = my_xreg[,iter_vars]
    )


  }


