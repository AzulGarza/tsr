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

    fixed_no_errors <- NULL
    last_fixed <- NULL
    # ARMA PART
    if(intercept){
      fixed <- c(ar1 = NA, intercept = NA)
    } else {
      fixed <- c(ar1 = NA, intercept = 0)
    }
    terms <- list()
    terms[["ar"]][[1]] <- NA
    terms[["intercept"]][[1]] <- NA

    for(ar in 1:max_ar){
      #print(fixed)
      model <-
        arima_error(
          df,
          c(ar, 0, 0),
          fixed = fixed,
          xreg = NULL
        )

      if(model == "error"){
        # En caso de error el modelo es el último que no posee errores
        # Actualizar los terms
        terms[["ar"]][[ar]] <- 0
      } else {
        fixed_no_errors <- fixed
        tidy_model <- tidy(coeftest(model))
        drop_model <- drop_no_sign(tidy_model)
        terms <- actual_fixed_terms(drop_model, size = list("ar" = ar, "ma" = 0))
      }

      # Update model only if another iteration
      if(ar < max_ar){
        #New fixed
        fixed <- update_fixed_terms(terms, add = list("ar" = 1, "ma" = 0))
        last_fixed <- fixed
      }
    }

    #CONTROLS PART
    controls <- names(xreg_df)

    max_terms <- list("ar" = max_ar, "ma" = 0)

    my_xreg <- lag_data(xreg_df, max_cont)

    iter_vars <- NULL
    add_list <- list("ar" = 0, "ma" = 0)
    for(control in controls){
      max_terms[[control]] <- 1
      add_list[[control]] <- 1
      fixed <- update_fixed_terms(terms, terms = names(max_terms), add = add_list)
      for(ind in 1:max_cont){
        iter_vars <- c(iter_vars, paste0(control, ind))
        model <-
          arima_error(
            df,
            c(max_ar, 0, 0),
            fixed = fixed,
            xreg = my_xreg[,iter_vars] # Seleccionar variables
          )
        #print(fixed)
        #print(model)
        if(model == "error"){
          # En caso de error el modelo es el último que no posee errores
          # Actualizar los terms
          terms[[control]][[ind]] <- 0
          max_terms[[control]] <- ind
        } else {
          fixed_no_errors <- fixed
          # Update model only if another iteration
          tidy_model <- tidy(coeftest(model))

          drop_model <- drop_no_sign(tidy_model)
          max_terms[[control]] <- ind

          terms <- actual_fixed_terms(drop_model, size = max_terms, terms = names(max_terms))
        }

        if(ind < max_cont){
          fixed <- update_fixed_terms(terms, terms = names(max_terms), add = add_list)
        } else {
          add_list[[control]] <- 0
          fixed <- update_fixed_terms(terms, terms = names(max_terms), add = add_list)
        }
      }
    }

    #print("-----MA PART----")
    #MA PART!
    add_list[["ma"]] <- 1
    #terms <- actual_fixed_terms(drop_model, size = max_terms, terms = names(max_terms))
    fixed <- update_fixed_terms(terms, terms = names(max_terms), add = add_list)
    fixed_no_errors <- fixed

    for(ma in 1:max_ma){
      #print(fixed)
      model <-
        arima_error(
          df,
          c(max_ar, 0, ma),
          fixed = fixed,
          xreg = my_xreg[,iter_vars]
        )
      #print(fixed)
      #print(model)

      if(model == "error"){
        # En caso de error el modelo es el último que no posee errores
        # Actualizar los terms
        terms[["ma"]][[ma]] <- 0
        max_terms[["ma"]] <- ma
      } else {
        fixed_no_errors <- fixed
        tidy_model <- tidy(coeftest(model))
        drop_model <- drop_no_sign(tidy_model)
        max_terms[["ma"]] <- ma
        terms <- actual_fixed_terms(drop_model, size = max_terms, terms = names(max_terms))
      }

      # Update model only if another iteration

      if(ma < max_ma){
        fixed <- update_fixed_terms(terms, terms = names(max_terms),add = add_list)
        last_fixed <- fixed
      } else {
        add_list[["ma"]] <- 0
        #fixed <- update_fixed_terms(terms, terms = names(max_terms),add = add_list)
      }

    }
    #print(fixed_no_errors)
    fixed_no_errors <- update_fixed_terms(terms, terms = names(max_terms),add = add_list)
    #print(fixed_no_errors)
    arima_error(
      df,
      c(max_ar, 0, max_ma),
      fixed = fixed_no_errors,
      xreg = my_xreg[,iter_vars]
    )


  }

