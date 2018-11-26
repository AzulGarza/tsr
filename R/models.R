library(lmtest)
library(forecast)


# Models ------------------------------------------------------------------

ma_model <- function(formula, data, ma_order = NULL, complete_data = F,
                     include_mean = T){
  # complete_data if the data is complete respect to the formula
  # include_mean is the intercept

  #FORMULA
  info_formula <- info_formula(formula)

  # DATA
  new_data <- data
  if(!complete_data){
    # If data is not complete
    new_data <- model.frame(info_formula[["formula"]], data)
  }

  #VARS
  dependent <- info_formula[["dep_var"]]
  independent <- info_formula[["indep_var"]]


  #PARAMETERS TO ESTIMATE
  pars <- ma_order
  order_ma <- rep(0, 3)
  if(is.vector(ma_order)){
    if(include_mean){
      pars <- c(ma_order, NA, rep(NA, length(independent)))
    } else {
      pars <- c(ma_order, rep(NA, length(independent)))
    }

    order_ma[3] <- length(ma_order)
  }



  # FITTING MODEL
  arima(
    x = new_data[, dependent],
    order = order_ma,
    fixed = pars,
    xreg = new_data[, independent],
    include.mean = include_mean
  )

}


dynamic_model <- function(formula, data, order = c(1, 0, 0),
                          include_mean = T, type_model = "mdg",
                          verbose = F, n_forecast = NULL){
  # FORMULA
  dynamic_spec <- get_dynamic_model(formula, order = order[1:2], type_model)
  inertia_control_formula <- dynamic_spec[["formula"]]
  inertia_control <- c(dynamic_spec[["inertia"]], dynamic_spec[["control"]])
  independent <- dynamic_spec[["indep_var"]]
  dependent <- dynamic_spec[["dep_var"]]

  # ANCHOR
  anchor <-
    paste0(
      dependent,
      " ~ ",
      paste(
        independent,
        collapse = " + "
      )
    )

  # DATA
  new_data <- model.frame(inertia_control_formula, data)

  ### INERTIA AND CONTROL PART
  # FITTING MODEL
  main_model <- anchor

  for(term in inertia_control){
    # MODEL UPDATED WITH THE THERM
    new_model <- paste(main_model, term, sep = " + ")

    #FITTING MODEL
    fit_model <-
      ma_model(
        new_model,
        new_data,
        #ma_order = rep(NA, order[3]),
        complete_data = T,
        include_mean = include_mean
      )

    if(verbose){
      print(new_model)
    }

    #IMPROVING THE MODEL PART
    is_significant <- check_significance_var(fit_model, term, verbose = verbose)

    if(is_significant){
      # IF the last term is significant
      fix_model <-
        fix_significance_model(
          fit_model,
          new_model,
          independent,
          new_data,
          include_mean,
          NULL,
          verbose = verbose
        )

      new_model <- fix_model[["model"]]

    } else {
      # IF the last term is not significant simply we remove it
      # from the formula
      new_model <- drop_term_formula(new_model, term)
    }


    main_model <- new_model
  }

  ### MA PART
  ma_terms <- NULL
  for(ma_term in 1:order[3]){
    # DEFINING MA TERMS
    ma_chr <- paste0("ma", ma_term)
    ma_terms <- c(ma_terms, NA)

    # FITTING MODEL WITH MA PART
    fit_model <-
      ma_model(
        new_model,
        new_data,
        ma_order = ma_terms,
        complete_data = T,
        include_mean = include_mean
      )

    if(verbose){
      print(new_model)
      print(ma_terms)
    }

    #IMPROVING THE MODEL MA PART
    if(check_significance_var(fit_model, ma_chr, verbose = verbose)){
      # IF the last term is significant
      fix_model <-
        fix_significance_model(
          fit_model,
          main_model,
          independent,
          new_data,
          include_mean,
          ma_terms,
          verbose = verbose
        )

      new_model <- fix_model[["model"]]
      ma_terms <- fix_model[["ma_order"]]

    } else {
      # IF the last term is not significant simply we remove it
      # from the terms
      ma_terms[ma_term] <- 0
    }

    main_model <- new_model

  }

  final_model <-
    ma_model(
      main_model,
      new_data,
      ma_order = ma_terms,
      complete_data = T,
      include_mean = include_mean
    )

  # FIT MODEL
  attr(final_model, "fit_model") <-
    list(
      main_model = main_model,
      new_data = new_data,
      ma_order = ma_terms,
      complete_data = T,
      include_mean = include_mean
    )

  # FITTED VALUES FROM THE MODEL
  diff_rows <- nrow(data) - nrow(new_data)

  fitted <- new_data[[dependent]] - residuals(final_model)
  # if(type_model == "mce"){
  #   diff_rows <- diff_rows - 1
  #   original_dep <- dynamic_spec[["original_dep"]]
  #   #print(fitted)
  #   #print(length(fitted))
  #   fitted <- diffinv(fitted, xi = data[[original_dep]][diff_rows+1])
  # }
  fitted <- c(rep(NA, diff_rows), fitted)


  # LAST RESULTS
  orig_formula <- info_formula(formula)

  # RETURN
  print(final_model)
  attr(final_model, "fitted") <- fitted
  attr(final_model, "type") <- type_model
  attr(final_model, "dep_var") <- orig_formula[["dep_var"]]
  attr(final_model, "indep_var") <- orig_formula[["indep_var"]]
  attr(final_model, "residuals") <- residuals(final_model)

  #FORECASTING
  if(length(n_forecast)){
    assign("new_data", new_data, envir = .GlobalEnv)
    assign("dependent", dependent, envir = .GlobalEnv)
    assign("independent", independent, envir = .GlobalEnv)
    predict <- forecast(final_model, h = n_forecast)
    rm(list = lis(new_data, dependent, independet))
    attr(final_model, "predicted") <- predict
  }

  #RETURNS

  final_model

  # list(
  #   model = final_model,
  #   fitted = fitted,
  #   type = type_model,
  #   dep_var = orig_formula[["dep_var"]],
  #   indep_var = orig_formula[["indep_var"]]
  # )

}
