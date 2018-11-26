drop_no_sign <- function(tidy_model, significance = 0.01, intercept = T){
  # tidy_model: Tidy model get by broom::tidy
  # If intercept = T conserves the intercept no matter what
  # returns the significant terms of the tidy model
  #print(tidy_model)
  my_tidy_model <- as.data.table(tidy_model)
  if(intercept){
    my_tidy_model[p.value<significance | term == "intercept"]
  } else {
    my_tidy_model[p.value<significance]
  }
}

actual_fixed_terms <- function(drop_model, size, terms = c("ar", "ma")){
  # drop_model: result of drop_no_sign applied
  # returns the fixed argument of _arima_
  # with respect to the therm column of
  # the drop_model
  # size is a list of the form
  # size = list(ar = , ma = )
  #terms <- c("ar", "ma")

  terms_sign <-
    map(
      terms,
      ~drop_model[str_detect(term, .), parse_number(term)]
    )


  terms_sign <- set_names(terms_sign, terms)

  #ars <- drop_model[str_detect(term, "ar"), parse_number(term)]

  #mas <- drop_model[str_detect(term, "ma"), parse_number(term)]

  res <-
    map(
      terms,
      function(term){
        vec <- vector("integer", size[[term]])
        vec[terms_sign[[term]]] <- NA
        return(vec)
      }
    )

  res <- set_names(res, terms)

  if(nrow(drop_model[term == "intercept"])){
    # If the model has intercept...
    res[["intercept"]] <- NA
  }

  return(res)
}

update_fixed_terms <- function(fixed_terms, terms = c("ar", "ma"), add = list(ar = 1, ma = 1)){
  # fixed_terms: object resulting from actual_fixed_terms
  # add_ar: additional ars
  # add_ma: additional_mas

  # Auteregresive and ma part!
  terms_arma <- c("ar", "ma")

  size <-
    map(
      terms_arma,
      ~length(fixed_terms[[.]]) + add[[.]]
    )

  size <- set_names(size, terms_arma)

  names <-
    map(
      terms_arma,
      function(term){
        term_size <- size[[term]]
        if(term_size){
          str_c(term, 1:term_size)
        }
      }
    )

  names <- reduce(names, c)

  res <-
    map(
      terms_arma,
      ~ c(fixed_terms[[.]], rep(NA, add[[.]]))
    )

  if(!is.null(fixed_terms[["intercept"]])){
    # If the model has intercept...
    names <- c(names, "intercept")
    res[["intercept"]] <- NA
  }

  res <- set_names(reduce(res, c), names)

  # Controls part
  control_terms <- setdiff(terms, c(terms_arma, "intercept"))

  if(length(control_terms)){
    # Ih there are control variables
    size_control <-
      map(
        control_terms,
        ~length(fixed_terms[[.]]) + add[[.]]
      )

    size_control <- set_names(size_control, control_terms)

    names_control <-
      map(
        control_terms,
        function(term){
          term_size <- size_control[[term]]
          if(term_size){
            str_c(term, 1:term_size)
          }
        }
      )

    names_control <- reduce(names_control, c)

    res_control <-
      map(
        control_terms,
        ~ c(fixed_terms[[.]], rep(NA, add[[.]]))
      )

    res_control <- set_names(reduce(res_control, c), names_control)
  } else {
    res_control <- NULL
  }

  return(c(res, res_control))
}

# Control variables


lag_data <- function(df, lags = 4){
  # Adds the number of lags to the df
  df <- as.data.table(df)

  df_names <- names(df)

  new_df <-
    map(
      0:lags,
      function(lag){
        new_df <- map(df, ~dplyr::lag(., lag))
        new_df <- set_names(new_df, str_c(df_names, lag + 1))
      }
    )

  reduce(new_df, bind_cols)

}


# collapse ----------------------------------------------------------------

collapse_terms <- function(actual_ft){
  # Get a list from actual_fixed_terms
  collapsed_terms <-
    actual_ft %>%
      map(
        function(term){
          if(sum(is.na(term)) == 0){
            return(integer(0))
          }

          csum <- cumsum(is.na(term))
          max_pos <- which.max(csum)

          return(term[1:max_pos])

        }
      )

  collapsed_terms <-
    set_names(
      collapsed_terms,
      names(actual_ft)
    )

  collapsed_terms
}

get_max_terms <- function(collapsed_terms){
  max_terms <-
    collapsed_terms %>%
      map(~length(.))

  max_terms <- set_names(max_terms, names(collapsed_terms))
}

# ERRORS ------------------------------------------------------------------

arima_error <- function(df, order, fixed, xreg){
  tryCatch(
    {
      arima(
        df,
        order,
        fixed = fixed,
        transform.pars = F,
        xreg = xreg,
        method = "CSS-ML"
      )
    },
    error = function(e) print(e) #"error"
  )
}



