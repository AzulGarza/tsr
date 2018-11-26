library(modelr)
library(formula.tools)
library(purrr)
library(dplyr)
library(broom)
library(stringr)
library(readr)

# working with specific terms ---------------------------------------------

## WORKING WITH FORMULA

info_formula <- function(formula){
  formula <- as.formula(formula)

  # VARS info
  dependent <- lhs.vars(formula)
  if(!is.character(dependent)){
    dependent <- labels(dependent)
  }
  independent <- labels(terms(formula))

  # RETURN
  list(
    formula = formula,
    dep_var = dependent,
    indep_var = independent
  )

}

add_scape <- function(string){
  new_string <- str_replace_all(string, "\\(", "\\\\(")
  new_string <- str_replace_all(new_string, "\\)", "\\\\)")

  new_string
}

drop_term_formula <- function(formula, chr_term){
  #INIT INFO
  info_formula <- info_formula(formula)
  dependent <- info_formula[["dep_var"]]
  independent <- info_formula[["indep_var"]]

  #NEW TERM
  new_chr_term <-
    add_scape(
      str_remove_all(
        chr_term,
        "\\s"
      )
    )

  #DETECT TERM
  detect_drop <-
    str_detect(
      str_remove_all(independent, "\\s"),
      new_chr_term
    )

  #CONSTRUCT NEW FORMULA
  new_formula <-
    paste0(
      dependent,
      " ~ ",
      paste0(
        independent[!detect_drop],
        collapse = " + "
      )
    )

  #RETURN
  new_formula
}

# WORKING WITH TERMS
add_terms <- function(var, order, type_model = "mdg", inertia = F){
  # Initial case
  type <- "lag"
  init <- 1

  # Other models

  if(type_model %in% c("mdgd", "mce")){
    var <- paste0("c(NA, diff(", var, "))")
    order <- order - 1
    if(!inertia){
      init <- 0
    }
  }

  map_chr(init:order, ~paste0(type, "(", var, ",", .x, ")"))
}



get_dynamic_model <- function(formula, order = c(1, 1), type_model = "mdg"){
  ###
  # Update formula to include intertia orders and control orders
  # IN mdgd type model the order muste be at least 2

  info_formula <- info_formula(formula)
  formula <- info_formula[["formula"]]
  # VARIABLES

  dependent <- info_formula[["dep_var"]]
  independent <- info_formula[["indep_var"]]

  # ORDERS
  inertia_order <- order[1]
  control_order <- order[2]


  # INTERTIA part

  inertia_terms <- add_terms(dependent, inertia_order, type_model, inertia = T)
  inertia_formula <- paste0(inertia_terms, collapse = " + ")

  # CONTROLS part

  new_control_vars_l <-
    map(
      independent,
      add_terms,
      control_order,
      type_model
    )
  new_control_vars <- reduce(new_control_vars_l, c)
  controls_formula <- paste0(new_control_vars, collapse = " + ")

  # TOTAL FORMULA
  formula_chr <- as.character(formula)


  #DIFFERENT MODELS
  if(type_model == "mdgd"){
    # NEW VARS FOR THIS MODEL
    new_dependent <- dependent
    new_independent <- paste0("lag(", dependent, ",1)")
    new_independent <- c(independent, new_independent)

    # NEW FORMULA FOR MDGD
    new_formula_chr <-
      paste(
        formula_chr,
        paste0(
          new_independent,
          collapse = " + "
        ),
        inertia_formula,
        controls_formula,
        sep = " + "
      )

  } else if(type_model == "mce"){
    # NEW VARIABLES FOR MCE
    new_dependent <- str_glue("c(NA, diff({dependent}))")
    new_independent <-
      c(
        str_glue("lag({dependent}, 1)"),
        str_glue("lag({independent}, 1)")
      )
    new_independent_for <- paste(new_independent, collapse = " + ")

    #NEW FORMULA
    new_formula_chr <-
      paste(
        str_glue("{new_dependent} ~ {new_independent_for}"),
        inertia_formula,
        controls_formula,
        sep = " + "
      )
  } else {
    # SAME VARIABLES
    new_dependent <- dependent
    new_independent <- independent

    #NEW MODEL
    new_formula_chr <-
      paste(
        formula_chr,
        inertia_formula,
        controls_formula,
        sep = " + "
      )
  }

  # RETURN
  list(
    original_formula = formula,
    formula = as.formula(new_formula_chr),
    dep_var = new_dependent,
    indep_var = new_independent,
    inertia = inertia_terms,
    control = new_control_vars,
    original_dep = dependent
  )

}



# remove terms ------------------------------------------------------------

remove_no_significant <- function(fitted_model, anchor,
                                  significance = 0.05, verbose = F){
  # ANALYSIS ONLY IF IS NOT AN ANCHOR
  test <- coeftest(fitted_model)
  test <- tidy(test)
  test <- filter(test, !term %in% anchor)

  # ANALYSIS OF PVALUE
  all_terms_significant <- all(test[["p.value"]] < significance)
  drop_term <- NA

  if(!all_terms_significant){
    # IF ALL EXCEPT ANCHOR ARE SIGNIFICANT
    max_no_significant <- filter(test, p.value == max(p.value, na.rm = T))
    drop_term <- pull(max_no_significant, term)
  }

  if(verbose){
    print(test)
    print(drop_term)
  }

  #RETURN
  list(
    all_significant = all_terms_significant,
    drop_term = drop_term
  )
}



# improve model -----------------------------------------------------------

check_significance_var <- function(fitted_model, var, significance = 0.05, verbose = F){
  # TEST
  test <- coeftest(fitted_model)
  test <- tidy(test)
  var_result <-
    filter(
      test,
      str_remove_all(term, "\\s") == str_remove_all(var, "\\s")
    )

  if(verbose){
    print(test)
  }
  # RESULT
  if(pull(var_result, "p.value") < significance){
    return(T)
  }

  return(F)
}



fix_significance_model <- function(fitted_model, model, indep_var,
                                   data_complete, include_mean, ma_order,
                                   verbose = F){
  #DISCLAIMER
  if(verbose){
    print("---FIXING MODEL---")
  }

  # NULL HYPOTHESIS
  new_model <- model
  new_ma_order <- ma_order

  # CHECK OF SIGNIFICANCE OF THE MODEL
  improve_model <-
    remove_no_significant(
      fitted_model,
      c(indep_var, "intercept"),
      verbose = verbose
    )
  significant_model <- improve_model[["all_significant"]]
  drop_term <- improve_model[["drop_term"]]

  while(!significant_model){
    # MEANWHILE NO MODEL SIGNIFICANCE
    if(str_detect(drop_term, "^ma\\d$")){
      ma_term <- parse_number(drop_term)
      new_ma_order[ma_term] <- 0
    } else {
      new_model <- drop_term_formula(new_model, drop_term)
    }

    # REFIT MODEL
    fit_model <-
      ma_model(
        new_model,
        data_complete,
        ma_order = new_ma_order,
        complete_data = T,
        include_mean = include_mean
      )

    improve_model <-
      remove_no_significant(
        fit_model,
        c(indep_var, "intercept"),
        verbose = verbose
      )
    significant_model <- improve_model[["all_significant"]]
    drop_term <- improve_model[["drop_term"]]
  }

  if(verbose){
    print("---END FIXING MODEL---")
  }

  ## RETURN
  list(
    model = new_model,
    ma_order = new_ma_order
  )
}
