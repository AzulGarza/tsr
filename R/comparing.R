#UPDATE THE MODEL
update_dyn_model <- function(model, new_dfr){
  # new_dfr same size of model data
  # ELEMENTOS FOR FITTING
  fit_model_elements <- attr(model, "fit_model")

  # FIRST ELEMENTS
  model_data <- fit_model_elements[["new_data"]]
  n_row_dfr <- nrow(new_dfr)
  n_row_model <- nrow(model_data)
  diff_rows <- n_row_dfr - n_row_model
  # print(n_row_dfr)
  # print(n_row_model)
  # print(diff_rows)

  # NEW DATA
  sel_new_dfr <- slice(new_dfr, (diff_rows+1):n_row_dfr)
  #print(nrow(sel_new_dfr))
  new_data <- bind_cols(model_data, sel_new_dfr)
  new_main_model <-
    paste0(
      fit_model_elements[["main_model"]],
      " + ",
      paste(
        colnames(new_dfr),
        collapse = " + "
      )
    )

  # REFITTING MODEL
  ma_model(
    new_main_model,
    new_data,
    ma_order = fit_model_elements[["ma_terms"]],
    complete_data = fit_model_elements[["complete_data"]],
    include_mean = fit_model_elements[["include_mean"]]
  )

}

# COMPARING MODELS

compare_dyn_models <- function(list_models){
  # LIST MODELS MUSTE BE NAMED

  # INITIAL VARS
  n_models <- length(list_models)
  names_models <- names(list_models)
  errors <-
    map(
      list_models,
      function(model){
        errors_model <- attr(model, "residuals")
        #print(as.vector(errors_model))
        tibble(error = as.vector(errors_model))
      }
    )

  # MATRIX
  comparison_matrix <- matrix(0, n_models, n_models)
  colnames(comparison_matrix) <- names_models
  rownames(comparison_matrix) <- names_models

  for(model in names_models){
    other_models <- names_models[names_models != model]
    # STATISTICS FROM OTHER MODELS
    statistics <-
      map_dbl(
        other_models,
        function(other_model){
          updated_model <- update_dyn_model(list_models[[model]], errors[[other_model]])
          test <- coeftest(updated_model)
          test_dfr <- tidy(test)

          unlist(test_dfr[test_dfr$term == "error", "statistic"])
        }
      )

    # RESULT
    names(statistics) <- other_models
    comparison_matrix[model, other_models] <- statistics
  }

  comparisons <- combn(names_models, 2, simplify = F)

  results_comp <-
    map_chr(
      comparisons,
      function(comparison){
        comparison_data <-comparison_matrix[comparison, comparison]
        if(abs(comparison_data[1,2]) > abs(comparison_data[2, 1])){
          return(str_glue("{comparison[1]} > {comparison[2]}"))
        } else {
          return(str_glue("{comparison[2]} > {comparison[1]}"))
        }
      }
    )

  # RESULTS OF COMPARISON
  attr(comparison_matrix, "results") <- results_comp

  # FINAL RESULTS
  comparison_matrix

}
