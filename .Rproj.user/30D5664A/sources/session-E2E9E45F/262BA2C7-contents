#' AIC weight calculator
#'
#' @param aic_df A dataframe of AIC scores where each column is a model type and each row a simulation. AIC weights will be presented row wise.
#'
#' @return An AIC weight table.
#'
#' @export
calculate_aic_weights <- function(aic_df) {
  aic_df %>%
    rowwise() %>%
    mutate(across(everything(), ~ exp(-0.5 * (. - min(c_across(everything())))))) %>%
    mutate(across(everything(), ~ . / sum(c_across(everything()))))
}

