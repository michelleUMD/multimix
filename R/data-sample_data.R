#' Example long-format binary dataset used for model fitting
#'
#' This dataset contains longitudinal drug administration data for
#' a post-operative inflammatory medication, used as an example for multimix modeling.
#'
#' @format A data frame with X rows and Y variables:
#' \describe{
#'   \item{Subject_ID}{Integer patient identifier}
#'   \item{Time}{Numeric time of measurement (months post-op)}
#'   \item{Binary_outcome}{0/1 binary if patient is on drug}
#' }
"sample_data"
