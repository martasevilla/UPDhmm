#' HMM data with alelle frequencies in different states
#'
#'
#' It contains the emission probabilites, transmission probabilites ...
#'
#'
#' @usage data(hmm)
#' @format A list with 5 different elements
#' \describe{
#' \item{states}{5 different possible states}
#' \item{symbols}{code symbol use for the genotype's combination}
#' \item{startProbs}{the initial probabilites of every state}
#' \item{transProbs}{probability between changing from one state to another}
#' \item{emsissionProbs}{given a certain genotype's combination what are the odds of every possible state}
#' }
#' @source Created in house for calculating Uniparental disomy events
#' @examples data(hmm)
"hmm"
