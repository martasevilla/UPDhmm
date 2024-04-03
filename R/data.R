#' HMM data for predicting UPD events in trio genomic data
#'
#' @usage data(hmm)
#' @format A list with 5 different elements
#' @description
#' This dataset provides Hidden Markov Model (HMM) parameters for predicting
#' uniparental disomy (UPD) events in trio genomic data.
#'
<<<<<<< HEAD
#' \describe{
=======
#'\describe{
>>>>>>> upstream/devel
#'   \item{states}{Five different possible states.}
#'   \item{symbols}{Code symbols used for genotype combinations.}
#'   \item{startProbs}{The initial probabilities of each state.}
#'   \item{transProbs}{Probabilities of transitioning from one state
#'   to another.}
#'   \item{emissionProbs}{Given a certain genotype combination, the odds
#'   of each possible state.}
<<<<<<< HEAD
#' }
#' @source Created in-house based on basic Mendelian rules for calculating
#' UPD events.
#' @examples data(hmm)
=======
#'}
#'@source Created in-house based on basic Mendelian rules for calculating
#' UPD events.
#'@examples data(hmm)
>>>>>>> upstream/devel
"hmm"
