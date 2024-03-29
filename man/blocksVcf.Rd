% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blocksVcf.R
\name{blocksVcf}
\alias{blocksVcf}
\title{Function to simplify contiguous variants with the same state into
blocks.}
\usage{
blocksVcf(df)
}
\arguments{
\item{df}{data.frame resulting from the \code{as_df_vcf} function.}
}
\value{
data.frame containing information on the chromosome, start
#' position of the block, end position of the block, and predicted state.
}
\description{
Function to simplify contiguous variants with the same state into
blocks.
}
