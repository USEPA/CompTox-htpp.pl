% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{curvePlots_htppCatMah}
\alias{curvePlots_htppCatMah}
\title{plot category Mahalanobis distances and save the plots to the location specified by file_path.}
\usage{
curvePlots_htppCatMah(
  file_path,
  study_name,
  mongoUrl = "",
  use_db = T,
  json_collection_path = ""
)
}
\arguments{
\item{file_path}{character: file path to where category mah plots will be created}

\item{study_name}{character string: name of study to be used for plot titles, should follow a similar naming convention used in other functions such as "viability_controlPlot_htppWell.R"}

\item{mongoUrl}{character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl}

\item{use_db}{boolean: Determines whether mongoDB will be used or not; default is TRUE}

\item{json_collection_path}{character: Full file path to where JSON collections will be stored}
}
\value{
a summary table for debugging
}
\description{
plot category Mahalanobis distances and save the plots to the location specified by file_path.
}
