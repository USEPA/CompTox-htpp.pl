% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{viability_controlPlot_htppWellNorm}
\alias{viability_controlPlot_htppWellNorm}
\title{Function to plot normalized well values based on vehicle and viability controls and write the plots to the location specified in file_path}
\usage{
viability_controlPlot_htppWellNorm(
  mongoUrl = "",
  file_path,
  vehicle_control,
  viability_positive_control,
  study_name,
  use_db = T,
  json_collection_path = ""
)
}
\arguments{
\item{mongoUrl}{A string consisting of the databased url and required credentials, generated by the mongoURL function}

\item{file_path}{A string consisting of the path where the plots will be stored}

\item{vehicle_control}{A string containing the vehicle chemical (such as dmso) used in this experiment}

\item{viability_positive_control}{A string containing the viability positive control chemical}

\item{study_name}{A string containing the name of the study}

\item{use_db}{boolean: Determines whether mongoDB will be used or not; default is TRUE}

\item{json_collection_path}{character: Full file path to where JSON collections will be stored}
}
\description{
Function to plot normalized well values based on vehicle and viability controls and write the plots to the location specified in file_path
}
