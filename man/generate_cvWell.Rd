% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellViability.R
\name{generate_cvWell}
\alias{generate_cvWell}
\title{Create mongo collection for cell viability by well (cv_well) from well treated collection}
\usage{
generate_cvWell(
  file_path,
  mongoUrl = "",
  rerun = FALSE,
  use_db = T,
  json_collection_path = ""
)
}
\arguments{
\item{file_path}{character string: file path to the top level directory of cell viability Harmony files for an HTPP dataset (i.e., the directory above plate-level directories)}

\item{mongoUrl}{character string: URL to connect to MongoDB for HTPP dataset; can be created using the mongoURL function in htpp.pl}

\item{rerun}{boolean: rerun = TRUE will drop existing cv_well collection and reinsert; FALSE by default}

\item{use_db}{boolean: Determines whether mongoDB will be used or not; default is TRUE}

\item{json_collection_path}{character: Full file path to where JSON collections will be stored}
}
\description{
Create mongo collection for cell viability by well (cv_well) from well treated collection
}
