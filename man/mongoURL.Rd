% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunctions.R
\name{mongoURL}
\alias{mongoURL}
\title{Combines credentials into a MongoURL to access a database}
\usage{
mongoURL(
  host,
  user,
  passwd,
  db,
  authSource = "admin",
  authMechanism = "SCRAM-SHA-256"
)
}
\arguments{
\item{host}{The database host}

\item{user}{The MongoDB username}

\item{passwd}{The password for that database}

\item{db}{The name of the database}

\item{authSource}{The authentication source}

\item{authMechanism}{The authentication mechanism}
}
\value{
a MongoDB URL granting read access to the database
}
\description{
Combines credentials into a MongoURL to access a database
}
\examples{
mongoURL(host="test", user="readonly", passwd="passwd", db="test_db")
}
