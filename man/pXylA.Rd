% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{pXylA}
\alias{pXylA}
\title{GFP driven by pXylA}
\format{A data frame with 1134 rows and 5 variables:
\describe{
   \item{sample}{sample number (1--6)}
   \item{time}{time in hours ()}
   \item{measurement}{biomass/fluorescence measurements}
   \item{meas_type}{whether or not the measurement is of biomass or fluorescence ("biomass", "fluorescence")}
   \item{induced}{whether or not the sample is induced or control (1 (induced), 0 (uninduced))}
}}
\usage{
pXylA
}
\description{
A dataset containing biomass and fluorescence data for E. coli strains with plasmid-based
expression of GFP driven by a xylose-inducible pXylA promoter. 3 replicates of induced
and 3 of uninduced cultures are included. Note: RAW DATA - does not include background
correction described in publication.
}
\keyword{datasets}
