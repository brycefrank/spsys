% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsample.R
\name{transform_coords}
\alias{transform_coords}
\title{Transforms a set of population unit coordinates to an integer index.}
\usage{
transform_coords(coords, d_x = NA, d_y = NA)
}
\arguments{
\item{coords}{A dataframe of x,y coordinates of population unit locations.}

\item{d_x}{The distance between points in the x dimension}

\item{d_y}{The distance between points in the y dimension}
}
\value{
A dataframe with two columns, r and c, that correspond to row and column indices
}
\description{
Note that input coordinates must be spatially aligned vertically and horizontally. This
may not always be the case, e.g. for some spatial projections that curve geographic space.
}
\keyword{internal}
