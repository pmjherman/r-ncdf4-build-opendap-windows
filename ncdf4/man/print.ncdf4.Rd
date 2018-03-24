\name{print.ncdf4}
\alias{print.ncdf4}
\title{Print Information About a netCDF File}
\description{
 Prints information about a netCDF file, including the variables and dimensions
 it contains.
}
\usage{
 \method{print}{ncdf4}( x, ... )
}
\arguments{
 \item{x}{An object of class "ncdf4".}
 \item{...}{Extra arguments are passed to the generic print function.}
}
\references{
 http://dwpierce.com/software
}
\details{
 NetCDF files contain variables, which themselves have dimensions.  
 This routine prints out useful information about a netCDF file's
 variables and dimensions.  It is overloaded on the regular print
 function, so if "nc" is an object of class "ncdf4", then just
 calling \code{print(nc)} will suffice.  Objects of class
 "ncdf4" are returned from \code{\link[ncdf4]{nc_open}} and \code{\link[ncdf4]{nc_create}}.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncvar_def}}
}
\examples{
# Open a netCDF file, print information about it
nc <- nc_open( "salinity.nc" )
print(nc)
}
\keyword{utilities}
