\name{nc_redef}
\alias{nc_redef}
\title{Puts a netCDF file back into define mode}
\description{
 Puts a netCDF that is not currently in define mode back into define mode.
}
\usage{
 nc_redef( nc )
}
\arguments{
 \item{nc}{An object of class \code{ncdf4} (as returned by either 
 function \code{\link[ncdf4]{nc_open}}(..., write=TRUE)
 or function \code{\link[ncdf4]{nc_create}}, indicating what file to operate on.}
}
\references{
 http://dwpierce.com/software
}
\details{
 Typically, users will never need this function.

 NetCDF files can be in "define mode", at which time dimensions and variables can be
 defined, or new attributes added to a file, or in "data mode", at which time data
 can be read from the file.
 This call puts a file that is currently in data mode back into define mode.
 This functionality is handled transparently by the ncdf4 library, so users
 will never need to call this unless they are doing advanced manipulations
 of netcdf files.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{nc_enddef}}.
}
\examples{
# This function is for advanced useage only, and will never
# be needed by the typical users R code.
}
\note{
 The typical user will never need this call, nor will ever have to worry about
 "define mode" or "data mode".  THIS CALL IS PROVIDED FOR ADVANCED USERS ONLY!
 If the user goes through this package's standard functional interface,
 the file will always automatically be set to whatever mode it needs to be in
 without the user having to do anything.
 For example, the call to write an attribute (\code{\link[ncdf4]{ncatt_put}})
 handles this automatically.
}
\keyword{utilities}
