\name{nc_enddef}
\alias{nc_enddef}
\title{Takes a netCDF file out of define mode}
\description{
 Changes a netCDF that is currently in define mode back into data mode.
}
\usage{
 nc_enddef( nc )
}
\arguments{
 \item{nc}{An object of class \code{ncdf4} (as returned by either 
 function \code{\link[ncdf4]{nc_open}}
 or function \code{\link[ncdf4]{nc_create}}, indicating what file to operate upon.}
}
\references{
 http://dwpierce.com/software
}
\details{
 NOTE: typical users will never need to use this function.

 NetCDF files can be in "define mode", at which time dimensions and variables can be
 defined, or new attributes added to a file, or in "data mode", at which time data
 can be read from the file.
 This call puts a file that is currently in define mode back into data mode.
 The ncdf4 package manages this process transparently, so normally, an end user
 will not need to call this explicitly.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{nc_redef}}.
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
