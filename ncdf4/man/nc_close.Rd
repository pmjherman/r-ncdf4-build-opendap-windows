\name{nc_close}
\alias{nc_close}
\title{Close a netCDF File}
\description{
 Closes an open netCDF file, which flushes any unwritten data to disk.
 Always close a netCDF file when you are done with it!
 You are risking data loss otherwise.
}
\usage{
 nc_close( nc )
}
\arguments{
 \item{nc}{An object of class \code{ncdf4} (as returned by either function
 \code{\link[ncdf4]{nc_open}} or function \code{\link[ncdf4]{nc_create}}.}
}
\references{
 http://dwpierce.com/software
}
\details{
 Data written to a netCDF file is cached in memory, for better performance.
 This data is only written out to disk when the file is closed.  Therefore, always
 remember to close a netCDF file when done with it.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{nc_sync}}.
}
\examples{
\dontrun{nc <- nc_open("salinity.nc")}
\dontrun{data <- ncvar_get( nc )	# Read the "only" var in the file}
\dontrun{nc_close(nc)}
}
\keyword{utilities}
