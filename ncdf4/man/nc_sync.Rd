\name{nc_sync}
\alias{nc_sync}
\title{Synchronize (flush to disk) a netCDF File}
\description{
 Flushes any pending operations on a netCDF file to disk.
}
\usage{
 nc_sync( nc )
}
\arguments{
 \item{nc}{An object of class \code{ncdf4} that is opened for writing (as returned by either 
 function \code{\link[ncdf4]{nc_open}}\code{(..., write=TRUE)}
 or function \code{\link[ncdf4]{nc_create}}, indicating what file is being written to.}
}
\references{
 http://dwpierce.com/software
}
\details{
 Data in a netCDF file is cached in memory, for better performance.
 An example of when this might be bad is if a long-running job writes one
 timestep of the output file at a time; if the job crashes near the end, 
 the results of many timesteps might be lost.  In such an event, the user
 can manually force any cached data to be written to disk using this call.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\examples{
# The time you would use the sync.ncdf function is when you have an unlimited
# dimension and are writing to the file timestep-by-timestep. Make a netCDF file 
# that has an unlimited dimension for illustration.
nx <- 5
ny <- 8
dimx <- ncdim_def( "X", "meters", 1:nx )
dimy <- ncdim_def( "Y", "meters", 1:ny )
dimt <- ncdim_def( "Time", "days since 1900-01-01", 0, unlim=TRUE )

vartemp <- ncvar_def( "Temperature", "degC", list(dimx,dimy,dimt), 1.e30 )
nc  <- nc_create( "temperature.nc", vartemp )

nt <- 10  # Imagine this is actually some very large number of timesteps
for( i in 1:nt ) {
	# Long, slow computation to get the data ... for illustration, we just
	# use the following:
	data <- runif(nx*ny)

	# Write the data to this timestep
	ncvar_put( nc, vartemp, data, start=c(1,1,i), count=c(nx,ny,1) )

	# Write the time value for this timestep as well
	timeval <- i*10
	ncvar_put( nc, dimt, timeval, start=i, count=1 )

	# Flush this timestep's data to the file so we dont lose it
	# if there is a crash or other problem
	nc_sync( nc )
	}

# Always remember to close the file when done!!
nc_close(nc)
}
\keyword{utilities}
