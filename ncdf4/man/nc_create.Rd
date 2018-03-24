\name{nc_create}
\alias{nc_create}
\title{Create a netCDF File}
\description{
 Creates a new netCDF file on disk, given the variables the new file is to contain. 
}
\usage{
 nc_create( filename, vars, force_v4=FALSE, verbose=FALSE )
}
\arguments{
 \item{filename}{Name of the  netCDF file to be created.}
 \item{vars}{Either an object of class \code{ncvar4} describing the variable to be created, or a vector (or list) of such objects to be created.}
 \item{force_v4}{If TRUE, then the created output file will always be in netcdf-4 format (which
 supports more features, but
 cannot be read by version 3 of the netcdf library).  If FALSE, then the file is created
 in netcdf version 3 format UNLESS the user has requested features that require version 4.
 Deafult is FALSE.}
 \item{verbose}{If TRUE, then information is printed while the file is being created.}
}
\value{
 An object of class \code{ncdf4}, which has the fields described in \code{\link[ncdf4]{nc_open}}.
}
\references{
 http://dwpierce.com/software
}
\details{
 This routine creates a new netCDF file on disk.  The routine must be called with the variables
 that will be created in the file.
 Keep in mind that the new file may not actually be written to disk until 
 \code{\link[ncdf4]{nc_close}} is called.  Always call
 \code{\link[ncdf4]{nc_close}} when you are done with your file, or
 before exiting R!
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncdim_def}}, \code{\link[ncdf4]{ncvar_def}}. 
}
\examples{
# Define an integer dimension 
dimState <- ncdim_def( "StateNo", "count", 1:50 )

# Make an integer variable.  Note that an integer variable can have
# a double precision dimension, or vice versa; there is no fixed
# relationship between the precision of the dimension and that of the
# associated variable.  We just make an integer variable here for
# illustration purposes.
varPop <- ncvar_def("Pop", "count", dimState, -1, 
		longname="Population", prec="integer")

# Create a netCDF file with this variable
ncnew <- nc_create( "states_population.nc", varPop )

# Write some values to this variable on disk.
popAlabama <- 4447100
ncvar_put( ncnew, varPop, popAlabama, start=1, count=1 )

nc_close(ncnew)
}
\keyword{utilities}
