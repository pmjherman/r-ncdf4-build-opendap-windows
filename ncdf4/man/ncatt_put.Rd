\name{ncatt_put}
\alias{ncatt_put}
\title{Put an attribute into a netCDF file}
\description{
 Writes an attribute to a netCDF file.
}
\usage{
 ncatt_put( nc, varid, attname, attval, prec=NA, verbose=FALSE, 
 	definemode=FALSE )
}
\arguments{
 \item{nc}{An object of class \code{ncdf4} (as returned from \code{\link[ncdf4]{nc_open}}),
 indicating what file to write to.}
 \item{varid}{The variable whose attribute is to be written.  Can be a
 character string with the variable's name, an object of class \code{ncvar4},
 or an id contained in the "id" field of a \code{ncvar} object. As a
 special case, if varid==0, then a global attribute is written
 instead of a variable's attribute.}
 \item{attname}{Name of the attribute to write.}
 \item{attval}{Attribute to write.}
 \item{prec}{Precision to write the attribute.  If not specified, the written precision is the same as the variable whose attribute this is.  This can be overridden by specifying this argument with a value of "short", "float", "double", or "text".}
 \item{verbose}{Can be set to TRUE if additional information is desired while the attribute is being created.}
 \item{definemode}{If FALSE (the default), it is assumed that the file is NOT already in define mode.  Since the file must be in define mode for this call to work, the file will be put in define mode, the attribute defined, and then the file taken out of define mode.  If this argument is set to TRUE, it is assumed the file is already in define mode, and the file is also left in define mode.  If you don't know what any of this means, just leave
 this at the default value.}
}
\references{
 http://dwpierce.com/software
}
\details{
 This function write an attribute to a netCDF variable (or a global attribute
 to a netCDF file, if the passed argument "varid" is zero).
 The type of the written variable can be controlled by the "prec" argument, if
 the default behavior (the precision follows that of the associated variable) is
 not wanted.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncatt_get}}.
}
\examples{
# Make a simple netCDF file
filename <- "atttest_types.nc"
dim <- ncdim_def( "X", "inches", 1:12 )
var <- ncvar_def( "Data", "unitless", dim, -1 ) 
ncnew <- nc_create( filename, var )

# Define some attributes of various types
attvaldbl <- 3.1415926536
ncatt_put( ncnew, var, "testatt_dbl", attvaldbl, prec="double" )
attvalfloat <- c(1.0,4.0,9.0,16.0)
ncatt_put( ncnew, var, "testatt_float", attvalfloat )
# varid=0 means it is a global attribute
ncatt_put( ncnew, 0, "globalatt_int", 32000, prec="int" ) 
ncatt_put( ncnew, 0, "globalatt_short", 7, prec="short" )
ncatt_put( ncnew, 0, "description", 
	"this is a test file with attributes of various types")
nc_close(ncnew)
}
\keyword{utilities}
