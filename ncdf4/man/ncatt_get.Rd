\name{ncatt_get}
\alias{ncatt_get}
\title{Get attribute from netCDF file}
\description{
 Reads an attribute from a netCDF file.
}
\usage{
 ncatt_get( nc, varid, attname=NA, verbose=FALSE )
}
\arguments{
 \item{nc}{An object of class \code{ncdf4} (as returned from \code{\link[ncdf4]{nc_open}}),
 indicating what file to read from.}
 \item{varid}{The variable whose attribute is to be read.  Can be a
 character string with the variable's name or an object of class \code{ncvar4}.
 As a special case, if varid==0, then a global (file) 
 attribute will be read rather than a particular variable's attribute.}
 \item{attname}{Name of the attribute to read; if not specified, a list
 containg ALL attributes of the selected variable or file is returned.}
 \item{verbose}{If TRUE, then debugging information is printed.}
}
\value{
 If an attribute name is supplied (i.e., argument
 attname is given), this returns 
 a list with two components, "hasatt" and "value".
 "hasatt" is TRUE if the named attribute was found, and FALSE otherwise.
 "value" is the (possibly vector) value of the attribute.
 If the on-disk type of the attribute is short or integer, then an integer
 value is returned.  If the on-disk type is float or double, than a 
 double value is returned.  If the on-disk type is character, than a
 character string is returned.

 If no attribute name is supplied, then this returns a list containing
 ALL the attributes for the specified variable
 along with their associated values.  For example,
 if attlist is the list returned by this call, then names(attlist) shows
 all the attributes defined for the variable, and attlist[[N]] is the
 value of the N'th attribute.
}
\references{
 http://dwpierce.com/software
}
\details{
 This function gets an attribute from a netCDF variable (or a global attribute
 from a netCDF file, if the passed argument "varid" is zero).
 Multiple attributes are returned in a vector.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncatt_put}}.
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

# Now illustrate the use of the ncatt_get function by reading them back in
doitfor <- function( nc, var, attname ) {
	av <- ncatt_get( nc, var, attname )
	if( av$hasatt ) {
		print(paste("File",nc$filename,", var",var,"DOES have attribute",
				attname))
		print(paste("Storage mode:",storage.mode(av$value)))
		print("Attribute value:")
		print(av$value)
	} else {
		print(paste("File",nc$filename,", var",var,"does NOT have",
			"attribute", attname))
		}
}
	
nc <- nc_open( filename )
var <- "Data"
doitfor( nc, var, "testatt_dbl" )
doitfor( nc, var, "testatt_float" )
doitfor( nc, var, "testatt_wacko" )
doitfor( nc, 0,   "globalatt_int" )
doitfor( nc, 0,   "globalatt_short" )
doitfor( nc, 0,   "description" )
}
\keyword{utilities}
