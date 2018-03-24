\name{ncvar_change_missval}
\alias{ncvar_change_missval}
\title{Change the Missing Value For a netCDF Variable}
\description{
 Changes the missing_value attribute for a netCDF variable.
}
\usage{
 ncvar_change_missval( nc, varid, missval ) 
}
\arguments{
 \item{nc}{An object of class \code{ncdf4}, as returned by \code{\link[ncdf4]{nc_open}}\code{(...,write=TRUE)} 
 or \code{\link[ncdf4]{nc_create}}.}
 \item{varid}{Either the name of the variable or an \code{ncvar} object indicating 
 whose missing value will be changed.}
 \item{missval}{The missing value to change to.}
}
\references{
 http://dwpierce.com/software
}
\details{
 Note: this specialty function is only used to change a variable's missing value
 after it has already been defined,
 which is rare.  The proper way to set
 a variable's missing value in the first place is by setting the missing value argument to
 routine \code{\link[ncdf4]{ncvar_def}} appropriately.

 Missing values are special values in netCDF files whose value is to be taken
 as indicating the data is "missing".
 This is a convention, and is indicated by the netCDF variable having an
 attribute named "missing_value" that holds this number.
 This function sets the "missing_value" attribute for a variable.

 R uses a similar concept to indicate missing values, the "NA" value.
 When the ncdf library reads in data set from a pre-existing file, all data values
 that equal that variable's missing value attribute appear to the R code
 as being "NA" values.
 When the R code writes values to a netCDF variable, any "NA" values are
 set to that variable's missing value before being written out.
 This makes the mapping between netCDF's "missing_value" attribute and
 R's "NA" values transparent to the user.

 For this to work, though, the user still has to specify a missing value
 for a variable.  Usually this is specified when the variable is created,
 as a required argument to \code{\link[ncdf4]{ncvar_def}}.
 However, sometimes it is useful to add (or change) a missing value for variable
 that already exists in a disk file.  This function enables that.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncvar_def}}.
}
\examples{
# Make an example netCDF file with a given missing value.  We will
# then change the missing value in the file using ncvar_change_missval

origMissVal <- -1.
dimX   <- ncdim_def( "X", "meters", 1:7 )
varAlt <- ncvar_def( "Altitude", "km", dimX, origMissVal )
ncnew  <- nc_create( "transect.nc", varAlt )
data <- c(10.,2.,NA,1.,7.,NA,8.)
ncvar_put( ncnew, varAlt, data )
nc_close(ncnew)

# At this point, the actual data values in the netCDF
# file will be: 10 2 -1 1 7 -1 8
# because the "NA" values were filled with the missing
# value, -1.  Also, the missing_value attribute of variable
# "varAlt" will be equal to -1.

# Now change the missing value to something else.  Remember
# we have to open the file as writable to be able to change
# the missing value on disk!

newMissVal <- 999.9
nc <- nc_open( "transect.nc", write=TRUE )
varname <- "Altitude"
data <- ncvar_get( nc, varname )  # data now has: 10., 2., NA, 1., 7., NA, 8.
print(data)
ncvar_change_missval( nc, varname, newMissVal )
ncvar_put( nc, varname, data )
nc_close(nc)

# Now, the actual data values in the netCDF file will be:
# 10 2 999.9 1 7 999.9 8
# and the variables "missing_value" attributre will be 999.9

# **NOTE** that we had to explicitly read in the data and write
# it out again in order for the on-disk missing values in the
# data array to change!  The on-disk missing_value attribute for
# the variable is set automatically by this function, but it is 
# up to you whether or not you want to read in all the existing
# data and change the values to the new missing value.
}
\keyword{utilities}
