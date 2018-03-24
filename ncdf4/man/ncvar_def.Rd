\name{ncvar_def}
\alias{ncvar_def}
\title{Define a netCDF Variable}
\description{
 Defines a netCDF variable.  This variable initially only exists in memory.
 It is later written to disk using \code{nc_create()}.
}
\usage{
 ncvar_def( name, units, dim, missval, longname=name, prec="float", 
 shuffle=FALSE, compression=NA, chunksizes=NA, verbose=FALSE )
}
\arguments{
 \item{name}{Name of the variable to be created (character string).  
 The name can 
 optionally have forward slashes in it, in which case the variable will be defined 
 in the indicated group. For example, a variable
 named \code{model3/run1/Temperature} will define a group named \code{model3}, with a
 subgroup named \code{run1}, which will hold a variable named \code{Temperature}. Using
 groups forces a netcdf version 4 file to be written. Note that older software might not
 be able to read netcdf version 4 files.}
 \item{units}{The variable's units (character string). Or, pass a zero length string ('') to have no units attribute.}
 \item{dim}{The variable's dimension(s) (one or a list of "ncdim4" class objects, as returned by \code{ncdim_def}).  To create a variable with NO dimensions, pass an empty list ("list()"). }
 \item{missval}{The variable's missing value. If NO missing value is desired, pass a NULL, or omit this argument entirely. If a NaN missing value is desired, pass an NA.}
 \item{longname}{Optional longer name for the variable, which is assigned to the variable's "long_name" attribute. For example, a variable named "TS" might have the longname "Surface Temperature"}
 \item{prec}{Precision of the created variable.  Valid options: 'short' 'integer' 'float' 'double' 'char'.}
 \item{shuffle}{Turns on (if TRUE) or off (if FALSE, the default) the shuffle filter.  According to 
 netcdf docs, turning the shuffle filter on can improve compression for integer variables.
 Turning the shuffle filter on forces the created file to be in netcdf version 4 format,
 which will not be compatible with older software that only reads netcdf version 3 files.}
 \item{compression}{If set to an integer between 1 (least compression) and 9 (most compression), this
 enables compression for the variable as it is written to the file.
 Turning compression on forces the created file to be in netcdf version 4 format,
 which will not be compatible with older software that only reads netcdf version 3 files.}
 \item{chunksizes}{If set, this must be a vector of integers with a length equal to the number
 of dimensions in the variable. When data from this variable is written to the file, it
 will be buffered in blocks as indicated by the chunksize. The order of dimensions
 in this vector is the standard R ordering of XYZT. In some instances, setting a 
 chunksize that reflects how the variable's data will be read or written can greatly
 reduce read or write times. See the netcdf documentation for more detail on how to set this
 parameter. Enabling this feature forces the created file to be in netcdf version 4 format,
  which will not be compatible with older software that only reads netcdf version 3 files.}
 \item{verbose}{Print debugging information.}
}
\value{
 An object of class \code{ncvar4} that can later be passed to 
 \code{nc_create()}.
}
\references{
 http://dwpierce.com/software
}
\details{
 This routine creates a netCDF variable in memory.  The variable can then
 be passed to the routine \code{nc_create} when writing a file to disk.

 Note that this interface to the netCDF library includes more than the
 minimum required by the netCDF standard.  I.e., the netCDF standard allows
 variables with no units or missing values.  This call requires units and 
 a missing value,
 as it is useful to ensure that all variables have units and
 missing values, and considerably easier to include them in this call than it is
 to add them later.  The units and missing value are implemented through attributes
 to the variable, named "units" and "missing_value", respectively.
 This is standard practice in netCDF files.

 After a variable is defined with this call, and created on disk using
 \code{\link[ncdf4]{nc_create}}, then data values for the variable can be written
 to disk using \code{\link[ncdf4]{ncvar_put}}.

 This function returns a \code{ncvar} object, which describes the newly-created
 variable.  However, the \code{ncvar} object is used for more than just
 creating new variables.  The function \code{\link[ncdf4]{nc_open}} returns a \code{ncdf4}
 class object that itself contains a list of \code{ncvar4} objects that describe
 the variables in an existing, on-disk netCDF file.  (Note that coordinate
 variables are NOT included in this list.  Attributes of the coordinate variables
 are kept in the \code{ncdim4} class object instead.)

 The \code{ncvar4} class object has the following fields, which are all read-only:
 1) name, which is a character string containing the name of the variable;
 2) units, which is a character string containing the contents of the
 variable's "units" attribute; 3) missval, which contains the contents of the
 variable's "missing_value" attribute; 4) longname, which is the 
 contents of the variable's "long_name" attribute, or defaults to the name
 of the variable if there is no "long_name" attribute; 5) ndims, which is the
 number of dimensions this variable has; 6) dim, which is a list of objects of
 class "ncdim4" (see \code{\link[ncdf4]{ncdim_def}}), and describe this
 variable's dimensions; 7) unlim, which is TRUE if this variable has an unlimited
 dimension and FALSE otherwise; 8) varsize, which is a convenience array
 that gives the shape of the variable (in XYZT ordering).

 Note that the missval attribute does not need to be used much in R, because
 R's special value NA is fully supported.  I.e., when data is read in from an
 existing file, any values equal to the "missing" value are set to NA.  When
 data is written out, any NAs are set equal to the missing value.  If not
 explicitly set by the user, a default value of 1.e30 is used for the missing
 value.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncdim_def}}, \code{\link[ncdf4]{nc_create}}, 
 \code{\link[ncdf4]{ncvar_put}}.
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

# Add source info metadata to file
ncatt_put( ncnew, 0, "source", "Census 2000 from census bureau web site")

nc_close(ncnew)

# Now illustrate some manipulations of the var.ncdf object
filename <- "states_population.nc"
nc <- nc_open(filename)
print(paste("File",nc$filename,"contains",nc$nvars,"variables"))
for( i in 1:nc$nvars ) {
	v <- nc$var[[i]]
	print(paste("Here is information on variable number",i))
	print(paste("   Name: ",v$name))
	print(paste("   Units:",v$units))
	print(paste("   Missing value:",v$missval))
	print(paste("   # dimensions :",v$ndims))
	print(paste("   Variable size:",v$varsize))
	}

# Illustrate creating variables of various types.  You will find
# that the type of the missing_value attribute automatically follows
# the type of the variable.
dimt <- ncdim_def( "Time", "days", 1:3 ) 
missval <- -1
varShort <- ncvar_def( "varShort", "meters", dimt, missval, prec="short")
varInt   <- ncvar_def( "varInt",   "meters", dimt, missval, prec="integer")
varFloat <- ncvar_def( "varFloat", "meters", dimt, missval, prec="single")
varDouble<- ncvar_def( "varDouble","meters", dimt, missval, prec="double")
nctypes <- nc_create("vartypes.nc", list(varShort,varInt,varFloat,varDouble) )
nc_close(nctypes)
}
\keyword{utilities}
