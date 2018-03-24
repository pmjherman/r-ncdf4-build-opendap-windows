\name{ncdim_def}
\alias{ncdim_def}
\title{Define a netCDF Dimension}
\description{
 Defines a netCDF dimension.  This dimension initially only exists in memory.
 The dimension is later added to a netCDF variable using a call to \code{ncvar_def()}, 
 and written to disk using \code{nc_create()}.
}
\usage{
 ncdim_def( name, units, vals, unlim=FALSE, 
 	create_dimvar=TRUE, calendar=NA, longname=name )
}
\arguments{
 \item{name}{Name of the dimension to be created (character string).
 The dimension name can 
 optionally have forward slashes in it, in which case the dimension will be defined 
 in the indicated group. For example, a dimension
 named \code{model3/run1/Longitude} will define a group named \code{model3}, with a
 subgroup named \code{run1}, which will hold a dimension named \code{Longitude}. Using
 groups forces a netcdf version 4 file to be written. Note that older software might not
 be able to read netcdf version 4 files.}
 \item{units}{The dimension's units (character string).}
 \item{vals}{The dimension's values (vector of numeric type). If integers are passed, the associated dimensional variable will be integer type; otherwise, it will be double precision.}
 \item{unlim}{If TRUE, this dimension is unlimited.  Unlimited dimensions are convenient for storing, for example, data that extends over time; 
 the time dimension can be made unlimited, and extended as needed.  Or, an unlimited dimension could be the number of stations, and extended
 as more stations come on-line.  Note that in netCDF version 4, multiple dimensions can be unlimited.  In netCDF version 3, there could
 only be one unlimited dimension, typically the time dimension.}
 \item{create_dimvar}{If TRUE, a dimensional variable (aka coordinate variable) will be created for this dimension. Note: if this is set to FALSE, then 'units' must be an empty string.  It is good practice to always leave this as TRUE.}
 \item{calendar}{If set, the specified string will be added as an attribute named "calendar" to the dimension variable.  Used almost exclusively with unlimited time dimensions.  Useful values include "standard" (or "gregorian"), "noleap" (or "365_day"), and "360_day"). }
 \item{longname}{If set, AND create_dimvar is TRUE, 
 then the created dimvar will have a long_name attribute with this value.}
}
\value{
 An object of class \code{ncdim4} that can later be passed to 
 \code{ncvar_def()}.
}
\references{
 http://dwpierce.com/software
}
\details{
 This routine creates a netCDF dimension in memory.  The created dimension can then
 later be passed to the routine \code{ncvar_def()} when defining a variable.

 Note that this interface to the netCDF library by default includes that more than the
 minimum required by the netCDF standard.  I.e., the netCDF standard allows
 dimensions with no units or values.  This call encourages creating dimensions
 that have units and values,
 as it is useful to ensure that all dimensions have units and
 values, and considerably easier to include them in this call than it is
 to add them later.  The units and values are implemented through "dimensional
 variables," which are variables with the same name as the dimension.  By default, these
 dimensional variables are created automatically -- there is no need for the
 user to create them explicitly.  Dimensional variables are standard practice
 in netCDF files.  To suppress the creation of the dimensional variable for
 the dimension, set passed parameter create_dimvar to FALSE.  As a
 check, if create_dimvar is FALSE, you must ALSO pass an empty string ('') as the unit, and the values
 must be simple integers from 1 to the length of the dimension (e.g., 1:10 to
 make a dimension of length 10).
 This empahsizes that without a dimensional variable, a netCDF file cannot
 store a dimension's units or values.

 The dimensional variable is usually created as a double precision floating
 point.  The other possibility is to pass integer values (using \code{as.integer},
 for example), in which case the dimensional variable with be integer.

 The return value of this function is an object of class \code{ncdim4}, which
 describes the newly created dimension.
 The \code{ncdim} object is used for more than just creating a new 
 dimension, however.
 When opening an existing file, function \code{\link[ncdf4]{nc_open}} returns a 
 \code{ncdf4} class object, which itself has a list of \code{ncdim} objects
 that describe all the dimensions in that existing file.

 The \code{ncdim} object has the following fields, which are all read only:
 1) name, which is a character string containing the name of the dimension;
 2) units, which is a character string containing the units for the dimension,
 if there are any (technically speaking, this is the "units" attribute of the
 associated coordinate variable); 3) vals, which is a vector containing the
 dimension's values (i.e., the values of the associated coordinate variable,
 or, if there is none, an integer sequence from 1 to the length of the dimension);
 3) len, which is the length of this dimension; 4) unlim, which is a boolean
 indicating whether or not this is an unlimited dimension; 5) (optional) calendar,
 which is set if and only if the on-disk dimvar had an attribute named
 "calendar" (in which case, it is set to the value of that attribute).
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncvar_def}}, \code{\link[ncdf4]{nc_create}}
}
\examples{
# Define some straightforward dimensions
x <- ncdim_def( "Lon", "degreesE", 0.5:359.5)
y <- ncdim_def( "Lat", "degreesN", as.double(-89:89))
t <- ncdim_def( "Time", "days since 1900-01-01", 1:10, unlim=TRUE)

# Make a variable with those dimensions.  Note order: time is LAST
salinity <- ncvar_def("Salinity",    "ppt",  list(x,y,t), 1.e30 )

# Create a netCDF file with this variable
ncnew <- nc_create( "salinity.nc", salinity )

nc_close(ncnew)

# Now, illustrate some manipulations of the ncdim object.
filename <- "salinity.nc"
nc <- nc_open( filename )
print(paste("File",filename,"contains",nc$ndims,"dimensions"))
for( i in 1:nc$ndims ) {
	print(paste("Here is information about dimension number",i,":"))
	d <- nc$dim[[i]]
	print(paste("    Name  :",d$name))
	print(paste("    Units :",d$units))
	print(paste("    Length:",d$len))
	print("    Values:")
	print(d$vals)
	print(paste("    Unlimited:",d$unlim))
	}
}
\note{It is good practice, but not necessary, to pass the dimension's values 
 to this routine when the dimension is created.  It is also possible to 
 write them later with a call to 'ncvar_put', using as the dimension
 name as the 'varid' in the call.  This is useful when creating large variables
 with long unlimited dimensions; it can take a long time to write out the
 unlimited dimension's values.  In this case, it can be more efficient to
 step through the file, writing one timestep at a time, and write that 
 timestep's dimensional value at the same time.
}
\keyword{utilities}
