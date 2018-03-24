\name{nc_open}
\alias{nc_open}
\title{Open a netCDF File}
\description{
 Opens an existing netCDF file for reading (or, optionally, writing).
}
\usage{
 nc_open( filename, write=FALSE, readunlim=TRUE, verbose=FALSE, 
 	auto_GMT=TRUE, suppress_dimvals=FALSE )
}
\arguments{
 \item{filename}{Name of the existing netCDF file to be opened.}
 \item{write}{If FALSE (default), then the file is opened read-only.  If TRUE, then writing to the file is allowed.}
 \item{readunlim}{When invoked, this function reads in the values of all dimensions from 
 the associated variables.  This can be slow for a large file with a long unlimited 
 dimension.  If set to FALSE, the values for the unlimited dimension are not 
 automatically read in (they can be read in later, manually, using \code{ncvar_get()}).}
 \item{verbose}{If TRUE, then messages are printed out during execution of this function.}
 \item{auto_GMT}{If TRUE, then GMT files are automatically detected. Does not yet do anything.}
 \item{suppress_dimvals}{If TRUE, then NO dimensional values are automatically read in from the file. (Use
 this if there are so many dimensional values that a out-of-memory error is generated).}
}
\value{
 An object of class \code{ncdf4} that has the fields described above.
}
\references{
 http://dwpierce.com/software
}
\details{
 This routine opens an existing netCDF file for reading (or, if write=TRUE, for writing).
 To create a new netCDF file, use \code{\link[ncdf4]{nc_create}} instead.

 In addition to simply opening the file, information about the file and its contents
 is read in and stored in the returned object, which is of class \code{ncdf4}.
 This class has the following user-accessible fields, all of which are read-only: 1) filename,
 which is a character string holding the name of the file; 2) ndims, which is an
 integer holding the number of dimensions in the file; 3) nvars, which is an integer
 holding the number of the variables in the file that are NOT coordinate variables
 (aka dimensional variables); 4) natts, which is an integer holding the number of
 global attributes; 5) unlimdimid, which is an integer holding the dimension id of
 the unlimited dimension, or -1 if there is none; 6) dim, which is a list of
 objects of class \code{ncdim4}; 7) var, which is a list of objects of class
 \code{ncvar4}; 8) writable, which is TRUE or FALSE, depending on whether the file
 was opened with write=TRUE or write=FALSE.
 
 The concept behind the R interface to a netCDF file is that the \code{ncdf4} object
 returned by this function, as well as the list of \code{ncdim4} objects contained
 in the ncdf object's "dim" list and the \code{ncvar4} objects contained in the
 ncdf object's "var" list, completely describe the netCDF file. I.e., they hold
 the entire contents of the file's metadata.  Therefore, there are no R interfaces
 to the explicit netCDF query functions, such as "nc_inq_nvars" or "nc_inq_natts".
 The upshot is, look in the ncdf4 object or its children to get information about
 the netCDF file.  (Note: the \code{ncdim4} object is described in the help
 file for \code{\link[ncdf4]{ncdim_def}}; the \code{ncvar4} object is described
 in the help file for \code{\link[ncdf4]{ncvar_def}}).

 Missing values: R uses "NA" as a missing value. Netcdf files have various 
 standards for indicating a missing value. The most common is that a variable
 will have an attribute named "_FillValue" indicating the value that should
 be interpreted as a missing value. (For example, the _FillValue attribute might
 have the value of 1.e30, indicating that any data in the netcdf file with
 a value of 1.e30 should be interpreted as a missing value.) 
 If the "_FillValue" attribute is found, then the ncdf4
 package transparently maps all the netcdf file's missing values to NA's; this
 is the most common case.
 (The attribute "missing_value" is also recognized if there is no "_FillValue
 attribute.) If the netcdf file does not have a missing value, then the ncdf4
 package assigns a default missing value of 1.e30 to the netcdcf file so that R NA's,
 which are always possible in the R environment, can be sensibly handled
 in the netcdf file. On rare occasions this can cause problems with non-compliant
 or incorrect netcdf files that implicitly use some particular value,
 for example 9.96921e+36, to indicate a missing value but without 
 setting a proper _FillValue attribute. The best way to fix such netcdf files is to explicitly
 put in the correct _FillValue attribute using an \code{\link[ncdf4]{ncatt_put}} call.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncdim_def}}, \code{\link[ncdf4]{ncvar_def}}, \code{\link[ncdf4]{ncatt_put}}. 
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

# Now open the file and read its data
ncold <- nc_open("states_population.nc")
data <- ncvar_get(ncold)
print("here is the data in the file:")
print(data)
nc_close( ncold )
}
\keyword{utilities}
