\name{ncvar_put}
\alias{ncvar_put}
\title{Write data to a netCDF file}
\description{
 Writes data to an existing netCDF file.  The variable to be written to
 must already exist on disk (i.e., you must call either
 \code{\link[ncdf4]{nc_create}} or \code{\link[ncdf4]{nc_open}}
 before calling this function).
}
\usage{
 ncvar_put( nc, varid, vals, start=NA, count=NA, verbose=FALSE ) 
}
\arguments{
 \item{nc}{An object of class \code{ncdf4} (as returned by either function
 \code{\link[ncdf4]{nc_open}}
 or \code{\link[ncdf4]{nc_create}}), indicating what file to write to.}
 \item{varid}{What variable to write the data to.  Can be a string with the name
 of the variable or an object of class \code{ncvar4}, as
 returned by \code{\link[ncdf4]{ncvar_def}} or  \code{\link[ncdf4]{nc_open}}.}
 \item{vals}{The values to be written.}
 \item{start}{A vector of indices indicating where to start writing the passed values
 (starting at 1).  
 The length of this vector must equal the number of dimensions the variable has.
 Order is X-Y-Z-T (i.e., the time dimension is last).  If not specified, writing 
 starts at the beginning of the file (1,1,1,...).}
 \item{count}{A vector of integers indicating the count of values to write
 along each dimension (order is X-Y-Z-T).  
 The length of this vector must equal the number of dimensions the variable has.
 If not specified and the variable does
 NOT have an unlimited dimension, the entire variable is written.  If the variable
 has an unlimited dimension, this argument must be specified.  As a special case, the value
 "-1" indicates that all entries along that dimension should be written.}
 \item{verbose}{If true, prints information while executing.}
}
\references{
 http://dwpierce.com/software
}
\details{
 This routine writes data values to a variable in a netCDF file.  The file
 should have either been created with \code{\link[ncdf4]{nc_create}}, or 
 opened with \code{\link[ncdf4]{nc_open}} called with parameter \code{write=TRUE}.

 Note that the data type (i.e., precision) of the values written to the file is determined when the
 variable is created; in particular, it does not matter what type you pass
 to this function to be written.  In other words, if the variable was created
 with type 'integer', passing double precision values to this routine will 
 still result in integer values being written to disk.

 Values of "NA" are supported; they are converted to the netCDF variable's
 missing value attribute before being written.  See \code{\link[ncdf4]{ncvar_change_missval}}
 for more information.

 Data in a netCDF file is conceived as being a multi-dimensional array.
 The number and length of dimensions is determined when the variable is
 created.  The 'start' and 'count' indices that this routine takes indicate
 where the writing starts along each dimension, and the count of values
 along each dimension to write.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncdim_def}}, 
 \code{\link[ncdf4]{nc_create}}, 
 \code{\link[ncdf4]{nc_open}}, 
 \code{\link[ncdf4]{ncvar_get}}.
}
\examples{
#--------------------------------
# Make a few dimensions we can use
#--------------------------------
nx <- 3
ny <- 4
nt <- 5
xvals <- (1:nx)*100.
dimX <- ncdim_def( "X", "meters", xvals )
dimY <- ncdim_def( "Y", "meters", (1:ny)*100. )
dimT <- ncdim_def( "Time", "seconds", (1:nt)/100., unlim=TRUE )

#-------------------------------------------------------------------
# Make varables of various dimensionality, for illustration purposes
#-------------------------------------------------------------------
mv <- 1.e30		# missing value to use
var1d <- ncvar_def( "var1d", "units", dimX, mv )
var2d <- ncvar_def( "var2d", "units", list(dimX,dimY), mv )
var3d <- ncvar_def( "var3d", "units", list(dimX,dimY,dimT), mv )

#---------------------
# Create the test file
#---------------------
nc <- nc_create( "writevals.nc", list(var1d,var2d,var3d) )

#----------------------------
# Write some data to the file
#----------------------------
data1d <- runif(nx)
ncvar_put( nc, var1d, data1d )	# no start or count: write all values
ncvar_put( nc, var1d, 27.5, start=3, count=1 ) # Write a value to the third slot

data2d <- runif(nx*ny)
ncvar_put( nc, var2d, data2d )	# no start or count: write all values

#--------------------------------
# Write a 1-d slice to the 2d var
#--------------------------------
ncvar_put( nc, var2d, data1d, start=c(1,2), count=c(nx,1) )	

#--------------------------------------------------------------
# Note how "-1" in the count means "the whole dimension length", 
# which equals nx in this case
#--------------------------------------------------------------
ncvar_put( nc, var2d, data1d, start=c(1,3), count=c(-1,1) )	

#-------------------------------------------------------------------------------
# The 3-d variable has an unlimited dimension.  We will loop over the timesteps,
# writing one 2-d slice per timestep.
#-------------------------------------------------------------------------------
for( i in 1:nt) 
	ncvar_put( nc, var3d, data2d, start=c(1,1,i), count=c(-1,-1,1) )

nc_close(nc)

#----------------------------------------------------------------------
# Illustrate creating a character type variable
#----------------------------------------------------------------------
cnames   <- c("red", "orange", "green", "yellow", "puce", "colorwithverylongname" )
nstrings <- length(cnames)

#--------------------------------------------------------------
# Make dimensions. Setting "dimnchar" to have a length of 12
# means that the maximum color name
# length can be 12.  Longer names will be truncated to this.
# We don't need dimvars for this example.
#--------------------------------------------------------------
dimnchar   <- ncdim_def("nchar",   "", 1:12, create_dimvar=FALSE )
dimcolorno <- ncdim_def("colorno", "", 1:nstrings, create_dimvar=FALSE )

#------------------------------------------------------------------------
# NOTE in the following call that units is set to the empty string (""),
# which suppresses creation of a units attribute, and the missing value
# is entirely omitted, which suppresses creation of the missing value att
#------------------------------------------------------------------------
varcolors  <- ncvar_def("colors", "", list(dimnchar, dimcolorno), 
			prec="char" )

ncid <- nc_create( "colornames.nc", list(varcolors) )

ncvar_put( ncid, "colors", cnames, verbose=TRUE )

nc_close( ncid )
}
\keyword{utilities}
