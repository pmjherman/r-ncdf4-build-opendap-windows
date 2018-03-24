\name{ncvar_add}
\alias{ncvar_add}
\title{Add New netCDF Variable to Existing File}
\description{
 Special purpose routine for adding a new variable to a netCDF file
 that already exists on disk.
}
\usage{
 ncvar_add( nc, v, verbose=FALSE, indefine=FALSE )
}
\arguments{
 \item{nc}{The already-existing netCDF file we want to add a new variable to. This must
 be a value of class "ncdf4" returned by a call to \code{\link[ncdf4]{nc_open}}\code{(...,write=TRUE).}}
 \item{v}{The variable to be added to the file.  This must be 
 a value of class "ncvar4" returned by a call to \code{\link[ncdf4]{ncvar_def}}.}
 \item{verbose}{If true, prints diagnostic messages.}
 \item{indefine}{If true, the file is assumed to already be in define mode.}
}
\value{
 A handle to the netCDF file that describes the newly modified file.
 This is an object of type 'ncdf', the same as returned by \code{\link[ncdf4]{nc_open}}
 or \code{\link[ncdf4]{nc_create}}.
}
\references{
 http://dwpierce.com/software
}
\details{
 There are two cases in which you might want to add a variable to a netCDF
 file.  The first, and most common way, is when you are creating a 
 new netCDF file.
 Usually when you create a netCDF file, you specify what variables you want
 the file to contain.
 This is the method
 most users will use to make netCDF files.  To do this, do NOT use this
 routine; instead, pass a list of the variables you wish to have created in the output
 file to routine \code{\link[ncdf4]{nc_create}}.

 The second, less common, case is when you already have an existing netCDF
 file on disk and wish to add a new variable to it.  In that case, use
 this routine.  First define the variable you want to add to the existing
 file using routine \code{\link[ncdf4]{ncvar_def}}; then add it to the already-existing
 and opened (for writing) netCDF file using this routine.  (This routine
 automatically creates any additional dimensions that are needed in the
 output file to handle the new variable.)

 NOTE that the return value of this routine should replace the old netCDF
 file handle that you were using.  This newly returned value reflects the
 modifications to the file that were accomplished by calling this routine.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncdim_def}}, \code{\link[ncdf4]{nc_create}}, 
 \code{\link[ncdf4]{ncvar_def}}.
}
\examples{
#===========================================================================
# PART 1.  MAKE A TEST NETCDF FILE THAT WE WILL ADD A VARIABLE TO IN PART 2.
#===========================================================================

#----------------
# Make dimensions
#----------------
xvals <- 1:360
yvals <- -90:90

nx <- length(xvals)
ny <- length(yvals)

xdim <- ncdim_def( 'Lon', 'degreesE', xvals )
ydim <- ncdim_def( 'Lat', 'degreesE', yvals )
tdim <- ncdim_def( 'Time', 'days since 1900-01-01', 0, unlim=TRUE )

#---------
# Make var
#---------
mv <- 1.e30     # missing value
var_temp <- ncvar_def( 'Temperature', 'K', list(xdim,ydim,tdim), mv )

#---------------------
# Make new output file
#---------------------
output_fname <- 'test_real3d.nc'
ncid_new <- nc_create( output_fname, list(var_temp))

#-------------------------------
# Put some test data in the file
#-------------------------------
data_temp <- array(0.,dim=c(nx,ny,1))
for( j in 1:ny )
for( i in 1:nx )
        data_temp[i,j,1] <- sin(i/10)*sin(j/10)

ncvar_put( ncid_new, var_temp, data_temp, start=c(1,1,1), count=c(nx,ny,1))

#--------------------------
# Close our new output file
#--------------------------
nc_close( ncid_new )

#===========================================================================
# PART 2.  ADD A NEW VARIABLE TO THE FILE
#===========================================================================

#---------------------------------------------------
# Open the existing file we're going to add a var to
#---------------------------------------------------
ncid_old <- nc_open( output_fname, write=TRUE )

#------------------------------------------------------------------
# Make a NEW variable to put into the file.  Have this new variable
# use the same dimensions already in the file
#------------------------------------------------------------------
xdim2 <- ncid_old$dim[['Lon']]
ydim2 <- ncid_old$dim[['Lat']]
tdim2 <- ncid_old$dim[['Time']]
mv2 <- 1.e30
var_q <- ncvar_def( 'Humidity', 'g/kg', list(xdim2,ydim2,tdim2), mv2 )

ncid_old <- ncvar_add( ncid_old, var_q )	# NOTE this returns a modified netcdf file handle 

#-------------------------------------------------------------------
# Make a DIFFERENT new var that will be added to the file.  This var
# uses a dim that does NOT already exist in the file.
#-------------------------------------------------------------------
zdim   <- ncdim_def( 'Level', 'hPa', seq(1000,100,by=-100))
var_cf <- ncvar_def( 'CloudFraction', 'percent', list(xdim2,ydim2,zdim,tdim2), mv2 )

ncid_old <- ncvar_add( ncid_old, var_cf )

print(ncid_old)

nc_close( ncid_old )
}
\keyword{utilities}
