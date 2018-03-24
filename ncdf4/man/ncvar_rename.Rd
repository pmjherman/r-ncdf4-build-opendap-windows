\name{ncvar_rename}
\alias{ncvar_rename}
\title{Rename an Existing Variable in a netCDF File}
\description{
 Renames an existing variable that currently is part of a netCDF file
 that is on disk.
}
\usage{
 ncvar_rename( nc, old_varname, new_varname, verbose=FALSE )
}
\arguments{
 \item{nc}{The already-existing netCDF file that we want to manipulate. This must
 be a value of class "ncdf4" returned by a call to \code{\link[ncdf4]{nc_open}}\code{(...,write=TRUE).}}
 \item{old_varname}{The variable in the file that is to be renamed.
 This can be a string with the name of the variable to be renamed,
 or a value of class "ncvar4" returned by a call to \code{ncvar_def()}.}
 \item{new_varname}{A string containing the new name of the variable.}
 \item{verbose}{If true, run verbosely.}
}
\value{
 The updated value of nc that contains the new name.
 This needs to replace the old value of nc in the code.
 I..e, \code{ncid <- ncvar_rename( ncid, ... )}.
}
\references{
 http://dwpierce.com/software
}
\details{
 This call allows you to rename a variable that already exists in a netCDF 
 file.

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
#===================================================================
# PART 1.  MAKE A TEST NETCDF FILE THAT WE WILL MANIPULATE IN PART 2
#===================================================================

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
# PART 2.  RENAME A NEW VARIABLE TO THE FILE
#===========================================================================

#-------------------------------------------------
# Open the existing file we're going to manipulate
#-------------------------------------------------
ncid_old <- nc_open( output_fname, write=TRUE )

old_varname <- 'Temperature'
new_varname <- 'T'

ncid_old <- ncvar_rename( ncid_old, old_varname, new_varname )

print(ncid_old)

nc_close( ncid_old )
}
\keyword{utilities}
