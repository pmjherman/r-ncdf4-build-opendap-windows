\name{ncvar_get}
\alias{ncvar_get}
\title{Read data from a netCDF file}
\description{
 Reads data from an existing netCDF file.  
}
\usage{
 ncvar_get(nc, varid=NA, start=NA, count=NA, verbose=FALSE,
 signedbyte=TRUE, collapse_degen=TRUE)
}
\arguments{
 \item{nc}{An object of class \code{ncdf4} (as returned by either 
 function \code{\link[ncdf4]{nc_open}}
 or function \code{\link[ncdf4]{nc_create}}), indicating what file to read from.}
 \item{varid}{What variable to read the data from.  Can be a string with the name
 of the variable or an object of class \code{ncvar4}
 If left unspecified, the function will determine if there
 is only one variable in the file and, if so, read from that.  If left unspecified
 and there are multiple variables in the file, an error is generated.
 This argument can also, optionally, specify the name of a dimension (usually the
 unlimited dimension) in order to read values from a coordinate variable.
 Note this is not usual practice, because the \code{ncdim} object already
 contains all the dimension's values in the field named "vals".  However, it can
 sometimes be faster to turn off this automatic reading of the unlimited dimension's
 values by using \code{nc_open(filename, readunlim=FALSE)}, then read the
 dimension values in later with this function.}
 \item{start}{A vector of indices indicating where to start reading the passed values
 (beginning at 1).  
 The length of this vector must equal the number of dimensions the variable has.
 Order is X-Y-Z-T (i.e., the time dimension is last).  If not specified, reading 
 starts at the beginning of the file (1,1,1,...).}
 \item{count}{A vector of integers indicating the count of values to read
 along each dimension (order is X-Y-Z-T).  
 The length of this vector must equal the number of dimensions the variable has.
 If not specified and the variable does
 NOT have an unlimited dimension, the entire variable is read.  
 As a special case, the value
 "-1" indicates that all entries along that dimension should be read.}
 \item{verbose}{If TRUE, then progress information is printed.}
 \item{signedbyte}{If TRUE (default), then on-disk byte variables are interpreted
 as signed.  This is in accord with the netCDF standard.  If FALSE, then
 on-disk byte variables are interpreted as unsigned.}
 \item{collapse_degen}{If TRUE (the default), then degenerate (length==1) dimensions
 in the returned array are removed.}
}
\references{
 http://dwpierce.com/software
}
\details{
 This routine reads data values from a variable in an existing netCDF file.  The file
 must already have been opened with a call to \code{\link[ncdf4]{nc_open}}.

 Returned values will be in ordinary R double precision if the netCDF 
 variable type is float or double.
 Returned values will be in R's integer storage mode if the netCDF 
 variable type is short or int.
 Returned values will be of character type if the netCDF variable is 
 of character type.

 Values of "NA" are supported; values in the data file that match the
 variable's missing value attribute are automatically converted to "NA" before
 being returned to the user.
 See \code{\link[ncdf4]{ncvar_change_missval}} for more information.

 Data in a netCDF file is conceived as being a multi-dimensional array.
 The number and length of dimensions is determined when the variable is
 created.  The 'start' and 'count' indices that this routine takes indicate
 where the writing starts along each dimension, and the count of values
 along each dimension to write.  Note that the special count value "-1"
 means "all the values along that dimension".

 If the variable in the netCDF file has a scale and/or offset attribute defined, 
 the returned data are automatically and silently scaled and/or offset as requested.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\seealso{ 
 \code{\link[ncdf4]{ncvar_put}}.
}
\examples{
# Start with the simplest example.  If the file only has one variable in it,
# you can read the data as easily as this:
#
nc <- nc_open("salinity.nc")
# NOTE how not specifying varid reads the "only" var in the file
data <- ncvar_get( nc )	
nc_close(nc)

# In this next example we read values from file "writevals.nc", which is created by
# the R code in the example section for function "ncvar_put".  We open the
# file with readunlim=FALSE for potentially faster access, and to illustrate
# (below) how to read in the unlimited dimension values.
#
nc <- nc_open( "writevals.nc", readunlim=FALSE )

print(paste("The file has",nc$nvars,"variables"))

# This illustrates how to read all the data from a variable
v1 <- nc$var[[1]]
data1 <- ncvar_get( nc, v1 )	# by default, reads ALL the data
print(paste("Data for var ",v1$name,":",sep=""))
print(data1)

# This shows how the shape of the read data is preserved
v2 <- nc$var[[2]]
data2 <- ncvar_get( nc, v2 )
print(paste("Var 2 has name",v2$name,"and is of shape",dim(data2),
	". Here are the values:"))
print(data2)

# This illustrates how to read data one timestep at a time.  In this
# example we will elaborately show how to deal with a variable whose
# shape is completely unknown (i.e., how many dimensions, and what their
# sizes are).  We will also, for illustration of a common case, show how
# to read in the values of the time dimension at each timestep.
v3      <- nc$var[[3]]
varsize <- v3$varsize
ndims   <- v3$ndims
nt      <- varsize[ndims]  # Remember timelike dim is always the LAST dimension!
for( i in 1:nt ) {
	# Initialize start and count to read one timestep of the variable.
	start <- rep(1,ndims)	# begin with start=(1,1,1,...,1)
	start[ndims] <- i	# change to start=(1,1,1,...,i) to read timestep i
	count <- varsize	# begin w/count=(nx,ny,nz,...,nt), reads entire var
	count[ndims] <- 1	# change to count=(nx,ny,nz,...,1) to read 1 tstep
	data3 <- ncvar_get( nc, v3, start=start, count=count )

	# Now read in the value of the timelike dimension
	timeval <- ncvar_get( nc, v3$dim[[ndims]]$name, start=i, count=1 )

	print(paste("Data for variable",v3$name,"at timestep",i,
		" (time value=",timeval,v3$dim[[ndims]]$units,"):"))
	print(data3)
	}

nc_close(nc)
}
\keyword{utilities}
