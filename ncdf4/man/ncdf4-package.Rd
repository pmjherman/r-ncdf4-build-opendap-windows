\name{ncdf4-package}
\alias{ncdf4-package}
\title{Read, write, and create netCDF files (including version 4 format)}
\description{
 Read from or write to existing netCDF format files, or create new ones.
}
\references{
 http://dwpierce.com/software
}
\details{
 More information on this package, including detailed installation instructions,
 can be found at http://dwpierce.com/software.

 The netCDF data file format from Unidata is a platform-independent, binary
 file that also contains metadata describing the contents and format of the data
 in the file.  Version 4 of the netcdf library stores data in HDF5 format files;
 earlier versions stored data in a custom format.
 The R package ncdf4 can read either format.

 NetCDF files contain one or more variables, which are usually
 structured as regular N-dimensional arrays.  For example, you might
 have a variable named "Temperature" that is a function of longitude,
 latitude, and height. NetCDF files also contain dimensions, which 
 describe the extent of the variables' arrays.  In our Temperature example,
 the dimensions are "longitude", "latitude", and "height". Data can be read from
 or written to variables in arbitrary hyperslabs (for example, you can 
 read or write all the Temperature values at a given height, or at a given
 latitude).  
 
 The R package 'ncdf4' allows
 reading from, writing to, and creation of netCDF files, either netCDF version
 3 or (optionally) netCDF version 4.  If you choose to create version 4
 output files, be aware that older netcdf software might only be
 able to read version 3 files.

 Note that both the
 netCDF library and the HDF5 library must already be installed on your machine for 
 this R interface to the library to work.

 If you are new to netCDF files, they can be a little overwhelming,
 so here is a brief sketch of what documentation you need to read next.  
 
 If you want to READ data 
 from an already-existing netCDF file, first call \code{\link[ncdf4]{nc_open}} 
 to open the file, then call \code{\link[ncdf4]{ncvar_get}} to read the data
 from a variable in the file.  
 
 If you want to WRITE data to a new netCDF file, the procedure is
 to first define the dimensions your data array has, then define the
 variable, then create the file.  So, first call
 \code{\link[ncdf4]{ncdim_def}} to define the dimensions that your data
 exists along (for example, latitude, longitude, and time). Then call
 \code{\link[ncdf4]{ncvar_def}} to define a variable 
 that uses those dimensions, and will hold your data.
 Then call \code{\link[ncdf4]{nc_create}}
 to create the netCDF file.
 Finally, call \code{\link[ncdf4]{ncvar_put}}
 to write your data to the newly created netCDF file, 
 and \code{\link[ncdf4]{nc_close}} when you are done.

 Not all features of netcdf-4 are supported yet.
 This version supports compression, chunking, groups, vlen strings, and multiple
 unlimited dimensions.
 User-defined types and vlen data arrays (variable-length arrays) are not supported yet.
}
\seealso{ 
 \code{\link[ncdf4]{print.ncdf4}},
 \code{\link[ncdf4]{nc_open}},
 \code{\link[ncdf4]{nc_close}}, 
 \code{\link[ncdf4]{nc_create}},
 \code{\link[ncdf4]{ncdim_def}},
 \code{\link[ncdf4]{ncvar_def}},
 \code{\link[ncdf4]{ncvar_get}},
 \code{\link[ncdf4]{ncvar_put}},
 \code{\link[ncdf4]{ncvar_change_missval}},
 \code{\link[ncdf4]{ncatt_get}},
 \code{\link[ncdf4]{ncatt_put}}, 
 \code{\link[ncdf4]{nc_sync}},
 \code{\link[ncdf4]{nc_redef}}.
}
\author{David W. Pierce \email{dpierce@ucsd.edu}}
\keyword{utilities}

