
#include <stdio.h>
#include <netcdf.h>
#include <string.h>

#include <Rdefines.h>

/* These same values are hard-coded into the R source. Don't change them! */
#define R_NC_TYPE_SHORT 1
#define R_NC_TYPE_INT   2
#define R_NC_TYPE_FLOAT 3
#define R_NC_TYPE_DOUBLE 4
#define R_NC_TYPE_TEXT  5

/*********************************************************************/
void R_nc4_enddef( int *ncid )
{
	int	err;
	err = nc_enddef(*ncid);
	if( err != NC_NOERR ) 
		Rprintf( "Error in R_nc4_enddef: %s\n", 
			nc_strerror(err) );
}

/*********************************************************************/
void R_nc4_sync( int *ncid )
{
	int	err;
	err = nc_sync(*ncid);
	if( err != NC_NOERR ) 
		Rprintf( "Error in R_nc4_sync: %s\n", 
			nc_strerror(err) );
}

/*********************************************************************/
void R_nc4_close( int *ncid )
{
	int	err;
	err = nc_close(*ncid);
	if( err != NC_NOERR ) 
		Rprintf( "Error in R_nc4_close: %s\n", 
			nc_strerror(err) );
}

