/*   #define DEBUG */

#include <stdio.h>
/*  #include <values.h> */

#include <R.h>
#include <Rinternals.h>
#include "netcdf.h"

/*********************************************************************/
void R_nc4_def_grp( int *parent_ncid, char **grpname, int *new_ncid, int *retval )
{
	*retval = nc_def_grp(*parent_ncid, grpname[0], new_ncid );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_def_grp: %s\n", 
			nc_strerror(*retval) );
		Rprintf( "Name of group that the error occurred on: \"%s\"\n", 
			grpname[0] );
		}
}

