
#include <stdio.h>
#include <netcdf.h>
#include <string.h>
#include <stdlib.h>

#include <Rdefines.h>

/* These same values are hard-coded into the R source. Don't change them! 
 * Note that they are not the same as values defined in the netcdf
 * library headers, since I don't want my R code to depend on those 
 */
#define R_NC_TYPE_SHORT  1
#define R_NC_TYPE_INT    2
#define R_NC_TYPE_FLOAT  3
#define R_NC_TYPE_DOUBLE 4
#define R_NC_TYPE_TEXT   5
#define R_NC_TYPE_BYTE   6

/* New types in netcdf version 4 */
#define R_NC_TYPE_UBYTE		7
#define R_NC_TYPE_USHORT	8
#define R_NC_TYPE_UINT		9
#define R_NC_TYPE_INT64		10
#define R_NC_TYPE_UINT64	11
#define R_NC_TYPE_STRING	12

void R_nc4_inq_varid_hier( int *ncid, char **varname, int *returned_grpid, int *returned_varid );

/*********************************************************************
 * Converts from type "nc_type" to an integer as defined in the beginning 
 * of this file.  We do NOT use the raw nc_type integers because then the
 * R code would have a dependency on the arbitrary values in the netcdf 
 * header files!
 */
int R_nc4_nctype_to_Rtypecode( nc_type nct )
{
	if( nct == NC_CHAR )
		return(R_NC_TYPE_TEXT);
	else if( nct == NC_SHORT )
		return(R_NC_TYPE_SHORT);
	else if( nct == NC_INT )
		return(R_NC_TYPE_INT);
	else if( nct == NC_FLOAT )
		return(R_NC_TYPE_FLOAT);
	else if( nct == NC_DOUBLE )
		return(R_NC_TYPE_DOUBLE);
	else if( nct == NC_BYTE )
		return(R_NC_TYPE_BYTE);

	else if( nct == NC_UBYTE )
		return(R_NC_TYPE_UBYTE);
	else if( nct == NC_USHORT )
		return(R_NC_TYPE_USHORT);
	else if( nct == NC_UINT )
		return(R_NC_TYPE_UINT);
	else if( nct == NC_INT64 )
		return(R_NC_TYPE_INT64);
	else if( nct == NC_UINT64 )
		return(R_NC_TYPE_UINT64);
	else if( nct == NC_STRING )
		return(R_NC_TYPE_STRING);
	else
		return(-1);
}

/*********************************************************************/
/* Returns a vector of dim sizes for the variable.
 * 'retval' is 0 for no error, or -1 for an error.
 */
void R_nc4_varsize( int *ncid, int *varid, int *ndims, int *varsize, int *retval )
{
	int 	i, err, dimid[NC_MAX_DIMS];
	size_t	dimlen;

	*retval = 0;

	/* Get ndims */
	err = nc_inq_varndims( *ncid, *varid, ndims );
	if( err != NC_NOERR ) {
		Rprintf( "Error in R_nc4_varsize on nc_inq_varndims call: %s\n", 
			nc_strerror(err) );
		*retval = -1;
		return;
		}

	/* Get dimids */
	err = nc_inq_vardimid( *ncid, *varid, dimid );
	if( err != NC_NOERR ) {
		Rprintf( "Error in R_nc4_varsize on nc_inq_vardimid call: %s\n",
			nc_strerror(err) );
		*retval = -1;
		return;
		}

	/* Get size of each dim */
	for( i=0; i<(*ndims); i++ ) {
		err = nc_inq_dimlen( *ncid, dimid[i], &dimlen );
		if( err != NC_NOERR ) {
			Rprintf( "Error in R_nc4_varsize on nc_inq_dimlen call: %s\n",
				nc_strerror(err) );
			*retval = -1;
			return;
			}
		varsize[i] = (int)dimlen;
		}
}

/*********************************************************************/
/* Returns 1 if this var has an unlimited dimension, 0 otherwise */
void R_nc4_inq_varunlim( int *ncid, int *varid, int *isunlim, int *retval )
{
	int	i, ndims, unlimdimid, dimids[MAX_NC_DIMS];

	/* Get the unlimited dim id */
	*retval = nc_inq_unlimdim( *ncid, &unlimdimid );
	if( *retval != NC_NOERR ) {
		Rprintf(  "Error in R_nc4_inq_varunlim while getting unlimdimid: %s\n", 
			nc_strerror(*retval) );
		return;
		}

	/* Get this var's ndims */
	*retval = nc_inq_varndims( *ncid, *varid, &ndims );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_inq_varunlim while getting ndims: %s\n", 
			nc_strerror(*retval) );
		Rprintf( "Using ncid=%d and varid=%d\n", 
			*ncid, *varid );
		return;
		}

	/* Get this var's dims */
	*retval = nc_inq_vardimid( *ncid, *varid, dimids );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_inq_varunlim while getting dimids: %s\n", 
			nc_strerror(*retval) );
		return;
		}

	*isunlim = 0;
	for( i=0; i<ndims; i++ )
		if( dimids[i] == unlimdimid ) {
			*isunlim = 1;
			break;
			}
}

/*********************************************************************/
/* Note that space for all returned values MUST already be declared! */
void R_nc4_inq_var( int *ncid, int *varid, char **varname, 
	int *type, int *ndims, int *dimids, int *natts, int *precint, int *retval )
{
	nc_type nct;

	*retval = nc_inq_var(*ncid, *varid, varname[0], &nct,  
		ndims, dimids, natts );
	*type = (int)nct;

	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_inq_var: %s\n", 
			nc_strerror(*retval) );

	*precint = R_nc4_nctype_to_Rtypecode(nct);
}

/*********************************************************************/
void R_nc4_inq_vartype( int *ncid, int *varid, int *precint, int *retval )
{
	nc_type nct;

	*retval = nc_inq_vartype( *ncid, *varid, &nct );

	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_inq_var: %s\n", 
			nc_strerror(*retval) );

	*precint = R_nc4_nctype_to_Rtypecode(nct);
}

/*********************************************************************/
/* Note that space for returned value MUST already be declared! */
void R_nc4_inq_varname( int *ncid, int *varid, char **varname, int *retval )
{
	*retval = nc_inq_varname(*ncid, *varid, varname[0] );
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_inq_varname: %s\n", 
			nc_strerror(*retval) );
}

/*********************************************************************/
void R_nc4_inq_varndims( int *ncid, int *varid, int *ndims, int *retval )
{
	*retval = nc_inq_varndims(*ncid, *varid, ndims );
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_inq_varndims: %s\n", 
			nc_strerror(*retval) );
}

/*********************************************************************/
/* Inputs:
 *	sx_fixmiss  : is 1 if we want to fix the missing values in this
 *  		routine, and 0 otherwise
 * 	sx_imvstate : 0=var has no missing value; 1=var has a NA for
 * 		the missing value; 2=var has a valid, non-NA missing value
 *
 * Returns a list with elements:
 *	$error	: 0 for success, -1 for error
 *	$data   : array of integer values read in from the netcdf file
 */
SEXP Rsx_nc4_get_vara_double( SEXP sx_ncid, SEXP sx_varid, SEXP sx_start, SEXP sx_count, 
	SEXP sx_fixmiss, SEXP sx_imvstate, SEXP sx_missval )
{
	SEXP	sx_retval, sx_retnames, sx_reterr, sx_retdata;
	int	ncid, varid, i, err, ndims, fixmiss, imvstate, scalar_var;
	double	*p_data, missval, mvtol;
	size_t	s_start[MAX_NC_DIMS], s_count[MAX_NC_DIMS], tot_size;
	char	vn[2048];

	/* Make space for our returned list, which will have
	 * two elements, named $error and $data
	 */
	PROTECT( sx_retval = allocVector( VECSXP, 2 ));	/* 2 elements in returned list */

	/* Set names for the returned list */
	PROTECT( sx_retnames = allocVector( STRSXP, 2 ));
	SET_STRING_ELT( sx_retnames, 0, mkChar("error") );
	SET_STRING_ELT( sx_retnames, 1, mkChar("data" ) );
	setAttrib( sx_retval, R_NamesSymbol, sx_retnames );
	UNPROTECT(1);     				/* done with sx_retnames */

	/* Set provisional 'no error' retval */
	PROTECT(sx_reterr = allocVector( INTSXP, 1 ));
	INTEGER(sx_reterr)[0] = 0;

	ncid  	   	= INTEGER(sx_ncid    )[0];
	varid 	   	= INTEGER(sx_varid   )[0];
	fixmiss		= INTEGER(sx_fixmiss )[0];
	imvstate	= INTEGER(sx_imvstate)[0];
	missval		= REAL   (sx_missval )[0];

	/* Get number of dimensions in this variable */
	err = nc_inq_varndims( ncid, varid, &ndims );
	if( err != NC_NOERR ) {
		Rprintf( "Error in R_nc4_get_vara_double while getting ndims: %s\n", 
			nc_strerror(err) );
		INTEGER(sx_reterr)[0] = -1;
		SET_VECTOR_ELT( sx_retval, 0, sx_reterr );
		UNPROTECT(2);
		return( sx_retval );
		}

	/* Sanity check -- number of start and count elements must
	   match number of dims in the var
	*/
	scalar_var = ((ndims==0) && (GET_LENGTH(sx_start)==1) && (INTEGER(sx_start)[0]==0) && (INTEGER(sx_count)[0]==1));
	if( (!scalar_var) && (ndims != GET_LENGTH(sx_start))) {
		Rprintf( "Error in R_nc4_get_vara_double: I think var has %d dimensions, but passed start array is length %d. They must be the same!\n",
			ndims, GET_LENGTH(sx_start) );
		INTEGER(sx_reterr)[0] = -1;
		SET_VECTOR_ELT( sx_retval, 0, sx_reterr );
		UNPROTECT(2);
		return( sx_retval );
		}

	if( (!scalar_var) && (ndims != GET_LENGTH(sx_count))) {
		Rprintf( "Error in R_nc4_get_vara_double: I think var has %d dimensions, but passed count array is length %d. They must be the same!\n",
			ndims, GET_LENGTH(sx_count) );
		INTEGER(sx_reterr)[0] = -1;
		SET_VECTOR_ELT( sx_retval, 0, sx_reterr );
		UNPROTECT(2);
		return( sx_retval );
		}

	/* Get total number of elements we will be reading so we can 
	 * allocate R space
	 */
	tot_size = 1L;
	for( i=0; i<ndims; i++ ) {
		s_start[i] = (size_t)(INTEGER(sx_start)[i]);
		s_count[i] = (size_t)(INTEGER(sx_count)[i]);
		tot_size *= s_count[i];
		}
		
	/* Make space for the returned data */
	PROTECT( sx_retdata = allocVector(REALSXP, tot_size));
	p_data = REAL( sx_retdata );

	/* Actually read in the data now */
	err = nc_get_vara_double( ncid, varid, s_start, s_count, p_data );
	if( err != NC_NOERR ) {
		nc_inq_varname( ncid, varid, vn );
		Rprintf( "Error in Rsx_nc4_get_vara_double: %s\n", 
			nc_strerror( err ) );
		Rprintf( "Var: %s  Ndims: %d   Start: ", vn, ndims );
		for( i=0; i<ndims; i++ ) {
			Rprintf( "%u", (unsigned int)s_start[i] );
			if( i < ndims-1 )
				Rprintf( "," );
			}
		Rprintf( " " );
		Rprintf( "Count: " );
		for( i=0; i<ndims; i++ ) {
			Rprintf( "%u", (unsigned int)s_count[i] );
			if( i < ndims-1 )
				Rprintf( "," );
			}
		Rprintf( "\n" );
		INTEGER(sx_reterr)[0] = -1;
		SET_VECTOR_ELT( sx_retval, 0, sx_reterr );
		UNPROTECT(3);
		return( sx_retval );
		}

	/* We only do this for imvstate==2, which is a valid, non-NA
	 * missing value. If the variable has NO missing value
	 * (imvstate==0) or a NA missing value (imvstate==1) then
	 * we skip this block.
	 */
	if( (fixmiss == 1) && (imvstate == 2)) {
		if( missval == 0.0 ) 
			mvtol = 1.e-10;	/* arbitrary -- what are you supposed to use?? */
		else
			mvtol = fabs( missval ) * 1.e-5;

		for( i=0L; i<tot_size; i++ ) {
			if( fabs( p_data[i] - missval ) < mvtol )
				p_data[i] = NA_REAL;
			}
		}

	/* OK we can successfully return now */
	SET_VECTOR_ELT( sx_retval, 0, sx_reterr  );
	SET_VECTOR_ELT( sx_retval, 1, sx_retdata );
	UNPROTECT(3);
	return( sx_retval );
}

/*********************************************************************/
/* Input value byte_style is 1 for signed, 2 for unsigned
 * 
 * Returns a list with elements:
 *	$error	: 0 for success, -1 for error
 *	$data   : array of integer values read in from the netcdf file
 */
SEXP Rsx_nc4_get_vara_int( SEXP sx_ncid, SEXP sx_varid, SEXP sx_start, 
	SEXP sx_count, SEXP sx_byte_style )
{
	SEXP	sx_retval, sx_retnames, sx_reterr, sx_retdata;
	int	ncid, varid, byte_style, i, err, ndims, *p_data, scalar_var;
	size_t	s_start[MAX_NC_DIMS], s_count[MAX_NC_DIMS], tot_size, k;
	char	vn[2048];
	nc_type	nct;

	/* Make space for our returned list, which will have
	 * two elements, named $error and $data
	 */
	PROTECT( sx_retval = allocVector( VECSXP, 2 ));	/* 2 elements in returned list */

	/* Set names for the returned list */
	PROTECT( sx_retnames = allocVector( STRSXP, 2 ));
	SET_STRING_ELT( sx_retnames, 0, mkChar("error") );
	SET_STRING_ELT( sx_retnames, 1, mkChar("data" ) );
	setAttrib( sx_retval, R_NamesSymbol, sx_retnames );
	UNPROTECT(1);     				/* done with sx_retnames */

	/* Set provisional 'no error' retval */
	PROTECT(sx_reterr = allocVector( INTSXP, 1 ));
	INTEGER(sx_reterr)[0] = 0;

	ncid  	   = INTEGER(sx_ncid      )[0];
	varid 	   = INTEGER(sx_varid     )[0];
	byte_style = INTEGER(sx_byte_style)[0];

	/* Get number of dimensions in this variable */
	err = nc_inq_varndims( ncid, varid, &ndims );
	if( err != NC_NOERR ) {
		Rprintf( "Error in R_nc4_get_vara_int while getting ndims: %s\n", 
			nc_strerror(err) );
		INTEGER(sx_reterr)[0] = -1;
		SET_VECTOR_ELT( sx_retval, 0, sx_reterr );
		UNPROTECT(2);
		return( sx_retval );
		}

	/* Sanity check -- number of start and count elements must
	   match number of dims in the var
	*/
	scalar_var = ((ndims==0) && (GET_LENGTH(sx_start)==1) && (INTEGER(sx_start)[0]==0) && (INTEGER(sx_count)[0]==1));
	if( (!scalar_var) && (ndims != GET_LENGTH(sx_start))) {
		Rprintf( "Error in R_nc4_get_vara_int: I think var has %d dimensions, but passed start array is length %d. They must be the same!\n",
			ndims, GET_LENGTH(sx_start) );
		INTEGER(sx_reterr)[0] = -1;
		SET_VECTOR_ELT( sx_retval, 0, sx_reterr );
		UNPROTECT(2);
		return( sx_retval );
		}

	if( (!scalar_var) && (ndims != GET_LENGTH(sx_count))) {
		Rprintf( "Error in R_nc4_get_vara_int: I think var has %d dimensions, but passed count array is length %d. They must be the same!\n",
			ndims, GET_LENGTH(sx_count) );
		INTEGER(sx_reterr)[0] = -1;
		SET_VECTOR_ELT( sx_retval, 0, sx_reterr );
		UNPROTECT(2);
		return( sx_retval );
		}

	/* Get total number of elements we will be reading so we can 
	 * allocate R space
	 */
	tot_size = 1L;
	for( i=0; i<ndims; i++ ) {
		s_start[i] = (size_t)(INTEGER(sx_start)[i]);
		s_count[i] = (size_t)(INTEGER(sx_count)[i]);
		tot_size *= s_count[i];
		}
		
	/* Make space for the returned data */
	PROTECT( sx_retdata = allocVector(INTSXP, tot_size));
	p_data = INTEGER( sx_retdata );

	/* Actually read in the data now */
	err = nc_get_vara_int( ncid, varid, s_start, s_count, p_data );
	if( err != NC_NOERR ) {
		nc_inq_varname( ncid, varid, vn );
		Rprintf( "Error in Rsx_nc4_get_vara_int: %s\n", 
			nc_strerror( err ) );
		Rprintf( "Var: %s  Ndims: %d   Start: ", vn, ndims );
		for( i=0; i<ndims; i++ ) {
			Rprintf( "%u", (unsigned int)s_start[i] );
			if( i < ndims-1 )
				Rprintf( "," );
			}
		Rprintf( " Count: " );
		for( i=0; i<ndims; i++ ) {
			Rprintf( "%u", (unsigned int)s_count[i] );
			if( i < ndims-1 )
				Rprintf( "," );
			}
		Rprintf( "\n" );
		INTEGER(sx_reterr)[0] = -1;
		SET_VECTOR_ELT( sx_retval, 0, sx_reterr );
		UNPROTECT(3);
		return( sx_retval );
		}

	err = nc_inq_vartype( ncid, varid, &nct );
	if( nct == NC_BYTE ) {
		/* netcdf library for reading from byte to int, as we do here,
		 * is SIGNED.  So, if user requests signed, we don't have to
		 * do anything; only adjust if user asks for unsigned.
		 */
		if( byte_style == 2 ) {	/* unsigned */
			for( k=0L; k<tot_size; k++ ) {
				if( p_data[k] < 0 )
					p_data[k] += 256;
				}
			}
		}

	/* OK we can successfully return now */
	SET_VECTOR_ELT( sx_retval, 0, sx_reterr  );
	SET_VECTOR_ELT( sx_retval, 1, sx_retdata );
	UNPROTECT(3);
	return( sx_retval );
}

/*********************************************************************/
void R_nc4_get_vara_text( int *ncid, int *varid, int *start, 
	int *count, char **tempstore, char **data, int *retval )
{
	int	i, err, ndims;
	size_t	s_start[MAX_NC_DIMS], s_count[MAX_NC_DIMS], nstr, slen;
	char	vn[2048];

	err = nc_inq_varndims(*ncid, *varid, &ndims );
	if( err != NC_NOERR ) 
		Rprintf( "Error in R_nc4_get_vara_text while getting ndims: %s\n", 
			nc_strerror(*retval) );

	nstr = 1L;
	for( i=0; i<ndims; i++ ) {
		s_start[i] = (size_t)start[i];
		s_count[i] = (size_t)count[i];
		if( i < (ndims-1) ) 
			nstr *= s_count[i];
		}
	slen = s_count[ndims-1];

	*retval = nc_get_vara_text(*ncid, *varid, s_start, s_count, tempstore[0] );

	if( *retval != NC_NOERR ) {
		nc_inq_varname( *ncid, *varid, vn );
		Rprintf( "Error in R_nc4_get_vara_text: %s\n", 
			nc_strerror(*retval) );
		Rprintf( "Var: %s  Ndims: %d   Start: ", vn, ndims );
		for( i=0; i<ndims; i++ ) {
			Rprintf( "%u", (unsigned int)s_start[i] );
			if( i < ndims-1 )
				Rprintf( "," );
			}
		Rprintf( " Count: " );
		for( i=0; i<ndims; i++ ) {
			Rprintf( "%u", (unsigned int)s_count[i] );
			if( i < ndims-1 )
				Rprintf( "," );
			}
		Rprintf( "\n" );
		}

	/* Now copy each string over to the final data array */
	for( i=0; i<nstr; i++ ) {
		strncpy( data[i], tempstore[0]+i*slen, slen );
		data[i][slen] = '\0';
		}
}

/*********************************************************************/
/* Returns -1 if the dim is not found in the file. Otherwise,
 * returns the dimension's length
 */
void R_nc4_inq_dimlen( int *ncid, char **dimname, int *dimlen )
{
	int err, dimid;
	err = nc_inq_dimid(*ncid, dimname[0], &dimid );
	if( err != NC_NOERR ) {
		*dimlen = -1;
		return;
		}

	size_t st_dimlen;
	err = nc_inq_dimlen( *ncid, dimid, &st_dimlen );
	*dimlen = (int)st_dimlen;
}

/*********************************************************************/
/* Returns -1 if the dim is not found in the file */
void R_nc4_inq_dimid( int *ncid, char **dimname, int *dimid )
{
	int err;
	err = nc_inq_dimid(*ncid, dimname[0], dimid );
	if( err != NC_NOERR ) 
		*dimid = -1;
}

/*********************************************************************/
/* See also the version for hierarchial groups, R_nc4_inq_varid_hier below. */ 
/* Returns -1 if the var is not found in the file */
void R_nc4_inq_varid( int *ncid, char **varname, int *varid )
{
	int err;
	err = nc_inq_varid(*ncid, varname[0], varid );
	if( err != NC_NOERR ) 
		*varid = -1;
}

/*********************************************************************/
void R_nc4_inq_unlimdim( int *ncid, int *unlimdimid, int *retval )
{
	*retval = nc_inq_unlimdim(*ncid, unlimdimid );
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_inq_unlimdim: %s\n", 
			nc_strerror(*retval) );
}

/*********************************************************************/
/* NOTE that 'dimids' must ALREADY be allocated to ndims
*/
void R_nc4_inq_dimids( int *ncid, int *dimids, int *retval )
{
	int dummy_ndims, include_parents;

	include_parents = 0;
        *retval = nc_inq_dimids( *ncid, &dummy_ndims, dimids, include_parents );
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_inq_dimids call with ncid=%d: %s\n", 
			*ncid, nc_strerror(*retval) );
}

/*********************************************************************/
void R_nc4_inq_dim( int *ncid, int *dimid, char **dimname, int *dimlen, int *unlim,
				int *retval )
{
	char name[NC_MAX_NAME];
	size_t len;
	int	i, n_unlimdims, *unlimids;

	*retval = nc_inq_dim(*ncid, *dimid, name, &len);
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_inq_dim on nc_inq_dim call with ncid=%d and dimid=%d: %s\n", 
			*ncid, *dimid, nc_strerror(*retval) );
		return;
		}

	*dimlen = (int)len;
	/* NOTE NOTE NOTE!! This assumes that the calling process
	 * allocated storage of at least NC_MAX_NAME!
	 */
	strcpy(dimname[0], name);

	/* Now set *unlim to 1 if this is an unlimited dimension,
	 * and to 0 otherwise.  Start by getting NUMBER of
	 * unlimited dimensions.
	 */
	*retval = nc_inq_unlimdims( *ncid, &n_unlimdims, NULL );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_inq_dim on nc_inq_unlimdims call (1): %s\n", 
			nc_strerror(*retval) );
		return;
		}

	/* There ARE no unlmited dims, so this can't be one! */
	if( n_unlimdims == 0 ) {
		*unlim = 0;
		return;
		}

	/* Get IDs of the unlimited dims */
	unlimids = (int *)malloc( sizeof(int)*n_unlimdims );
	if( unlimids == NULL ) {
		Rprintf( "Error in R_nc4_inq_dim: failed to allocate %d ints\n",
			n_unlimdims );
		*retval = -1;
		return;
		}

	*retval = nc_inq_unlimdims( *ncid, NULL, unlimids );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_inq_dim on nc_inq_unlimdims call (2): %s\n", 
			nc_strerror(*retval) );
		return;
		}

	/* See if there is a match between this dimid and the list of unlimited dim ids */
	*unlim = 0;
	for( i=0; i<n_unlimdims; i++ )
		if( unlimids[i] == *dimid ) {
			*unlim = 1;
			break;
			}
	free( unlimids );
}

/*********************************************************************/
void R_nc4_inq( int *ncid, int *ndims, int *nvars, int *natts,
	int *retval )
{
	*retval = nc_inq(*ncid, ndims, nvars, natts, NULL );
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_inq: %s\n", 
			nc_strerror(*retval) );
}

/*********************************************************************/
/* cmode is 0 for read only, 1 for write access.
 */
void R_nc4_open( char **filename, int *cmode, int *ncid, int *retval )
{
	int	nc_mode;

	if( *cmode == 0 )
		nc_mode = 0;
	else if( *cmode == 1 )
		nc_mode = NC_WRITE;
	else
		{
		Rprintf( "Error in R_nc4_open: bad mode passed.  Must be 0 (read) or 1 (write)\n");
		*retval = -1;
		return;
		}

	*retval = nc_open(R_ExpandFileName(filename[0]), nc_mode, ncid);
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_open: %s\n", 
			nc_strerror(*retval) );
}

/*********************************************************************/
void R_nc4_create( char **filename, int *cmode, int *ncid, int *retval )
{
	int	nc_cmode, flag_NC_NOCLOBBER, flag_NC_SHARE, flag_NC_64BIT_OFFSET,
		flag_NC_NETCDF4;

	flag_NC_NOCLOBBER       = 1;
	flag_NC_SHARE           = 2;
	flag_NC_64BIT_OFFSET    = 4;
	flag_NC_NETCDF4         = 8;
	
	/* cmode is passed in our own R bit values, not the actual
	   netcdf file values.  Convert to netcdf values.
	 */
	nc_cmode = 0;
	if( *cmode & flag_NC_NOCLOBBER )
		nc_cmode += NC_NOCLOBBER;
	if( *cmode & flag_NC_SHARE )
		nc_cmode += NC_SHARE;
	if( *cmode & flag_NC_64BIT_OFFSET )
		nc_cmode += NC_64BIT_OFFSET;
	if( *cmode & flag_NC_NETCDF4 )
		nc_cmode += NC_NETCDF4;

	*retval = nc_create(R_ExpandFileName(filename[0]), nc_cmode, ncid);
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_create: %s (creation mode was %d)\n", 
			nc_strerror(*retval), nc_cmode );
}

/*********************************************************************/
nc_type R_nc4_ttc_to_nctype( int type_to_create )
{
	/* These values are from the R code, they must match exactly */
	if( type_to_create == 1 )
		return( NC_SHORT );
	if( type_to_create == 2 )
		return( NC_INT );
	if( type_to_create == 3 )
		return( NC_FLOAT );
	if( type_to_create == 4 )
		return( NC_DOUBLE );
	if( type_to_create == 5 )
		return( NC_CHAR );
	if( type_to_create == 6 )
		return( NC_BYTE );

	error("Error, R_nc4_ttc_to_nctype passed unknown value: %d\n",
		type_to_create );
	return( NC_BYTE ); /* -Wall */
}

/*********************************************************************/
/* NOTE that NA's are handled through this vector.  According to the 
 * docs as of 2010-11-02, NA's have the value MIN_INT */
void R_nc4_put_att_logical( int *ncid, int *varid, char **attname, 
		int *type_to_create, int *natts, int *attribute, int *retval )
{
	int	R_NA_val;
	float	C_NA_val_f;
	double 	C_NA_val_d;

	/* Rprintf( "in R_nc4_put_att_logical with val=%d\n", *attribute ); */
	nc_type ttc;
	ttc = R_nc4_ttc_to_nctype( *type_to_create );

	R_NA_val = -2147483648;		/* From R docs */

	if( *attribute == R_NA_val ) {
		/* Rprintf( "PUTTING a NA -- float \n" ); */
		/* Put a NA */
		if( ttc == NC_FLOAT ) {
			C_NA_val_f = 0./0.;
			*retval = nc_put_att_float(*ncid, *varid, attname[0], 
				ttc, *natts, &C_NA_val_f );
			if( *retval != NC_NOERR ) 
				Rprintf( "Error in R_nc4_put_att_logical: %s\n", 
					nc_strerror(*retval) );
			}

		else if( ttc == NC_DOUBLE ) {
			/* Rprintf( "PUTTING a NA -- double \n" ); */
			C_NA_val_d = 0./0.;
			*retval = nc_put_att_double(*ncid, *varid, attname[0], 
				ttc, *natts, &C_NA_val_d );
			if( *retval != NC_NOERR ) 
				Rprintf( "Error in R_nc4_put_att_logical: %s\n", 
					nc_strerror(*retval) );
			}
		else
			{
			Rprintf( "Error in R_nc4_put_att_logical: asked to put a NA value, but the variable's type is not a double or float, which are the only two types that have a defined NaN value\n" );
			*retval = -1;
			return;
			}

		}
	else
		*retval = nc_put_att_int(*ncid, *varid, attname[0], 
			ttc, *natts, attribute );

	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_put_att_logical: %s\n", 
			nc_strerror(*retval) );
}

/*********************************************************************/
void R_nc4_put_att_int( int *ncid, int *varid, char **attname, 
		int *type_to_create, int *natts, int *attribute, int *retval )
{
	nc_type ttc;
	ttc = R_nc4_ttc_to_nctype( *type_to_create );
	*retval = nc_put_att_int(*ncid, *varid, attname[0], 
		ttc, *natts, attribute );
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_put_att_int: %s\n", 
			nc_strerror(*retval) );
}

/* NOTE: if you think a R_nc4_put_att_float function should be here, it means
 * you are confused.  The type attribute on the end of these function calls
 * (such as R_nc4_put_att_XXX) means "XXX" is the R storage mode, NOT the
 * type of the attribute!  R storage modes can only be integer, text, or
 * double (for our purposes); there is no "float" storage mode in R. Hence,
 * there is no R_nc4_put_att_float function.  The type of attribute to be
 * created can, of course, be float, and that is indicated by type_to_create==3.
 */

/*********************************************************************/
void R_nc4_put_att_double( int *ncid, int *varid, char **attname, 
		int *type_to_create, int *natts, double *attribute, int *retval )
{
	nc_type ttc;
	ttc = R_nc4_ttc_to_nctype( *type_to_create );
	*retval = nc_put_att_double(*ncid, *varid, attname[0], 
		ttc, *natts, attribute );
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_put_att_double: %s\n", 
			nc_strerror(*retval) );
}

/*********************************************************************/
void R_nc4_put_att_text( int *ncid, int *varid, char **attname, 
		int *type_to_create, int *natts, char **attribute, int *retval )
{
	nc_type ttc;
	size_t attlen;

	ttc = R_nc4_ttc_to_nctype( *type_to_create );
	/* For some reason the C interface does not include the nc_type for this call */

	attlen = strlen(attribute[0]);
	*retval = nc_put_att_text(*ncid, *varid, attname[0], 
		attlen, attribute[0] );
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_put_att_text: %s\n", 
			nc_strerror(*retval) );
}

/*********************************************************************/
/* NOTE that space for the attribute must already be allocated! */
void R_nc4_get_att_int( int *ncid, int *varid, char **attname, 
			int *attribute, int *retval )
{
	*retval = nc_get_att_int(*ncid, *varid, attname[0], attribute);
}

/*********************************************************************/
/* NOTE that space for the attribute must already be allocated! */
void R_nc4_get_att_double( int *ncid, int *varid, char **attname, 
			double *attribute, int *retval )
{
	*retval = nc_get_att_double(*ncid, *varid, attname[0], attribute);
}

/*********************************************************************/
/* If returned value 'retval' is 0, then returned value 'type' will hold
 * the type of the named attribute, and returned value 'attlen' will
 * hold the length of the named attribute.  If returned value 'retval'
 * is NOT 0, then the specified variable did not have an attribute
 * named 'attname'.
 */
void R_nc4_inq_att( int *ncid, int *varid, char **attname, int *type,
			int *attlen, int *retval )
{
	size_t  s_attlen;
	nc_type nctype;

	*retval = nc_inq_att(*ncid, *varid, attname[0], &nctype, &s_attlen );
	if( (*retval != 0) && (*retval != NC_ENOTATT)) 
		Rprintf( "Error in R_nc4_inq_att: while looking for attribute %s, got error %s\n",
			attname[0], nc_strerror(*retval) );

	if( *retval == 0 ) {
		*type = R_nc4_nctype_to_Rtypecode(nctype);
		if( *type == -1 ) {
			if( nctype == NC_BYTE )
				Rprintf( "Error in R_nc4_inq_att: not set up to handle attributes of type \"BYTE\"!  Netcdf type code: %d Attribute name: %s\n", nctype, attname[0] );
			else
				{
				Rprintf( "Error in R_nc4_inq_att: not set up to handle attributes of this type!  Netcdf type code: %d Attribute name: %s\n", nctype, attname[0] );
				*retval = -1;
				}
			}

		*attlen = (int)s_attlen;
		}
}

/*********************************************************************/
/* NOTE that space for the attribute must already be allocated! */
void R_nc4_get_att_text( int *ncid, int *varid, char **attname, 
			char **attribute, int *retval )
{
	int	err;
	size_t	attlen;

	*retval = nc_get_att_text(*ncid, *varid, attname[0], 
		attribute[0] );
	/* append a NULL */
	if( *retval != NC_NOERR ) {
		Rprintf( "R_nc4_get_att_text: error encountered on call to nc_get_att_text: %s\n", nc_strerror(*retval) );
		strcpy( attribute[0], "\0" );
		return;
		}
	err = nc_inq_attlen( *ncid, *varid, attname[0], &attlen );
	if( err != NC_NOERR ) {
		strcpy( attribute[0], "\0" );
		return;
		}
	attribute[0][attlen] = '\0';
}

/*********************************************************************/
/* With strings, we explicitly allocate space for them here.
 * Note that the string array is returned as the function return value; the
 * error code is passed in and out as an argument.
 */
SEXP R_nc4_get_att_string( SEXP sx_ncid, SEXP sx_varid, SEXP sx_attname, SEXP sx_attlen, SEXP sx_ierr_returned )
{
	int	ncid, varid, attlen, ierr, i, slen;
	char	**strings;
	const char *attname = CHAR(STRING_ELT(sx_attname,0));
	char	*tstr;
	SEXP	sx_string, sx_ret_str_array;

	INTEGER(sx_ierr_returned)[0] = 0;	/* initialize to "no error" */

	ncid   = INTEGER(sx_ncid  )[0];
	varid  = INTEGER(sx_varid )[0];
	attlen = INTEGER(sx_attlen)[0];	/* for strings, this is the NUMBER of strings, not any string's length */

	if( attlen < 1 ) {
		Rf_error( "Error, in call to R_nc4_get_att_string, number of strings (attlen) must be >= 1\n" );
		INTEGER(sx_ierr_returned)[0] = -1;
		return( R_NilValue );
		}

	/* Make space in C to hold array of pointers to strings */
	strings = (char **)R_alloc( attlen, sizeof(char *) );

	/* Read the strings from the netcdf file. Note that the netcdf library
	 * call allocates space for the strings themselves
	 */
	if( (ierr = nc_get_att_string( ncid, varid, attname, strings )) != NC_NOERR ) {
		Rf_error( "Error, in call to R_nc4_get_att_string, failed to get the strings:\n" );
		Rf_error( nc_strerror( ierr ) );
		INTEGER(sx_ierr_returned)[0] = -1;
		return( R_NilValue );
		}

	/* Initialize our return value
	 */
	sx_ret_str_array = PROTECT( NEW_CHARACTER( attlen ));	/* make array with "attlen" strings */

	/* Copy from our C array to an R array that we make. We do this so 
	 * that we can free the C netcdf-library allocated space, which
	 * otherwise would never be freed
	 */
	for( i=0; i<attlen; i++ ) {
		slen = strlen( strings[i] );
		tstr = R_alloc( slen+1, sizeof(char) );
		strncpy(tstr, strings[i], slen);
		tstr[slen] = '\0';
		sx_string = PROTECT( mkChar(tstr) );

		SET_STRING_ELT( sx_ret_str_array, i, sx_string );
		}

	/* Free C memory allocated by netcdf library
	 */
	nc_free_string( attlen, strings );

	UNPROTECT(attlen+1);

	return( sx_ret_str_array );
}

/*********************************************************************/
SEXP Rsx_nc4_put_vara_double( SEXP sx_ncid, SEXP sx_varid, SEXP sx_start,
	SEXP sx_count, SEXP sx_data )
{
	int 	ncid, varid, i, ndims, err, verbose, scalar_var;
	size_t  s_start[MAX_NC_DIMS], s_count[MAX_NC_DIMS];
	char    varname[MAX_NC_NAME];
	SEXP	sx_retval;

	ncid  = INTEGER(sx_ncid )[0];
	varid = INTEGER(sx_varid)[0];

	verbose = 0;

	if( verbose ) {
		err = nc_inq_varname( ncid, varid, varname );
		Rprintf( "Rsx_nc4_put_vara_double: entering with ncid=%d, varid=%d  (varname=%s)\n", 
				ncid, varid, varname );
		}

	/* Set provisional 'no error' retval */
	PROTECT(sx_retval = NEW_NUMERIC(1));
	NUMERIC_POINTER(sx_retval)[0] = 0;

	/* Get # of dims for this var */
	err = nc_inq_varndims( ncid, varid, &ndims );
	if( err != NC_NOERR ) {
		Rprintf( "Error on nc_inq_varndims call in Rsx_nc4_put_vara_double: %s\n", 
			nc_strerror(err) );
		NUMERIC_POINTER(sx_retval)[0] = -1;
		UNPROTECT(1);
		return( sx_retval );
		}
	if( verbose ) 
		Rprintf( "Rsx_nc4_put_vara_double: for this var ndims=%d\n", ndims );

	/* Sanity check -- number of start and count elements must
	   match nubmer of dims in the var
	*/
	scalar_var = ((ndims==0) && (GET_LENGTH(sx_start)==1) && (INTEGER(sx_start)[0]==0) && (INTEGER(sx_count)[0]==1));
	if( (!scalar_var) && (ndims != GET_LENGTH(sx_start))) {
		Rprintf( "Error in Rsx_nc4_put_vara_double: I think var has %d dimensions, but passed start array is length %d. They must be the same!\n",
			ndims, GET_LENGTH(sx_start) );
		NUMERIC_POINTER(sx_retval)[0] = -1;
		UNPROTECT(1);
		return( sx_retval );
		}

	if( (!scalar_var) && (ndims != GET_LENGTH(sx_count))) {
		Rprintf( "Error in Rsx_nc4_put_vara_double: I think var has %d dimensions, but passed count array is length %d. They must be the same!\n",
			ndims, GET_LENGTH(sx_count) );
		NUMERIC_POINTER(sx_retval)[0] = -1;
		UNPROTECT(1);
		return( sx_retval );
		}

	/* Copy over from ints to size_t */
	for( i=0; i<ndims; i++ ) {
		s_start[i] = (size_t)(INTEGER(sx_start)[i]);
		s_count[i] = (size_t)(INTEGER(sx_count)[i]);
		}

	if( verbose ) {
		Rprintf( "Rsx_nc4_put_vara_double: about to write with start=" );
		for( i=0; i<ndims; i++ ) 
			Rprintf("%d ", s_start[i] );
		Rprintf( "   count=" );
		for( i=0; i<ndims; i++ ) 
			Rprintf("%d ", s_count[i] );
		Rprintf( "\n" );
		}

	err = nc_put_vara_double( ncid, varid, s_start, s_count, REAL(sx_data) );
	if( err != NC_NOERR ) {
		Rprintf( "Error in Rsx_nc4_put_vara_double: %s\n", nc_strerror(err) );
		NUMERIC_POINTER(sx_retval)[0] = -1;
		UNPROTECT(1);
		return( sx_retval );
		}

	if( verbose ) 
		Rprintf( "Rsx_nc4_put_vara_double: successful return\n" );

	UNPROTECT(1);
	return( sx_retval );
}

/*********************************************************************/
SEXP Rsx_nc4_put_vara_int( SEXP sx_ncid, SEXP sx_varid, SEXP sx_start,
	SEXP sx_count, SEXP sx_data )
{
	int 	ncid, varid, i, ndims, err, scalar_var;
	size_t s_start[MAX_NC_DIMS], s_count[MAX_NC_DIMS];
	SEXP	sx_retval;

	/* Set provisional 'no error' retval */
	PROTECT(sx_retval = NEW_NUMERIC(1));
	NUMERIC_POINTER(sx_retval)[0] = 0;

	ncid  = INTEGER(sx_ncid )[0];
	varid = INTEGER(sx_varid)[0];

	/* Get # of dims for this var */
	err = nc_inq_varndims( ncid, varid, &ndims );
	if( err != NC_NOERR ) {
		Rprintf( "Error on nc_inq_varndims call in Rsx_nc4_put_vara_int: %s\n", 
			nc_strerror(err) );
		NUMERIC_POINTER(sx_retval)[0] = -1;
		UNPROTECT(1);
		return( sx_retval );
		}

	/* Sanity check -- number of start and count elements must
	   match nubmer of dims in the var
	*/
	scalar_var = ((ndims==0) && (GET_LENGTH(sx_start)==1) && (INTEGER(sx_start)[0]==0) && (INTEGER(sx_count)[0]==1));
	if( (!scalar_var) && (ndims != GET_LENGTH(sx_start))) {
		Rprintf( "Error in Rsx_nc4_put_vara_int: I think var has %d dimensions, but passed start array is length %d. They must be the same!\n",
			ndims, GET_LENGTH(sx_start) );
		NUMERIC_POINTER(sx_retval)[0] = -1;
		UNPROTECT(1);
		return( sx_retval );
		}

	if( (!scalar_var) && (ndims != GET_LENGTH(sx_count))) {
		Rprintf( "Error in Rsx_nc4_put_vara_int: I think var has %d dimensions, but passed count array is length %d. They must be the same!\n",
			ndims, GET_LENGTH(sx_count) );
		NUMERIC_POINTER(sx_retval)[0] = -1;
		UNPROTECT(1);
		return( sx_retval );
		}

	/* Copy over from ints to size_t */
	for( i=0; i<ndims; i++ ) {
		s_start[i] = (size_t)(INTEGER(sx_start)[i]);
		s_count[i] = (size_t)(INTEGER(sx_count)[i]);
		}

	/* Actually write the data now */
	err = nc_put_vara_int(ncid, varid, s_start, s_count, INTEGER(sx_data) );
	if( err != NC_NOERR ) {
		Rprintf( "Error in Rsx_nc4_put_vara_int: %s\n", 
			nc_strerror(err) );
		NUMERIC_POINTER(sx_retval)[0] = -1;
		UNPROTECT(1);
		return( sx_retval );
		}

	UNPROTECT(1);
	return( sx_retval );
}

/**************************************************************************************************************/
void R_nc4_put_vara_text( int *ncid, int *varid, int *start,
	int *count, char **data, int *retval )
{
	int	ndims, err, idx_string, idx_char, idx_j, idx_k;
	size_t 	s_start[MAX_NC_DIMS], s_count[MAX_NC_DIMS], slen, slen2use;
	long	i, j, k, stridx, nstrings, nj, nk;

	*retval = 0;

	/* Get # of dims for this var */
	err = nc_inq_varndims( *ncid, *varid, &ndims );
	if( err != NC_NOERR )
		Rprintf( "Error (loc #1) on nc_inq_ndims call in R_nc_put_vara_int: %s\n", 
			nc_strerror(*retval) );

	/* Copy over from ints to size_t.  Remember things are in C style at
	 * this point, so the rightmost dim is the number of characters
	 */
	for( i=0; i<ndims; i++ ) {
		s_start[i] = (size_t)start[i];
		s_count[i] = (size_t)count[i];
		}

	/* Chars are an unusually difficult because R seems to store
	 * them as an array of character pointers, while netcdf stores
	 * them as a monolithic block (like any other var type).  We
	 * must convert between these representations.
	 */
	slen = s_count[ndims-1];
	if( ndims == 1 ) {
		*retval = nc_put_vara_text(*ncid, *varid, s_start, s_count, data[0] );
		if( *retval != NC_NOERR ) 
			Rprintf( "Error (loc #2) in R_nc_put_vara_int: %s\n", 
				nc_strerror(*retval) );
		}

	else if( ndims == 2 ) {
		idx_string = 0;
		idx_char   = 1;
		nstrings = s_count[idx_string];  /* number of character strings */
/* Rprintf( "nstrings=%ld\n", nstrings ); */
		for( i=0L; i<nstrings; i++ ) {
			slen2use = ((slen < strlen(data[i])) ? slen : strlen(data[i]));
/* Rprintf( "slen2use=%ld\n", slen2use ); */
			s_count[idx_string] = 1L;
			s_count[idx_char  ] = slen2use;
			s_start[idx_string] = i + start[idx_string];
			s_start[idx_char  ] = 0L;
			*retval = nc_put_vara_text(*ncid, *varid, s_start, s_count, data[i] );
			if( *retval != NC_NOERR ) {
				Rprintf( "Error (loc #3) in R_nc_put_vara_text: %s\n", 
					nc_strerror(*retval) );
				Rprintf( "Here was C-style start I tried: %ld %ld\n", 
					s_start[0], s_start[1] );
				Rprintf( "Here was C-style count I tried: %ld %ld\n", 
					s_count[0], s_count[1] );
				return;
				}
			}
		}

	else if( ndims == 3 ) {
		idx_j      = 0;
		idx_string = 1;
		idx_char   = 2;
		nj       = s_count[idx_j];
		nstrings = s_count[idx_string];  /* number of character strings PER j */
/* Rprintf( "nstrings=%ld  jn=%ld\n", nstrings, nj ); */
		stridx = 0L;
		for( j=0L; j<nj; j++ ) 
		for( i=0L; i<nstrings; i++ ) {
			slen2use = ((slen < strlen(data[stridx])) ? slen : strlen(data[stridx]));
/* Rprintf( "slen2use=%ld\n", slen2use ); */
			s_count[idx_j     ] = 1L;
			s_count[idx_string] = 1L;
			s_count[idx_char  ] = slen2use;
			s_start[idx_j     ] = j + start[idx_j];
			s_start[idx_string] = i + start[idx_string];
			s_start[idx_char  ] = 0L;
/* Rprintf( "writing following string: >%s<\n", data[stridx] ); */
			*retval = nc_put_vara_text(*ncid, *varid, s_start, s_count, data[stridx++] );
			if( *retval != NC_NOERR ) {
				Rprintf( "Error (loc #4) in R_nc_put_vara_text: %s\n", 
					nc_strerror(*retval) );
				Rprintf( "Here was C-style start I tried: %ld %ld %ld\n", 
					s_start[0], s_start[1], s_start[2] );
				Rprintf( "Here was C-style count I tried: %ld %ld %ld\n", 
					s_count[0], s_count[1], s_count[2] );
				return;
				}
			}
		}

	else if( ndims == 4 ) {
		idx_k      = 0;
		idx_j      = 1;
		idx_string = 2;
		idx_char   = 3;
		nk       = s_count[idx_k];
		nj       = s_count[idx_j];
		nstrings = s_count[idx_string];  /* number of character strings PER j,k */
/* Rprintf( "nstrings=%ld  nj=%ld  nk=%ld\n", nstrings, nj, nk ); */
		stridx = 0L;
		for( k=0L; k<nk; k++ ) 
		for( j=0L; j<nj; j++ ) 
		for( i=0L; i<nstrings; i++ ) {
			slen2use = ((slen < strlen(data[stridx])) ? slen : strlen(data[stridx]));
/* Rprintf( "slen2use=%ld\n", slen2use ); */
			s_count[idx_k     ] = 1L;
			s_count[idx_j     ] = 1L;
			s_count[idx_string] = 1L;
			s_count[idx_char  ] = slen2use;
			s_start[idx_k     ] = k + start[idx_k];
			s_start[idx_j     ] = j + start[idx_j];
			s_start[idx_string] = i + start[idx_string];
			s_start[idx_char  ] = 0L;
/* Rprintf( "writing following string: >%s<\n", data[stridx] ); */
			*retval = nc_put_vara_text(*ncid, *varid, s_start, s_count, data[stridx++] );
			if( *retval != NC_NOERR ) {
				Rprintf( "Error (loc #5) in R_nc_put_vara_text: %s\n", 
					nc_strerror(*retval) );
				Rprintf( "Here was C-style start I tried: %ld %ld %ld %ld\n", 
					s_start[0], s_start[1], s_start[2], s_start[3] );
				Rprintf( "Here was C-style count I tried: %ld %ld %ld %ld\n", 
					s_count[0], s_count[1], s_count[2], s_count[3] );
				return;
				}
			}
		}

	else
		{
		*retval = -1;
		Rprintf("Error in R_nc_put_vara_text: unhandled case.  I only handle char dims with # of dims up to 4.  Was passed # dims = %d\n", ndims );
		return;
		}
}

/*********************************************************************/
void R_nc4_def_var_byte( int *ncid, char **varname, int *ndims, int *dimids, 
	int *varid, int *retval )
{
	*retval = nc_def_var(*ncid, varname[0], 
		NC_BYTE, *ndims, dimids, varid );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_def_var_byte: %s\n", 
			nc_strerror(*retval) );
		Rprintf( "Name of variable that the error occurred on: \"%s\"\n", 
			varname[0] );
		if( *retval == NC_ENAMEINUSE )
			Rprintf( "I.e., you are trying to add a variable with that name to the file, but it ALREADY has a variable with that name!\n");
		}
}

/*********************************************************************/
void R_nc4_def_var_int( int *ncid, char **varname, int *ndims, int *dimids, 
	int *varid, int *retval )
{
	*retval = nc_def_var(*ncid, varname[0], 
		NC_INT, *ndims, dimids, varid );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_def_var_int: %s\n", 
			nc_strerror(*retval) );
		Rprintf( "Name of variable that the error occurred on: \"%s\"\n", 
			varname[0] );
		if( *retval == NC_ENAMEINUSE )
			Rprintf( "I.e., you are trying to add a variable with that name to the file, but it ALREADY has a variable with that name!\n");
		}
}

/*********************************************************************/
void R_nc4_def_var_short( int *ncid, char **varname, int *ndims, int *dimids, 
	int *varid, int *retval )
{
	*retval = nc_def_var(*ncid, varname[0], 
		NC_SHORT, *ndims, dimids, varid );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_def_var_short: %s\n", 
			nc_strerror(*retval) );
		Rprintf( "Name of variable that the error occurred on: \"%s\"\n", 
			varname[0] );
		if( *retval == NC_ENAMEINUSE )
			Rprintf( "I.e., you are trying to add a variable with that name to the file, but it ALREADY has a variable with that name!\n");
		}
}

/*********************************************************************/
void R_nc4_def_var_float( int *ncid, char **varname, int *ndims, int *dimids, 
	int *varid, int *retval )
{
	*retval = nc_def_var(*ncid, varname[0], 
		NC_FLOAT, *ndims, dimids, varid );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_def_var_float: %s\n", 
			nc_strerror(*retval) );
		Rprintf( "Name of variable that the error occurred on: \"%s\"\n", 
			varname[0] );
		if( *retval == NC_ENAMEINUSE )
			Rprintf( "I.e., you are trying to add a variable with that name to the file, but it ALREADY has a variable with that name!\n");
		}
}

/*********************************************************************/
void R_nc4_def_var_double( int *ncid, char **varname, int *ndims, int *dimids, 
	int *varid, int *retval )
{
	*retval = nc_def_var(*ncid, varname[0], 
		NC_DOUBLE, *ndims, dimids, varid );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_def_var_double: %s\n", 
			nc_strerror(*retval) );
		Rprintf( "Name of variable that the error occurred on: \"%s\"\n", 
			varname[0] );
		if( *retval == NC_ENAMEINUSE )
			Rprintf( "I.e., you are trying to add a variable with that name to the file, but it ALREADY has a variable with that name!\n");
		}
}

/*********************************************************************/
void R_nc4_def_var_char( int *ncid, char **varname, int *ndims, int *dimids, 
	int *varid, int *retval )
{
	*retval = nc_def_var(*ncid, varname[0], 
		NC_CHAR, *ndims, dimids, varid );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_def_var_char: %s\n", 
			nc_strerror(*retval) );
		Rprintf( "Name of variable that the error occurred on: \"%s\"\n", 
			varname[0] );
		if( *retval == NC_ENAMEINUSE )
			Rprintf( "I.e., you are trying to add a variable with that name to the file, but it ALREADY has a variable with that name!\n");
		}
}

/*********************************************************************/
void R_nc4_def_dim( int *ncid, char **dimname, int *size, int *dimid, 
	int *retval )
{
	*retval = nc_def_dim(*ncid, dimname[0], 
		*size, dimid );
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_def_dim: %s\n", 
			nc_strerror(*retval) );
}

/*********************************************************************/
void R_nc4_redef( int *ncid )
{
	int	err;
	err = nc_redef(*ncid);
	if( err != NC_NOERR ) 
		Rprintf( "Error in R_nc4_redef: %s\n", 
			nc_strerror(err) );
}

/*********************************************************************/
void R_nc4_rename_var( int *ncid, int *varid, char **newname, int *retval )
{
	*retval = nc_rename_var( *ncid, *varid, newname[0] );
	if( *retval != NC_NOERR ) 
		Rprintf( "Error in R_nc4_rename_var: %s\n", 
			nc_strerror(*retval) );
}

/*********************************************************************/
void R_nc4_inq_attname( int *ncid, int *varid, int *attnum, char **attname, int *retval )
{
	*retval = nc_inq_attname( *ncid, *varid, *attnum, attname[0] );
	if( *retval != NC_NOERR ) {
		Rprintf( "Error in R_nc4_inq_attname: %s\n",
			nc_strerror(*retval) );
		}
}

/*********************************************************************/
/* Inputs:
	root_id:	netcdfID of the root group
	varid:		var ID to set the deflate params for
	ndims:		number of dims in the variable
	storage:	either 1 for NC_CONTIGUOUS or 2 for NC_CHUNKED
	chunksizesp:	pointer to array indicating chunk sizes along each dimension
   Outputs:
	ierr:		0 on success, otherwise an error was encountered
*/
void R_nc4_def_var_chunking( int *root_id, int *varid, int *ndims, int *storage, int *chunksizesp, int *ierr )
{
	int	stor_param, i;
	size_t	sizet_chunkparam[MAX_NC_DIMS];

	if( (*storage != 1) && (*storage != 2)) {
		Rprintf( "R_nc4_def_var_chunking: bad value of storage parameter passed: %d. Must be 1 (NC_CONTIGUOUS) or 2 (NC_CHUNKED).\n",
			*storage );
		*ierr = -1;
		return;
		}

	if( *storage == 1 )
		stor_param = NC_CONTIGUOUS;
	else
		stor_param = NC_CHUNKED;

	if( *ndims > MAX_NC_DIMS ) {
		Rprintf( "R_nc4_def_var_chunking: Error: ndims too large, max is %d but passed value is %d\n",
			MAX_NC_DIMS, *ndims );
		*ierr = -1;
		return;
		}

	for( i=0; i<(*ndims); i++ )
		sizet_chunkparam[i] = (size_t)(*(chunksizesp+i));

	*ierr = nc_def_var_chunking( *root_id, *varid, stor_param, sizet_chunkparam );
	if( *ierr != NC_NOERR ) {
		Rprintf( "Error in R_nc4_def_var_chunking: %s\n", nc_strerror(*ierr) );
		}
}

/*********************************************************************/
/* Inputs:
	root_id:	netcdfID of the root group
	varid:		var ID to set the deflate params for
	ndims:		number of dims in the variable
   Outputs:
	storage:	either 1 for NC_CONTIGUOUS or 2 for NC_CHUNKED
	chunksizesp:	pointer to array indicating chunk sizes along each dimension. NOTE:
			THIS MUST BE ALLOCATED BEFORE CALLING THIS ROUTINE!!
			Returned values are in C order (which is opposite R order)
	ierr:		0 on success, otherwise an error was encountered
*/
void R_nc4_inq_var_chunking( int *root_id, int *varid, int *ndims, int *storage, int *chunksizesp, int *ierr )
{
	int	stor_param, i;
	size_t	sizet_chunkparam[MAX_NC_DIMS];

	*ierr = nc_inq_var_chunking( *root_id, *varid, &stor_param, sizet_chunkparam );
	if( *ierr != NC_NOERR ) {
		Rprintf( "Error in R_nc4_inq_var_chunking: %s\n", nc_strerror(*ierr) );
		return;
		}

	if( stor_param == NC_CONTIGUOUS )
		*storage = 1;
	else if( stor_param == NC_CHUNKED )
		*storage = 2;
	else
		{
		Rprintf( "Error in R_nc4_inq_var_chunking: obtained value of storage is neither NC_CONTIGUOUS nor NC_CHUNKED!  Value=%d\n",
			stor_param );
		*ierr = -1;
		return;
		}

	for( i=0; i<(*ndims); i++ )
		chunksizesp[i] = (int)(sizet_chunkparam[i]);
}

/*********************************************************************/
/* Inputs:
	root_id:	netcdfID of the root group
	varid:		var ID to inquire about the deflate params for
   Outputs:
	shuffle:	if 1, shuffle filter is turned on for this var; otherwise, this is set to 0
	deflate:	if 1, deflate (compression) filter is turned on for this var; otherwise, this is set to 0
	deflate_level:	If deflate == 1, this is set to the deflate level for this var (0-9)
	ierr:		0 on success, otherwise an error was encountered
*/
void R_nc4_inq_var_deflate( int *root_id, int *varid, int *shuffle, int *deflate, int *deflate_level, int *ierr )
{
	*ierr = nc_inq_var_deflate( *root_id, *varid, shuffle, deflate, deflate_level );
	if( *ierr != NC_NOERR ) {
		Rprintf( "Error in R_nc4_inq_var_deflate: %s\n", nc_strerror(*ierr) );
		}
}

/*********************************************************************/
/* Inputs:
	root_id:	netcdfID of the root group
	varid:		var ID to set the deflate params for
	shuffle:	if 1, turn on shuffle filter; otherwise, must be 0
	deflate:	if 1, turn on deflate (compression) filter; otherwise, must be 0
	deflate_level:	must be bewteen 0 and 9.  If deflate != 0, set deflate to this level.
   Outputs:
	ierr:		0 on success, otherwise an error was encountered
*/
void R_nc4_def_var_deflate( int *root_id, int *varid, int *shuffle, int *deflate, int *deflate_level, int *ierr )
{
	if( (*shuffle != 0) && (*shuffle != 1)) {
		Rprintf( "R_nc4_def_var_deflate: bad value of shuffle parameter passed: %d. Must be 0 or 1.\n",
			*shuffle );
		*ierr = -1;
		return;
		}

	if( (*deflate != 0) && (*deflate != 1)) {
		Rprintf( "R_nc4_def_var_deflate: bad value of deflate parameter passed: %d. Must be 0 or 1.\n",
			*deflate );
		*ierr = -1;
		return;
		}

	if( (*deflate_level < 0) || (*deflate_level > 9)) {
		Rprintf( "R_nc4_def_var_deflate: bad value of deflate_level parameter passed: %d. Must be between 0 and 9, inclusive.\n",
			*deflate_level );
		*ierr = -1;
		return;
		}

	*ierr = nc_def_var_deflate( *root_id, *varid, *shuffle, *deflate, *deflate_level );
	if( *ierr != NC_NOERR ) {
		Rprintf( "Error in R_nc4_def_var_deflate: %s\n", nc_strerror(*ierr) );
		}
}

/*********************************************************************/
/* Inputs:
	root_id:	netcdfID of the root group
   Outputs:
   	ngroups:	How many groups are BELOW the specified root group
	ierr:		0 on success, otherwise an error was encountered
*/
void R_nc4_inq_ngroups( int *root_id, int *ngroups, int *ierr )
{
	*ierr = nc_inq_grps( *root_id, ngroups, NULL );
	if( *ierr != NC_NOERR ) {
		Rprintf( "Error in R_nc4_inq_ngroups: %s\n",
			nc_strerror(*ierr) );
		}
}

/*********************************************************************/
/* Inputs:
	root_id:	netcdfID of the root group
   Outputs:
   	gids:		IDs of the groups below the root group, must already be allocated!
	ierr:		0 on success, otherwise an error was encountered
*/
void R_nc4_inq_groupids( int *root_id, int *gids, int *ierr )
{
	*ierr = nc_inq_grps( *root_id, NULL, gids );
	if( *ierr != NC_NOERR ) {
		Rprintf( "Error in R_nc4_inq_groupids: %s\n",
			nc_strerror(*ierr) );
		}
}

/*********************************************************************/
/* Code kindly donated by Martin Morgan of the Fred Hutchinson
 * Cancer Research Center
 */
SEXP R_nc4_blankstring(SEXP size)
{
	int	i, len;
	char	*str;
	SEXP	blank, string;

	if (!IS_INTEGER(size) || LENGTH(size)!=1)
		Rf_error("'size' must be 'integer(1)'");

	len = INTEGER(size)[0];
	if (len < 0)
		Rf_error("'size' must be >= 0");

	str = R_alloc(len + 1, sizeof(char));
	for (i = 0; i < len; ++i)
		str[i] = ' ';
	str[len] = '\0';

	blank  = PROTECT(NEW_CHARACTER(1));
	string = PROTECT(mkChar(str));
	SET_STRING_ELT(blank, 0, string);
	UNPROTECT(2);

	return(blank); 
}

/*********************************************************************/
/* Input: 
	sx_root_id: 	integer group id
   Outputs:
   	sx_ierr:	set to 0 on success, -1 on failure
	returned value: a string witht he group name
*/
SEXP R_nc4_grpname(SEXP sx_root_id, SEXP sx_ierr_retval)
{
	int	root_id, ierr;
	size_t	nchar;
	char	*str;
	SEXP	sx_retval, sx_string;

	INTEGER(sx_ierr_retval)[0] = 0;		/* start with 'no error' condition */

	/* Make sure our passed args are the right type and size */
	if (!IS_INTEGER(sx_root_id) || LENGTH(sx_root_id)!=1)
		Rf_error("Passed argument sx_root_id must be 'integer(1)'");

	if (!IS_INTEGER(sx_ierr_retval) || LENGTH(sx_ierr_retval)!=1)
		Rf_error("Passed argument sx_ierr_retval must be 'integer(1)'");

	/* Turn R sx_value into ordinary integer */
	root_id = INTEGER(sx_root_id)[0];

	/* Get # of chars in group name */
	ierr = nc_inq_grpname_len( root_id, &nchar );
	if( ierr != 0 ) {
		Rprintf( "Error in R_nc_grpname: on call to nc_inq_grpname_len: %s\n", 
				nc_strerror(ierr) );
		INTEGER(sx_ierr_retval)[0] = -1;	/* Indicate error condition */
		return( NULL );
		}

	/* Make space for C string to be returned */
	str = R_alloc(nchar + 1, sizeof(char));
	ierr = nc_inq_grpname( root_id, str );
	if( ierr != 0 ) {
		Rprintf( "Error in R_nc_grpname: on call to nc_inq_grpname: %s\n", 
				nc_strerror(ierr) );
		INTEGER(sx_ierr_retval)[0] = -1;	/* Indicate error condition */
		return( NULL );
		}
	str[nchar] = '\0';

	sx_retval = PROTECT(NEW_CHARACTER(1));
	sx_string = PROTECT(mkChar(str));
	SET_STRING_ELT(sx_retval, 0, sx_string);
	UNPROTECT(2);

	return(sx_retval); 
}

/*********************************************************************/
/* Input: 
	sx_root_id: 	integer group id
   Outputs:
   	sx_ierr:	set to 0 on success, -1 on failure
	returned value: an integer as follows:
		NC_FORMAT_CLASSIC:	1
		NC_FORMAT_64BIT:	2
		NC_FORMAT_NETCDF4:	3
		NC_FORMAT_NETCDF4_CLASSIC: 4
*/
SEXP R_nc4_inq_format(SEXP sx_root_id, SEXP sx_ierr_retval)
{
	int	root_id, ierr;
	int	iformat, iretval;
	SEXP	sx_retval;

	INTEGER(sx_ierr_retval)[0] = 0;		/* start with 'no error' condition */

	/* Make sure our passed args are the right type and size */
	if (!IS_INTEGER(sx_root_id) || LENGTH(sx_root_id)!=1)
		Rf_error("Passed argument sx_root_id must be 'integer(1)'");

	if (!IS_INTEGER(sx_ierr_retval) || LENGTH(sx_ierr_retval)!=1)
		Rf_error("Passed argument sx_ierr_retval must be 'integer(1)'");

	/* Turn R sx_value into ordinary integer */
	root_id = INTEGER(sx_root_id)[0];

	/* Get format */
	ierr = nc_inq_format( root_id, &iformat );
	if( ierr != 0 ) {
		Rprintf( "Error in R_nc4_inq_format: on call to nc_inq_format: %s\n", 
				nc_strerror(ierr) );
		INTEGER(sx_ierr_retval)[0] = -1;	/* Indicate error condition */
		return( NULL );
		}

	switch(iformat) {
		
		case NC_FORMAT_CLASSIC:
			iretval = 1;
			break;

		case NC_FORMAT_64BIT:
			iretval = 2;
			break;

		case NC_FORMAT_NETCDF4:
			iretval = 3;
			break;

		case NC_FORMAT_NETCDF4_CLASSIC:
			iretval = 4;
			break;

		default:
			Rprintf( "Error in R_nc4_inq_format: unrecognized format integer returned: %d\n", iformat );
			INTEGER(sx_ierr_retval)[0] = -1;	/* Indicate error condition */
			return( NULL );
		}

	PROTECT( sx_retval = allocVector(INTSXP, 1));
	INTEGER(sx_retval)[0] = iretval;

	UNPROTECT(1);

	return(sx_retval); 
}

/*********************************************************************/
/* This goes through an input array, and replaces all NA's in that array
 * with the passed value.
 */
SEXP R_nc4_set_NA_to_val_double(SEXP sx_dat, SEXP sx_val )
{
	int		nels;
	size_t		i;
	double		val = REAL(sx_val)[0];

	nels = GET_LENGTH( sx_dat );

	for( i=0L; i<nels; i++ ) {
		if( ISNA(REAL(sx_dat)[i]) ) 
			NUMERIC_POINTER(sx_dat)[i] = val;
		}

	return( R_NilValue );
}

/************************************************************************************************/
/* Returned index_first_slash is 0-based counting 
 */
int R_nc4_util_nslashes( char *s, int *idx_first_slash ) 
{
	int	i, nslashes;

	nslashes = 0;
	*idx_first_slash = -1;
	for( i=0; i<strlen(s); i++ ) {
		if( s[i] == '/' ) {
			nslashes++;
			if( *idx_first_slash == -1 ) 
				*idx_first_slash = i;
			}
		}

	return( nslashes );
}

/****************************************************************************************/
void R_nc4_inq_varid_hier_inner( int *ncid, char *varname, int *returned_grpid, int *returned_varid )
{
	int ierr, nslashes, idx_first_slash, gid;
	char gname[MAX_NC_NAME];

/* Rprintf("R_nc4_inq_varid_hier_inner: Entering with varname=>%s<\n", varname ); */

	/* Passed var name must not start with a slash */
	if( varname[0] == '/' ) {
		Rprintf( "Error in R_nc4_varid_hier: passed varname must not start with a slash!\n" );
		*returned_varid = -1;
		*returned_grpid = -1;
		return;
		}

	/* If there are no forward slashes in the name, then it is a simple var name */
	nslashes = R_nc4_util_nslashes( varname, &idx_first_slash );
/* Rprintf("R_nc4_inq_varid_hier_inner: varname >%s< has nslashes=%d and idx first slash=%d\n", varname, nslashes, idx_first_slash ); */
	if( nslashes == 0 ) {
		*returned_grpid = *ncid;
		ierr = nc_inq_varid(*ncid, varname, returned_varid );
		if( ierr != NC_NOERR ) {
			*returned_varid = -1;
			*returned_grpid = -1;
			}
/* Rprintf("R_nc4_inq_varid_hier_inner: Found id for var >%s< varid:%d gid:%d\n", varname, *returned_varid, *returned_grpid ); */
		return;
		}

	/* If we get here then the varname must have at least one slash, which indicates
	 * a group name under the current group name. Get the ID of that group
	 */
	strncpy( gname, varname, idx_first_slash );
	gname[ idx_first_slash ] = '\0';
	ierr = nc_inq_grp_ncid( *ncid, gname, &gid );
	if( ierr != NC_NOERR ) {
		Rprintf( "Error in R_nc4_varid_hier: extracted groupname not found: >%s< (coding error?)\n", gname );
		*returned_varid = -1;
		*returned_grpid = -1;
		return;
		}

	/* Call ourselves recursively with the subgroup id and the varname stripped
	 * of its leading group name
	 */
/* Rprintf("R_nc4_inq_varid_hier_inner: calling myself recursively with varname=>%s<\n", varname+idx_first_slash+1 ); */
	R_nc4_inq_varid_hier_inner( &gid, varname+idx_first_slash+1, returned_grpid, returned_varid );
}

/****************************************************************************************/
/* This version of inq_varid works with hierarchial var names such as "model1/TS" and
 * returns both the group ID and the var ID in that group, given the root group id
 * (netcdf file ID).
 *
 * If the passed var is not found in the file, then -1 is returned.
 *
 * The passed slash-ful varname must not start or end with a slash.
 */
void R_nc4_inq_varid_hier( int *ncid, char **varname, int *returned_grpid, int *returned_varid )
{
/* Rprintf("R_nc4_inq_varid_hier: entering for var >%s<\n", varname[0] ); */
	R_nc4_inq_varid_hier_inner( ncid, varname[0], returned_grpid, returned_varid );
}

