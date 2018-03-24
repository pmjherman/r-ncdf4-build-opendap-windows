/*   #define DEBUG */

#include <stdio.h>
/*  #include <values.h> */

#include <R.h>
#include <Rinternals.h>
#include "netcdf.h"

SEXP R_ncu4_getListElement(SEXP list, char *str);
int R_ncu4_get_varsize( int ncid, int varid, int ndims, size_t *varsize );
SEXP R_nc4_get_vara_numvarid( SEXP sx_nc, SEXP sx_varid, SEXP sx_start, SEXP sx_count ) ;

/*********************************************************************************
 * Return the varid of the only var in the file, or -1 if more than one.
 */
int R_ncu4_varid_onlyvar( int ncid )
{
	int	ierr, nvars, varid, i, dimid;
	char	varname[MAX_NC_NAME];

	varid = -1;
	ierr = nc_inq_nvars( ncid, &nvars );
	if( ierr != NC_NOERR )
		error("Error reading from netcdf file!");

	for( i=0; i<nvars; i++ ) {
		ierr = nc_inq_varname( ncid, i, varname );
		if( ierr != NC_NOERR )
			error("Error reading from netcdf file!");
		ierr = nc_inq_dimid( ncid, varname, &dimid );
		if( ierr != NC_NOERR ) {
			/* Did NOT find this var as a dim, means it is NOT a dimvar */
			if( varid != -1 )
				/* More than one non-dimvar dim in this file */
				return( -1 );
			varid = i;
			}
		}
	return( varid );
}

/*================================================================================================================
 * Given as inputs start_arg and count_arg (which can be -1's for example) this calculates
 * the actual start and count to use, and returns them.
 */
void R_ncu4_calc_start_count( int ncid, int varid, int *start_arg, int len_start, 
	int *count_arg, int len_count, size_t *varsize,
	int ndims, size_t *start, size_t *count )
{
	int i, j, tmp[MAX_NC_DIMS], n_nondegen_dims;	

	/*
	Rprintf( "Passed count: [" );
	for( i=0; i<len_count; i++ ) {
		Rprintf( "%d", count_arg[i] );
		if( i < (len_count-1))
			Rprintf( "," );
		}
	Rprintf( "]\n" );
	*/

	/*---------------------------------------------------------------- 
	 * If start is a single '-1', then figure out actual start to use.
	 * Note: 'start_arg' is what was passed to this routine.
	 * 'start' is the actual start to use, which may be somewhat 
	 * different.
	 *---------------------------------------------------------------*/
	if( (len_start == 1) && (start_arg[0] == -1)) {
		/*----------------------------------------------------
		 * User did not specify a start arg -- just start at 0
		 *---------------------------------------------------*/
		for( i=0; i<ndims; i++ )
			start[i] = 0L;
		}
	else
		{
		/*-------------------------------------------------------
		 * User specified start ... switch from R to C convention.
		 * R convention is (xyzt) order with 1-based counting. C
		 * convention is (tzyx) order with 0-based counting.
		 *------------------------------------------------------*/
		for( i=0; i<len_start; i++ )
			tmp[i] = start_arg[len_start-i-1] - 1;
		for( i=0; i<len_start; i++ )
			start_arg[i] = tmp[i];

		/*-------------------------------------------
		 * Make sure passed start arg has enough dims 
		 *------------------------------------------*/
		if( len_start != ndims ) {
			/*-------------------------------------------------------------- 
			 * Well, allow 1 special case ... user can specify values just
			 * for non-degenerate dims (ones that have length>1).  To figure
			 * this out we need to get the vector of dimlengths for this var.
			 *-------------------------------------------------------------*/
			if( R_ncu4_get_varsize( ncid, varid, ndims, varsize ) == -1 )
				error("read of netcdf file failed when getting variable size");
			n_nondegen_dims = 0;
			for( i=0; i<ndims; i++ )
				if( varsize[i] > 1 )
					n_nondegen_dims++;
			if( len_start != n_nondegen_dims ) {
				error( "Error, passed argument 'start' has length %d, but must either have a length equal to the number of dimensions (%d) OR the number of non-degenerate dimensions (%d)\n",
					len_start, ndims, n_nondegen_dims );
				}
			/*-----------------------------------------------------------------
			 * If we get here, user has specified non-degen dimensions only ... 
			 * translate this to a start string we can use.
			 *----------------------------------------------------------------*/
			j = 0;
			for( i=0; i<ndims; i++ ) {
				if( varsize[i] == 1 )
					start[i] = 0L;
				else
					start[i] = start_arg[j++];
				}
			}
		else
			{
			/*---------------------------------------------------
			 * user specified enough dims ... just copy them over
			 *--------------------------------------------------*/
			for( i=0; i<ndims; i++ )
				start[i] = (size_t)start_arg[i];
			}
		}

	/*---------------------------------------------------------------- 
	 * If count is a single '-1', then figure out actual count to use.
	 * Note: 'count_arg' is what was passed to this routine.
	 * 'count' is the actual count to use, which may be somewhat 
	 * different.
	 *---------------------------------------------------------------*/
	if( (len_count == 1) && (count_arg[0] == -1)) {
		/*----------------------------------------------------
		 * User did not specify a count arg -- do entire var, 
		 * taking start into account.
		 *---------------------------------------------------*/
		for( i=0; i<ndims; i++ )
			count[i] = varsize[i] - start[i];
		}
	else
		{
		/*-------------------------------------------------------
		 * User specified count ... switch from R to C convention.
		 * R convention is (xyzt) order, C convention is (tzyx).
		 *------------------------------------------------------*/
		for( i=0; i<len_count; i++ )
			tmp[i] = count_arg[len_count-i-1];
		for( i=0; i<len_count; i++ )
			count_arg[i] = tmp[i];

		/*-------------------------------------------
		 * Make sure passed count arg has enough dims 
		 *------------------------------------------*/
		if( len_count != ndims ) {
			/*-------------------------------------------------------------- 
			 * Well, allow 1 special case ... user can specify values just
			 * for non-degenerate dims (ones that have length>1).  To figure
			 * this out we need to get the vector of dimlengths for this var.
			 *-------------------------------------------------------------*/
			n_nondegen_dims = 0;
			for( i=0; i<ndims; i++ )
				if( varsize[i] > 1 )
					n_nondegen_dims++;
			if( len_count != n_nondegen_dims ) {
				error( "Error, passed argument 'count' has length %d, but must either have a length equal to the number of dimensions (%d) OR the number of non-degenerate dimensions (%d)\n",
					len_count, ndims, n_nondegen_dims );
				}
			/*-----------------------------------------------------------------
			 * If we get here, user has specified non-degen dimensions only ... 
			 * translate this to a count string we can use.
			 *----------------------------------------------------------------*/
			j = 0;
			for( i=0; i<ndims; i++ ) {
				if( varsize[i] == 1 )
					count[i] = 1L;
				else
					count[i] = count_arg[j++];
				}
			}
		else
			{
			/*---------------------------------------------------
			 * user specified enough dims ... just copy them over
			 *--------------------------------------------------*/
			for( i=0; i<ndims; i++ )
				if( count_arg[i] == -1 )
					count[i] = varsize[i] - start[i];
				else
					count[i] = (size_t)count_arg[i];
			}
		}

	/*
	Rprintf( "Final count to use: [" );
	for( i=0; i<len_count; i++ ) {
		Rprintf( "%ld", count[i] );
		if( i < (len_count-1))
			Rprintf( "," );
		}
	Rprintf( "]\n" );
	*/
}


/***************************************************************************************
 * Given a characer variable name, this reads the data for that var from the file.
 * Return value: a list with two elements:
 *	[[1]]: an integer error code. '0' means no error.
 *	[[2]]: a real or integer array that contains the desired data.
 *
 */
SEXP R_nc4_get_vara_charvarid( SEXP sx_nc, SEXP sx_varid, SEXP sx_start, SEXP sx_count ) 
{
	char	errmess[1024];
	int	varid, ierr, ncid;
	SEXP	retval, sx_numvarid;
	const char *varname = CHAR(STRING_ELT(sx_varid,0)); 

	ncid    = INTEGER(R_ncu4_getListElement( sx_nc, "id" ))[0];

	ierr = nc_inq_varid( ncid, varname, &varid );
	if( ierr != NC_NOERR ) {
		sprintf( errmess, "the passed variable name [%s] does not exist in the file!",
			varname );
		error( errmess );
		}

	PROTECT( sx_numvarid = allocVector( INTSXP, 1 ));
	INTEGER(sx_numvarid)[0] = varid+1;	/* +1 to conform to R standards */

	retval = R_nc4_get_vara_numvarid( sx_nc, sx_numvarid, sx_start, sx_count );

	UNPROTECT(1);
	return( retval );
}

/***************************************************************************************
 * Given a numeric varid, this reads the data from the file.
 * Does not return on errors.
 */
SEXP R_nc4_get_vara_numvarid( SEXP sx_nc, SEXP sx_varid, SEXP sx_start, SEXP sx_count ) 
{
	int 	varid, ncid, ndims, len_start, len_count, i, j, ierr,
		start_arg[MAX_NC_DIMS], count_arg[MAX_NC_DIMS],
		*data_addr_i, missval_i, ndims_cgt1;
	SEXP 	rv_data = R_NilValue /* -Wall */, sx_ncdf_var, sx_dim;
	size_t	start[MAX_NC_DIMS], count[MAX_NC_DIMS], varsize[MAX_NC_DIMS], tot_var_size,
		i_szt;
	double	*data_addr_d, missval_d, missval_tol;
	nc_type	vartype;

	/*--------------------------------------------------------------------------- 
	 * On entry, the following are guaranteed to be integers:
	 *	varid
	 *	*start
	 *	*count
	 *
	 * Note that varid, start, and/or count could be a single '-1' if the user
	 * has not specified the start and count to use.
	 * 'sx_nc' is guaranteed to be the full object of class 'ncdf'.
	 *----------------------------------------------------------------------------*/


	varid = INTEGER(sx_varid)[0];
	ncid  = INTEGER(R_ncu4_getListElement( sx_nc, "id" ))[0];
	sx_ncdf_var = R_ncu4_getListElement( sx_nc, "var" );

	/*-----------------------------------------------------------------------
	 * Copy passed start and count to local vars so we can modify them safely
	 *----------------------------------------------------------------------*/
	len_start = length(sx_start);
	for( i=0; i<len_start; i++ )
		start_arg[i] = INTEGER(sx_start)[i];
	len_count = length(sx_count);
	for( i=0; i<len_count; i++ )
		count_arg[i] = INTEGER(sx_count)[i];
	
	/*-----------------------------------------
	 * Get varid to use, if passed value is -1.
	 *----------------------------------------*/
	if( varid == -1 ) {
		/*----------------------------------------------------
		 * Get how many vars are in this file ... if only one,
		 * use that one.  Otherwise, signal error.
		 *---------------------------------------------------*/
		varid = R_ncu4_varid_onlyvar( ncid );
		if( varid == -1 ) 
			error( "Error: no var specified, and the file has more than one valid var!" );
		}
	else
		varid--;	/* go from R to C indexing */
	
	/*--------------------------------------------------------
	 * Get # of dims for this var, as a check to make sure any
	 * passed 'start' and 'count' are correct.
	 *-------------------------------------------------------*/
	ierr = nc_inq_varndims( ncid, varid, &ndims );
	if( ierr != NC_NOERR )
		error( "Internal error in ncdf package, routine R_nc4_get_vara_numvarid: failed to get ndims for var!\n" );

	/*------------------------------------------------------
	 * Get our variable's size, and the start & count to use
	 *-----------------------------------------------------*/
	R_ncu4_get_varsize( ncid, varid, ndims, varsize );
	R_ncu4_calc_start_count( ncid, varid, start_arg, len_start, count_arg, len_count, 
			varsize, ndims, start, count );

	/*------------------------------------------------------------
	 * Allocate space for data, depending on the type of var it is
	 *-----------------------------------------------------------*/
	ierr = nc_inq_vartype( ncid, varid, &vartype );
	if( ierr != NC_NOERR )
		error( "Internal error in ncdf package, routine R_nc4_get_vara_numvarid: failed to get type for var!\n" );

	tot_var_size = 1L;
	for( i=0; i<ndims; i++ ) {
		tot_var_size *= count[i];
		}

	switch( vartype ) {
		case NC_CHAR:
			error( "chars not handled yet, use old interface" );
			break;

		case NC_BYTE:
		case NC_SHORT:
		case NC_INT:
			/*---------------
			 * Allocate space
			 *--------------*/
			PROTECT( rv_data = allocVector( INTSXP, tot_var_size ));
			data_addr_i = &(INTEGER(rv_data)[0]);	/* Is this guaranteed to work?  Dunno. */

			/*--------------
			 * Read the data
			 *-------------*/
			ierr        = nc_get_vara_int( ncid, varid, start, count, data_addr_i );
			if( ierr != NC_NOERR )
				error( "Error while trying to read int data from file!" );

			/*---------------------
			 * Handle missing value
			 *--------------------*/
			ierr = nc_get_att_int( ncid, varid, "missing_value", &missval_i );
			if( ierr != NC_NOERR )
				/* No missing value attribute found, use default value */
				missval_i = NC_FILL_INT;
			for( i_szt=0L; i_szt<tot_var_size; i_szt++ ) 
				if( data_addr_i[i_szt] == missval_i )
					data_addr_i[i_szt] = NA_INTEGER;
			break;

		case NC_FLOAT:
		case NC_DOUBLE:
			/*---------------
			 * Allocate space
			 *--------------*/
			PROTECT( rv_data = allocVector( REALSXP, tot_var_size ));
			data_addr_d = &(REAL(rv_data)[0]);	/* Is this guaranteed to work?  Dunno. */

			/*--------------
			 * Read the data
			 *-------------*/
			ierr        = nc_get_vara_double( ncid, varid, start, count, data_addr_d );
			if( ierr != NC_NOERR )
				error( "Error while trying to read real data from file!" );

			/*---------------------
			 * Handle missing value
			 *--------------------*/
			ierr = nc_get_att_double( ncid, varid, "missing_value", &missval_d );
			if( ierr != NC_NOERR )
				/* No missing value attribute found, use default value */
				missval_d = 1.e30;
			missval_tol = 1.e-5*fabs(missval_d);
			for( i_szt=0L; i_szt<tot_var_size; i_szt++ ) 
				if( fabs(data_addr_d[i_szt] - missval_d) < missval_tol )
					data_addr_d[i_szt] = NA_REAL;
			break;

		default:
			error( "unhandled var type when allocating var space in R_nc4_get_vara_numvarid");
		}

	/*-----------------------------------------
	 * Set our dims (note: non-degenerate only)
	 *----------------------------------------*/
	ndims_cgt1 = 0;  
	for( i=0; i<ndims; i++ )
		if( count[i] > 1 )
			ndims_cgt1++;
	if( ndims_cgt1 == 0 ) {
		PROTECT( sx_dim = allocVector( INTSXP, 1 ));
		INTEGER(sx_dim)[0] = 1;
		}
	else
		{
		PROTECT( sx_dim = allocVector( INTSXP, ndims_cgt1 ));
		j = 0;
		for( i=0; i<ndims; i++ )
			if( count[i] > 1 ) {
				INTEGER(sx_dim)[ndims_cgt1-j-1] = count[i];
				j++;
				}
		}
	setAttrib( rv_data, R_DimSymbol, sx_dim );

	UNPROTECT(2);
	return(rv_data);
}

/*******************************************************************************
 * Internal utility function to get the vector of the variable's dim sizes.
 * Returns 0 on success, -1 if an error is encountered.
 */
int R_ncu4_get_varsize( int ncid, int varid, int ndims, size_t *varsize )
{
	int ierr, i, dimids[MAX_NC_DIMS];
	size_t len;

	ierr = nc_inq_vardimid( ncid, varid, dimids );
	if( ierr != NC_NOERR ) {
		error( "Internal error in ncdf package, routine R_ncu4_get_varsize: error while reading file to get var's dimids!\n" );
		return(-1);
		}

	for( i=0; i<ndims; i++ ) {
		ierr = nc_inq_dimlen( ncid, dimids[i], &len );
		if( ierr != NC_NOERR ) {
			error( "Internal error in ncdf package, routine R_ncu4_get_varsize: error while reading file to get dim's length!\n" );
			return(-1);
			}
		varsize[i] = len;
		}

	return(0);
}

/*******************************************************************************
 * Internal utility function that returns '1' if the passed var name is
 * the name of a dimvar, '0' if it is NOT the name of a dimvar, and '-1' on error.
 */
int R_ncu4_isdimvar( int ncid, char *name ) 
{
	int 	i, ndims, ierr;
	char 	dimname[MAX_NC_NAME];

	ierr = nc_inq_ndims( ncid, &ndims );
	if( ierr != NC_NOERR ) {
		error( "Internal error in ncdf package, routine R_ncu4_isdimvar: error while reading file to get ndims!\n" );
		return( -1 );
		}

	for( i=0; i<ndims; i++ ) {
		ierr = nc_inq_dimname( ncid, i, dimname );
		if( ierr != NC_NOERR ) {
			error( "Internal error in ncdf package, routine R_ncu4_isdimvar: error while reading file to get dim name!\n" );
			return( -1 );
			}
		if( strcmp( name, dimname ) == 0 )
			return(1);
		}

	return(0);
}

/*********************************************************************************
 * Internal utility function to get an element from a list and return it,
 * given the name.
 */
SEXP R_ncu4_getListElement(SEXP list, char *str)
{
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;

/*Rprintf( "list has %d elements; looking for one named %s\n", length(list), str);*/
	for (i = 0; i < length(list); i++) {
/*Rprintf( "element %d is named %s\n", i, CHAR(STRING_ELT(names, i)) );*/
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			return(VECTOR_ELT(list, i));
			}
		}

Rprintf( "warning, no match found for element %s\n", str );
	return elmt;
}

/*********************************************************************************
 * Read vlen strings given the numeric varid, start, and count to use
 */
SEXP R_nc4_get_vara_string( SEXP sx_nc, SEXP sx_varid, SEXP sx_start, SEXP sx_count ) 
{
	SEXP	sx_retval, sx_retnames, sx_retstrings, sx_reterror;
	int	i, ierr, nchar, varid, ncid, ndims, count_int[MAX_NC_DIMS], start_int[MAX_NC_DIMS], len_count, len_start; 
	size_t	count[MAX_NC_DIMS], start[MAX_NC_DIMS], tot_count, isz;
	char 	**ss;

	/* Convert passed parameters (which are in R format) into C format */
	ncid  = INTEGER(sx_nc   )[0];
	varid = INTEGER(sx_varid)[0];

	len_start = length(sx_start);
	for( i=0; i<len_start; i++ ) {
		start_int[i] = INTEGER(sx_start)[i];
		start    [i] = (size_t)(start_int[i]);
		}

	len_count = length(sx_count);
	for( i=0; i<len_count; i++ ) {
		count_int[i] = INTEGER(sx_count)[i];
		count    [i] = (size_t)(count_int[i]);
		}

	PROTECT( sx_retval   = allocVector( VECSXP, 2 ));       /* 2 elements in the returned list: $error, $strings */

	/* Set the names for the returned list */
	PROTECT( sx_retnames = allocVector( STRSXP, 2 ));       /* 2 elements in the returned list */
	SET_STRING_ELT( sx_retnames, 0, mkChar("error") );
	SET_STRING_ELT( sx_retnames, 1, mkChar("data") );
	setAttrib( sx_retval, R_NamesSymbol, sx_retnames );
	UNPROTECT(1); 

	/* Set the return error code */
	PROTECT( sx_reterror = allocVector( INTSXP, 1 ));
	INTEGER( sx_reterror)[0] = 0;

	/* Get number of dims in the var */
	ierr = nc_inq_varndims( ncid, varid, &ndims );

	/*--------------------------------------------------------------
	 * At this point we have all the C values we need:
	 *  	ncid
	 *	varid (numeric)
	 *	ndims
	 *	start
	 *	count
	 *---------------------------------------------------------------*/
	tot_count = 1L;
	for( i=0; i<ndims; i++ ) 
	 	tot_count *= count[i];

	ss = (char **)malloc( sizeof( char *) * tot_count );
	if( ss == NULL ) {
		INTEGER( sx_reterror)[0] = -1;
		error("ncdf4 library: routine R_nc4_get_vara_string: Error trying to allocate space to read the vlen strings: total count of strings requested: %ld\n", tot_count );
		}

	if( (ierr = nc_get_vara_string( ncid, varid, start, count, ss )) != 0 ) {
		INTEGER( sx_reterror)[0] = -2;
		error("ncdf4 library: routine R_nc4_get_vara_string: Error reading vlen strings: %s\n",
			nc_strerror(ierr));
		}

	PROTECT( sx_retstrings = allocVector( STRSXP, tot_count ));
	for( isz=0L; isz<tot_count; isz++ ) {
		nchar = strlen( ss[isz] );
		SET_STRING_ELT( sx_retstrings, isz, mkChar(ss[isz]) );
		}

	SET_VECTOR_ELT( sx_retval, 0, sx_reterror   );
	SET_VECTOR_ELT( sx_retval, 1, sx_retstrings );

	UNPROTECT(3);	

	/* Free netcdf string storage */
	nc_free_string( tot_count, ss );

	return( sx_retval );
}

