AC_INIT(src/Makevars.in)

echo "configure.ac: starting"

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
        echo "Could not determine R_HOME"
	exit 1
fi

#---------------------------------------------------------------------------
# The way R's package installation mechanism works, the C compiler and flags
# are set by R, and not changable at the package level. So this is really
# just for information purposes. Changing these values will NOT result in 
# a different compiler being used!!
#---------------------------------------------------------------------------
CC=`${R_HOME}/bin/R CMD config CC`
CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`

# for ncdf4, I've chosen to REQUIRE the netcdf v4 interface
AC_PATH_NETCDF
if test $NETCDF_V4 = "no"; then
	echo "-------------------------------------------------------------------"
	echo "ncdf4 REQUIRES the netcdf library be version 4 or above,"
	echo "AND installed with HDF-5 support (i.e., the netcdf library must be"
	echo "compiled with the --enable-netcdf-4 flag). If you don't want to install"
	echo "the full version of netcdf-4 with HDF-5 support, then please install"
	echo "the old, deprecated ncdf package instead."
	echo "-------------------------------------------------------------------"
	exit -1
fi

HAVE_NETCDF=TRUE

NETCDF_LDFLAGS="$NETCDF_RPATH $NETCDF_LDFLAGS $LDFLAGS"

AC_SUBST(NETCDF_CPPFLAGS)
AC_SUBST(NETCDF_LDFLAGS)
AC_SUBST(NETCDF_RPATH)
AC_SUBST(HAVE_NETCDF)

AC_OUTPUT(src/Makevars)

echo " "
echo "**********************  Results of ncdf4 package configure *******************"
echo " "
echo "netCDF v4 CPP flags      = $NETCDF_CPPFLAGS"
echo "netCDF v4 LD flags       = $NETCDF_LDFLAGS"
echo "netCDF v4 runtime path   = $NETCDF_RPATH"
echo " "
echo "netCDF C compiler used   = $NETCDF_CC"
echo "R      C compiler used   = $CC $CFLAGS"
echo " "
echo "******************************************************************************"
echo " "

#---------------------------------------------------------------------
# Check to see if a different compiler is being used for netcdf and R
#---------------------------------------------------------------------
set $NETCDF_CC
NETCDF_CC_ARG1=$1
set $CC
R_CC_ARG1=$1
TEST_EMIT_WARNING=0
if test x$NETCDF_CC_ARG1 != x$R_CC_ARG1 -o $TEST_EMIT_WARNING -eq 1; then

	#------------------------------------------------------------------------------
	# By default, we emit a warning if the compilers of R and netcdf are different.
	# This happens with some regularity to naive users.
	#------------------------------------------------------------------------------
	emit_warning=1

	#-------------------------------------------------------------------------------------------------
	# Skip warning if one is gcc and the other is clang -- seems to bother some overly sensitive folks
	#-------------------------------------------------------------------------------------------------
	NETCDF_IS_GCC=0
	if echo "$NETCDF_CC_ARG1" | grep -q "gcc"; then
		NETCDF_IS_GCC=1
	else
		AC_PROG_CC([$NETCDF_CC_ARG1])
		if test "z$GCC" = zyes; then
			NETCDF_IS_GCC=1
		fi
	fi
	NETCDF_IS_CLANG=0
	if echo "$NETCDF_CC_ARG1" | grep -q "clang"; then
		NETCDF_IS_CLANG=1
	fi

	R_IS_GCC=0
	if echo "$R_CC_ARG1" | grep -q "gcc"; then
		R_IS_GCC=1
	else
		AC_PROG_CC([$R_CC_ARG1])
		if test "z$GCC" = zyes; then
			R_IS_GCC=1
		fi
	fi
	R_IS_CLANG=0
	if echo "$R_CC_ARG1" | grep -q "clang"; then
		R_IS_CLANG=1
	fi

	#-------------------------------------------------------------
	# We don't bother emitting the warning if the two compilers in
	# use are clang and gcc -- supposedly they are compatible
	#-------------------------------------------------------------
	if test $NETCDF_IS_CLANG -eq 1 -a $R_IS_GCC   -eq 1; then
		emit_warning=0
	fi
	if test $NETCDF_IS_GCC   -eq 1 -a $R_IS_CLANG -eq 1; then
		emit_warning=0
	fi

	#-----------------------------------------------------
	# Both same, but different exact names of compilers...
	# not sure if this actually guarantees success??
	#-----------------------------------------------------
	if test $NETCDF_IS_CLANG -eq 1 -a $R_IS_CLANG -eq 1; then
		emit_warning=0
	fi
	if test $NETCDF_IS_GCC -eq 1 -a $R_IS_GCC -eq 1; then
		emit_warning=0
	fi

	echo "netcdf cc arg1: $NETCDF_CC_ARG1  R cc arg1: $R_CC_ARG1 R_IS_GCC: $R_IS_GCC R_IS_CLANG: $R_IS_CLANG NETCDF_IS_GCC: $NETCDF_IS_GCC NETCDF_IS_CLANG: $NETCDF_IS_CLANG  emit_warning: $emit_warning"

	if test $emit_warning -eq 1 -o $TEST_EMIT_WARNING -eq 1; then
		echo " "
		echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
		echo "Dear ncdf4 user: The software configuration tool thinks that perhaps different compilers are being"
		echo "used for R and the netcdf library. However, we all know that software configuration tools are"
		echo "not that smart. If ncdf4 compiles fine for you, then please ignore this message. However, if"
		echo "installation fails, it might be that you are using different compilers for R and the netcdf library."
		echo "If so, then you must create a new version of the netcdf library that was compiled with the same"
		echo "compiler that was used to install R (see information above). Then, try installing the ncdf4"
		echo "package again using that new netcdf library. If necessary, you can specify the location of"
		echo "the new nc-config (the one that uses the same compiler used by R) like this:"
		echo " "
		echo "R CMD INSTALL --configure-args=\"--with-nc-config=/home/joe/bin/nc-config\""
		echo " "
		echo "(obviously, put in the correct path to YOUR nc-config, don't just use my example location!)"
		echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
		echo " "
	fi

fi

