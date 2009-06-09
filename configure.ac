#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

m4_ifdef([AC_OPENMP], [], [

# _AC_LANG_OPENMP
# ---------------
# Expands to some language dependent source code for testing the presence of
# OpenMP.
AC_DEFUN([_AC_LANG_OPENMP],
[_AC_LANG_DISPATCH([$0], _AC_LANG, $@)])

# _AC_LANG_OPENMP(C)
# ------------------
m4_define([_AC_LANG_OPENMP(C)],
[
#ifndef _OPENMP
 choke me
#endif
#include <omp.h>
int main () { return omp_get_num_threads (); }
])

# _AC_LANG_OPENMP(C++)
# --------------------
m4_copy([_AC_LANG_OPENMP(C)], [_AC_LANG_OPENMP(C++)])

# _AC_LANG_OPENMP(Fortran 77)
# ---------------------------
m4_define([_AC_LANG_OPENMP(Fortran 77)],
[AC_LANG_FUNC_LINK_TRY([omp_get_num_threads])])

# _AC_LANG_OPENMP(Fortran)
# ---------------------------
m4_copy([_AC_LANG_OPENMP(Fortran 77)], [_AC_LANG_OPENMP(Fortran)])

# AC_OPENMP
# ---------
# Check which options need to be passed to the C compiler to support OpenMP.
# Set the OPENMP_CFLAGS / OPENMP_CXXFLAGS / OPENMP_FFLAGS variable to these
# options.
# The options are necessary at compile time (so the #pragmas are understood)
# and at link time (so the appropriate library is linked with).
# This macro takes care to not produce redundant options if $CC $CFLAGS already
# supports OpenMP. It also is careful to not pass options to compilers that
# misinterpret them; for example, most compilers accept "-openmp" and create
# an output file called 'penmp' rather than activating OpenMP support.
AC_DEFUN([AC_OPENMP],
[
  OPENMP_[]_AC_LANG_PREFIX[]FLAGS=
  AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--disable-openmp], [do not enable OpenMP support])])
  if test "$enable_openmp" != no; then
    AC_CACHE_CHECK([for $CC option to support OpenMP],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
      [AC_LINK_IFELSE([_AC_LANG_OPENMP],
	 [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'],
	 [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
	  dnl Try these flags:
	  dnl   GCC >= 4.2           -fopenmp
	  dnl   SunPRO C             -xopenmp
	  dnl   Intel C              -openmp
	  dnl   SGI C, PGI C         -mp
	  dnl   Tru64 Compaq C       -omp
	  dnl   IBM C (AIX, Linux)   -qsmp=omp
	  dnl If in this loop a compiler is passed an option that it doesn't
	  dnl understand or that it misinterprets, the AC_LINK_IFELSE test
	  dnl will fail (since we know that it failed without the option),
	  dnl therefore the loop will continue searching for an option, and
	  dnl no output file called 'penmp' or 'mp' is created.
	  for ac_option in -fopenmp -xopenmp -openmp -mp -omp -qsmp=omp; do
	    ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
	    _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
	    AC_LINK_IFELSE([_AC_LANG_OPENMP],
	      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option])
	    _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
	    if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != unsupported; then
	      break
	    fi
	  done])])
    case $ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp in #(
      "none needed" | unsupported)
        ;; #(
      *)
        OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp ;;
    esac
  fi
  AC_SUBST([OPENMP_]_AC_LANG_PREFIX[FLAGS])
])

])


dnl AC_PATH_HEADER(variable,header,action-if-not-found,default-locations)
dnl Note that the above may set variable to an empty value if the header is
dnl already in the include path+6511
AC_DEFUN([AC_PATH_HEADER], [{
AC_MSG_CHECKING([for $2])
AC_CACHE_VAL(ac_cv_pheader_$1,
[
  ac_found=
  ac_dummy="ifelse([$4], , :/usr/local/include, [$4])"
  save_CPPFLAGS="$CPPFLAGS"
  IFS="${IFS= 	}"; ac_save_ifs="$IFS"; IFS=":"
  for ac_dir in $ac_dummy; do
    IFS="$ac_save_ifs"
    if test -z "$ac_dir"
    then
      CPPFLAGS="$save_CPPFLAGS"
    else
      CPPFLAGS="-I$ac_dir $save_CPPFLAGS"
    fi
    AC_TRY_COMPILE([#include <stdio.h> #include <math.h>
		    #include "$2"],[IDL_VPTR drift_dec;],ac_found=1;ac_cv_pheader_$1="$ac_dir";break)
  done
  CPPFLAGS="$save_CPPFLAGS"
  ifelse([$3],,,[if test -z "$ac_found"
    then
      $3
    fi
  ])
])
$1="$ac_cv_pheader_$1"
if test -n "$ac_found" -o -n "[$]$1"
then
  AC_MSG_RESULT([$]$1)
  IDL_HEADER="$ac_cv_pheader_$1"
else
  AC_MSG_ERROR(idl_export.h not found)
fi
AC_SUBST($1)
}])


AC_INIT(cleantools, 1.0.3, jdowell@astro.indiana.edu)
AC_CONFIG_SRCDIR([lib_version.c])
#AC_CONFIG_HEADER([config.h])


PACKAGE_DATE2=$(grep clean lib_version.c)
PACKAGE_DATE=$(echo $PACKAGE_DATE2 | sed 's/.*\([[0-9]]\{4\}\)\([[0-9]]\{2\}\)\([[0-9]]\{2\}\).*/\1\/\2\/\3/')
AC_SUBST(PACKAGE_DATE)
AC_MSG_NOTICE(configuring clean tools v$PACKAGE_VERSION; released $PACKAGE_DATE)

# Checks for programs.
AC_PROG_CC([icc gcc])
if test $CC = gcc; then
  AC_MSG_CHECKING([whether gcc version is >= 4.3])
  AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[
	#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3)
	 int maj_ver;
	#else
	 gcc is too old, libgomp cannot by dlopen()ed
	#endif]]), [
		AC_MSG_RESULT(yes)
		AC_OPENMP([C])],
  [if test "$enable_openmp" != no; then
     AC_MSG_RESULT(no)
     AC_MSG_WARN(gcc < 4.3; disabling OpenMP support)
   else
     AC_MSG_RESULT(no)
   fi])
else
  AC_OPENMP([C])
fi

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([m], [exp])
AC_CHECK_LIB([m], [pow])
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdio.h math.h string.h])

AC_ARG_WITH([idl-inc],
        AS_HELP_STRING([--with-idl-inc=PATH], [specify the location of the IDL include directory]),
        [
	 AC_PATH_HEADER([IDL_HEADER], [idl_export.h], 
		[AC_MSG_ERROR([cannot find idl_export.h; try --with-idl-inc=PATH])], 
		["$withval":"/usr/local/itt/idl/external/include":"/usr/local/rsi/idl/external/include"])], 
	[
	 AC_PATH_HEADER([IDL_HEADER], [idl_export.h], 
		[AC_MSG_ERROR([cannot find idl_export.h; try --with-idl-inc=PATH])], 
		["/usr/local/itt/idl/external/include":"/usr/local/rsi/idl/external/include"])])

# Checks for typedefs, structures, and compiler characteristics.


# Checks for library functions.


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
