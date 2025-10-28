dnl This file is part of rseg
dnl
dnl Copyright (C) 2025: Andrew D. Smith
dnl
dnl Authors: Andrew D. Smith
dnl
dnl This is free software: you can redistribute it and/or modify it under the
dnl terms of the GNU General Public License as published by the Free Software
dnl Foundation, either version 3 of the License, or (at your option) any later
dnl version.
dnl
dnl This software is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
dnl for more details.

AC_INIT([rseg], [2.0.0], [andrewds@usc.edu],
        [rseg], [https://github.com/smithlabcode/rseg])
dnl the config.h is #included in the sources for version info
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects foreign])

AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
AC_PROG_CXX

# Set CXXFLAGS only if the user hasn't set them
AS_IF([test "x$CXXFLAGS_set" != "xset"], [
  CXXFLAGS="-O3 -DNDEBUG"
])
AC_MSG_NOTICE([CXXFLAGS is $CXXFLAGS])

AX_CXX_COMPILE_STDCXX(23, [noext], [mandatory])
AC_PROG_RANLIB

dnl recursively configure smithlab_cpp
AC_CONFIG_SUBDIRS([src/smithlab_cpp])

dnl check for HTSLib if requested
hts_fail_msg="
Failed to locate HTSLib on your system. Please use the LDFLAGS and CPPFLAGS
variables to specify the directories where the HTSLib library and headers can
be found.
"

gsl_fail_msg="
Failed to locate the GNU Scientific Library (GSL) on your system. Please use
the LDFLAGS and CPPFLAGS variables to specify the directories where the GSL
library and headers can be found.
"

dnl arg for using libdeflate, which might happen by default anyway
AC_ARG_WITH([libdeflate],
  [AS_HELP_STRING([--with-libdeflate], [use libdeflate for BAM output])],
  [with_libdeflate=yes], [with_libdeflate=no])

dnl check for required libraries
AC_SEARCH_LIBS([pthread_create], [pthread], [], [AC_MSG_FAILURE(["pthread library not found"])])
AC_SEARCH_LIBS([gzopen], [z], [], [AC_MSG_FAILURE(["Zlib library not found"])])
AS_IF([test "x$with_libdeflate" = "xyes"],
  [
    AC_SEARCH_LIBS([libdeflate_deflate_compress], [deflate], [],
                   [AC_MSG_ERROR([--with-libdeflate specified but libdeflate not found])])
  ])
AC_SEARCH_LIBS([hts_version], [hts], [], [AC_MSG_FAILURE([$hts_fail_msg])])
AC_SEARCH_LIBS([cblas_dgemm], [gslcblas], [], [AC_MSG_FAILURE([$gsl_fail_msg])])
AC_SEARCH_LIBS([gsl_blas_dgemm], [gsl], [], [AC_MSG_FAILURE([$gsl_fail_msg])])

AC_CONFIG_FILES([
Makefile
])

AC_OUTPUT
