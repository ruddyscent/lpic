# Process this file with autoconf to produce a configure script.
AC_INIT([lpic], [1.3.1], [pic@lichters.net, kemp@physics.unr.edu], [lpic])
AC_PREREQ(2.50)
AC_CONFIG_SRCDIR([./])
AC_CONFIG_AUX_DIR([config])
AM_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE([foreign])

AC_COPYRIGHT([WARANTY NOTICE AND COPYRIGHT

Copyright (C) 1994-1997 Roland Lichters

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software 
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

])

# Help the compiler find external software

AC_ARG_WITH([pvm-include],
            AC_HELP_STRING([--with-pvm-include=INCLUDE_PATH],
                           [Supply the location of pvm header directory]),
            [lpic_pvm_include_path=$withval],
            [lpic_pvm_include_path=""])
if test [ -n "$lpic_pvm_include_path" ] ; then
   AC_SUBST([CPPFLAGS],["${CPPFLAGS} -I${lpic_pvm_include_path}"])
fi
AC_ARG_WITH([pvm-lib],
            AC_HELP_STRING([--with-pvm-lib=LIB_PATH],
                           [Supply the location of the pvm library]),
            [lpic_pvm_lib_path=$withval],
            [lpic_pvm_lib_path=""])
if test [ -n "$lpic_pvm_lib_path" ] ; then
   AC_SUBST([LDFLAGS],["${LDFLAGS} -L${lpic_pvm_lib_path}"])
fi


# Language
AC_LANG([C++])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_AWK

# Default compilation variables
#AC_SUBST([CXXFLAGS],[${CXXFLAGS='-g -O2 -Wall'}])
AC_SUBST([CXXFLAGS],[${CXXFLAGS='-g -O2 -Wno-deprecated'}])
AC_SUBST([CPPFLAGS],[${CPPFLAGS='-g -O2 -Wno-deprecated'}])

# Check for mandatory features and header files

# Checks for header files.
#AC_HEADER_STDC
#AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
#AC_C_CONST
#AC_C_INLINE
#AC_TYPE_SIZE_T
#AC_STRUCT_TM

# Checks for library functions.
#AC_FUNC_ERROR_AT_LINE
#AC_FUNC_VPRINTF
#AC_CHECK_FUNCS([floor pow sqrt strstr strtol])

# Check for tools needed for building documentation
AC_PATH_PROG([DOXYGEN], [doxygen])
AC_PATH_PROG([SED], [sed])
AC_PATH_PROG([LATEX], [latex])
AC_PATH_PROG([PDFLATEX], [pdflatex])

# site configuration

AC_ARG_ENABLE([pvm],
              AC_HELP_STRING([--enable-pvm],
                             [If enabled, a parallel version of lpic using pvm
                              will be built.]),
              [lpic_pvm=$enableval],
              [lpic_pvm=no])
if test "$lpic_pvm" = "yes" ; then
   # check pvm header and library by compiling a one-liner  
   # see acinclude.m4
   LPIC_CHECK_PVM
   LPIC_CHECK_FUNC([pvm_parent], [], [pvm3.h])

   AC_DEFINE([LPIC_PARALLEL],[1],
             [Define this if a parallelized version of lpic should
              be built.)])
   AC_DEFINE([LPIC_PVM],[1],
             [Define this if a parallelized version of lpic using
              PVM should be built.)])
fi
AM_CONDITIONAL(PVM, test x$lpic_pvm = xyes)

AC_ARG_ENABLE([mpi],
              AC_HELP_STRING([--enable-mpi],
                             [If enabled, a parallel version of lpic using mpi
                              will be built.]),
              [lpic_mpi=$enableval],
              [lpic_mpi=no])
if test "$lpic_mpi" = "yes" ; then
   # check pvm header and library by compiling a one-liner  
   # see acinclude.m4
   LPIC_CHECK_MPI
   # todo:
   #LPIC_CHECK_FUNC([pvm_parent], [], [pvm3.h])

   AC_DEFINE([LPIC_PARALLEL],[1],
             [Define this if a parallelized version of lpic should
              be built.)])
   AC_DEFINE([LPIC_MPI],[1],
             [Define this if a parallelized version of lpic using
              MPI should be built.)])
fi
AM_CONDITIONAL(MPI, test x$lpic_mpi = xyes)

AC_ARG_ENABLE([tracing],
              AC_HELP_STRING([--enable-tracing],
                             [If enabled, tracing messages might be emitted
                              by the library depending on run-time
                              settings. Enabling this option can degrade
                              performance.]),
              [rl_tracing=$enableval],
              [rl_tracing=no])
if test "$lpic_tracing" = "yes" ; then
   AC_DEFINE([LPIC_TRACING],[1],
             [Define this if tracing messages should allowed (whether they
              are actually emitted will depend on run-time settings.)])
fi

# manual configurations for specific hosts
case $host in
  powerpc-apple-darwin*)
    AC_DEFINE([LPIC_PATCH_DARWIN],[1],
              [Define if running on a Mac OS X machine.])
    ;;
  sparc-sun-solaris*)
    AC_DEFINE([LPIC_PATCH_SOLARIS],[1],
              [Define if running on a Sun Solaris machine.])
    ;;
esac

# done, output the configured files
AC_CONFIG_FILES([
    Makefile
    config/Makefile
    lpic/Makefile
    lpic/src/Makefile
    post/Makefile
    post/src/Makefile
    fresnel/Makefile])
AC_OUTPUT
