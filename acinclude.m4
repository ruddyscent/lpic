# LPIC_TO_UPPER(STRING)
# ---------------------
# Convert the passed string to uppercase
AC_DEFUN([LPIC_TO_UPPER],
[translit([$1],[abcdefghijklmnopqrstuvwxyz.],[ABCDEFGHIJKLMNOPQRSTUVWXYZ_])])

# Check whether the pvm headers are available
AC_DEFUN([LPIC_CHECK_PVM],
[AC_MSG_CHECKING([for pvm header])
 AC_TRY_COMPILE(
    [@%:@include <pvm3.h>],
    [],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([pvm header not found])])
])

# Check whether the mpi headers are available
AC_DEFUN([LPIC_CHECK_MPI],
[AC_MSG_CHECKING([for mpi header])
 AC_TRY_COMPILE(
    [@%:@include <mpi.h>],
    [],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([mpi header not found])])
])

# Check whether FUNCTION (including namespace qualifier) is available
# when compiling after including HEADER. 
AC_DEFUN([LPIC_CHECK_FUNC],
[AC_MSG_CHECKING([for compilation with $1])
 AC_TRY_COMPILE(
    [@%:@include <$3>],
    [$1($2);],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([function $1 not found])])
])
