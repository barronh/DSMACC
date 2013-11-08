AC_INIT(DSMACC, 0.1)
test -z $FCFLAGS && FCFLAGS="-O3"
AC_PROG_LEX()
AC_PROG_CC()
AC_MSG_CHECKING([System])
UNAME_SYSTEM=`(uname -s) 2>/dev/null` || UNAME_SYSTEM=unknown
AC_MSG_RESULT([${UNAME_SYSTEM}])

# Adding compiler specific flags that have no macros in autoconf
case ${UNAME_SYSTEM} in
Darwin)
  AC_MSG_CHECKING([Adding include path for Mac where malloc.h is in /usr/include/sys])
  CFLAGS="$CFLAGS -I/usr/include/sys"
  ;;
esac

AC_MSG_CHECKING([Adding 32-bit option for KPP compilation])
case $CC in 
gcc)
  CFLAGS="$CFLAGS -m32"
  ;;
icc)
  CFLAGS="$CFLAGS -mia32"
  ;;
pgc*)
  CFLAGS="$CFLAGS -pc 32"
  ;;
esac

AC_PROG_FC()

# Adding compiler specific flags that have no macros in autoconf
AC_MSG_CHECKING([Adding no automatic for $FC])
case $FC in 
gfortran)
  FCFLAGS="$FCFLAGS -fno-automatic"
  ;;
ifort)
  FCFLAGS="$FCFLAGS -noautomatic"
  ;;
pgf*)
  FCFLAGS="$FCFLAGS -save"
  ;;
esac

AC_CHECK_HEADERS([malloc.h])
AC_FC_PP_SRCEXT(f90)
AC_FC_CHECK_BOUNDS()
AC_FC_IMPLICIT_NONE()
AC_OUTPUT(Makefile.defs)
