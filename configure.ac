AC_INIT(DSMACC, 0.1)
test -z $FCFLAGS && FCFLAGS="-O3"
AC_PROG_LEX()
AC_PROG_CC()
AC_MSG_CHECKING([System])
UNAME_SYSTEM=`(uname -s) 2>/dev/null` || UNAME_SYSTEM=unknown
AC_MSG_RESULT([${UNAME_SYSTEM}])
case ${UNAME_SYSTEM} in
Darwin)
  CFLAGS="$CFLAGS -I/usr/include/sys"
  ;;
esac
case $CC in 
gcc)
  CFLAGS="$CFLAGS -m32"
  ;;
icc)
  CFLAGS="$CFLAGS -mia32"
  ;;
esac

AC_PROG_FC()

case $FC in 
gfortran)
  FCFLAGS="$FCFLAGS -fno-automatic"
  ;;
ifort)
  FCFLAGS="$FCFLAGS -noautomatic"
  ;;
esac
AC_CHECK_HEADERS([malloc.h])
AC_FC_PP_SRCEXT(f90)
AC_FC_CHECK_BOUNDS()
AC_FC_IMPLICIT_NONE()
AC_OUTPUT(Makefile.defs)
