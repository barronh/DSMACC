AC_INIT(DSMACC, 0.1)
test -z $FCFLAGS && FCFLAGS="-O3"
test -z $CFLAGS && CFLAGS="-I/usr/include/sys"
AC_PROG_LEX()
AC_PROG_CC()
case $CC in 
gcc)
  CFLAGS="$CFLAGS -m32"
  ;;
icc)
  CFLAGS="$CFLAGS -mia32"
  ;;
esac
AC_PROG_CPP()
AC_PROG_FC()
AC_CHECK_HEADERS([malloc.h])
AC_FC_PP_SRCEXT(f95)
AC_FC_CHECK_BOUNDS()
AC_FC_IMPLICIT_NONE()
AC_OUTPUT(Makefile.defs)
