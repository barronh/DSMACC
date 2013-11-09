AC_INIT(DSMACC, 0.1)
if test -z $FCFLAGS; then
     FCFLAGS="-O3"
     fflags_update=1
else
     fflags_update=0
fi
if test -z $CFLAGS; then
     cflags_update=1
else
     cflags_update=0
fi
AC_PROG_LEX()
AC_PROG_CC()
AC_MSG_CHECKING([System])
UNAME_SYSTEM=`(uname -s) 2>/dev/null` || UNAME_SYSTEM=unknown
AC_MSG_RESULT([${UNAME_SYSTEM}])

# Adding compiler specific flags that have no macros in autoconf
case ${UNAME_SYSTEM} in
Darwin)
  AC_MSG_NOTICE([Adding include path for Mac where malloc.h is in /usr/include/sys])
  test $cflags_update -gt 0 && CFLAGS="$CFLAGS -I/usr/include/sys"
  ;;
esac
if test $cflags_update -gt 0; then
AC_MSG_NOTICE([Adding 32-bit option for KPP compilation; to disable use CFLAGS=" "])
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
fi
AC_PROG_FC()

# Adding compiler specific flags that have no macros in autoconf
if test $fflags_update -gt 0; then
AC_MSG_NOTICE([Adding no automatic for $FC; to disable use FCFLAGS=" "])
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

fi
AC_CHECK_HEADERS([malloc.h])
AC_FC_PP_SRCEXT(f90)
AC_FC_CHECK_BOUNDS()
AC_FC_IMPLICIT_NONE()
AC_OUTPUT(Makefile.defs)
