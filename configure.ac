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
AC_PROG_YACC
if test -z $FCFLAGS; then
    AC_MSG_NOTICE([bison/byacc/yacc not found. If KPP will not install, you may need to install bison])
fi
AC_PROG_LEX
AC_PROG_CC
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lfl':
AC_CHECK_LIB([fl], [yywrap],,AC_MSG_ERROR([Flex library `libfl` is required by kpp and not found]), )

AC_MSG_CHECKING([System])
UNAME_SYSTEM=`(uname -s) 2>/dev/null` || UNAME_SYSTEM=unknown
AC_MSG_RESULT([${UNAME_SYSTEM}])

if test $cflags_update -gt 0; then
AC_MSG_NOTICE([Updating CFLAGS automatically KPP compilation; to disable use CFLAGS=" "])
# Adding compiler specific flags that have no macros in autoconf
case ${UNAME_SYSTEM}:${CC} in
Darwin:gcc)
  AC_MSG_NOTICE([Adding include path for Mac where malloc.h is in /usr/include/malloc])
  CFLAGS="$CFLAGS -arch i386 -I/usr/include/malloc"
  ;;
Darwin:*)
  AC_MSG_NOTICE([Adding 32-bit option for KPP compilation; to disable use CFLAGS=" "])
  CFLAGS="$CFLAGS -I/usr/include/malloc"
  ;;
*:gcc)
  AC_MSG_NOTICE([Adding 32-bit option for KPP compilation; to disable use CFLAGS=" "])
  CFLAGS="$CFLAGS -m32"
  ;;
esac
AC_MSG_NOTICE([${CFLAGS}])

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

AC_CHECK_HEADERS([malloc.h],,AC_MSG_ERROR([`malloc.h` is required and not found]), )
AC_CHECK_HEADERS([stdio.h stdlib.h string.h unistd.h])
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FC_PP_SRCEXT(f90)
AC_FC_CHECK_BOUNDS()
AC_FC_IMPLICIT_NONE()
AC_OUTPUT(Makefile.defs)
