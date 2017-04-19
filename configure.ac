AC_INIT(DSMACC, 0.1)

if test -z $FCFLAGS; then
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
AC_CHECK_PROG([FOUND_YACC],[$YACC],[yes],[AC_MSG_NOTICE([bison/byacc/yacc not found. If KPP will not install, you may need to install bison])])
AC_PROG_LEX
AC_CHECK_PROG([FOUND_LEX],[$LEX],[yes],[AC_MSG_NOTICE([flex/lex not found. If KPP will not install, you may need to install lex])])
AC_PROG_CC
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lfl':
AC_CHECK_LIB([fl], [yywrap],,AC_MSG_ERROR([Flex library `libfl` is required by kpp and not found]), )

AC_MSG_CHECKING([System])
UNAME_SYSTEM=`(uname -s) 2>/dev/null` || UNAME_SYSTEM=unknown
AC_MSG_RESULT([${UNAME_SYSTEM}])

if test $cflags_update -eq 0; then
  AC_MSG_NOTICE([You have disabled automatic FCFLAGS by setting your own. Unset FCFLAGS to allow configure to automatically manage flags.])
else
AC_MSG_NOTICE([Updating CFLAGS automatically KPP compilation; to disable use CFLAGS=" "])
# Adding compiler specific flags that have no macros in autoconf
case ${UNAME_SYSTEM}:${CC} in
Darwin:gcc)
  AC_MSG_NOTICE([Adding -arch i386 option to prevent failure of test])
  CFLAGS="$CFLAGS -arch i386"
  AC_MSG_NOTICE([Adding include path for Mac where malloc.h is in /usr/include/malloc])
  CFLAGS="$CFLAGS -I/usr/include/malloc"
  ;;
Darwin:*)
  AC_MSG_NOTICE([Adding malloc subfolder for mac support.])
  CFLAGS="$CFLAGS -I/usr/include/malloc"
  ;;
esac
AC_MSG_NOTICE([${CFLAGS}])

fi
AC_PROG_FC()

# Adding compiler spec	ic flags that have no macros in autoconf
if test $fflags_update -eq 0; then
  AC_MSG_NOTICE([You have disabled automatic FCFLAGS by setting your own. Unset FCFLAGS to allow configure to automatically manage flags.])
else
  AC_MSG_NOTICE([Updating FCFLAGS automatically; to disable set your own FCFLAGS])
  AC_MSG_NOTICE([Adding no automatic for $FC; to disable set your own FCFLAGS])
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

# Checks for header files.
AC_FUNC_ALLOCA
AC_CHECK_HEADERS([malloc.h],,AC_MSG_ERROR([`malloc.h` is required and not found]), )
AC_CHECK_HEADERS([stdio.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor ftime getcwd memmove memset pow sqrt strchr])

AC_FC_PP_SRCEXT(f90)
AC_FC_CHECK_BOUNDS()
AC_FC_IMPLICIT_NONE()
AC_OUTPUT(Makefile.defs)
