/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     JACOBIAN = 258,
     DOUBLE = 259,
     FUNCTION = 260,
     DEFVAR = 261,
     DEFRAD = 262,
     DEFFIX = 263,
     SETVAR = 264,
     SETRAD = 265,
     SETFIX = 266,
     HESSIAN = 267,
     STOICMAT = 268,
     STOCHASTIC = 269,
     DECLARE = 270,
     INITVALUES = 271,
     EQUATIONS = 272,
     LUMP = 273,
     INIEQUAL = 274,
     EQNEQUAL = 275,
     EQNCOLON = 276,
     LMPCOLON = 277,
     LMPPLUS = 278,
     SPCPLUS = 279,
     SPCEQUAL = 280,
     ATOMDECL = 281,
     CHECK = 282,
     CHECKALL = 283,
     REORDER = 284,
     MEX = 285,
     DUMMYINDEX = 286,
     EQNTAGS = 287,
     LOOKAT = 288,
     LOOKATALL = 289,
     TRANSPORT = 290,
     TRANSPORTALL = 291,
     MONITOR = 292,
     USES = 293,
     SPARSEDATA = 294,
     WRITE_ATM = 295,
     WRITE_SPC = 296,
     WRITE_MAT = 297,
     WRITE_OPT = 298,
     INITIALIZE = 299,
     XGRID = 300,
     YGRID = 301,
     ZGRID = 302,
     USE = 303,
     LANGUAGE = 304,
     INTFILE = 305,
     DRIVER = 306,
     RUN = 307,
     INLINE = 308,
     ENDINLINE = 309,
     PARAMETER = 310,
     SPCSPC = 311,
     INISPC = 312,
     INIVALUE = 313,
     EQNSPC = 314,
     EQNSIGN = 315,
     EQNCOEF = 316,
     RATE = 317,
     LMPSPC = 318,
     SPCNR = 319,
     ATOMID = 320,
     LKTID = 321,
     MNIID = 322,
     INLCTX = 323,
     INCODE = 324,
     SSPID = 325,
     EQNLESS = 326,
     EQNTAG = 327,
     EQNGREATER = 328,
     TPTID = 329,
     USEID = 330
   };
#endif
/* Tokens.  */
#define JACOBIAN 258
#define DOUBLE 259
#define FUNCTION 260
#define DEFVAR 261
#define DEFRAD 262
#define DEFFIX 263
#define SETVAR 264
#define SETRAD 265
#define SETFIX 266
#define HESSIAN 267
#define STOICMAT 268
#define STOCHASTIC 269
#define DECLARE 270
#define INITVALUES 271
#define EQUATIONS 272
#define LUMP 273
#define INIEQUAL 274
#define EQNEQUAL 275
#define EQNCOLON 276
#define LMPCOLON 277
#define LMPPLUS 278
#define SPCPLUS 279
#define SPCEQUAL 280
#define ATOMDECL 281
#define CHECK 282
#define CHECKALL 283
#define REORDER 284
#define MEX 285
#define DUMMYINDEX 286
#define EQNTAGS 287
#define LOOKAT 288
#define LOOKATALL 289
#define TRANSPORT 290
#define TRANSPORTALL 291
#define MONITOR 292
#define USES 293
#define SPARSEDATA 294
#define WRITE_ATM 295
#define WRITE_SPC 296
#define WRITE_MAT 297
#define WRITE_OPT 298
#define INITIALIZE 299
#define XGRID 300
#define YGRID 301
#define ZGRID 302
#define USE 303
#define LANGUAGE 304
#define INTFILE 305
#define DRIVER 306
#define RUN 307
#define INLINE 308
#define ENDINLINE 309
#define PARAMETER 310
#define SPCSPC 311
#define INISPC 312
#define INIVALUE 313
#define EQNSPC 314
#define EQNSIGN 315
#define EQNCOEF 316
#define RATE 317
#define LMPSPC 318
#define SPCNR 319
#define ATOMID 320
#define LKTID 321
#define MNIID 322
#define INLCTX 323
#define INCODE 324
#define SSPID 325
#define EQNLESS 326
#define EQNTAG 327
#define EQNGREATER 328
#define TPTID 329
#define USEID 330




/* Copy the first part of user declarations.  */
#line 37 "scan.y"

  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>
  #include <unistd.h>
  #include "scan.h"

  #define __YYSCLASS

  #define YYDEBUG 1
  extern char yytext[];
  extern FILE * yyin;
  
  int nError   = 0;
  int nWarning = 0;

  int crt_section;
  int eqState;
  int isPhoto = 0;

  char crt_term[ 30 ];
  char crt_coef[ 30 ];

  char * InlineBuf;
  int InlineLen;

  void SemicolonError();
  int yyerrflag=0;

  void ParserErrorMessage();
  void yyerror(char *);



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 71 "scan.y"
{
  char str[200];
}
/* Line 193 of yacc.c.  */
#line 284 "y.tab.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 297 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  124
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   192

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  77
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  35
/* YYNRULES -- Number of rules.  */
#define YYNRULES  111
/* YYNRULES -- Number of states.  */
#define YYNSTATES  202

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   330

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    76,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     8,    11,    14,    17,    20,    23,
      26,    29,    32,    35,    38,    41,    44,    47,    50,    53,
      56,    59,    62,    65,    68,    71,    74,    77,    80,    83,
      85,    87,    89,    91,    93,    95,    97,   100,   103,   106,
     109,   112,   115,   120,   123,   126,   129,   132,   135,   138,
     141,   143,   147,   150,   153,   155,   159,   162,   165,   167,
     171,   174,   177,   179,   183,   186,   189,   191,   195,   198,
     201,   203,   207,   210,   213,   215,   219,   222,   225,   227,
     229,   233,   235,   239,   241,   244,   246,   250,   253,   256,
     260,   264,   267,   270,   275,   279,   282,   284,   288,   291,
     294,   298,   301,   303,   306,   308,   312,   315,   318,   322,
     326,   329
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      78,     0,    -1,    79,    -1,    79,    78,    -1,     3,    55,
      -1,    12,    55,    -1,    15,    55,    -1,    13,    55,    -1,
       4,    55,    -1,    29,    55,    -1,    30,    55,    -1,    31,
      55,    -1,    32,    55,    -1,     5,    55,    -1,    14,    55,
      -1,    26,    81,    -1,    27,    81,    -1,     6,    93,    -1,
       7,    93,    -1,     8,    93,    -1,     9,    91,    -1,    10,
      91,    -1,    11,    91,    -1,    16,    99,    -1,    17,   101,
      -1,    18,   109,    -1,    33,    83,    -1,    37,    85,    -1,
      35,    87,    -1,    28,    -1,    34,    -1,    36,    -1,    40,
      -1,    41,    -1,    42,    -1,    43,    -1,    48,    55,    -1,
      49,    55,    -1,    44,    55,    -1,    45,    55,    -1,    46,
      55,    -1,    47,    55,    -1,    53,    68,   111,    54,    -1,
      53,     1,    -1,    50,    55,    -1,    51,    55,    -1,    52,
      55,    -1,    38,    89,    -1,    39,    55,    -1,    80,    76,
      -1,    76,    -1,    81,    82,    80,    -1,    82,    80,    -1,
       1,    80,    -1,    65,    -1,    83,    84,    80,    -1,    84,
      80,    -1,     1,    80,    -1,    66,    -1,    85,    86,    80,
      -1,    86,    80,    -1,     1,    80,    -1,    67,    -1,    87,
      88,    80,    -1,    88,    80,    -1,     1,    80,    -1,    74,
      -1,    89,    90,    80,    -1,    90,    80,    -1,     1,    80,
      -1,    75,    -1,    91,    92,    80,    -1,    92,    80,    -1,
       1,    80,    -1,    70,    -1,    93,    94,    80,    -1,    94,
      80,    -1,     1,    80,    -1,    95,    -1,    96,    -1,    56,
      25,    97,    -1,    56,    -1,    97,    24,    98,    -1,    98,
      -1,    64,    56,    -1,    56,    -1,    99,   100,    80,    -1,
     100,    80,    -1,     1,    80,    -1,    57,    19,    58,    -1,
     101,   102,    80,    -1,   102,    80,    -1,     1,    80,    -1,
     104,   105,   106,   103,    -1,   105,   106,   103,    -1,    62,
     103,    -1,    62,    -1,    71,    72,    73,    -1,   107,    20,
      -1,   107,    21,    -1,   107,    60,   108,    -1,    60,   108,
      -1,   108,    -1,    61,    59,    -1,    59,    -1,   109,   110,
      80,    -1,   110,    80,    -1,     1,    80,    -1,    63,    23,
     110,    -1,    63,    22,    63,    -1,   111,    69,    -1,    69,
      -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    95,    95,    96,    98,   101,   104,   107,   110,   113,
     116,   119,   122,   125,   128,   131,   133,   135,   137,   139,
     141,   143,   145,   147,   149,   151,   153,   155,   157,   159,
     161,   163,   165,   167,   169,   171,   173,   175,   177,   179,
     181,   183,   185,   190,   192,   194,   196,   198,   200,   204,
     207,   209,   210,   211,   214,   221,   222,   223,   226,   230,
     231,   232,   235,   239,   240,   241,   244,   248,   249,   250,
     253,   257,   258,   259,   262,   270,   271,   272,   275,   276,
     278,   286,   294,   295,   297,   300,   304,   305,   306,   309,
     312,   313,   314,   319,   324,   329,   333,   337,   341,   344,
     347,   350,   353,   357,   361,   366,   367,   368,   371,   374,
     379,   383
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "JACOBIAN", "DOUBLE", "FUNCTION",
  "DEFVAR", "DEFRAD", "DEFFIX", "SETVAR", "SETRAD", "SETFIX", "HESSIAN",
  "STOICMAT", "STOCHASTIC", "DECLARE", "INITVALUES", "EQUATIONS", "LUMP",
  "INIEQUAL", "EQNEQUAL", "EQNCOLON", "LMPCOLON", "LMPPLUS", "SPCPLUS",
  "SPCEQUAL", "ATOMDECL", "CHECK", "CHECKALL", "REORDER", "MEX",
  "DUMMYINDEX", "EQNTAGS", "LOOKAT", "LOOKATALL", "TRANSPORT",
  "TRANSPORTALL", "MONITOR", "USES", "SPARSEDATA", "WRITE_ATM",
  "WRITE_SPC", "WRITE_MAT", "WRITE_OPT", "INITIALIZE", "XGRID", "YGRID",
  "ZGRID", "USE", "LANGUAGE", "INTFILE", "DRIVER", "RUN", "INLINE",
  "ENDINLINE", "PARAMETER", "SPCSPC", "INISPC", "INIVALUE", "EQNSPC",
  "EQNSIGN", "EQNCOEF", "RATE", "LMPSPC", "SPCNR", "ATOMID", "LKTID",
  "MNIID", "INLCTX", "INCODE", "SSPID", "EQNLESS", "EQNTAG", "EQNGREATER",
  "TPTID", "USEID", "';'", "$accept", "program", "section", "semicolon",
  "atomlist", "atomdef", "lookatlist", "lookatspc", "monitorlist",
  "monitorspc", "translist", "transspc", "uselist", "usefile",
  "setspclist", "setspcspc", "species", "spc", "spcname", "spcdef",
  "atoms", "atom", "initvalues", "assignment", "equations", "equation",
  "rate", "eqntag", "lefths", "righths", "expresion", "term", "lumps",
  "lump", "inlinecode", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,    59
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    77,    78,    78,    79,    79,    79,    79,    79,    79,
      79,    79,    79,    79,    79,    79,    79,    79,    79,    79,
      79,    79,    79,    79,    79,    79,    79,    79,    79,    79,
      79,    79,    79,    79,    79,    79,    79,    79,    79,    79,
      79,    79,    79,    79,    79,    79,    79,    79,    79,    80,
      80,    81,    81,    81,    82,    83,    83,    83,    84,    85,
      85,    85,    86,    87,    87,    87,    88,    89,    89,    89,
      90,    91,    91,    91,    92,    93,    93,    93,    94,    94,
      95,    96,    97,    97,    98,    98,    99,    99,    99,   100,
     101,   101,   101,   102,   102,   103,   103,   104,   105,   106,
     107,   107,   107,   108,   108,   109,   109,   109,   110,   110,
     111,   111
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     1,
       1,     1,     1,     1,     1,     1,     2,     2,     2,     2,
       2,     2,     4,     2,     2,     2,     2,     2,     2,     2,
       1,     3,     2,     2,     1,     3,     2,     2,     1,     3,
       2,     2,     1,     3,     2,     2,     1,     3,     2,     2,
       1,     3,     2,     2,     1,     3,     2,     2,     1,     1,
       3,     1,     3,     1,     2,     1,     3,     2,     2,     3,
       3,     2,     2,     4,     3,     2,     1,     3,     2,     2,
       3,     2,     1,     2,     1,     3,     2,     2,     3,     3,
       2,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    29,
       0,     0,     0,     0,     0,    30,     0,    31,     0,     0,
       0,    32,    33,    34,    35,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     2,     4,     8,    13,
       0,    81,    17,     0,    78,    79,    18,    19,     0,    74,
      20,     0,    21,    22,     5,     7,    14,     6,     0,     0,
      23,     0,     0,   104,     0,     0,     0,    24,     0,     0,
       0,     0,   102,     0,     0,    25,     0,     0,    54,    15,
       0,    16,     9,    10,    11,    12,     0,    58,    26,     0,
       0,    66,    28,     0,     0,    62,    27,     0,     0,    70,
      47,     0,    48,    38,    39,    40,    41,    36,    37,    44,
      45,    46,    43,     0,     1,     3,    50,    77,     0,     0,
      76,    73,     0,    72,    88,     0,     0,    87,    92,   101,
     103,     0,     0,    91,     0,     0,     0,    98,     0,   107,
       0,     0,     0,   106,    53,     0,    52,    57,     0,    56,
      65,     0,    64,    61,     0,    60,    69,     0,    68,   111,
       0,    49,    85,     0,    80,    83,    75,    71,    89,    86,
      97,    90,     0,    96,    94,    99,   100,   109,   108,   105,
      51,    55,    63,    59,    67,    42,   110,    84,     0,    93,
      95,    82
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,    45,    46,   127,    89,    90,    98,    99,   106,   107,
     102,   103,   110,   111,    60,    61,    52,    53,    54,    55,
     174,   175,    70,    71,    77,    78,   184,    79,    80,   145,
      81,    82,    85,    86,   170
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -91
static const yytype_int16 yypact[] =
{
     112,    -6,    -3,     4,    11,    11,    11,     3,     3,     3,
      30,    33,    39,    40,     8,     1,    21,     9,     9,   -91,
      41,    42,    43,    45,     5,   -91,     6,   -91,    15,     0,
      46,   -91,   -91,   -91,   -91,    48,    49,    51,    52,    54,
      55,    57,    58,    76,    10,    63,   112,   -91,   -91,   -91,
      56,   108,    14,    56,   -91,   -91,    14,    14,    56,   -91,
      64,    56,    64,    64,   -91,   -91,   -91,   -91,    56,   116,
      79,    56,    56,   -91,   -38,    78,    94,   -33,    56,   -20,
     -20,    -7,   -91,    56,    34,   104,    56,    56,   -91,   103,
      56,   103,   -91,   -91,   -91,   -91,    56,   -91,   105,    56,
      56,   -91,    95,    56,    56,   -91,   106,    56,    56,   -91,
      97,    56,   -91,   -91,   -91,   -91,   -91,   -91,   -91,   -91,
     -91,   -91,   -91,   101,   -91,   -91,   -91,    98,   -32,    56,
      98,    98,    56,    98,    98,   117,    56,    98,    98,   -91,
     -91,   107,    56,    98,   -20,   114,    27,   -91,   -38,    98,
     115,   104,    56,    98,    98,    56,    98,    98,    56,    98,
      98,    56,    98,    98,    56,    98,    98,    56,    98,   -91,
     -40,   -91,   -91,   121,   155,   -91,    98,    98,   -91,    98,
     -91,    98,   114,   114,   -91,   -91,   -91,   -91,   -91,    98,
      98,    98,    98,    98,    98,   -91,   -91,   -91,   -32,   -91,
     -91,   -91
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -91,   135,   -91,   -53,   164,   -25,   -91,    86,   -91,    77,
     -91,    83,   -91,    80,    60,   -18,    85,   -21,   -91,   -91,
     -91,   -12,   -91,   118,   -91,   110,   -90,   -91,   113,    47,
     -63,   -71,   -91,   -65,   -91
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
     130,   108,    72,   139,    58,   131,    96,   100,   133,    68,
      87,   122,    50,   147,   195,   134,   104,   146,   137,   138,
     152,    73,    83,    75,   172,   143,    73,    74,    75,   196,
     149,   129,   173,   153,   154,   129,   129,   156,    76,    73,
      74,    75,   132,   157,   132,   132,   159,   160,   185,    47,
     162,   163,    48,   148,   165,   166,   150,   151,   168,    49,
      73,    74,    75,   124,   155,    69,   155,    51,    62,    63,
      51,    97,    76,    59,    88,   109,   176,   186,   123,   177,
     101,   146,   105,   179,    84,    64,   188,   148,    65,   181,
      56,    57,   199,   200,    66,    67,    92,    93,    94,   189,
      95,   112,   190,   113,   114,   191,   115,   116,   192,   117,
     118,   193,   119,   120,   194,     1,     2,     3,     4,     5,
       6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,   121,   126,   128,    59,   135,    69,   140,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,   141,    84,    88,   101,
     169,    97,   109,   105,   171,   178,   183,   197,   187,   198,
     180,   125,    91,   164,   158,   161,   201,   142,   136,     0,
     167,   182,   144
};

static const yytype_int16 yycheck[] =
{
      53,     1,     1,    74,     1,    58,     1,     1,    61,     1,
       1,     1,     1,    20,    54,    68,     1,    80,    71,    72,
      85,    59,     1,    61,    56,    78,    59,    60,    61,    69,
      83,    52,    64,    86,    87,    56,    57,    90,    71,    59,
      60,    61,    60,    96,    62,    63,    99,   100,    21,    55,
     103,   104,    55,    60,   107,   108,    22,    23,   111,    55,
      59,    60,    61,     0,    89,    57,    91,    56,     8,     9,
      56,    66,    71,    70,    65,    75,   129,   148,    68,   132,
      74,   144,    67,   136,    63,    55,   151,    60,    55,   142,
       5,     6,   182,   183,    55,    55,    55,    55,    55,   152,
      55,    55,   155,    55,    55,   158,    55,    55,   161,    55,
      55,   164,    55,    55,   167,     3,     4,     5,     6,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    55,    76,    25,    70,    19,    57,    59,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    50,    51,    52,    53,    72,    63,    65,    74,
      69,    66,    75,    67,    76,    58,    62,    56,    63,    24,
      73,    46,    18,   106,    98,   102,   198,    77,    70,    -1,
     110,   144,    79
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,    45,    46,    47,    48,
      49,    50,    51,    52,    53,    78,    79,    55,    55,    55,
       1,    56,    93,    94,    95,    96,    93,    93,     1,    70,
      91,    92,    91,    91,    55,    55,    55,    55,     1,    57,
      99,   100,     1,    59,    60,    61,    71,   101,   102,   104,
     105,   107,   108,     1,    63,   109,   110,     1,    65,    81,
      82,    81,    55,    55,    55,    55,     1,    66,    83,    84,
       1,    74,    87,    88,     1,    67,    85,    86,     1,    75,
      89,    90,    55,    55,    55,    55,    55,    55,    55,    55,
      55,    55,     1,    68,     0,    78,    76,    80,    25,    94,
      80,    80,    92,    80,    80,    19,   100,    80,    80,   108,
      59,    72,   102,    80,   105,   106,   107,    20,    60,    80,
      22,    23,   110,    80,    80,    82,    80,    80,    84,    80,
      80,    88,    80,    80,    86,    80,    80,    90,    80,    69,
     111,    76,    56,    64,    97,    98,    80,    80,    58,    80,
      73,    80,   106,    62,   103,    21,   108,    63,   110,    80,
      80,    80,    80,    80,    80,    54,    69,    56,    24,   103,
     103,    98
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 4:
#line 99 "scan.y"
    { CmdJacobian( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 5:
#line 102 "scan.y"
    { CmdHessian( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 6:
#line 105 "scan.y"
    { CmdDeclareValues( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 7:
#line 108 "scan.y"
    { CmdStoicmat( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 8:
#line 111 "scan.y"
    { CmdDouble( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 9:
#line 114 "scan.y"
    { CmdReorder( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 10:
#line 117 "scan.y"
    { CmdMex( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 11:
#line 120 "scan.y"
    { CmdDummyindex( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 12:
#line 123 "scan.y"
    { CmdEqntags( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 13:
#line 126 "scan.y"
    { CmdFunction( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 14:
#line 129 "scan.y"
    { CmdStochastic( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 15:
#line 132 "scan.y"
    {;}
    break;

  case 16:
#line 134 "scan.y"
    {;}
    break;

  case 17:
#line 136 "scan.y"
    {;}
    break;

  case 18:
#line 138 "scan.y"
    {;}
    break;

  case 19:
#line 140 "scan.y"
    {;}
    break;

  case 20:
#line 142 "scan.y"
    {;}
    break;

  case 21:
#line 144 "scan.y"
    {;}
    break;

  case 22:
#line 146 "scan.y"
    {;}
    break;

  case 23:
#line 148 "scan.y"
    {;}
    break;

  case 24:
#line 150 "scan.y"
    {;}
    break;

  case 25:
#line 152 "scan.y"
    {;}
    break;

  case 26:
#line 154 "scan.y"
    {;}
    break;

  case 27:
#line 156 "scan.y"
    {;}
    break;

  case 28:
#line 158 "scan.y"
    {;}
    break;

  case 29:
#line 160 "scan.y"
    { CheckAll(); ;}
    break;

  case 30:
#line 162 "scan.y"
    { LookAtAll(); ;}
    break;

  case 31:
#line 164 "scan.y"
    { TransportAll(); ;}
    break;

  case 32:
#line 166 "scan.y"
    { WriteAtoms(); ;}
    break;

  case 33:
#line 168 "scan.y"
    { WriteSpecies(); ;}
    break;

  case 34:
#line 170 "scan.y"
    { WriteMatrices(); ;}
    break;

  case 35:
#line 172 "scan.y"
    { WriteOptions(); ;}
    break;

  case 36:
#line 174 "scan.y"
    { CmdUse( (yyvsp[(2) - (2)].str) ); ;}
    break;

  case 37:
#line 176 "scan.y"
    { CmdLanguage( (yyvsp[(2) - (2)].str) ); ;}
    break;

  case 38:
#line 178 "scan.y"
    { DefineInitializeNbr( (yyvsp[(2) - (2)].str) ); ;}
    break;

  case 39:
#line 180 "scan.y"
    { DefineXGrid( (yyvsp[(2) - (2)].str) ); ;}
    break;

  case 40:
#line 182 "scan.y"
    { DefineYGrid( (yyvsp[(2) - (2)].str) ); ;}
    break;

  case 41:
#line 184 "scan.y"
    { DefineZGrid( (yyvsp[(2) - (2)].str) ); ;}
    break;

  case 42:
#line 186 "scan.y"
    { 
		    AddInlineCode( (yyvsp[(2) - (4)].str), InlineBuf );
                    free( InlineBuf );
		  ;}
    break;

  case 43:
#line 191 "scan.y"
    { ParserErrorMessage(); ;}
    break;

  case 44:
#line 193 "scan.y"
    { CmdIntegrator( (yyvsp[(2) - (2)].str) ); ;}
    break;

  case 45:
#line 195 "scan.y"
    { CmdDriver( (yyvsp[(2) - (2)].str) ); ;}
    break;

  case 46:
#line 197 "scan.y"
    { CmdRun( (yyvsp[(2) - (2)].str) ); ;}
    break;

  case 47:
#line 199 "scan.y"
    {;}
    break;

  case 48:
#line 201 "scan.y"
    { SparseData( (yyvsp[(2) - (2)].str) );
                  ;}
    break;

  case 49:
#line 205 "scan.y"
    { ScanWarning("Unnecessary ';'");
                  ;}
    break;

  case 53:
#line 212 "scan.y"
    { ParserErrorMessage(); ;}
    break;

  case 54:
#line 215 "scan.y"
    { switch( crt_section ) {
                      case ATOMDECL: DeclareAtom( (yyvsp[(1) - (1)].str) ); break;
                      case CHECK:    SetAtomType( (yyvsp[(1) - (1)].str), DO_CHECK ); break;
                    }
                  ;}
    break;

  case 57:
#line 224 "scan.y"
    { ParserErrorMessage(); ;}
    break;

  case 58:
#line 227 "scan.y"
    { AddLookAt( (yyvsp[(1) - (1)].str) );
                  ;}
    break;

  case 61:
#line 233 "scan.y"
    { ParserErrorMessage(); ;}
    break;

  case 62:
#line 236 "scan.y"
    { AddMonitor( (yyvsp[(1) - (1)].str) );
                  ;}
    break;

  case 65:
#line 242 "scan.y"
    { ParserErrorMessage(); ;}
    break;

  case 66:
#line 245 "scan.y"
    { AddTransport( (yyvsp[(1) - (1)].str) );
                  ;}
    break;

  case 69:
#line 251 "scan.y"
    { ParserErrorMessage(); ;}
    break;

  case 70:
#line 254 "scan.y"
    { AddUseFile( (yyvsp[(1) - (1)].str) );
                  ;}
    break;

  case 73:
#line 260 "scan.y"
    { ParserErrorMessage(); ;}
    break;

  case 74:
#line 263 "scan.y"
    { switch( crt_section ) {
                      case SETVAR: SetSpcType( VAR_SPC, (yyvsp[(1) - (1)].str) ); break;
                      case SETRAD: SetSpcType( RAD_SPC, (yyvsp[(1) - (1)].str) ); break;
                      case SETFIX: SetSpcType( FIX_SPC, (yyvsp[(1) - (1)].str) ); break;
                    }
                  ;}
    break;

  case 77:
#line 273 "scan.y"
    { ParserErrorMessage(); ;}
    break;

  case 80:
#line 279 "scan.y"
    { switch( crt_section ) {
                      case DEFVAR: DeclareSpecies( VAR_SPC, (yyvsp[(1) - (3)].str) ); break;
                      case DEFRAD: DeclareSpecies( RAD_SPC, (yyvsp[(1) - (3)].str) ); break;
                      case DEFFIX: DeclareSpecies( FIX_SPC, (yyvsp[(1) - (3)].str) ); break;
                    } 
                  ;}
    break;

  case 81:
#line 287 "scan.y"
    { switch( crt_section ) {
                      case DEFVAR: DeclareSpecies( VAR_SPC, (yyvsp[(1) - (1)].str) ); break;
                      case DEFRAD: DeclareSpecies( RAD_SPC, (yyvsp[(1) - (1)].str) ); break;
                      case DEFFIX: DeclareSpecies( FIX_SPC, (yyvsp[(1) - (1)].str) ); break;
                    } 
                  ;}
    break;

  case 84:
#line 298 "scan.y"
    { AddAtom( (yyvsp[(2) - (2)].str), (yyvsp[(1) - (2)].str) );
                  ;}
    break;

  case 85:
#line 301 "scan.y"
    { AddAtom( (yyvsp[(1) - (1)].str), "1" );
                  ;}
    break;

  case 88:
#line 307 "scan.y"
    { ParserErrorMessage(); ;}
    break;

  case 89:
#line 310 "scan.y"
    { AssignInitialValue( (yyvsp[(1) - (3)].str), (yyvsp[(3) - (3)].str) ); ;}
    break;

  case 92:
#line 315 "scan.y"
    { ParserErrorMessage();
                    eqState = LHS; 
                  ;}
    break;

  case 93:
#line 320 "scan.y"
    { eqState = LHS;
                    StoreEquationRate( (yyvsp[(4) - (4)].str), (yyvsp[(1) - (4)].str) ); 
                    CheckEquation();
                  ;}
    break;

  case 94:
#line 325 "scan.y"
    { eqState = LHS;
                    StoreEquationRate( (yyvsp[(3) - (3)].str), "          " ); 
                    CheckEquation();
                  ;}
    break;

  case 95:
#line 330 "scan.y"
    { strcpy( (yyval.str), (yyvsp[(1) - (2)].str) );
                    strcat( (yyval.str), (yyvsp[(2) - (2)].str) ); 
                  ;}
    break;

  case 96:
#line 334 "scan.y"
    { strcpy( (yyval.str), (yyvsp[(1) - (1)].str) );
                  ;}
    break;

  case 97:
#line 338 "scan.y"
    { strcpy( (yyval.str), (yyvsp[(2) - (3)].str) );
                  ;}
    break;

  case 98:
#line 342 "scan.y"
    { eqState = RHS; ;}
    break;

  case 99:
#line 345 "scan.y"
    { eqState = RAT; ;}
    break;

  case 100:
#line 348 "scan.y"
    { ProcessTerm( eqState, (yyvsp[(2) - (3)].str), crt_coef, crt_term ); 
                  ;}
    break;

  case 101:
#line 351 "scan.y"
    { ProcessTerm( eqState, (yyvsp[(1) - (2)].str), crt_coef, crt_term );
                  ;}
    break;

  case 102:
#line 354 "scan.y"
    { ProcessTerm( eqState, "+", crt_coef, crt_term );
                  ;}
    break;

  case 103:
#line 358 "scan.y"
    { strcpy( crt_term, (yyvsp[(2) - (2)].str) );
                    strcpy( crt_coef, (yyvsp[(1) - (2)].str) );  
                  ;}
    break;

  case 104:
#line 362 "scan.y"
    { strcpy( crt_term, (yyvsp[(1) - (1)].str) );         
                    strcpy( crt_coef, "1" ); 
                  ;}
    break;

  case 107:
#line 369 "scan.y"
    { ParserErrorMessage(); ;}
    break;

  case 108:
#line 372 "scan.y"
    { AddLumpSpecies( (yyvsp[(1) - (3)].str) );
                  ;}
    break;

  case 109:
#line 375 "scan.y"
    {
                    AddLumpSpecies( (yyvsp[(1) - (3)].str) );
                    CheckLump( (yyvsp[(3) - (3)].str) );  
                  ;}
    break;

  case 110:
#line 380 "scan.y"
    {
		    InlineBuf = AppendString( InlineBuf, (yyvsp[(2) - (2)].str), &InlineLen, MAX_INLINE );
		  ;}
    break;

  case 111:
#line 384 "scan.y"
    {
		    InlineBuf = malloc( MAX_INLINE ); 
                    InlineLen = MAX_INLINE;
		    strcpy( InlineBuf, (yyvsp[(1) - (1)].str));
		  ;}
    break;


/* Line 1267 of yacc.c.  */
#line 2174 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 390 "scan.y"


void yyerror( char * str )
{
}

void ParserErrorMessage()
{
  /* yyerrok; */
/*
  Message("[%d,%s] -> [%d,%s]", crtTokType, crtToken, nextTokType, nextToken );  
*/
  if( crtToken[0] == ';' ) {
    ParserError("Misplaced ';'");
    return;
  }
  switch( crtTokType ) {
    case ATOMID:
      ParserError("Missing ';' after '%s'", crtToken );
      break; 

    case SPCSPC: 
      ParserError("Missing ';' or '+' after '%s'", crtToken );
      break; 
    case SPCNR:
      ParserError("Missing species after '%s'", crtToken );
      break; 
    case SPCPLUS:
      ParserError("Missing atom after '%s'", crtToken );
      break; 
    case SPCEQUAL:
      ParserError("Invalid '=' after '%s'", crtToken );
      break; 

    case INISPC: 
      ParserError("Missing '=' after '%s'", crtToken );
      break; 
    case INIEQUAL: 
      ParserError("Missing value after '%s'", crtToken );
      break; 
    case INIVALUE: 
      ParserError("Missing ';' after '%s'", crtToken );
      break; 

    case EQNSPC: 
      ParserError("Missing '+' or '=' after '%s'", crtToken );
      break; 
    case EQNEQUAL: 
      ParserError("Invalid right hand side of equation");
      break; 
    case EQNCOLON: 
      ParserError("Missing rate after '%s'", crtToken );
      break; 
    case EQNSIGN: 
      ParserError("Missing coeficient after '%s'", crtToken );
      break; 
    case EQNCOEF: 
      ParserError("Missing species after '%s'", crtToken );
      break; 
    case RATE: 
      ParserError("Missing ';' after '%s'", crtToken );
      break; 

    case LMPSPC: 
      ParserError("Missing '+' or ':' or ';' after '%s'", crtToken );
      break; 
    case LMPPLUS: 
      ParserError("Missing species after '%s'", crtToken );
      break; 
    case LMPCOLON: 
      ParserError("Missing species after '%s'", crtToken );
      break; 
    case INLINE:
      ParserError("Missing inline option after '%s'", crtToken );
      break;

    default:
      ParserError("Syntax error after '%s'", crtToken ); 
  }
}


int Parser( char * filename )
{
extern int yydebug;
FILE *f;

  crt_filename = filename;

  f = fopen( crt_filename, "r" );
  if( f == 0 ) {
    FatalError(7,"%s: File not found", crt_filename);
  } 
  
  yyin = f;
  nError   = 0;
  nWarning = 0;
  yydebug = 0;

  yyparse();

  fclose( f );

  return nError;
}          


