/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

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




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 72 "scan.y"
{
  char str[200];
}
/* Line 1529 of yacc.c.  */
#line 203 "y.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

