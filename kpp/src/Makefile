##############################################################################
# -*- sh -*-
##############################################################################
#
#  KPP - The Kinetic PreProcessor
#        Builds simulation code for chemical kinetic systems
#
#  Copyright (C) 1995-1997 Valeriu Damian and Adrian Sandu
#  Copyright (C) 1997-2005 Adrian Sandu
#
#  KPP is free software; you can redistribute it and/or modify it under the
#  terms of the GNU General Public License as published by the Free Software
#  Foundation (http://www.gnu.org/copyleft/gpl.html); either version 2 of the
#  License, or (at your option) any later version.
#
#  KPP is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
#  details.
#
#  You should have received a copy of the GNU General Public License along
##  with this program; if not, consult http://www.gnu.org/copyleft/gpl.html or
#  write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
#  Boston, MA  02111-1307,  USA.
#
#  Adrian Sandu
#  Computer Science Department
#  Virginia Polytechnic Institute and State University
#  Blacksburg, VA 24060
#  E-mail: sandu@cs.vt.edu
#
##############################################################################

# In order to compile KPP you have to provide the following information:

# 1. The name of the compiler you want to use. Normaly this 
#    is either GNU C compiler (gcc) or the native compiler (cc)
#    You can use the complete pathname if the compiler is not in $PATH 
#    Note that for SUN machines is better to use gcc.
#    For GNU C compiler use:
#      CC=gcc
#    For the native compiler use:
#      CC=cc

CC=gcc

# 2. Platform independent C compiler flags. By default "-O" is used which 
#    turns on optimisation. If you are experiencing problems you may try 
#    "-g" to include debuging informations.

CFLAGS=-g -O2
CPPFLAGS=

# 3. The name of your lexical analizer. KPP requires FLEX to be used.
#    FLEX is a public domain lexical analizer and you can download it from
#    http://www.gnu.org/software/flex/ or any other mirror site. If flex
#    directory is not included in your path use the complete pathname.

FLEX=flex

# 4. Linker flags to ensure access to the FLEX library (libfl.a).
#    If libfl.a is in a standard location, this can be blank.
#    Otherwise LDFLAGS=-L<DIRPATH> where DIRPATH is the path that
#    contains libfl.a
LDFLAGS=

# 4. Flex library. Can either be -lfl or the full path to a libfl.a
FLEX_LIB=-lfl

# 5. The name of the fortran compiler you want to use. Normaly this 
#    is either GNU Fortran compiler (gfortran), the Intel Fortran
#    compiler (ifort), or the Portland Group Compiler.
#    For GNU C compiler use:
#      FC=gfortran
#    For the Intel compiler use:
#      FC=ifort

FC=gfortran

# 6. Platform independent C compiler flags. By default "-O" is used which 
#    turns on optimisation. If you are experiencing problems you may try 
#    "-g" to include debuging informations.

F90FLAGS=-cpp -g -O2 -fno-automatic -fcheck=bounds -fimplicit-none


BISON=bison -d 

all: kpp

.c.o:
	@echo "  "$(CC) $(CC_FLAGS) $(CFLAGS) -c $*.c
	@$(CC) $(CC_FLAGS) $(CFLAGS) -c $*.c

OBJS  = \
	y.tab.o \
	lex.yy.o \
	scanner.o \
	scanutil.o \
	kpp.o \
	gen.o \
	code.o \
	code_c.o \
	code_f77.o \
	code_f90.o \
	code_matlab.o \
	debug.o

kpp:    $(OBJS)
	@echo "  "$(CC) $(CC_FLAGS) $(CFLAGS) $(OBJS) $(LDFLAGS) $(FLEX_LIB) -o kpp
	@$(CC) $(CC_FLAGS) $(CFLAGS) $(OBJS) $(LDFLAGS) $(FLEX_LIB) -o kpp
	@mkdir -p ../bin
	@mv kpp ../bin/

clean:  
	@rm -f *~ *.o

maintainer-clean: clean
	@rm -f lex.yy.c y.tab.c y.tab.h

lex.yy.c: scan.l scan.h
	@echo "  "$(FLEX) scan.l
	@$(FLEX) -olex.yy.c scan.l 

y.tab.c: scan.y scan.h
	@echo "  "$(BISON) scan.y
	@$(BISON) -o y.tab.c scan.y

flex: lex.yy.c y.tab.c

scanner.o: scan.h gdata.h
scanutil.o: scan.h
kpp.o: gdata.h
gen.o: gdata.h code.h
debug.o: gdata.h
code.o: gdata.h code.h

code_c.o: gdata.h code.h
code_f.o: gdata.h code.h
