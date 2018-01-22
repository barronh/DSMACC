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

CC=@CC@

# 2. Platform independent C compiler flags. By default "-O" is used which 
#    turns on optimisation. If you are experiencing problems you may try 
#    "-g" to include debuging informations.

CFLAGS=@CFLAGS@
CPPFLAGS=@CPPFLAGS@

# 3. The name of your lexical analizer. KPP requires FLEX to be used.
#    FLEX is a public domain lexical analizer and you can download it from
#    http://www.gnu.org/software/flex/ or any other mirror site. If flex
#    directory is not included in your path use the complete pathname.

FLEX=@LEX@

# 4. Linker flags to ensure access to the FLEX library (libfl.a).
#    If libfl.a is in a standard location, this can be blank.
#    Otherwise LDFLAGS=-L<DIRPATH> where DIRPATH is the path that
#    contains libfl.a
LDFLAGS=@LDFLAGS@

# 4. Flex library. Can either be -lfl or the full path to a libfl.a
FLEX_LIB=@LEXLIB@

# 5. The name of the fortran compiler you want to use. Normaly this 
#    is either GNU Fortran compiler (gfortran), the Intel Fortran
#    compiler (ifort), or the Portland Group Compiler.
#    For GNU C compiler use:
#      FC=gfortran
#    For the Intel compiler use:
#      FC=ifort

FC=@FC@

# 6. Platform independent C compiler flags. By default "-O" is used which 
#    turns on optimisation. If you are experiencing problems you may try 
#    "-g" to include debuging informations.

F90FLAGS=@FCFLAGS_f90@ @FCFLAGS@

export KPP_HOME=$(PWD)/../kpp

##############################################################################
all:
	make source
#	make pykpp
#	make test

# link TUV inputs for the run
inputs:
	ln -fs ../tuv_new/INPUTS ../tuv_new/DATA* .

# build the standard cri model including f90 and objects
source: cri_mcmj.kpp ../src/global.inc ../src/mcmratedef.inc ../driver.f90 ../src/tuv.inc ../src/io.inc
	cd ../src && make distclean && $(KPP_HOME)/bin/kpp ../pysrc/cri_mcmj.kpp dsmacc && make libtuv libfastjx #&& rm ../bin/dsmacc


# Run the model and make a plot
test:
	python driver.py
	python plot.py

# complete list of all f90 source files
SRCS1 = $(wildcard ../src/dsmacc_*.f90)

# the object files are the same as the source files but with suffix ".o"
OBJS1 := $(SRCS1:.f90=.o) 
OBJS := $(OBJS1)
LIBS := ../tuv_new/libtuv.a ../UCI_fastJX72e/libfastjx.a

pykpp: kpp.pyf
	$(F2PY) --f90exec=gfortran -m kpp kpp.pyf pyint.f90 pyglob.f90 pyrate.f90 pymon.f90 ../src/dsmacc_Parameters.f90 ../src/dsmacc_Precision.f90
	$(F2PY) --f90exec=gfortran -c -m kpp kpp.pyf pyint.f90 pyglob.f90 pyrate.f90 pymon.f90 -I../src/ ../src/dsmacc_Parameters.f90  ../src/dsmacc_Precision.f90 $(OBJS) $(LIBS)

kpp.pyf: ../src/dsmacc_Precision.f90
	python makekpp_pyf.py

clean:
	rm -f kpp.pyf photolysis.txt *.so *.tsv *.png INPUTS DATA* *.dat