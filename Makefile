##############################################################################
# -*- sh -*-
##############################################################################
#
#  DSMACC - The Dynamically Simple Model of Atmospheric Chemical Complexity
#        uses KPP to build simulation code for chemical kinetic systems
#
##############################################################################

# Basic compilers and compiler options are selected using configure

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

# 3. KPP requires a lexical analizer like FLEX to be used.
#    FLEX is a public domain lexical analizer and you can download it from
#    http://www.gnu.org/software/flex/ or any other mirror site. If flex
#    directory is not included in your path use the complete pathname.

FLEX=flex

# 4. Flex library. Can either be -lfl or the full path to a libfl.a
FLEX_LIB=-lfl


# 5. Linker flags to ensure access to the libraries including libfl.a
#    If libfl.a is in a standard location, this can be blank.
#    Otherwise LDFLAGS=-L<DIRPATH> where DIRPATH is the path that
#    contains libfl.a
LDFLAGS=

# 6. The name of the fortran compiler you want to use. Normaly this 
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

TOPDIR=$(realpath .)

export KPP_HOME=$(PWD)/kpp

all: testrun

testrun: kpp/bin/kpp
	cd cri && make

kpp/bin/kpp:
	make -C kpp

clean:
	cd src && make -i distclean
	cd cri && make -i clean
	cd tuv_new && make -i clean
	cd UCI_fastJX72e && make -i clean
	cd kpp && make -i clean

distclean: clean
	rm -f bin/dsmacc
	cd kpp && make -i distclean
	rm -rf autom4te.cache config.status config.log tuvlog.txt
