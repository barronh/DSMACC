# Makefile for TUV 4.2
# Use with disort (discrete ordinate), or ps2str (2 stream approximation, 
# pseudo-spherical correction)
#----------
# EXC      : name of executable
# INCLUDES : required include files
# USE_INCL : object files referencing include file (params)
# FOBJS    : all required object files that do not use the include file
#

EXC = tuv

INCLUDES = params 

USE_INCL = TUV.o \
           grids.o \
           rdinp.o rdetfl.o rdxs.o \
           swphys.o swbiol.o swchem.o rxn.o qys.o \
           wshift.o \
	   vpair.o vptmp.o vpo3.o \
	   odrl.o odo3.o \
           setaer.o setalb.o setcld.o setsnw.o \
           setno2.o seto2.o setso2.o \
           sphers.o  \
	   la_srb.o \
           rtrans.o \
	   savout.o

FOBJS = numer.o functs.o orbit.o

#----------
# FC   : FORTRAN compiler 
#        Linux users:  try FC = g77
#        Cray users :  try FC = f90
FC = gfortran

# FFLAGS : command line options to compiler call (if not set, default is
#          probably some basic optimization level)
# FFLAGS = 
#FFLAGS =  -x f77-cpp-input -DACETONEQY=BLITZ2004QY -fbounds-check -fno-automatic -ffixed-line-length-72
#FFLAGS   = -m32 -x f77-cpp-input -DACETONEQY=BLITZ2004QY -I/usr/local/include -fcheck=all -fno-automatic 
TUVF90FLAGS   = -m32 -I/Users/barronh/include -fbounds-check -O3
# LIBS  : libraries required
# LIBS =  

#----------

$(EXC):		$(FOBJS) $(USE_INCL) main.o
		$(FC) -m32 -DACETONEQY=BLITZ2004QY -I/usr/local/include $(FOBJS) main.o  $(USE_INCL) $(LIBS) -o $@

lib:	$(FOBJS) $(USE_INCL)

$(USE_INCL) $(LIBS) -o $@

$(USE_INCL):	$(INCLUDES)

.f.o:
	$(FC) $(FFLAGS) -c $*.f

swchem.o: swchem.f
	$(FC) $(FFLAGS) -DVERBOSE=1 -c $< 

clean:		
		rm -f core $(EXC) $(USE_INCL) $(FOBJS) main.o