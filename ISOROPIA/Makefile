#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: Makefile (in the ISOROPIA/ subdirectory)
#
# !DESCRIPTION: This makefile compiles the ISOROPIA code.  
#  Object files (*.o) are bundled into the libIsoropia.a library 
#  (located in the LIB directory).  Module files (*.mod) are copied to 
#  the MOD directory. 
#\\
#\\
# !REMARKS:
# To build the programs, call "make" with the following syntax:
#                                                                             .
#   make -jN TARGET REQUIRED-FLAGS [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% NOTE: Normally you will not have to call this Makefile directly,     %%%
# %%% it will be called automatically from the Makefile in the directory   %%%
# %%% just above this one!                                                 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                                                             .
# Makefile uses the following variables:
#                                                                             .
# Variable   Description
# --------   -----------
# SHELL      Specifies the shell for "make" to use (usually SHELL=/bin/sh)
# ROOTDIR    Specifies the root-level directory of the GEOS-Chem code
# HDR        Specifies the directory where GEOS-Chem include files are found
# LIB        Specifies the directory where library files (*.a) are stored
# MOD        Specifies the directory where module files (*.mod) are stored
# AR         Sys var w/ name of library creator program (i.e., "ar", "ranlib")
# MAKE       Sys var w/ name of Make command (i.e, "make" or "gmake")
# R8         Specifies the command to treat "REAL" as "REAL*8"
#
# !REVISION HISTORY: 
#  21 Dec 2009 - C. Carouge  - Initial version
#  22 Aug 2011 - R. Yantosca - Add "-fp-model source" flag for IFORT compiler,
#                              which prevents random numerical noise
#  25 Aug 2011 - R. Yantosca - Remove -fp-model source flag here, as this is
#                              now added to FFLAGS in Makefile_header.mk
#  19 Mar 2014 - R. Yantosca - Add more visible comment section dividers
#  21 Nov 2014 - R. Yantosca - Compile w/o includes for ESMF, MAPL, FVdycore
#  04 Jun 2015 - R. Yantosca - Also remove *.mod, *.a files with "make clean"
#  04 Jun 2015 - R. Yantosca - Add debug target, remove help
#EOP
#------------------------------------------------------------------------------
#BOC

###############################################################################
###                                                                         ###
###  Initialization section                                                 ###
###                                                                         ###
###############################################################################

PRECISION=8
R8                 := -fdefault-real-8 -fdefault-double-8
FFLAGS             :=-cpp -w -std=legacy -fautomatic -fno-align-commons
FFLAGS             += -fconvert=big-endian -fPIC
FFLAGS             += -fno-range-check -DUSE_REAL8
F90ISO=gfortran $(FFLAGS)
# Directories
ROOT    :=..
LIB     :=.
MOD     :=.

# List of source files
SOURCES :=$(wildcard *.f) #$(wildcard *.F90)

# List of object files (replace .F and .F90 extensions with *.o)
TMP     :=$(SOURCES:.f=.o)
OBJECTS :=$(TMP:.F90=.o)

# List of module files.  Convert to lowercase, then prefix directory name.
MODULES :=$(OBJECTS:.o=.mod)
MODULES :=$(shell echo $(MODULES) | tr A-Z a-z)
MODULES :=$(foreach I,$(MODULES),$(MOD)/$(I))

# Library file
LIBRARY :=libIsoropia.a

###############################################################################
###                                                                         ###
###  Makefile targets: type "make help" for a complete listing!             ###
###                                                                         ###
###############################################################################

.PHONY: clean debug

lib: $(OBJECTS)
	$(AR) crs $(LIBRARY) $(OBJECTS)
#	mv $(LIBRARY) $(LIB)

clean:
	@echo "===> Making clean in directory: ISOROPIA <==="
	@rm -f *.o *.mod *.a *.x

slowclean:
	@echo "===> Making slowclean in directory: ISOROPIA <==="
	@rm -f $(OBJECTS) $(MODULES) $(LIBRARY) $(LIB)/$(LIBRARY)

debug:
	@echo "Targets : $(MAKECMDGOALS)"
	@echo "ROOT    : $(ROOT)"
	@echo "LIB     : $(LIB)"
	@echo "MOD     : $(MOD)"
	@echo "F90     : $(F90)"
	@echo "OBJECTS : $(OBJECTS)"
	@echo "MODULES : $(MODULES)"
	@echo "LIBRARY : $(LIBRARY)"

###############################################################################
###                                                                         ###
###  Dependencies listing                                                   ###
###  (grep "USE " to get the list of module references!)                    ###
###                                                                         ###
###  From this list of dependencies, the "make" utility will figure out     ###
###  correct order of compilation (so we don't have to do that ourselves).  ###
###  This also allows us to compile on multiple processors with "make -j".  ###
###                                                                         ###
###  NOTES:                                                                 ###
###  (1) Only specify object-file dependencies that are within this         ###
###       directory.  Object files in other directories will be referenced  ### 
###       at link-time.                                                     ###
###  (2) For "make -jN" (i.e. compile N files simultaneously), all files    ###
###       in this directory must have a listed dependency.                  ###
###                                                                         ###
###############################################################################

%.o: %.f isrpia.inc
ifeq ($(PRECISION),8)
	$(F90ISO) $(R8) -c $<
endif
ifeq ($(PRECISION),4)
	$(F90ISO) $(R8) -c $<
endif

isocom.o:  isocom.f isrpia.inc
isofwd.o:  isofwd.f isrpia.inc
isorev.o:  isorev.f isrpia.inc
#EOC

