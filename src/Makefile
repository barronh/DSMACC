F90     = gfortran -x f95-cpp-input -arch i386 -m32 -O3
FC     = gfortran -x f95-cpp-input -arch i386 -m32 -O3
LD     = gfortran -arch i386 -m32 
F90FLAGS   = -I/usr/local/include -fbounds-check -fimplicit-none -fno-automatic 
##############################################################################

PROG = ../bin/dsmacc 

# complete list of all f90 source files
SRCS1 = $(wildcard dsmacc_*.f90)
SRCS2 = $(wildcard ../tuv_old/*.f)
SRCS3 = constants.f90 

# the object files are the same as the source files but with suffix ".o"
OBJS1 := $(SRCS1:.f90=.o) 
OBJS2 := $(SRCS2:.f=.o)
OBJS3 := $(SRCS3:.f90=.o) 
OBJS := $(OBJS1) $(OBJS2) $(OBJS3)
MAKEFILE_INC = depend.mk

# If you don't have the perl script sfmakedepend, get it from:
# http://www.arsc.edu/~kate/Perl
F_makedepend = ../sfmakedepend --file=$(MAKEFILE_INC)

# the dependencies depend on the link
# the executable depends on depend and also on all objects
# the executable is created by linking all objects
$(PROG): dsmacc_Main.f90 $(MAKEFILE_INC) $(OBJS) dsmacc_Main.f90
	@mkdir -p ../bin
	$(LD) $(F90FLAGS) $(OBJS) -o $@

# update file dependencies
depend $(MAKEFILE_INC): $(SRCS1) $(SRCS2) $(SRCS3) dsmacc_Main.f90
	$(F_makedepend) dsmacc_*.f90 constants.f90 tuv_old/*.f

clean:
	rm -f $(OBJS) *.mod *.log *~ depend.mk.old $(SRCS1) Makefile_dsmacc dsmacc.map tuvlog.txt ../bin/dsmacc test/*_1.dat.diurnal test/*_1.dat.timeseries

distclean: clean
	rm -f $(PROG)
	rm -f depend.mk* 
	rm -f *.nc
	rm -f *.dat

# all object files *.o depend on their source files *.f90
# the object files are created with the "-c" compiler option
%.o: %.f90
	$(F90) $(F90FLAGS) $(LINCLUDES) -c $<

tuv_old/%.o: %.f
	$(F90) $(F90FLAGS) $(LINCLUDES) -c $<

# list of dependencies (via USE statements)
include depend.mk

