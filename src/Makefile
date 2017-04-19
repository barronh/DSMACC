-include ../Makefile.defs

##############################################################################

PROG = ../bin/dsmacc 

# complete list of all f90 source files
SRCS1 = $(wildcard dsmacc_*.f90)
SRCS2 = $(wildcard ../tuv_new/*.f)
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
	$(FC) $(F90FLAGS) $(OBJS) -o $@

# update file dependencies
depend $(MAKEFILE_INC): $(SRCS1) $(SRCS2) $(SRCS3)
	$(F_makedepend) dsmacc_*.f90 constants.f90 tuv_new/*.f

clean:
	rm -f $(OBJS) *.mod *.log *~ depend.mk.old

distclean: clean
	rm -f $(PROG)
	rm -f depend.mk* 
	rm -f *.nc
	rm -f *.dat
	rm -f $(SRCS1) Makefile_dsmacc dsmacc.map tuvlog.txt ../bin/dsmacc
	rm -f fflags

# all object files *.o depend on their source files *.f90
# the object files are created with the "-c" compiler option
%.o: %.f90
	$(FC) $(F90FLAGS) $(LINCLUDES) -c $<

../tuv_new/%.o: ../tuv_new/%.f
	$(FC) $(TUVF90FLAGS) $(LINCLUDES) -o $@ -c $<

# list of dependencies (via USE statements)
include depend.mk
