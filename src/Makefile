-include ../Makefile.defs

##############################################################################

PROG = ../bin/dsmacc 

all: libtuv libfastjx $(PROG)

libtuv:
	cd ../tuv_new; make

libfastjx:
	cd ../UCI_fastJX72e/; make

# complete list of all f90 source files
SRCS1 = $(wildcard ./dsmacc_*.f90)
SRCS3 = constants.f90 

# the object files are the same as the source files but with suffix ".o"
OBJS1 := $(SRCS1:.f90=.o) 
OBJS := $(OBJS1)
LIBS = -L../tuv_new/ -ltuv -L../UCI_fastJX72e/ -lfastjx

# the dependencies depend on the link
# the executable depends on depend and also on all objects
# the executable is created by linking all objects
$(PROG): dsmacc_Main.f90 $(OBJS) dsmacc_Main.f90
	@mkdir -p ../bin
	$(FC) $(F90FLAGS) $(OBJS) $(LIBS) -o $@

# update file dependencies
depend depend.mk: $(SRCS1) $(SRCS2) $(SRCS3)
	../sfmakedepend --file=depend.mk --srcdir=. --moddir=.

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

../UCI_fastJX72e/%.o: ../UCI_fastJX72e/%.F90
	$(FC) $(F90FLAGS) -J../UCI_fastJX72e $(LINCLUDES) -o $@ -c $<

# list of dependencies (via USE statements)
#include depend.mk
./dsmacc_Function.o: ./dsmacc_Parameters.o
./dsmacc_Global.o: ../tuv_new/params
./dsmacc_Global.o: ./dsmacc_Parameters.o
./dsmacc_Initialize.o: ./dsmacc_Global.o ./dsmacc_Parameters.o ./dsmacc_Util.o
./dsmacc_Integrator.o: ./dsmacc_Function.o ./dsmacc_Global.o
./dsmacc_Integrator.o: ./dsmacc_Jacobian.o ./dsmacc_LinearAlgebra.o
./dsmacc_Integrator.o: ./dsmacc_Parameters.o ./dsmacc_Rates.o
./dsmacc_Jacobian.o: ./dsmacc_JacobianSP.o ./dsmacc_Parameters.o
./dsmacc_LinearAlgebra.o: ./dsmacc_JacobianSP.o ./dsmacc_Parameters.o
./dsmacc_Main.o: ./dsmacc_Global.o ./dsmacc_Integrator.o ./dsmacc_Monitor.o
./dsmacc_Main.o: ./dsmacc_Parameters.o ./dsmacc_Rates.o ./dsmacc_Util.o
./dsmacc_Model.o: ./dsmacc_Function.o ./dsmacc_Global.o ./dsmacc_Integrator.o
./dsmacc_Model.o: ./dsmacc_Jacobian.o ./dsmacc_LinearAlgebra.o
./dsmacc_Model.o: ./dsmacc_Monitor.o ./dsmacc_Parameters.o ./dsmacc_Precision.o
./dsmacc_Model.o: ./dsmacc_Rates.o ./dsmacc_Stoichiom.o ./dsmacc_Util.o
./dsmacc_Parameters.o: ./dsmacc_Precision.o
./dsmacc_Rates.o: ./dsmacc_Global.o ./dsmacc_Parameters.o
./dsmacc_Stoichiom.o: ./dsmacc_Parameters.o ./dsmacc_StoichiomSP.o
./dsmacc_StoichiomSP.o: ./dsmacc_Precision.o
./dsmacc_Util.o: ./dsmacc_Global.o ./dsmacc_Integrator.o ./dsmacc_Monitor.o
./dsmacc_Util.o: ./dsmacc_Parameters.o
