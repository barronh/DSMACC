-include ../Makefile.defs

MODEL=dsmacc

##############################################################################

PROG = ../bin/$(MODEL) 

all: libtuv libfastjx $(PROG)

libtuv:
	cd ../tuv_new; make

libfastjx:
	cd ../UCI_fastJX72e/; make

# complete list of all f90 source files
SRCS1 = $(wildcard ./$(MODEL)_*.f90)
SRCS3 = constants.f90 

# the object files are the same as the source files but with suffix ".o"
OBJS1 := $(SRCS1:.f90=.o) 
OBJS := $(OBJS1)
LIBS = -L../tuv_new/ -ltuv -L../UCI_fastJX72e/ -lfastjx

# the dependencies depend on the link
# the executable depends on depend and also on all objects
# the executable is created by linking all objects
$(PROG): $(MODEL)_Main.f90 $(OBJS) $(MODEL)_Main.f90
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
	rm -f $(SRCS1) Makefile_$(MODEL) $(MODEL).map tuvlog.txt ../bin/$(MODEL)
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
./$(MODEL)_Function.o: ./$(MODEL)_Parameters.o
./$(MODEL)_Global.o: ../tuv_new/params
./$(MODEL)_Global.o: ./$(MODEL)_Parameters.o
./$(MODEL)_Initialize.o: ./$(MODEL)_Global.o ./$(MODEL)_Parameters.o ./$(MODEL)_Util.o
./$(MODEL)_Integrator.o: ./$(MODEL)_Function.o ./$(MODEL)_Global.o
./$(MODEL)_Integrator.o: ./$(MODEL)_Jacobian.o ./$(MODEL)_LinearAlgebra.o
./$(MODEL)_Integrator.o: ./$(MODEL)_Parameters.o ./$(MODEL)_Rates.o
./$(MODEL)_Jacobian.o: ./$(MODEL)_JacobianSP.o ./$(MODEL)_Parameters.o
./$(MODEL)_LinearAlgebra.o: ./$(MODEL)_JacobianSP.o ./$(MODEL)_Parameters.o
./$(MODEL)_Main.o: ./$(MODEL)_Global.o ./$(MODEL)_Integrator.o ./$(MODEL)_Monitor.o
./$(MODEL)_Main.o: ./$(MODEL)_Parameters.o ./$(MODEL)_Rates.o ./$(MODEL)_Util.o
./$(MODEL)_Model.o: ./$(MODEL)_Function.o ./$(MODEL)_Global.o ./$(MODEL)_Integrator.o
./$(MODEL)_Model.o: ./$(MODEL)_Jacobian.o ./$(MODEL)_LinearAlgebra.o
./$(MODEL)_Model.o: ./$(MODEL)_Monitor.o ./$(MODEL)_Parameters.o ./$(MODEL)_Precision.o
./$(MODEL)_Model.o: ./$(MODEL)_Rates.o ./$(MODEL)_Stoichiom.o ./$(MODEL)_Util.o
./$(MODEL)_Parameters.o: ./$(MODEL)_Precision.o
./$(MODEL)_Rates.o: ./$(MODEL)_Global.o ./$(MODEL)_Parameters.o
./$(MODEL)_Stoichiom.o: ./$(MODEL)_Parameters.o ./$(MODEL)_StoichiomSP.o
./$(MODEL)_StoichiomSP.o: ./$(MODEL)_Precision.o
./$(MODEL)_Util.o: ./$(MODEL)_Global.o ./$(MODEL)_Integrator.o ./$(MODEL)_Monitor.o
./$(MODEL)_Util.o: ./$(MODEL)_Parameters.o
