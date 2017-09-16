-include ../Makefile.defs

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
	cd ../src && make distclean && $(KPP_HOME)/bin/kpp ../pysrc/cri_mcmj.kpp dsmacc #&& make && rm ../bin/dsmacc


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