-include ../Makefile.defs

##############################################################################
all:
	make inputs
	make source
	make pykpp
	make test

# link TUV inputs for the run
inputs:
	ln -fs ../tuv_new/INPUTS ../tuv_new/DATA* .

# build the standard cri_tuv model including f90 and objects
source:
	cd ../cri_tuv; make cri

# Run the model and make a plot
test:
	python driver.py
	python -c "import pandas as pd; data = pd.read_csv('cri_isop.tsv', delimiter = '\t'); data.eval('JDAY=time/3600./24.', inplace = True); ax = data.plot(x = 'JDAY', y = ['O3', 'C5H8']); ax.figure.savefig('cri_isop.png')"

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
	rm -f kpp.pyf photolysis.txt *.so *.tsv *.png INPUTS DATA*