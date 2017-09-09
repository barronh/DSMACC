export KPP_HOME=$(PWD)/../kpp

all: cri cri_conc_1.dat figs

help:
	@echo "$ make all"
	@echo "will build and run a timeseries box model"
	@echo "using the CRI_subset.kpp and tuv photolysis"
	
.SECONDARY: cri_rate_1.dat

figs: cri_conc_1.dat_O3.png

cri_conc_1.dat_O3.png:
	python ../plot.py -v "O3;label='Ozone',color='red'" cri_conc_1.dat

INPUTS:
	ln -s ../tuv_new/INPUTS .
	ln -s ../tuv_new/DATA* .

%_conc_1.dat: % %_init.dat INPUTS
	./$< cri_init.dat $*_conc_ $*_rate_  > $*.log 2>&1

cri: cri.kpp ../src/global.inc ../src/mcmratedef.inc ../driver.f90 ../src/tuv.inc ../src/io.inc
	cd ../src && $(KPP_HOME)/bin/kpp ../cri/cri.kpp dsmacc && make && mv ../bin/dsmacc ../cri/cri

clean:
	rm -f cri cri_conc_*.dat cri_rate_*.dat cri.log *.png DATA* INPUTS photolysis.txt
