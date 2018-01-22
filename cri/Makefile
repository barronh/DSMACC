export KPP_HOME=$(PWD)/../kpp

all: cri cri_conc_1.dat figs

help:
	@echo "$ make all"
	@echo "will build and run a timeseries box model"
	@echo "using the CRI_subset.kpp and tuv photolysis"
	
.SECONDARY: cri_rate_1.dat

figs: cri_conc_1.dat_O3.png

%_O3.png: %
	python ../plot.py -v "O3;label='Ozone',color='red'" $<

inputs:
	ln -s ../tuv_new/INPUTS .
	ln -s ../tuv_new/DATA* .
	ln -s ../data/photolysis.txt .
	ln -fs ../UCI_fastJX72e/CTM_GrdCld.dat ../UCI_fastJX72e/FJX_scat-UMa.dat ../UCI_fastJX72e/FJX_scat-aer.dat ../UCI_fastJX72e/FJX_scat-cld.dat ../UCI_fastJX72e/atmos_std.dat .
	ln -fs ../UCI_fastJX72e/FJX_spec_mcm.dat FJX_spec.dat
	ln -fs ../UCI_fastJX72e/FJX_j2j_mcm.dat FJX_j2j.dat

%_conc_1.dat: % %_init.dat inputs
	./$< cri_init.dat $*_conc_ $*_rate_  > $*.log 2>&1

src/Makefile:
	mkdir -p src && cd src && ln -s ../../src/* .

src/dsmacc_Main.f90: src/Makefile cri.kpp ../src/global.inc ../src/mcmratedef.inc ../driver.f90 ../src/tuv.inc ../src/io.inc
	make -C src distclean
	cd src && $(KPP_HOME)/bin/kpp ../cri.kpp dsmacc

cri: src/dsmacc_Main.f90
	make -C src PROG=../cri TOPDIR=../../

clean:
	rm -f cri cri_conc_*.dat cri_rate_*.dat cri.log *.png DATA* INPUTS photolysis.txt CTM_GrdCld.dat FJX_scat-UMa.dat FJX_scat-aer.dat FJX_scat-cld.dat atmos_std.dat FJX_spec.dat FJX_j2j.dat
	rm -rf src
