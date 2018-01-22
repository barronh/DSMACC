export KPP_HOME=$(PWD)/../kpp

all: geoschem geoschem_conc_1.dat figs

help:
	@echo "$ make all"
	@echo "will build and run a timeseries box model"
	@echo "using the geoschem_subset.kpp and tuv photolysis"
	
.SECONDARY: geoschem_rate_1.dat

figs: geoschem_conc_1.dat_O3.png

%_O3.png: %
	python ../plot.py -v "O3;label='Ozone',color='red'" $<

inputs:
	ln -s ../tuv_new/INPUTS .
	ln -s ../tuv_new/DATA* .
	ln -s ../data/photolysis.txt .
	ln -fs ../UCI_fastJX72e/CTM_GrdCld.dat ../UCI_fastJX72e/FJX_scat-UMa.dat ../UCI_fastJX72e/FJX_scat-aer.dat ../UCI_fastJX72e/FJX_scat-cld.dat ../UCI_fastJX72e/atmos_std.dat .
	ln -fs ../UCI_fastJX72e/FJX_spec_gc.dat FJX_spec.dat
	ln -fs ../UCI_fastJX72e/FJX_j2j_gc.dat FJX_j2j.dat

%_conc_1.dat: % %_init.dat inputs
	./$< geoschem_init.dat $*_conc_ $*_rate_  > $*.log 2>&1

src/Makefile:
	mkdir -p src && cd src && ln -s ../../src/* .

src/dsmacc_Main.f90: src/Makefile gckpp.kpp geoschem.inc ../src/global.inc ../driver.f90 ../src/tuv.inc ../src/io.inc
	make -C src distclean
	cd src && $(KPP_HOME)/bin/kpp ../gckpp.kpp dsmacc ; test $${?} -ne 4

geoschem: src/dsmacc_Main.f90
	make -C src PROG=../geoschem TOPDIR=../../

clean:
	rm -f geoschem geoschem_conc_*.dat geoschem_rate_*.dat geoschem.log *.png DATA* INPUTS photolysis.txt CTM_GrdCld.dat FJX_scat-UMa.dat FJX_scat-aer.dat FJX_scat-cld.dat atmos_std.dat FJX_spec.dat FJX_j2j.dat
	rm -rf src
