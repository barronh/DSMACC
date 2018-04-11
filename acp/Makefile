export KPP_HOME=$(PWD)/../kpp

all: acp acp_conc_1.dat figs

help:
	@echo "$ make all"
	@echo "will build and run a timeseries box model"
	@echo "using the acp.kpp and tuv photolysis"
	
.SECONDARY: acp_rate_1.dat

figs: acp_conc_1.dat_O3.png

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
	./$< acp_init.dat $*_conc_ $*_rate_  > $*.log 2>&1

src/Makefile:
	mkdir -p src && cd src && ln -s ../../src/* .

src/dsmacc_Main.f90: src/Makefile acp.kpp ../src/global.inc ../driver.f90 ../src/tuv.inc ../src/io.inc
	make -C src distclean
	cd src && $(KPP_HOME)/bin/kpp ../acp.kpp dsmacc

acp: src/dsmacc_Main.f90
	make -C src PROG=../acp TOPDIR=../../

clean:
	rm -f acp acp_conc_*.dat acp_rate_*.dat acp.log *.png DATA* INPUTS photolysis.txt CTM_GrdCld.dat FJX_scat-UMa.dat FJX_scat-aer.dat FJX_scat-cld.dat atmos_std.dat FJX_spec.dat FJX_j2j.dat
	rm -rf src
