F90     = gfortran -x f95-cpp-input -arch i386 -m32 -O3
FC     = gfortran -x f95-cpp-input -arch i386 -m32 -O3
LD     = gfortran -arch i386 -m32 
F90FLAGS   = -I/usr/local/include -fbounds-check -fimplicit-none -fno-automatic 
##############################################################################

export KPP_HOME=$(PWD)/kpp

all: check bin/dsmacc

source src/dsmacc_Main.f90: dsmacc.kpp global.inc rate.inc util.inc driver.f90 photolysis.inc kpp/bin/kpp
	cd src && rm -f depend.mk && ../kpp/bin/kpp ../dsmacc.kpp dsmacc

bin/dsmacc: src/dsmacc_Main.f90
	cd src && make

kpp/bin/kpp:
	cd kpp && make

check:
	cd test && make

clean:
	cd src && make -i distclean
	cd test && make -i clean
	cd kpp && make clean

distclean: clean
	rm -f bin/dsmacc
	cd kpp && make -i distclean