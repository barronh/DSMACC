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