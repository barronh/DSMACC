##############################################################################
-include Makefile.defs
export KPP_HOME=$(PWD)/kpp

all: check bin/dsmacc

source src/dsmacc_Main.f90: dsmacc.kpp global.inc rate.inc util.inc driver.f90 photolysis.inc kpp/bin/kpp
	cd src && rm -f depend.mk && ../kpp/bin/kpp ../dsmacc.kpp dsmacc

bin/dsmacc: src/dsmacc_Main.f90 Makefile.defs
	cd src && make

kpp/bin/kpp: Makefile.defs
	cd kpp && make

check: kpp/bin/kpp
	cd test && make

clean:
	cd src && make -i distclean
	cd test && make -i clean
	cd kpp && make clean

distclean: clean
	@rm -f bin/dsmacc
	cd kpp && make -i distclean
	@rm -rf autom4te.cache config.status config.log Makefile.deps Makefile.defs