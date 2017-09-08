##############################################################################
-include Makefile.defs
export KPP_HOME=$(PWD)/kpp

all: testrun

testrun: kpp/bin/kpp
	cd cri_tuv && make

kpp/bin/kpp: Makefile.defs
	cd kpp && make

clean:
	cd src && make -i distclean
	cd cri_tuv && make -i clean
	cd tuv_new && make -i clean
	cd UCI_fastJX72e && make -i clean
	cd kpp && make -i clean

distclean: clean
	@rm -f bin/dsmacc
	cd kpp && make -i distclean
	@rm -rf autom4te.cache config.status config.log Makefile.deps Makefile.defs tuvlog.txt