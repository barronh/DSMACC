export KPP_HOME=$(PWD)/../kpp
check:
	make cri.check
	make geos.check

%.check:
	make $*.check.timeseries.log
	make $*.check.diurnal.log
	@echo "Timeseries"
	@tail -n 3 $*.check.timeseries.log
	@echo "Diurnal"
	@tail -n 3 $*.check.diurnal.log

.SECONDARY: cri geos geos.Spec_1.dat.diurnal geos.Spec_1.dat.timeseries cri.Spec_1.dat.diurnal cri.Spec_1.dat.timeseries
%.check.diurnal.log: %.Spec_1.dat.diurnal check.py
	@python check.py check.$*.Spec_1.dat.diurnal $*.Spec_1.dat.diurnal > $@

%.check.timeseries.log: %.Spec_1.dat.timeseries check.py
	@python check.py check.$*.Spec_1.dat.timeseries $*.Spec_1.dat.timeseries > $@

%.Spec_1.dat.diurnal: % %.Init_cons.dat.diurnal
	ln -fs $*.Init_cons.dat.diurnal Init_cons.dat && \
	./$* >& $*.diurnal.log && mv Spec_1.dat $*.Spec_1.dat.diurnal && mv Rate_1.dat $*.Rate_1.dat.diurnal

%.Spec_1.dat.timeseries: % %.Init_cons.dat.timeseries
	ln -fs $*.Init_cons.dat.timeseries Init_cons.dat && \
	./$* >& $*.timeseries.log && mv Spec_1.dat $*.Spec_1.dat.timeseries && mv Rate_1.dat $*.Rate_1.dat.timeseries

cri: cri.kpp ../global.inc ../rate.inc ../util.inc ../driver.f90 ../photolysis.inc ../kpp/bin/kpp
	cd ../src && rm -f depend.mk && ../kpp/bin/kpp ../test/cri.kpp dsmacc && make && mv ../bin/dsmacc ../test/cri

geos: geos.kpp ../global.inc ../rate.inc ../util.inc ../driver.f90 ../photolysis.inc ../kpp/bin/kpp
	cd ../src && rm -f depend.mk && ../kpp/bin/kpp ../test/geos.kpp dsmacc && make && mv ../bin/dsmacc ../test/geos

clean:
	rm -f cri.Spec_1.dat.diurnal cri.Rate_1.dat.diurnal geos.Spec_1.dat.diurnal geos.Rate_1.dat.diurnal *.diurnal.log *.timeseries.log cri geos