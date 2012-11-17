check: ../bin/dsmacc
	make check_timeseries
	make check_diurnal
	@echo "Timeseries"
	@tail -n 1 timeseries.log
	@echo "Diurnal"
	@tail -n 1 diurnal.log

check_diurnal: ../bin/dsmacc Spec_1.dat.diurnal
	@python check.py Spec_1.dat.diurnal.check Spec_1.dat.diurnal >> diurnal.log

check_timeseries: Spec_1.dat.timeseries
	@python check.py Spec_1.dat.timeseries.check Spec_1.dat.timeseries >> timeseries.log

Spec_1.dat.diurnal: ../bin/dsmacc Init_cons.dat.diurnal
	ln -fs Init_cons.dat.diurnal Init_cons.dat && \
	../bin/dsmacc >& diurnal.log && mv Spec_1.dat Spec_1.dat.diurnal && mv Rate_1.dat Rate_1.dat.diurnal

Spec_1.dat.timeseries: ../bin/dsmacc Init_cons.dat.timeseries
	ln -fs Init_cons.dat.timeseries Init_cons.dat && \
	../bin/dsmacc >& timeseries.log && mv Spec_1.dat Spec_1.dat.timeseries && mv Rate_1.dat Rate_1.dat.timeseries
	
clean:
	rm -f Spec_1.dat.diurnal Rate_1.dat.diurnal Spec_1.dat.timeseries Rate_1.dat.timeseries timeseries.log diurnal.log