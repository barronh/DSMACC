export KPP_HOME=$(PWD)/../kpp
MODEL=dsmacc
help:
	@echo "`make all`"
	@echo "will run a timeseries box model and"
	@echo "a diural steady-state."
	@echo ""
	@echo "By default, the model is dsmacc"
	@echo "from the dsmacc.kpp file"
	@echo ""
	@echo "To build another model from a "
	@echo "different .kpp file, type `make -e MODEL=newname`"

all: check

show: check.timeseries.log check.diurnal.log
	@echo "$(MODEL) Timeseries Check"
	@tail -n 3 check.timeseries.log
	@echo "$(MODEL) Diurnal Check"
	@tail -n 3 check.diurnal.log

check: check.timeseries.log check.diurnal.log show

.SECONDARY: $(MODEL) Spec_1.dat.diurnal Spec_1.dat.timeseries Rate_1.dat.diurnal Rate_1.dat.timeseries

check.diurnal.log: Spec_1.dat.diurnal check.py
	@python check.py check.Spec_1.dat.diurnal Spec_1.dat.diurnal > $@

check.timeseries.log: Spec_1.dat.timeseries check.py
	@python check.py check.Spec_1.dat.timeseries Spec_1.dat.timeseries > $@

Spec_1.dat.timeseries: $(MODEL) Init_cons.dat.timeseries
	ln -fs Init_cons.dat.timeseries Init_cons.dat
	./$< > timeseries.log 2>&1
	mv Spec_1.dat Spec_1.dat.timeseries
	mv Rate_1.dat Rate_1.dat.timeseries

Spec_1.dat.diurnal: $(MODEL) Init_cons.dat.diurnal
	ln -fs Init_cons.dat.diurnal Init_cons.dat
	./$< 1> diurnal.log 2>&1
	mv Spec_1.dat Spec_1.dat.diurnal
	mv Rate_1.dat Rate_1.dat.diurnal

$(MODEL): $(MODEL).kpp ../global.inc ../rate.inc ../util.inc ../driver.f90 ../photolysis.inc ../fastjx.inc ../kpp/bin/kpp
	cd ../src && rm -f depend.mk && ../kpp/bin/kpp ../working/$(MODEL).kpp dsmacc && make && mv ../bin/dsmacc ../working/$(MODEL)

clean:
	rm -f $(MODEL) Spec_1.dat.* Rate_1.dat.*
