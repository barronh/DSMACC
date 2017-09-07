export KPP_HOME=$(PWD)/../kpp
MODEL=dsmacc

help:
	@echo "$ make all"
	@echo "will run a timeseries box model and"
	@echo "a diural steady-state."
	@echo ""
	@echo "By default, the model is dsmacc"
	@echo "from the dsmacc.kpp file"
	@echo ""
	@echo "To build another model from a "
	@echo "different .kpp file, type make -e MODEL=newname"

all: check

show: check.timeseries.log check.diurnal.log
	@echo "$(MODEL) Timeseries Check"
	@tail -n 3 check.timeseries.log
	@echo "$(MODEL) Diurnal Check"
	@tail -n 3 check.diurnal.log

check: check.timeseries.log check.diurnal.log show

INPUTS:
	ln -s ../tuv_new/INPUTS .
	ln -s ../tuv_new/DATA* .
.SECONDARY: $(MODEL) Spec_diurnal_1.dat Spec_timeseries_1.dat Rate_diurnal_1.dat Rate_timeseries_1.dat

check.diurnal.log: Spec_diurnal_1.dat ../test/check.py
	@python check.py check.Spec_1.dat.diurnal Spec_1.dat.diurnal > $@

check.timeseries.log: Spec_timeseries_1.dat ../test/check.py
	@python check.py check.Spec_1.dat.timeseries Spec_1.dat.timeseries > $@

Spec_timeseries_1.dat: $(MODEL) Init_cons.dat.timeseries
	ln -fs Init_cons.dat.timeseries Init_cons.dat
	./$< Init_cons.dat.timeseries Spec_timeseries_ Rate_timeseries_  > timeseries.log 2>&1

Spec_diurnal_1.dat: $(MODEL) Init_cons.dat.diurnal
	ln -fs Init_cons.dat.diurnal Init_cons.dat
	./$< Init_cons.dat.diurnal Spec_diurnal_ Rate_diurnal_ 1> diurnal.log 2>&1

$(MODEL): $(MODEL).kpp ../global.inc ../rate.inc ../util.inc ../driver.f90 ../photolysis.inc ../fastjx.inc ../kpp/bin/kpp
	cd ../src && rm -f depend.mk && ../kpp/bin/kpp ../working/$(MODEL).kpp dsmacc && make && mv ../bin/dsmacc ../working/$(MODEL)

clean:
	rm -f $(MODEL) Spec_1.dat.* Rate_1.dat.*
