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

all: figs

figs: timeseriesfigs


INPUTS:
	ln -s ../tuv_new/INPUTS .
	ln -s ../tuv_new/DATA* .

.SECONDARY: $(MODEL) Spec_diurnal_1.dat Spec_timeseries_1.dat Rate_diurnal_1.dat Rate_timeseries_1.dat

timeseriesfigs: ../plot.py Spec_timeseries_1.dat
	@python $^

Spec_timeseries_1.dat: $(MODEL) Init_cons.dat.timeseries
	ln -fs Init_cons.dat.timeseries Init_cons.dat
	./$< Init_cons.dat.timeseries Spec_timeseries_ Rate_timeseries_  > timeseries.log 2>&1

Spec_diurnal_1.dat: $(MODEL) Init_cons.dat.diurnal
	ln -fs Init_cons.dat.diurnal Init_cons.dat
	./$< Init_cons.dat.diurnal Spec_diurnal_ Rate_diurnal_ 1> diurnal.log 2>&1

$(MODEL): $(MODEL).kpp ../src/global.inc ../src/mcmratedef.inc ../src/io.inc ../driver.f90 ../src/tuv.inc ../src/fastjx.inc ../kpp/bin/kpp
	cd ../src && rm -f depend.mk && ../kpp/bin/kpp ../working/$(MODEL).kpp dsmacc && make && mv ../bin/dsmacc ../working/$(MODEL)

clean:
	rm -f $(MODEL) Spec*_1.dat.* Rate_1.dat.*
	rm INPUTS DATA*
