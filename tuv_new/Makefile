include ../Makefile.defs
PROG=libtuv.a
SRCS =	main.f90 functs.f grids.f la_srb.f numer.f odo3.f odrl.f orbit.f \
	qys.f rdetfl.f rdinp.f rdxs.f rtrans.f rxn.f savout.f setaer.f \
	setalb.f setcld.f setno2.f seto2.f setsnw.f setso2.f sphers.f \
	swbiol.f swchem.f swphys.f TUV.f vpair.f vpo3.f vptmp.f wshift.f

OBJS =	functs.o grids.o la_srb.o numer.o odo3.o odrl.o orbit.o qys.o \
	rdetfl.o rdinp.o rdxs.o rtrans.o rxn.o savout.o setaer.o setalb.o \
	setcld.o setno2.o seto2.o setsnw.o setso2.o sphers.o swbiol.o \
	swchem.o swphys.o TUV.o vpair.o vpo3.o vptmp.o wshift.o

LIBS =	


all: $(PROG)

$(PROG): $(OBJS)
	ar -rcs $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f.o:
	$(FC) $(TUVF90FLAGS) -c $<

main.o: params
grids.o: params
la_srb.o: params
odo3.o: params
odrl.o: params
rdetfl.o: params
rdinp.o: params
rdxs.o: params
rtrans.o: params
rxn.o: params
savout.o: params
setaer.o: params
setalb.o: params
setcld.o: params
setno2.o: params
seto2.o: params
setsnw.o: params
setso2.o: params
sphers.o: params
swbiol.o: params
swchem.o: params
swphys.o: params
TUV.o: params
vpair.o: params
vpo3.o: params
vptmp.o: params
wshift.o: params
