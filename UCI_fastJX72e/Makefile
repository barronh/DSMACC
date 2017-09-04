include ../Makefile.defs
PROG =	libfastjx.a

SRCS =	dsmacc_fjx72.F90 cld_sub_mod.F90 fjx_init_mod.F90 fjx_sub_mod.F90  \
	cmn_fjx_mod.F90 

OBJS =	cld_sub_mod.o cmn_fjx_mod.o fjx_dsmacc.o fjx_init_mod.o \
	fjx_sub_mod.o


all: $(PROG)

$(PROG): $(OBJS)
	ar -rcs $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
	$(FC) $(F90FLAGS) -c $<

fjx_dsmacc.o:  cmn_fjx_mod.o fjx_sub_mod.o fjx_init_mod.o cld_sub_mod.o
cld_sub_mod.o:   cmn_fjx_mod.o fjx_sub_mod.o
fjx_init_mod.o:  cmn_fjx_mod.o fjx_sub_mod.o
fjx_sub_mod.o:   cmn_fjx_mod.o