# Makefile for JX
#------------------------------------------------------------------------------

OBJ   = cmn_fjx_mod.o fjx_sub_mod.o fjx_init_mod.o cld_sub_mod.o fjx72.o
DOBJ   = cmn_fjx_mod.o fjx_sub_mod.o fjx_init_mod.o cld_sub_mod.o dsmacc_fjx72.o dsmacc_test.o

FC    = gfortran
FCFLAGS = -O2 -fconvert=big-endian

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

JX70    : $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $^

DSMACC_test: $(DOBJ)
	$(FC) $(FCFLAGS) -o $@ $^
#------------------------ Cleaning -----------------------------------

clean :
	rm fort.* *.o *.s *.mod core

