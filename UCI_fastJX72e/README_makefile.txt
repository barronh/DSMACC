# Makefile for JX
#------------------------------------------------------------------------------

OBJ   = cmn_fjx_mod.o fjx_sub_mod.o fjx_init_mod.o cld_sub_mod.o fjx72.o

FC    = ifort
FLAGS = -O2 -ip -fpp2 -W0 -assume byterecl

%.o: %.f
	$(FC) $(FLAGS) -c $<

jx    : $(OBJ)
	$(FC) $(FLAGS) -o JX70 $(OBJ)

#------------------------ Cleaning -----------------------------------

clean :
	rm fort.* *.o *.s *.mod core




#PC compiler instructions (digital fortran)

#df -c cmn_fjx_mod.f90
#df -c fjx_sub_mod.f90
#df -c fjx_init_mod.f90
#df -c cld_sub_mod.f90
#          and if no other *.obj in directory:
#df  fjx72.f90 *.obj
