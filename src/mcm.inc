// global.inc is required
#INCLUDE global.inc

// rate includes definition of mcm generic rate
// constant variables (e.g., 
#INCLUDE mcmratedef.inc
// io handles input output in the 
// DSMACC fashion
#INCLUDE io.inc

#DOUBLE ON 

// computer language for code produced by kpp
#LANGUAGE FORTRAN90

// initial concentrations
#INITVALUES

// initial mixing ratios
ALL_SPEC = 0.;

// integrator should usually be rosenbrock
#INTEGRATOR rosenbrock
//#INTEGRATOR kpp_lsode
// #INTEGRATOR ros2_manual
// #INTEGRATOR radau5

// set indices of unused species to ind_*=0
#DUMMYINDEX ON

// do not create Stoichiom files
#STOICMAT ON

// do not create Hessian files
#HESSIAN OFF

// do not create Mex files
#MEX OFF

// equation IDs
#EQNTAGS on
