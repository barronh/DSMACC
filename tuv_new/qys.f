* This file contains subroutines used for calculation of quantum yields for 
* various photoreactions:
*     qyacet - q.y. for acetone, based on Blitz et al. (2004)

********************************************************************************

      SUBROUTINE qyacet(w, T, M, fco, fac)

* Compute acetone quantum yields according to the parameterization of:
* Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and M. P. Chipperfield 
*       (2004), Pressure and temperature-dependent quantum yields for the 
*       photodissociation of acetone between 279 and 327.5 nm, Geophys. 
*       Res. Lett., 31, L06111, doi:10.1029/2003GL018793.

      IMPLICIT NONE

* input:
* w = wavelength, nm
* T = temperature, K
* m = air number density, molec. cm-3

      REAL w, T, M

* internal:

      REAL a0, a1, a2, a3, a4
      REAL b0, b1, b2, b3, b4
      REAL c3
      REAL cA0, cA1, cA2, cA3, cA4

* output
* fco = quantum yield for product CO
* fac = quantum yield for product CH3CO (acetyl radical)

      REAL fco, fac

*** set out-of-range values:
* use low pressure limits for shorter wavelengths
* set to zero beyound 327.5

      IF(w .LT. 279.) THEN
         fco = 0.05
         fac = 0.95
         RETURN
      ENDIF

      IF(w .GT. 327.5 ) THEN
         fco = 0.
         fac = 0.
         RETURN
      ENDIF

*** CO (carbon monoxide) quantum yields:

      a0 = 0.350 * (T/295.)**(-1.28)
      b0 = 0.068 * (T/295.)**(-2.65)
      cA0 = exp(b0*(w - 248.)) * a0 / (1. - a0)

      fco = 1. / (1 + cA0)

*** CH3CO (acetyl radical) quantum yields:

      IF(w .GE. 279. .AND. w .LT. 302.) THEN

         a1 = 1.600E-19 * (T/295.)**(-2.38)
         b1 = 0.55E-3   * (T/295.)**(-3.19)
         cA1 = a1 * EXP(-b1*((1.e7/w) - 33113.))
 
         fac = (1. - fco) / (1 + cA1 * M)

      ENDIF

      IF(w .GE. 302. .AND. w .LT. 327.5) THEN

         a2 = 1.62E-17 * (T/295.)**(-10.03)
         b2 = 1.79E-3  * (T/295.)**(-1.364)
         cA2 = a2 * EXP(-b2*((1.e7/w) - 30488.))


         a3 = 26.29   * (T/295.)**(-6.59)
         b3 = 5.72E-7 * (T/295.)**(-2.93)
         c3 = 30006   * (T/295.)**(-0.064)
         ca3 = a3 * EXP(-b3*((1.e7/w) - c3)**2)


         a4 = 1.67E-15 * (T/295.)**(-7.25)
         b4 = 2.08E-3  * (T/295.)**(-1.16)
         cA4 = a4 * EXP(-b4*((1.e7/w) - 30488.))

         fac = (1. - fco) * (1. + cA3 + cA4 * M) /
     $        ((1. + cA3 + cA2 * M)*(1. + cA4 * M))

      ENDIF

      RETURN
      END

********************************************************************************
