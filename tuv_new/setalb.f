
*=============================================================================*

      SUBROUTINE setalb(albnew,nw,wl,albedo)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set the albedo of the surface.  The albedo is assumed to be Lambertian,  =*
*=  i.e., the reflected light is isotropic, and independent of direction     =*
*=  of incidence of light.  Albedo can be chosen to be wavelength dependent. =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  ALBEDO  - REAL, surface albedo at each specified wavelength           (O)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input: (wavelength working grid data)

      INTEGER nw
      REAL wl(kw)

      REAL albnew

* output:
      REAL albedo(kw)

* local:
      INTEGER iw
*_______________________________________________________________________

      DO 10, iw = 1, nw - 1
         albedo(iw) = albnew
   10 CONTINUE

* alternatively, can input wavelenght-dependent values if avaialble.
*_______________________________________________________________________

      RETURN
      END
