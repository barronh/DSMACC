*=============================================================================*

      SUBROUTINE odo3(nz,z,nw,wl,o3xs,c, dto3)

*-----------------------------------------------------------------------------*
*=  NAME:  Optical Depths of O3
*=  PURPOSE:                                                                 =*
*=  Compute ozone optical depths as a function of altitude and wavelength    +*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
*=           grid                                                            =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  O3XS   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (I)=*
*=           each specified wavelength and altitude                          =*
*=  C      - REAL, ozone vertical column increments, molec cm-2, for each (I)=*
*=           layer                                                           =*
*=  DTO3   - REAL, optical depth due to ozone absorption at each          (O)=*
*=           specified altitude at each specified wavelength                 =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

********
* input:
********

* grids:

      INTEGER nw
      INTEGER nz
      REAL wl(kw)
      REAL z(kz)

* ozone absorption cross section, functions of wavelength and altitude

      REAL o3xs(kz,kw)

* ozone vertical column increments

      REAL c(kz)

********
* output:
********

      REAL dto3(kz,kw)

********
* internal:
********

      INTEGER iw, iz

*_______________________________________________________________________

* calculate ozone optical depth for each layer, with temperature 
* correction.  Output, dto3(kz,kw)

      DO 20, iw = 1, nw-1
         DO 10, iz = 1, nz - 1
            dto3(iz,iw) = c(iz) * o3xs(iz,iw)
   10    CONTINUE
   20 CONTINUE

*_______________________________________________________________________

      RETURN
      END
