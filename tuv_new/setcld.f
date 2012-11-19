
*=============================================================================*

      SUBROUTINE setcld(taucld,zbase,ztop
     $     ,nz,z,nw,wl,dtcld,omcld,gcld)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set cloud properties for each specified altitude layer.  Properties      =*
*=  may be wavelength dependent.                                             =*
*=  Assumes horizontally infinite homogeneous cloud layers.
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  DTCLD   - REAL, optical depth due to absorption by clouds at each     (O)=*
*=            altitude and wavelength                                        =*
*=  OMCLD   - REAL, single scattering albedo due to clouds at each        (O)=*
*=            defined altitude and wavelength                                =*
*=  GCLD    - REAL, cloud asymmetry factor at each defined altitude and   (O)=*
*=            wavelength                                                     =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=51)

***** input

* (grids)
      REAL wl(kw)
      REAL z(kz)
      INTEGER nz
      INTEGER nw

* new total cloud optical depth:

      REAL taucld
      REAL zbase, ztop

***** Output: 

      REAL dtcld(kz,kw), omcld(kz,kw), gcld(kz,kw)

***** specified default data:

      REAL zd(kdata), cd(kdata), omd(kdata), gd(kdata)
      REAL womd(kdata), wgd(kdata)
      REAL cldold

* other:

      REAL cz(kz)
      REAL omz(kz)
      REAL gz(kz)
      INTEGER i, iw, n
      REAL scale

* External functions:
      REAL fsum
      EXTERNAL fsum
*_______________________________________________________________________

* Set up clouds:
* All clouds are assumed to be infinite homogeneous layers
* Can have different clouds at different altitudes.
*   If multiple cloud layers are specified, non-cloudy layers
*   between them (if any) must be assigned zero optical depth.
* Set cloud optical properties:
*   cd(i) = optical depth of i_th cloudy layer
*   omd(i) = singel scattering albedo of i_th  cloudy layer
*   gd(i) = asymmetry factorof i_th  cloudy layer
* Cloud top and bottom can be set to any height zd(i), but if they don't
* match the z-grid (see subroutine gridz.f), they will be interpolated to
* the z-grid.

* Example:  set two separate cloudy layers:
*  cloud 1:  
*     base = 4 km
*     top  = 7 km
*     optical depth = 20.  (6.67 per km)
*     single scattering albedo = 0.9999
*     asymmetry factor = 0.85
*  cloud 2:
*     base = 9 km
*     top  = 11 km
*     optical depth = 5.  (2.50 per km)
*     single scattering albedo = 0.99999
*     asymmetry factor = 0.85

      n = 2
      
* cloud 1

      zd(1) = zbase
      cd(1) = taucld
      omd(1) = .9999
      gd(1) = .85
      zd(2) = ztop

******************
* compute integrals and averages over grid layers:
* for g and omega, use averages weighted by optical depth

      DO 10, i = 1, n-1
         womd(i) = omd(i) * cd(i)
         wgd(i) = gd(i) * cd(i)
 10   CONTINUE
      CALL inter3(nz,z,cz,  n, zd,cd, 0)
      CALL inter3(nz,z,omz, n, zd,womd, 0)
      CALL inter3(nz,z,gz , n, zd,wgd, 0)

      DO 15, i = 1, nz-1
         IF (cz(i) .GT. 0.) THEN
            omz(i) = omz(i)/cz(i)
            gz(i)  = gz(i) /cz(i)
         ELSE
            omz(i) = 1.
            gz(i) = 0.
         ENDIF
   15 CONTINUE
      
* assign at all wavelengths
* (can move wavelength loop outside if want to vary with wavelength)

      DO 20, iw = 1, nw-1
         DO 25, i = 1, nz-1
            dtcld(i,iw) = cz(i)
            omcld(i,iw) = omz(i)
            gcld (i,iw) = gz(i)
 25      CONTINUE
 20   CONTINUE
*_______________________________________________________________________

      RETURN
      END
