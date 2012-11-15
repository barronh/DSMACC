
*=============================================================================*

      SUBROUTINE setno2(no2new,
     $     nz,z,nw,wl,
     $     no2xs, tlay, dcol,
     $     dtno2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of NO2 molecules, and corresponding absorption=*
*=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
*=  that allows scaling of the entire profile to a given overhead NO2        =*
*=  column amount.                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NO2NEW - REAL, overhead NO2 column amount (molec/cm^2) to which       (I)=*
*=           profile should be scaled.  If NO2NEW < 0, no scaling is done    =*
*=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
*=           grid                                                            =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  NO2XS  - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)=*
*=           each specified wavelength                                       =*
*=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
*=  DTNO2  - REAL, optical depth due to NO2 absorption at each            (O)=*
*=           specified altitude at each specified wavelength                 =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=51)

********
* input:
********

* grids:

      REAL wl(kw)
      REAL z(kz)
      INTEGER nw
      INTEGER nz
      REAL no2new

* mid-layer temperature, layer air column

      REAL tlay(kz), dcol(kz)

********
* output:
********

      REAL dtno2(kz,kw)

********
* local:
********

* absorption cross sections 

      REAL no2xs(kw)
      REAL cz(kz)

* nitrogen dioxide profile data:

      REAL zd(kdata), no2(kdata)
      REAL cd(kdata)
      REAL hscale
      REAL colold, scale
      REAL sno2

* other:

      INTEGER i, l, nd

********
* External functions:
********

      REAL fsum
      EXTERNAL fsum

*_______________________________________________________________________
* Data input:

* Example:  set to 1 ppb in lowest 1 km, set to zero above that.
* - do by specifying concentration at 3 altitudes.

      nd = 3
      zd(1) = 0.
      no2(1) = 1. * 2.69e10

      zd(2) = 1.
      no2(2) = 1. * 2.69e10

      zd(3) = zd(2)* 1.000001
      no2(3) = 10./largest

* compute column increments (alternatively, can specify these directly)

      DO 11, i = 1, nd - 1
         cd(i) = (no2(i+1)+no2(i)) * 1.E5 * (zd(i+1)-zd(i)) / 2. 
   11 CONTINUE

* Include exponential tail integral from top level to infinity.
* fold tail integral into top layer
* specify scale height near top of data (use ozone value)

      hscale = 4.50e5
      cd(nd-1) = cd(nd-1) + hscale * no2(nd)

***********
*********** end data input.

* Compute column increments and total column on standard z-grid.  

      CALL inter3(nz,z,cz, nd,zd,cd, 1)

**** Scaling of vertical profile by ratio of new to old column:
* If old column is near zero (less than 1 molec cm-2), 
* use constant mixing ratio profile (nominal 1 ppt before scaling) 
* to avoid numerical problems when scaling.

      IF(fsum(nz-1,cz) .LT. 1.) THEN
         DO i = 1, nz-1
            cz(i) = 1.E-12 * dcol(i)
         ENDDO
      ENDIF
      colold = fsum(nz-1, cz)
      scale =  2.687e16 * no2new / colold
      DO i = 1, nz-1
         cz(i) = cz(i) * scale
      ENDDO

************************************
* calculate optical depth for each layer.  Output: dtno2(kz,kw)

      DO 20, l = 1, nw-1
         sno2 = no2xs(l)
         DO 10, i = 1, nz-1
            dtno2(i,l) = cz(i)*sno2
   10    CONTINUE
   20 CONTINUE
*_______________________________________________________________________

      RETURN
      END
