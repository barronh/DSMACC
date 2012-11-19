*=============================================================================*

      SUBROUTINE setso2(ipbl, zpbl, xpbl,
     $     so2new, nz, z, nw, wl, so2xs, 
     $     tlay, dcol,
     $     dtso2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of SO2 molecules, and corresponding absorption=*
*=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
*=  that allows scaling of the entire profile to a given overhead SO2        =*
*=  column amount.                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  SO2NEW - REAL, overhead SO2 column amount (molec/cm^2) to which       (I)=*
*=           profile should be scaled.  If SO2NEW < 0, no scaling is done    =*
*=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
*=           grid                                                            =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  SO2XS  - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)=*
*=           each specified wavelength                                       =*
*=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
*=  DTSO2  - REAL, optical depth due to SO2 absorption at each            (O)=*
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
      REAL so2new

* mid-layer temperature and layer air column

      REAL tlay(kz), dcol(kz)

********
* output:
********

      REAL dtso2(kz,kw)

********
* local:
********

* absorption cross sections 

      REAL so2xs(kw)
      REAL cz(kz)

* sulfur dioxide profile data:

      REAL zd(kdata), so2(kdata)
      REAL cd(kdata)
      REAL hscale
      REAL colold, scale
      REAL sso2
      REAL zpbl, xpbl
      INTEGER ipbl

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
      so2(1) = 1. * 2.69e10

      zd(2) = 1.
      so2(2) = 1. * 2.69e10

      zd(3) = zd(2)* 1.000001
      so2(3) = 10./largest

* compute column increments (alternatively, can specify these directly)

      DO 11, i = 1, nd - 1
         cd(i) = (so2(i+1)+so2(i)) * 1.E5 * (zd(i+1)-zd(i)) / 2. 
   11 CONTINUE

* Include exponential tail integral from top level to infinity.
* fold tail integral into top layer
* specify scale height near top of data (use ozone value)

      hscale = 4.50e5
      cd(nd-1) = cd(nd-1) + hscale * so2(nd)

***********
*********** end data input.

* Compute column increments on standard z-grid.  

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
      colold = fsum(nz-1,cz)
      scale =  2.687e16 * so2new / colold
      DO i = 1, nz-1
         cz(i) = cz(i) * scale
      ENDDO

*! overwrite for specified pbl height, set concentration here

      IF(ipbl .GT. 0) THEN
         write(*,*) 'pbl SO2 = ', xpbl, ' ppb'

         DO i = 1, nz-1
            IF (i .LE. ipbl) THEN
               cz(i) = xpbl*1.E-9 * dcol(i)
            ELSE
               cz(i) = 0.
            ENDIF
         ENDDO
      ENDIF

************************************
* calculate sulfur optical depth for each layer, with optional temperature 
* correction.  Output, dtso2(kz,kw)

      DO 20, l = 1, nw-1
         sso2 = so2xs(l)
         DO 10, i = 1, nz - 1

c Leaving this part in in case i want to interpolate between 
c the 221K and 298K data.
c
c            IF ( wl(l) .GT. 240.5  .AND. wl(l+1) .LT. 350. ) THEN
c               IF (tlay(i) .LT. 263.) THEN
c                  sso2 = s221(l) + (s263(l)-s226(l)) / (263.-226.) *
c     $                 (tlay(i)-226.)
c               ELSE
c                  sso2 = s263(l) + (s298(l)-s263(l)) / (298.-263.) *
c     $              (tlay(i)-263.)
c               ENDIF
c            ENDIF

            dtso2(i,l) = cz(i)*sso2

   10    CONTINUE
   20 CONTINUE
*_______________________________________________________________________

      RETURN
      END
