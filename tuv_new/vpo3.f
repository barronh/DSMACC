*=============================================================================*

      SUBROUTINE vpo3(ipbl, zpbl, mr_pbl, 
     $     to3new, nz, z, aircol, col)

*-----------------------------------------------------------------------------*
*=  NAME:  Vertical Profiles of Ozone = vpo3                                 =*
*=  PURPOSE:                                                                 =*
*=  Computes O3 column increments, col(i), molec cm-2 for each layer i of    =* 
*=  the working grid z(i = 1, nz).                                           =*
*=  Normally, col(i) values are computed from input vertical profiles of     =*
*=  concentrations (molec cm-3), that are then interpolated and integrated   =*
*=  over each layer.  The default example here uses the US Standard          =*
*=  Atmosphere (1976) mid-latitude concentration profile.                    =*
*=  Users can substitute their own concentration profile, as long as         =*
*=  appropriate adjustments are made in Section 1 to input the data, and in  =*
*=  section 2 to interpolate to the working grid.                            =*
*=  A scale factor is provided to allow changing the total column amount,    =*
*=  but conserving the shape of the profile.                                 =*
*=  An option to insert PBL pollutants is provided                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  TO3NEW - REAL,  Dobson Units, new total column amount from surface       =* 
*=            to space, to be scaled.  If TO3NEW < 0, no scaling is done.    =*
*=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
*=           grid                                                            =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  AIRCOL(KZ) = REAL, air column increment (molec cm-2), provided here in   =*
*=  case it is necessary to convert from mixing ratio units (e.g. ppb).      =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

********
* inputs:
********

      INCLUDE 'params'

*** from calling program:

      INTEGER nz
      REAL z(kz)
      REAL to3new
      REAL aircol(kz)

      REAL zpbl, mr_pbl
      INTEGER ipbl

* from data file:  concentration data as a function of altitude

      INTEGER kdata
      PARAMETER(kdata=150)
      REAL zd(kdata), xd(kdata)


********
* internal
********

      INTEGER i, nd
      REAL hscale
      REAL to3old, scale
      
      REAL con(kz)

      REAL rfact

********
* output:
********

      REAL col(kz)

********
* External functions:
********
      REAL fsum
      EXTERNAL fsum

* The objective of this subroutine is to calculate the vertical increments 
*   in the O3 column, for each layer of the working grid z(i = 1, nz).
* The input O3 profiles can be specified in different ways, and each case
*   will require careful consideration of the interpolation scheme.  Some
*   examples of possible input data are:
*    altitude vs. O3 concentration (number density), molec cm-3
*    altitude vs. O3 mixing ratio (e.g. parts per billion) relative to air
*    altitude vs. O3 column increments, molec cm-2, between specific altitudes.
* Special caution is required with mixed inputs, e.g. ppb in boundary layer and 
*   molec cm-3 above the boundary layer.

* _________SECTION 1:  Read in vertical profile of concentration
* Default is US Standard Atmosphere (1976)
* If a different vertical concentration profile is specified, the code
* in this section (Section 1) should be replaced accordingly

      WRITE(kout,*) 'ozone profile: USSA, 1976'
      OPEN(kin,FILE='DATAE1/ATM/ussa.ozone',STATUS='old')
      DO i = 1, 7
        READ(kin,*)
      ENDDO
      nd = 39
      DO i = 1, nd
         READ(kin,*) zd(i), xd(i)
      ENDDO
      CLOSE(kin)

* Ussa data stop at 74 km.  Add values up to 121 km, 
* assuming exponential decay from 74 km up, with scale height of
*  4.5 km.

      hscale = 4.5
      rfact = EXP(-1./hscale)
 10   CONTINUE
      nd = nd + 1
      zd(nd) = zd(nd-1) + 1.
      xd(nd) = xd(nd-1) * rfact
      IF(zd(nd) .GE. 121.) GO TO 19
      GO TO 10
 19   CONTINUE

*********** end data input.

* ________SECTION 2:  Compute column increments on standard z-grid.  

* linear interpolation

      CALL inter1(nz,z,con, nd,zd,xd)

* compute column increments

      DO i = 1, nz-1
         col(i) = 0.5 * (con(i) + con(i+1)) * (z(i+1) - z(i)) * 1.E5
      ENDDO

* Add exponential tail integral at top of atmosphere:
*   this is folded into layer nz-1, making this layer "heavy'.  
*   The layer nz is not used. The radiative transfer 
*   calculation is based on nz-1 layers (not nz).

      col(nz-1) = col(nz-1) + 1.E5 * hscale * con(nz)

***** Scaling to new total ozone
* to3old = total o3 column, in Dobson Units, old value
* to3new = total o3 column, in Dobson Units, new value
*    (1 DU = 2.687e16)
* If to3new is not negative, scale to new total column value, to3new :
* (to3new = 0. is a possible input, to see effect of zero ozone)

      IF (to3new .GT. nzero) THEN

         to3old = fsum(nz-1, col)/2.687e16
         IF(to3old .LT. pzero) STOP 'in vpo3: to3old is too small'
         scale = to3new/to3old
         DO i = 1, nz-1
            col(i) = col(i) * scale
            con(i) = con(i) * scale
         ENDDO
         con(nz) = con(nz) * scale

      ENDIF

*! overwrite column increments for specified pbl height
* use mixing ratio in pbl

      IF(ipbl .GT. 0) THEN
         write(*,*) 'pbl O3 = ', mr_pbl, ' ppb'

         DO i = 1, nz-1
            IF (i .LE. ipbl) THEN
               col(i) = mr_pbl*1.E-9 * aircol(i)
            ENDIF
         ENDDO
      ENDIF

*_______________________________________________________________________

      RETURN
      END
