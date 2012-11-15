*=============================================================================*

      SUBROUTINE vpo3(colnew,nz,z, cair, c)

*-----------------------------------------------------------------------------*
*=  NAME:  Vertical Profiles of Ozone = vpo3
*=  PURPOSE:                                                                 =*
*=  Computes column increments, c(iz), molec cm-2 for each layer iz.         =*
*=  Normally, c(iz) values are computed from input vertical profiles of      =*
*=  concentrations (molec cm-3), that are then interpolated and integrated   =*
*=  over each layer.  The default example here uses the US Standard          =*
*=  Atmosphere (1976) mid-latitude concentration profile.                    =*
*=  Users can substitute their own profile, as long as appropriate are made  =*
*=  adjustments to the code in Section 1.                                    =*
*=  A scale factor is provided to allow changing the total column amount,
*=  but conserving the shape of the profile.
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  COLNEW - REAL, new total column amount from surface to space             =*
*=           should be scaled.  If O3NEW < 0, no scaling is done             =*
*=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
*=           grid                                                            =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  CAIR(KZ) = REAL, air column increment (molec cm-2), provided here in     =*
*=  case it is necessary to convert from mixing ratio units (e.g. ppb).      =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

********
* inputs:
********

      INCLUDE 'params'

* from calling program:

      INTEGER nz
      REAL z(kz)
      REAL colnew
      REAL cair(kz)

* from data file:  concentration data as a function of altitude

      INTEGER kdata
      PARAMETER(kdata=150)
      REAL zd(kdata), xd(kdata)

********
* internal
********

      INTEGER i, nd
      REAL hscale
      REAL cd(kdata)
      REAL colold, dobold, scale
      INTEGER iz

********
* output:
********

      REAL c(kz)

********
* External functions:
********
      REAL fsum
      EXTERNAL fsum


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

* compute column increments

      DO 11, i = 1, nd - 1
         cd(i) = (xd(i+1)+xd(i)) * 1.E5 * (zd(i+1)-zd(i)) / 2. 
   11 CONTINUE

* Include exponential tail integral from infinity to 50 km,
* fold tail integral into top layer
* specify scale height near top of data.

      hscale = 4.50e5
      cd(nd-1) = cd(nd-1) + hscale * xd(nd)

* alternative input ozone concentration data could include, e.g., 
* a read file here:

*********** end data input.

* ________SECTION 2:  Compute column increments on standard z-grid.  

      CALL inter3(nz,z,c, nd,zd,cd, 1)

* If colnew is not negative, scale to new total column value, colnew :

      IF (colnew .GT. nzero) THEN

* For ozone, colnew is in Dobson Units (1 DU = 2.687e16).

         colold = fsum(nz-1, c)/2.687e16

* If old column is near zero (less than 1 molec cm-2), 
* use constant mixing ratio profile (nominal 1 ppt before scaling) 
* to avoid numerical problems when scaling.


         IF(fsum(nz-1,c) .LT. 1.) THEN
            DO iz = 1, nz-1
               c(iz) = 1.E-12 * cair(iz)
            ENDDO
         ENDIF
         scale = colnew/colold
         DO iz = 1, nz-1
            c(iz) = c(iz) * scale
         ENDDO

      ENDIF




*_______________________________________________________________________

      RETURN
      END
