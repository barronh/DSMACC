*=============================================================================*

      SUBROUTINE vpair(psurf, nz, z,
     $     airden, c)

*-----------------------------------------------------------------------------*
*=  NAME:  Vertial Profile of AIR
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of air molecules.  Subroutine includes a      =*
*=  shape-conserving scaling method that allows scaling of the entire        =*
*=  profile to a given sea-level pressure.                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  PSURF   - REAL, surface pressure (mb) to which profile should be      (I)=*
*=            scaled.  If PSURF < 0, no scaling is done                      =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  AIRDEN  - REAL, air density (molec/cc) at each specified altitude     (O)=* 
*=  C       - REAL, number of air molecules per cm^2 at each specified    (O)=*
*=            altitude layer                                                 =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=150)

* input:
      REAL z(kz)
      INTEGER nz
      REAL  psurf

* specified data:
      REAL zd(kdata), air(kdata)
      REAL hscale
      REAL cd(kdata)

* output:
* airden(iz) = air density (molec cm-3) at level iz
* cz(iz) =  column amount (molec cm-2) for layer iz      
      REAL airden(kz)
      REAL c(kz)

* local:
      REAL scale
      REAL pold
      REAL pconv
      PARAMETER(pconv = 980.665 * 1.E-3 * 28.9644 / 6.022169E23)

* other:
      REAL deltaz
      INTEGER i, nd

* External functions:
      REAL fsum
      EXTERNAL fsum

*_______________________________________________________________________

* _________SECTION 1:  Read in vertical profile of concentration

      WRITE(kout,*) 'air density: USSA, 1976'

      OPEN(kin,FILE='DATAE1/ATM/ussa.dens',STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      nd = 1
 4    CONTINUE
        READ(kin,*,END=5) zd(nd), air(nd)
        nd = nd+1
        GOTO 4
 5    CONTINUE
      CLOSE(kin)
      nd = nd-1

* compute column increments (logarithmic integrals)

      DO 6, i = 1, nd - 1
         deltaz = 1.E5 * (zd(i+1)-zd(i)) 
         cd(i) =  (air(i+1)-air(i)) /ALOG(air(i+1)/air(i)) * deltaz
C         cd(i) = (air(i+1)+air(i)) * deltaz / 2. 
    6 CONTINUE

* Include exponential tail integral from infinity to 50 km,
* fold tail integral into top layer
* specify scale height near top of data.

      hscale = 8.05e5
      cd(nd-1) = cd(nd-1) + hscale * air(nd)

********************** end data input.

* ________SECTION 2:  Compute column increments on standard z-grid.  
* Compute column increments on standard z-grid.  

      CALL inter3(nz,z,c, nd,zd,cd, 1)
      
* Compute air density at each level.  This is required by several other
* subroutines (e.g. pressure correction of quantum yields)

      CALL inter1(nz,z,airden,nd,zd,air)

* scale by input surface pressure:
* min value = 1 molec cm-2

      pold =  pconv * MAX(fsum(nz-1,c),1.)

      IF(psurf .GT. 0.) THEN
         scale = psurf/pold
      ELSE
         scale = 1.
      ENDIF

      DO i = 1, nz - 1
         c(i) = c(i) * scale
         airden(i) = airden(i) * scale
      ENDDO
      airden(nz) = airden(nz) * scale

*_______________________________________________________________________

      RETURN
      END
