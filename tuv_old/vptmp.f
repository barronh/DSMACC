*=============================================================================*

      SUBROUTINE vptmp(nz,z,tlev,tlay)

*-----------------------------------------------------------------------------*
*   NAME: Vertical Profile of TeMPerature
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of temperatures.  Temperature values are      =*
*=  needed to compute some cross sections and quantum yields.  Distinguish   =*
*=  between temperature at levels and layers.                                =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  TLEV    - REAL, temperature (K) at each specified altitude level      (O)=*
*=  TLAY    - REAL, temperature (K) at each specified altitude layer      (O)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=150)

* input: (altitude working grid)
      REAL z(kz)
      INTEGER nz

* output:
      REAL tlev(kz), tlay(kz)

* local:
      REAL zd(kdata), td(kdata)
      INTEGER i, nd
*_______________________________________________________________________


* read in temperature profile

      WRITE(kout,*) 'air temperature: USSA, 1976'

      OPEN(kin,FILE='DATAE1/ATM/ussa.temp',STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO
      nd = 1
 4    CONTINUE
         READ(kin,*,END=5) zd(nd), td(nd) 
         nd = nd+1
         GOTO 4
 5    CONTINUE
      CLOSE(kin)
      nd = nd-1

* use constant temperature to infinity:  

      zd(nd) = 1.E10

* alternative input temperature data could include, e.g., a read file here:

***********
*********** end data input.

* interpolate onto z-grid

      CALL inter1(nz,z,tlev,nd,zd,td)

* compute layer-averages

      DO 20, i = 1, nz - 1
         tlay(i) = (tlev(i+1) + tlev(i))/2.
 20   CONTINUE
*_______________________________________________________________________
      
      RETURN
      END
