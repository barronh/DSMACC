*=============================================================================*

      SUBROUTINE vpair(psurf, nz, z,
     $     con, col)

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
*= outputs are on z-grid:
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  CON     - REAL, air density (molec/cc) at each specified altitude     (O)=* 
*=  COL     - REAL, number of air molecules per cm^2 in each specified    (O)=*
*=            altitude layer (column vertical increment                      =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=150)

********* input:
    
      REAL z(kz)
      INTEGER nz

      REAL  psurf

* specified air profile data:
      INTEGER nd
      REAL zd(kdata), air(kdata)
      REAL hscale
      REAL cd(kdata)

********* output:
* con(iz) = air density (molec cm-3) at level iz
* col(iz) =  column amount (molec cm-2) for layer iz      
      REAL con(kz)
      REAL col(kz)

* local:
      REAL scale
      REAL pold
      REAL pconv
      PARAMETER(pconv = 980.665 * 1.E-3 * 28.9644 / 6.022169E23)

* other:
      INTEGER i
      REAL airlog(kz), conlog(kz)


* External functions:
      REAL fsum
      EXTERNAL fsum

*_______________________________________________________________________

* The objective of this subroutine is to take the input air profile and interpolate
*   it to the working grid, z(i = 1, nz).  The desired outputs are con(i) and col(i).
* Input vertical profiles can be specified in various ways. For example: 
*   altitude vs. concentration (molecules cm-3)
*   altitude vs. pressure (mbar) and temperature (K)
*   altitude vs. column increments (molec cm-2)
* The interpolation scheme will depend on the specific type of input data.
* Here, the US Standard Atmosphere is given as altitude vs. concentration (also called
*   number density)

* _________SECTION 1:  Read in vertical profile of concentration

      WRITE(kout,*) 'air concentrations: USSA, 1976'

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
* add 1 meter to top, to avoid interpolation end-problem if z-grid is up to 120 km
      zd(nd) = zd(nd) + 0.001

* scale height, km, at top of atmosphere:
      hscale = 8.01

********************** end data input.

* ________SECTION 2:  Compute column increments on standard z-grid.  
* For air, this is best done using logarithms of concentration.
*   take logs
*   if z-grid extends beyond available data, stop (no extrapolation allowed)
*   interpolate log of air(nd) onto z grid 
*   re-exponentiate to get gridded concentrations

      DO i = 1, nd
         airlog(i) = ALOG(air(i))
      ENDDO

      IF(z(nz) .GT. zd(nd)) STOP 'in vpair: ztop < zdata'
      CALL inter1(nz,z,conlog, nd,zd,airlog)

      DO i = 1, nz
         con(i) = EXP(conlog(i))
      ENDDO

* Find gridded column increments in z-grid:
*   use log intergration

      DO i = 1, nz-1
         col(i) = 1.E5*(z(i+1)-z(i)) * (con(i+1)-con(i)) /
     $        ALOG(con(i+1)/con(i))
      ENDDO

* Add exponential tail integral at top of atmosphere:
*   this is folded into layer nz-1, making this layer "heavy'.  
*   The layer nz is not used. The radiative transfer 
*   calculation is based on nz-1 layers (not nz).

      col(nz-1) = col(nz-1) + 1.E5 * hscale * con(nz)
  
* Scale by input surface pressure:
* min value = 1 molec cm-2

      pold =  pconv * MAX(fsum(nz-1,col),1.)

      IF(psurf .GT. 0.) THEN
         scale = psurf/pold
      ELSE
         scale = 1.
      ENDIF

      DO i = 1, nz - 1
         col(i) = col(i) * scale
         con(i) = con(i) * scale
      ENDDO
      con(nz) = con(nz) * scale

*_______________________________________________________________________

      RETURN
      END
