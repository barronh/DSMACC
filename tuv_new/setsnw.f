      SUBROUTINE setsnw(nz,z,nw,wl,dtsnw,omsnw,gsnw)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set optical and physical properties for snowpack.                        =*
*=  Currently for wavelength-independent properties.                         =* 
*=  Subroutine outputs spectral quantities.                                  =*
*=  Lee-Taylor, J., and S. Madronich (2002), Calculation of actinic fluxes   =*
*=  with a coupled atmosphere-snow radiative transfer model, J. Geophys.     =*
*=  Res., 107(D24) 4796 (2002) doi:10.1029/2002JD002084                      =*
*-----------------------------------------------------------------------------*
*=  USER-DEFINED VARIABLES:                                                  =*
*=  zs      - height (km) of snow layer boundary above GROUND level          =*
*=  snwdens - density (g/cm3)                                                =*
*=  ksct    - mass-specific scattering coefficient (m2/kg)                   =*
*=  csoot   - soot content (ng Carbon / g snow)                              =*
*=  snow    - (=T/F) switch for presence of snow
*=                                                                           =*
*=  PARAMETERS:                                                              =*
*=  nz      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  z       - REAL, specified altitude working grid (km)                  (I)=*
*=  nw      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  wl      - REAL, vector of lower limits of wavelength intervals in     (I)=*
*=            working wavelength grid                                        =*
*=  dtsnw   - REAL, optical depth due to absorption by snow at each       (O)=*
*=            altitude and wavelength                                        =*
*=  omsnw   - REAL, single scattering albedo due to snow at each          (O)=*
*=            defined altitude and wavelength                                =*
*=  gsnw    - REAL, snow asymmetry factor at each defined altitude and    (O)=*
*=            wavelength                                                     =*
*=  rabs    - absorption coefficient of snow, wavelength-dependent           =*
*=  rsct    - scattering coefficient of snow, assume wavelength-independent  =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  10/00  adapted from setcld.f, Julia Lee-Taylor, ACD, NCAR                =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Jula Lee-Taylor, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  julial@ucar.edu                                           =*
*-----------------------------------------------------------------------------*
*= Copyright (C)2000 University Corporation for Atmospheric Research         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=51)

* input: (grids)
      REAL wl(kw)
      REAL z(kz)
      INTEGER nz
      INTEGER nw

* Output: 
      REAL dtsnw(kz,kw), omsnw(kz,kw), gsnw(kz,kw)

* local:

* specified data:
      REAL zs(kdata),dzs
      REAL cd(kdata), omd(kdata), gd
      REAL snwdens(kdata)             ! snwdens = snow density, g/cm3
      REAL csoot(kdata)               ! conc of elemental carbon, ng/g
      REAL r_ice(kw),rsoot
      REAL womd(kdata), wgd(kdata)
      REAL rsct(kdata),ksct(kdata),rabs(kdata)

* other:
      REAL cz(kz),omz(kz),gz(kz)
      INTEGER i,is,iw,iz,nsl

* External functions:
      REAL fsum
      EXTERNAL fsum
*--------------------------------------------------------------------------
* SNOW PROPERTIES: USER-DEFINED
*--------------------------------------------------------------------------
** define "number of snow layers + 1" (0 = no snow, 2 = single snow layer)
      nsl = 0

      IF(nsl.GE.2)THEN
** define snow grid, zs(ns), in km above GROUND level
* NOTE: to get good vertical resolution, subroutine gridz (in grids.f) should 
*       be modified to include small (1cm - 1mm) layers near snowpack top.

        zs(1) = 0.0
        zs(2) = 0.001

** define snow scattering coefficient, ksct, m2/kg snow
* melting midlatitude maritime (mountain) snow, ksct = 1-5 m2/kg_snow
* warmer polar coastal/maritime snow,           ksct = 6-13 m2/kg_snow
* cold dry polar/tundra snow,                   ksct = 20-30 m2/kg_snow
* Fisher, King and Lee-Taylor (2005), JGR 110(D21301) doi:10.1029/2005JD005963

        ksct(1) =  25.                        ! m2.kg-1 snow 

** define snow density, snwdens, g/cm3

        snwdens(1) = 0.4                      ! g/cm3

** define soot content, csoot, ng/g elemental carbon

        csoot(1) = 0.                         ! ng/g elemental carbon

*----------------------------------------------------------------------------
* SNOW PROPERTIES: FROM LITERATURE
*--------------------------------------------------------------------------
* read absorption coefficients 
        CALL rdice_acff(nw,wl,r_ice)          ! cm^-1 ice

* absorption due to soot, assume wavelength-independent
* rsoot ~ 10 m2/gC @500nm : Warren & Wiscombe, Nature 313,467-470 (1985)
        rsoot = 10.                           ! m2/gC 

* asymmetry factor : Wiscombe & Warren, J. Atmos. Sci, 37, 2712-2733 (1980)
        gd = 0.89
*----------------------------------------------------------------------------

* loop snow layers, assigning optical properties at each wavelength
        DO 17, iw = 1, nw-1
          DO 11 is = 1,nsl-1
            rsct(is)=ksct(is)*snwdens(is)*1.e+3        ! m-1 
            rsct(is)=rsct(is)*(zs(is+1)-zs(is))*1.e+3  ! no units

            rabs(is) = (r_ice(iw)/0.9177*1.e5 + rsoot*csoot(is)) 
     $             * snwdens(is)*(zs(is+1)-zs(is))   ! no units 
  
            cd(is) = rsct(is) + rabs(is)
            omd(is)= rsct(is) / cd(is) 
 
            if(iw.EQ.1)then
              print*,"Snowpack: is =",is,"; zs =",zs(is)
              PRINT*,"          ksct =", ksct(is)
              PRINT*,"          density =",snwdens(is)
              PRINT*,"          csoot =",csoot(is)
              PRINT*, 'cd = ',cd(is),'  omd = ',omd(is),'  gd = ',gd
              WRITE(kout,*)'snwdens = ',snwdens,' g/cm3'
              WRITE(kout,*)'ksct_snow = ',ksct(is),' m2.kg-1'
              WRITE(kout,*)'soot = ',csoot(is),' ng/g' 
              WRITE(kout,*)'cd = ',cd(is),'omd = ',omd(is),'gd = ',gd
            endif

* compute integrals and averages over snow layers:
* for g and omega, use averages weighted by optical depth
            womd(is) = omd(is) * cd(is)
            wgd(is) = gd * cd(is)
   11     CONTINUE

* interpolate snow layers onto TUV altitude grid (gridz)
          CALL inter3(nz,z,cz, nsl,zs,cd, 0)
          CALL inter3(nz,z,omz,nsl,zs,womd, 0)
          CALL inter3(nz,z,gz ,nsl,zs,wgd, 0)

          DO 15, iz = 1, nz-1
            IF (cz(iz) .GT. 0.) THEN
              omz(iz) = omz(iz)/cz(iz)
              gz(iz)  = gz(iz) /cz(iz)
            ELSE
              omz(iz) = 0.
              gz(iz) = 0.
            ENDIF
            dtsnw(iz,iw) = cz(iz)
            omsnw(iz,iw) = omz(iz)
            gsnw(iz,iw)  = gz(iz)
   15     CONTINUE
   17   CONTINUE

        PRINT*,"Snowpack top: zs =",zs(nsl)

      ELSE ! no snow
        DO 16, iz = 1, nz-1
          cz(iz) = 0.
          omz(iz) = 1.
          gz(iz) = 0.
          DO 18, iw = 1, nw-1
            dtsnw(iz,iw) = cz(iz)
            omsnw(iz,iw) = omz(iz)
            gsnw(iz,iw)  = gz(iz)
   18     CONTINUE
   16   CONTINUE
      ENDIF ! snow exists

      RETURN
      END
*******************************************************************************
      SUBROUTINE rdice_acff(nw,wl,rabs)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read ice absorption coefficient.  Re-grid data to match                  =*
*=  specified wavelength working grid.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  RABS_ice - REAL, absorption coefficient (cm^-1) of ice at             (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  10/00  Created routine by editing rdh2oxs.                               =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input: (altitude working grid)
      INTEGER nw
      REAL wl(kw)

* output:

      REAL rabs(kw)

* local:
      REAL x1(kdata)
      REAL y1(kdata),y2(kdata),y3(kdata)
      REAL yg(kw)
      REAL a1, a2, dum
      INTEGER ierr
      INTEGER i,l,m, n, idum
      CHARACTER*40 fil
*_______________________________________________________________________

************* absorption cross sections:
* ice absorption cross sections from 

      fil = 'DATA/ice'
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/ICE_Perov.acff',STATUS='old')
      m = 17       ! header lines
      n = 79       ! data lines
      !OPEN(UNIT=kin,FILE='DATAJ1/ABS/ICE_min.acff',STATUS='old')
      !m = 13       ! header lines
      !n = 52       ! data lines

      DO 11, i = 1,m
         read(kin,*)
   11 CONTINUE
      DO 12, i = 1, n
         READ(kin,*) x1(i), y1(i)
   12 CONTINUE
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF
      
      DO 13, l = 1, nw-1
         rabs(l) = yg(l)
   13 CONTINUE

*_______________________________________________________________________

      RETURN
      END
