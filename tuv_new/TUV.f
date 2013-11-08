      subroutine tuv(iboxyear, iboxmonth, iboxday, ro3col, ralt, rlat,
     &  rlon, box_temp, raird, sza, svj_tj, jlabel)
*-----------------------------------------------------------------------------*
*=    Tropospheric Ultraviolet-Visible (TUV) radiation model                 =*
*=    Version 4.5                                                            =*
*=    Sep 2007                                                               =*
*-----------------------------------------------------------------------------*
*= Developed by Sasha Madronich with important contributions from:           =*
*= Chris Fischer, Siri Flocke, Julia Lee-Taylor, Bernhard Meyer,             =*
*= Irina Petropavlovskikh,  Xuexi Tie, and Jun Zen.                          =*
*= Special thanks to Knut Stamnes and co-workers for the development of the  =*
*= Discrete Ordinates code, and to Warren Wiscombe and co-workers for the    =*
*= development of the solar zenith angle subroutine. Citations for the many  =*
*= data bases (e.g. extraterrestrial irradiances, molecular spectra) may be  =*
*= found in the data files headers and/or in the subroutines that read them. =*
*=              To contact the author, write to:                             =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu  or tuv@acd.ucar.edu                       =*
*-----------------------------------------------------------------------------*
*= This program is free software; you can redistribute it and/or modify      =*
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
*= Copyright (C) 1994,95,96,97,98,99,2000,01,02,03, 04, 05, 06, 07           =*
*= by the University Corporation for Atmospheric Research                    =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* Include parameter file

      INCLUDE 'params'

* Added to use TUV.f as a function instead of a program
* BHH: for use with the DSMACC model
      REAL ro3col, ralt, raird, box_temp, rlat, rlon
      integer iboxyear, iboxmonth, iboxday
      integer i,j

* Wavelength grid:

      INTEGER nw, iw, nwint
      REAL wl(kw), wc(kw), wu(kw)
      REAL wstart, wstop

* Altitude grid

      INTEGER nz, iz, izout
      REAL z(kz), zstart, zstop, zout

* Solar zenith angle and azimuth
* slant pathlengths in spherical geometry

      REAL sza(kt)
      REAL zen
      INTEGER nid(0:kz)
      REAL dsdh(0:kz,kz)

* Extra terrestrial solar flux and earth-Sun distance ^-2

      REAL f(kw), etf(kw)
      REAL esfact(kt)

* Ozone absorption cross section

      REAL o3xs(kz,kw)

* O2 absorption cross section

      REAL o2xs(kz,kw), o2xs1(kw)

* SO2 absorption cross section
     
      REAL so2xs(kw)

* NO2 absorption cross section
     
      REAL no2xs(kw)

* Atmospheric optical parameters

      REAL tlev(kz), tlay(kz)
      REAL aircon(kz), aircol(kz), vcol(kz), scol(kz)
      REAL dtrl(kz,kw)
      REAL co3(kz)
      REAL dto3(kz,kw), dto2(kz,kw), dtso2(kz,kw), dtno2(kz,kw)
      REAL dtcld(kz,kw), omcld(kz,kw), gcld(kz,kw)
      REAL dtaer(kz,kw), omaer(kz,kw), gaer(kz,kw)
      REAL dtsnw(kz,kw), omsnw(kz,kw), gsnw(kz,kw)
      REAL albedo(kw)

* Spectral irradiance and actinic flux (scalar irradiance)

      REAL edir(kz), edn(kz), eup(kz)
      REAL sirrad(kz,kw)
      REAL fdir(kz), fdn(kz), fup(kz)
      REAL saflux(kz,kw)

* Spectral weighting functions and weighted radiation

      INTEGER ns, is
      REAL sw(ks,kw), rate(ks,kz), dose(ks)
      REAL drdw
      CHARACTER*50 slabel(ks)

* Photolysis coefficients (j-values)

      INTEGER nj, ij
      REAL sj(kj,kz,kw), valj(kj,kz)
      REAL djdw
      CHARACTER*50 jlabel(kj)

**** Re-scaling factors (can be read from input file)
* New surface albedo and surface pressure (milli bar)
* Total columns of O3, SO2, NO2 (Dobson Units)
* Cloud optical depth, altitude of base and top
* Aerosol optical depth at 550 nm, single scattering albedo, Angstrom alpha

      REAL alsurf, psurf
      REAL o3_tc, so2_tc, no2_tc
      REAL taucld, zbase, ztop
      REAL tauaer, ssaaer, alpha

* Location: Lat and Lon (deg.), surface elev (km)
* Altitude, temperature and pressure for specific outputs

      REAL lat, lon
      REAL zaird, ztemp

* Time and/or solar zenith angle
      
      INTEGER iyear, imonth, iday
      INTEGER it, nt
      REAL t(kt), tstart, tstop
      REAL tmzone, ut
      LOGICAL lzenit

* number of radiation streams

      INTEGER nstr

* input/output control

      LOGICAL intrct
      CHARACTER*6 inpfil, outfil

      INTEGER iout

      REAL dirsun, difdn, difup

      CHARACTER*1 again

* Save arrays for output:

      LOGICAL lirrad, laflux, lrates, ljvals, lmmech
      INTEGER isfix, ijfix, itfix, izfix, iwfix
      INTEGER nms, ims(ks), nmj, imj(kj)

      REAL svj_zj(kz,kj), svj_tj(kt,kj), svj_zt(kz,kt)
      REAL svr_zs(kz,ks), svr_ts(kt,ks), svr_zt(kz,kt)
      REAL svf_zw(kz,kw), svf_tw(kt,kw), svf_zt(kz,kt)
      REAL svi_zw(kz,kw), svi_tw(kt,kw), svi_zt(kz,kt)

* Parameters for shifting wavelengths air <--> vacuum

      INTEGER mrefr
      LOGICAL lrefr
      REAL airout

* Planetary boundary layer height and pollutant concentrations

      INTEGER ipbl
      REAL zpbl
      REAL o3pbl, so2pbl, no2pbl, aod330

* Other user-defined variables here:

      real xtrl, xto2, xto3, xtso2, xtno2, xtcld, xtaer, xtsnw



* --- END OF DECLARATIONS ---------------------------------------------
! BHH: Noting which TUV variables are internal
      write(*,*) "Most inputs come from the TUV usrinp file"
      write(*,*) "*Using model values for following inputs:"
      write(*,*) " - o3_tc"
      write(*,*) " - zout"
      write(*,*) " - lat"
      write(*,*) " - lon"
      write(*,*) " - ztemp"
      write(*,*) " - zaird"
      write(*,*) " - tstart = 0"
      write(*,*) " - zstop = 90"
      write(*,*) " - nt = 19"
      write(*,*) " - lzenith = .TRUE."

* re-entry point

 1000 CONTINUE

* Open log file:

c      OPEN(UNIT=kout,FILE='tuvlog',STATUS='UNKNOWN')
      OPEN(UNIT=kout,FILE='../'//'tuvlog'//'.txt',STATUS='UNKNOWN')

* ___ SECTION 1: SIMPLE INPUT VARIABLES --------------------------------
******* Read simple input variables from a file:

* can read interactively (intrct = .TRUE.) 
* or in batch mode (intrct = .FALSE.)

c      intrct = .TRUE.
      intrct = .FALSE.
      IF ( .NOT. intrct) inpfil = 'usrinp'

      CALL rdinp(intrct, 
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3_tc,  so2_tc, no2_tc,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ims,    slabel, imj,    jlabel)

      IF(outfil .EQ. 'screen') THEN
         iout = 6
      ELSE
         iout = 30
      ENDIF         


************* Can overwrite basic inputs here manually:
* Input and output files:
*   inpfil = input file name
*   outfil = output file name
* Radiative transfer scheme:
*   nstr = number of streams
*          If nstr < 2, will use 2-stream Delta Eddington
*          If nstr > 1, will use nstr-stream discrete ordinates
* Location (geographic):
*   lat = LATITUDE (degrees, North = positive)
*   lon = LONGITUDE (degrees, East = positive)
*   tmzone = Local time zone difference (hrs) from Universal Time (ut):  
*            ut = timloc - tmzone
* Date:
*   iyear = year (1950 to 2050)
*   imonth = month (1 to 12)
*   iday = day of month
* Time of day grid:
*   tstart = starting time, local hours
*   tstop = stopping time, local hours
*   nt = number of time steps
*   lzenit = switch for solar zenith angle (sza) grid rather than time 
*             grid. If lzenit = .TRUE. then 
*                tstart = first sza in deg., 
*                tstop = last sza in deg., 
*                nt = number of sza steps. 
*                esfact = 1. (Earth-sun distance = 1.000 AU)
* Vertical grid:
*   zstart = surface elevation above sea level, km
*   zstop = top of the atmosphere (exospheric), km
*   nz = number of vertical levels, equally spaced
*        (nz will increase by +1 if zout does not match altitude grid)
* Wavlength grid:
*   wstart = starting wavelength, nm
*   wstop  = final wavelength, nm
*   nwint = number of wavelength intervals, equally spaced
*           if nwint < 0, the standard atmospheric wavelength grid, not
*           equally spaced, from 120 to 735 nm, will be used. In this
*           case, wstart and wstop values are ignored.
* Surface condition:
*   alsurf = surface albedo, wavelength independent
*   psurf = surface pressure, mbar.  Set to negative value to use
*           US Standard Atmosphere, 1976 (USSA76)
* Column amounts of absorbers (in Dobson Units, from surface to space):
*          Vertical profile for O3 from USSA76.  For SO2 and NO2, vertical
*          concentration profile is 2.69e10 molec cm-3 between 0 and 
*          1 km above sea level, very small residual (10/largest) above 1 km.
*   o3_tc = ozone (O3)
*   so2_tc = sulfur dioxide (SO2)
*   no2_tc = nitrogen dioxide (NO2)
* Cloud, assumed horizontally uniform, total coverage, single scattering
*         albedo = 0.9999, asymmetry factor = 0.85, indep. of wavelength,
*         and also uniform vertically between zbase and ztop:
*   taucld = vertical optical depth, independent of wavelength
*   zbase = altitude of base, km above sea level
*   ztop = altitude of top, km above sea level
* Aerosols, assumed vertical provile typical of continental regions from
*         Elterman (1968):
*   tauaer = aerosol vertical optical depth at 550 nm, from surface to space. 
*           If negative, will default to Elterman's values (ca. 0.235 
*           at 550 nm).
*   ssaaer = single scattering albedo of aerosols, wavelength-independent.
*   alpha = Angstrom coefficient = exponent for wavelength dependence of 
*           tauaer, so that  tauaer1/tauaer2  = (w2/w1)**alpha.
* Directional components of radiation, weighting factors:
*   dirsun = direct sun
*   difdn = down-welling diffuse
*   difup = up-welling diffuse
*        e.g. use:
*        dirsun = difdn = 1.0, difup = 0 for total down-welling irradiance
*        dirsun = difdn = difup = 1.0 for actinic flux from all directions
*        dirsun = difdn = 1.0, difup = -1 for net irradiance
* Output altitude:
*   zout = altitude, km, for desired output.
*        If not within 1 m of altitude grid, an additional
*        level will be inserted and nz will be increased by +1.
*   zaird = air density (molec. cm-3) at zout.  Set to negative value for
*        default USSA76 value interpolated to zout.
*   ztemp = air temperature (K) at zout.  Set to negative value for
*        default USSA76 value interpolated to zout.
* Output options, logical switches:
*   lirrad = output spectral irradiance
*   laflux = output spectral actinic flux
*   lmmech = output for NCAR Master Mechanism use
*   lrates = output dose rates (UVB, UVA, CIE/erythema, etc.)
* Output options, integer selections:
*   isfix:  if > 0, output dose rate for action spectrum is=isfix, tabulated
*           for different times and altitudes.
*   ijfix:  if > 0, output j-values for reaction ij=ijfix, tabulated
*           for different times and altitudes.
*   iwfix:  if > 0, output spectral irradiance and/or spectral actinic
*           flux at wavelength iw=iwfix, tabulated for different times
*           and altitudes.
*   itfix:  if > 0, output spectral irradiance and/or spectral actinic
*           flux at time it=itfix, tabulated for different altitudes
*           and wavelengths.
*   izfix:  if > 0, output spectral irradiance and/or spectral actinic
*           flux at altitude iz=izfix, tabulated for different times
*           and wavelengths.
*   nms:    number of dose rates that will be reported. Selections must be 
*           made interactively, or by editing input file.
*   nmj:    number of j-values that will be reported. Selections must be 
*           made interactively, or by editing input file.
* The following default settings are also found in the input file 'defin1':

c      inpfil = defin1
c      outfil = usrout
c      nstr = -2
      lat = rlat
      lon = rlon
c      tmzone = 0.
      iyear = iboxyear
      imonth = iboxmonth
      iday = iboxday
c      zstart = max(0., ralt - .1)
c      zstop = ralt + .1
c      nz = 3
c      wstart = 280.
c      wstop = 420.
c      nwint = 140
      tstart = 0
      tstop = 90.
      nt = 19
      lzenit = .TRUE.
c      alsurf = 0.1
c      psurf = -999.
       o3_tc = ro3col
c      so2_tc = 0.
c      no2_tc = 0.
c      tcloud = 0.
c      zbase = 4.
c      ztop = 5.
c      tauaer = 0.235
c      ssaaer = 0.99
c      alpha = 1.
c      dirsun = 1.
c      difdn = 1.
c      difup = 0.
       zout = ralt
       zaird = raird
       ztemp = box_temp

c      lirrad = .TRUE.
c      laflux = .FALSE.
c      lmmech = .FALSE.
c      lrates = .TRUE.
c      isfix = 0
*      nms cannot be set here
c      ljvals = .FALSE.
c      ijfix = 0
*      nmj cannot be set here
c      iwfix = 0
c      itfix = 0
c      izfix = 0

      IF(nstr .LT. 2) THEN
         WRITE(kout,*) 'Delta-Eddington 2-stream radiative transfer' 
      ELSE
         WRITE(kout,*) 'Discrete ordinates ', 
     $        nstr, '-stream radiative transfer' 
      ENDIF

      WRITE(*,*) 'calculating....'

* ___ SECTION 2: SET GRIDS _________________________________________________

* wavelengths (creates wavelength grid: lower, center, upper of each bin)
* NOTE:  Wavelengths are in vacuum.  To use wavelengths in air, see 
* Section 3 below, where you must set lrefr= .TRUE.  

      CALL gridw(wstart, wstop, nwint,
     $     nw, wl, wc, wu)

* altitudes (creates altitude grid, locates index for selected output, izout)

      CALL gridz(zstart, zstop, nz, z, zout, izout)
      if(izfix .gt. 0) izout = izfix

* time/zenith (creates time/zenith angle grid, starting at tstart)
      CALL gridt(lat, lon, tmzone,
     $     iyear, imonth, iday,
     $     lzenit, tstart, tstop,
     $     nt, t, sza, esfact)

* ___ SECTION 3: SET UP VERTICAL PROFILES OF TEMPERATURE, AIR DENSITY, and OZONE______

***** Temperature vertical profile, Kelvin 
*   can overwrite temperature at altitude z(izout)

      CALL vptmp(nz,z, tlev,tlay)

* Enabling box_temp to be used after a call to tuv
* BHH: For use as a diagnostic tool
      IF(ztemp .GT. nzero) then
        tlev(izout) = ztemp
      else
        box_temp = tlev(izout)
      endif

*****  Air density (molec cm-3) vertical profile 
*   can overwrite air density at altitude z(izout)

      CALL vpair(psurf, nz, z,
     $     aircon, aircol)
      IF(zaird .GT. nzero) aircon(izout) = zaird

***** Correction for air-vacuum wavelength shift:
* The TUV code assumes that all working wavelengths are strictly IN-VACUUM. This is assumed for ALL
* spectral data including extraterrestrial fluxes, ozone (and other) absorption cross sections,
* and various weighting functons (action spectra, photolysis cross sections, instrument spectral
* response functions).  

*  Occasionally, users may want their results to be given for wavelengths measured IN-AIR.
*   The shift between IN-VACUUM and IN-AIR wavelengths depends on the index of refraction
*   of air, which in turn depends on the local density of air, which in turn depends on
*   altitude, temperature, etc.
*  Here, we provide users with the option to use a wavelength grid IN-AIR, at the air density
*   corresponding to the output altitude, airout = aircon(izout), by setting the logical variable
*   lrefr = .TRUE.  (default is lrefr = .FALSE.).  The wavelengths specified in gridw.f will be assumed
*   to be IN-AIR, and will be shifted here to IN-VACUUM values to carry out the calculatons.  
*   The actual radiative transfer calculations will be done strictly with IN-VACUUM values.  
*   If this shift is applied (i.e., if lrefr = .TRUE.), the wavelength grid will be shifted back to air
*   values just before the output is written.
*  Note:  if this option is used (lref = .TRUE.), the wavelength values will be correct ONLY at the 
*   selected altitude, iz = iout.  The wavelength shift will be INCORRECT at all other altitudes.
*  Note:  This option cannot be changed interactively in the input table.  It must be changed here.

C      lrefr = .TRUE.
      lrefr = .FALSE.
      IF(lrefr) THEN
         airout = aircon(izout)
         write(*,*) 'applying air to vacuum wavelength shift', 
     $        izout, airout
         mrefr = 1
         CALL wshift(mrefr, nw,    wl, airout)
         CALL wshift(mrefr, nwint, wc, airout)
         CALL wshift(mrefr, nwint, wu, airout)
      ENDIF

*************************
*! PBL pollutants will be added if zpbl > 0.
* CAUTIONS:  
* 1. The top of the PBL, zpbl in km, should be on one of the z-grid altitudes.
* 2. Concentrations, column increments, and optical depths
*       will be overwritten between surface and zpbl.
* 3. Inserting PBL constituents may change their total column amount.
* 4. Above pbl, the following are used:
*       for O3:  USSA or other profile
*       for NO2 and SO2: set to zero.
*       for aerosols: Elterman
* Turning on pbl will affect subroutines:
* vpo3, setno2, setso2, and setaer. See there for details

      zpbl = -999.
C      zpbl = 3.

* locate z-index for top of pbl

      ipbl = 0
      IF(zpbl. GT. 0.) THEN
         DO iz = 1, nz-1
            IF(z(iz+1) .GT. z(1) + zpbl*1.00001) GO TO 19
         ENDDO
 19      CONTINUE
         ipbl = iz - 1
         write(*,*) 'top of PBL index, height (km) ', ipbl, z(ipbl)

* specify pbl concetrations, in parts per billion

         o3pbl = 100.
         so2pbl = 10.
         no2pbl = 50.

* PBL aerosol optical depth at 330 nm
* (to change ssa and g of pbl aerosols, go to subroutine setair.f)

         aod330 = 0.8

      ENDIF

***** Ozone vertical profile

      CALL vpo3(ipbl, zpbl, o3pbl,
     &       o3_tc, nz, z, aircol, co3)

* ___ SECTION 3: READ SPECTRAL DATA ____________________________

* read (and grid) extra terrestrial flux data:
      
      CALL rdetfl(nw,wl, f)

* read cross section data for 
*    O2 (will overwrite at Lyman-alpha and SRB wavelengths
*            see subroutine la_srb.f)
*    O3 (temperature-dependent)
*    SO2 
*    NO2

      CALL rdo2xs(nw,wl, o2xs1)
      CALL rdo3xs(nw,wl,nz,tlay, o3xs)
      CALL rdso2xs(nw,wl, so2xs)
      CALL rdno2xs(nw,wl, no2xs)


****** Spectral weighting functions 
* (Some of these depend on temperature T and pressure P, and therefore
*  on altitude z.  Therefore they are computed only after the T and P profiles
*  are set above with subroutines settmp and setair.)
* Photo-physical   set in swphys.f (transmission functions)
* Photo-biological set in swbiol.f (action spectra)
* Photo-chemical   set in swchem.f (cross sections x quantum yields)
* Physical and biological weigthing functions are assumed to depend
*   only on wavelength.
* Chemical weighting functions (product of cross-section x quantum yield)
*   for many photolysis reactions are known to depend on temperature
*   and/or pressure, and therefore are functions of wavelength and altitude.
* Output:
* from pphys & pbiol:  s(ks,kw) - for each weighting function slabel(ks)
* from pchem:  sj(kj,kz,kw) - for each reaction jlabel(kj)
* For pchem, need to know temperature and pressure profiles.

      CALL swphys(nw,wl,wc, ns,sw,slabel)
      CALL swbiol(nw,wl,wc, ns,sw,slabel)
      CALL swchem(nw,wl,nz,tlev,aircon,nj,sj,jlabel)

c      CALL swbiol2(nw,wl,wc, ns,sw,slabel)
c      CALL swbiol3(nw,wl,wc, ns,sw,slabel)

* The following lines serve to print a list of the spectral weighting
* functions.  If new functions (e.g. action spectra, photo-reactions)
* are added, this list should be used to replace the list in the
* default input files (defin1, defin2, etc.).  The true/false toggle
* will be set to F, and should be changed manually to select weighting
* functions for output. Note that if many more functions are added, it
* may be necessary to increase the parameters ks and kj in the include
* file 'params'
* The program will stop after writing this list.
* Comment out these lines when not generating a new list.

c       OPEN(UNIT=50,FILE='spectra.list',STATUS='NEW')
c       WRITE(50,500)
c  500  FORMAT(5('='),1X,'Available spectral weighting functions:')
c       DO is = 1, ns
c          WRITE(50,505) is, slabel(is)
c       ENDDO
c       WRITE(50,510)
c  510  FORMAT(5('='),1X,'Available photolysis reactions')
c       DO ij = 1, nj
c          WRITE(50,505) ij, jlabel(ij)
c       ENDDO
c  505  FORMAT('F',I3,1X,A50)
c       WRITE(50,520)
c  520  FORMAT(66('='))
c       CLOSE (50)
c       STOP

      DO ij = 1, nj
         WRITE(kout,505) ij, jlabel(ij)
      ENDDO
 505  FORMAT('F',I3,1X,A50)

******

* ___ SECTION 4: SET ATMOSPHERIC OPTICAL DEPTH INCREMENTS _____________________

* Rayleigh optical depth increments:

      CALL odrl(nz, z, nw, wl, aircol, dtrl)
      
* O2 vertical profile and O2 absorption optical depths
*   For now, O2 densitiy assumed as 20.95% of air density, can be changed
*   in subroutine.
*   Optical depths in Lyman-alpha and SRB will be over-written
*   in subroutine la_srb.f

      CALL seto2(nz,z,nw,wl,aircol,o2xs1, dto2)

* Ozone optical depths

      CALL odo3(nz,z,nw,wl,o3xs,co3, dto3)

* SO2 vertical profile and optical depths

      CALL setso2(ipbl, zpbl, so2pbl,
     $     so2_tc, nz, z, nw, wl, so2xs,
     $     tlay, aircol,
     $     dtso2)

* NO2 vertical profile and optical depths

      CALL setno2(ipbl, zpbl, no2pbl, 
     $     no2_tc, nz, z, nw, wl, no2xs,
     $     tlay, aircol,
     $     dtno2)

* Cloud vertical profile, optical depths, single scattering albedo, asymmetry factor

      CALL setcld(taucld,zbase,ztop,
     $     nz,z,nw,wl,
     $     dtcld,omcld,gcld)

* Aerosol vertical profile, optical depths, single scattering albedo, asymmetry factor

      CALL setaer(ipbl, zpbl, aod330,
     $     tauaer, ssaaer, alpha,
     $     nz, z, nw, wl,
     $     dtaer, omaer, gaer)

* Snowpack physical and optical depths, single scattering albedo, asymmetry factor

      CALL setsnw(
     $     nz,z,nw,wl,
     $     dtsnw,omsnw,gsnw)

* Surface albedo

      CALL setalb(alsurf,nw,wl,
     $     albedo)

* ___ SECTION 5: TIME/SZA LOOP  _____________________________________

* Initialize any time-integrated quantities here

      CALL zero1(dose,ks)

* Loop over time or solar zenith angle (zen):

      DO 20, it = 1, nt

         zen = sza(it)

!         WRITE(*,200) it, zen, esfact(it)
!         WRITE(kout,200) it, zen, esfact(it)
 200     FORMAT('step = ', I4,' sza = ', F9.3, 
     $        ' Earth-sun factor = ', F10.7)

* correction for earth-sun distance

         DO iw = 1, nw - 1
            etf(iw) = f(iw) * esfact(it)
         ENDDO

* ____ SECTION 6: CALCULATE ZENITH ANGLE-DEPENDENT QUANTITIES __________

* slant path lengths for spherical geometry

         CALL sphers(nz,z,zen, dsdh,nid)
         CALL airmas(nz, dsdh,nid, aircol,vcol,scol)

* Recalculate effective O2 optical depth and cross sections for Lyman-alpha
* and Schumann-Runge bands, must know zenith angle
* Then assign O2 cross section to sj(1,*,*)

         CALL la_srb(nz,z,tlev,nw,wl,vcol,scol,o2xs1,
     $        dto2,o2xs)
         CALL sjo2(nz,nw,o2xs,1, sj)

* ____ SECTION 7: WAVELENGTH LOOP ______________________________________

* initialize for wavelength integration

         CALL zero2(rate,ks,kz)
         CALL zero2(valj,kj,kz)

***** Main wavelength loop:

         DO 10, iw = 1, nw-1

** monochromatic radiative transfer. Outputs are:
*  normalized irradiances     edir(iz), edn(iz), eup(iz) 
*  normalized actinic fluxes  fdir(iz), fdn(zi), fup(iz)
*  where 
*  dir = direct beam, dn = down-welling diffuse, up = up-welling diffuse

            CALL rtlink(nstr, nz,
     $           iw, albedo(iw), zen,
     $           dsdh,nid,
     $           dtrl,
     $           dto3,
     $           dto2,
     $           dtso2,
     $           dtno2,
     $           dtcld, omcld, gcld,
     $           dtaer,omaer,gaer,
     $           dtsnw,omsnw,gsnw,
     $           edir, edn, eup, fdir, fdn, fup)

* Spectral irradiance, W m-2 nm-1, down-welling:

            DO iz = 1, nz
               sirrad(iz,iw) = etf(iw) * 
     $           (dirsun*edir(iz) + difdn*edn(iz) + difup*eup(iz))
            ENDDO

* Spectral actinic flux, quanta s-1 nm-1 cm-2, all directions:
*    units conversion:  1.e-4 * (wc*1e-9) / (hc = 6.62E-34 * 2.998E8) 

            DO iz = 1, nz
               saflux(iz,iw) = etf(iw)* 5.039e11 * wc(iw) *
     $              (dirsun*fdir(iz) + difdn*fdn(iz) + difup*fup(iz))
            ENDDO

*** Accumulate weighted integrals over wavelength, at all altitudes:

            DO 18 iz = 1, nz

* Weighted irradiances (dose rates) W m-2

               DO is = 1, ns
                  drdw = sirrad(iz,iw) * sw(is,iw) 
                  rate(is,iz) = rate(is,iz) + drdw * (wu(iw) - wl(iw))
               ENDDO

* Photolysis rate coefficients (J-values) s-1

               DO ij = 1, nj
                  djdw = saflux(iz,iw) * sj(ij,iz,iw)
                  valj(ij,iz) = valj(ij,iz) + djdw * (wu(iw) - wl(iw))
               ENDDO

 18         CONTINUE

* Save irradiances and actinic fluxes for output

            CALL saver1(it, itfix, iw, iwfix,  nz, izout,
     $           sirrad, saflux,
     $           svi_zw, svf_zw, svi_zt, svf_zt, svi_tw, svf_tw)

 10      CONTINUE

*^^^^^^^^^^^^^^^^ end wavelength loop
 
* Save dose rates and j-values for output

         CALL saver2(it,itfix, nz,izout, ns,isfix,ims, nj,ijfix,imj,
     $        rate, valj,
     $        svr_zs, svj_zj, svr_zt, svj_zt, svr_ts, svj_tj)
*BHH: I don't know why, but sza is reseting itself through the saver2 routine
*         sza(it) = zen

 20   CONTINUE

*^^^^^^^^^^^^^^^^ end time/zenith loop

** reset wavelength scale if needed:
      
      IF(lrefr) THEN
         write(*,*) 'applying vacuum to air wavelength shift'
     $        , izout, airout
         mrefr = -mrefr
         CALL wshift(mrefr, nw, wl, airout)
         CALL wshift(mrefr, nwint, wc, airout)
         CALL wshift(mrefr, nwint, wu, airout)
      ENDIF

** output 

C      call outpt1( outfil, iout, 
C     $     lirrad, laflux, lrates, ljvals, lmmech, lzenit,
C     $     nms, ims, nmj, imj,
C     $     nz, z, tlev, aircon, izout,
C     $     nw, wl, etf, iwfix,
C     $     nt, t, sza, itfix,
C     $     ns, slabel, isfix, nj, jlabel, ijfix,
C     $     svj_zj, svj_tj, svj_zt,
C     $     svr_zs, svr_ts, svr_zt,
C     $     svf_zw, svf_tw, svf_zt,
C     $     svi_zw, svi_tw, svi_zt )

*_______________________________________________________________________


      IF(intrct) THEN
         WRITE(*,*) 'do you want to do another calculation?'
         WRITE(*,*) 'y = yes'
         WRITE(*,*) 'any other key = no'
         READ(*,1001) again
 1001    FORMAT(A1)
         IF(again .EQ. 'y' .OR. again .EQ. 'Y') GO TO 1000
      ENDIF

      CLOSE(iout)
      CLOSE(kout)
      END SUBROUTINE tuv

      subroutine spline (n, x, y, b, c, d)
c======================================================================
c  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
c  for a cubic interpolating spline
c
c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
c
c    for  x(i) .le. x .le. x(i+1)
c
c  input..
c
c    n = the number of data points or knots (n.ge.2)
c    x = the abscissas of the knots in strictly increasing order
c    y = the ordinates of the knots
c
c  output..
c
c    b, c, d  = arrays of spline coefficients as defined above.
c
c  using  p  to denote differentiation,
c
c    y(i) = s(x(i))
c    b(i) = sp(x(i))
c    c(i) = spp(x(i))/2
c    d(i) = sppp(x(i))/6  (derivative from the right)
c
c  the accompanying function subprogram  seval  can be used
c  to evaluate the spline.
c======================================================================
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
      integer nm1, ib, i
      double precision t
c
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
c
c  set up tridiagonal system
c
c  b = diagonal, d = offdiagonal, c = right hand side.
c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
  10  continue
c
c  end conditions.  third derivatives at  x(1)  and  x(n)
c  obtained from divided differences
c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
c
c  forward elimination
c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
c
c  back substitution
c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
c
c  c(i) is now the sigma(i) of the text
c
c  compute polynomial coefficients
c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end

      double precision function seval(n, u, x, y, b, c, d, debug)
c======================================================================
c  this subroutine evaluates the cubic spline function
c
c    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
c
c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
c
c  if  u .lt. x(1) then  i = 1  is used.
c  if  u .ge. x(n) then  i = n  is used.
c
c  input..
c
c    n = the number of data points
c    u = the abscissa at which the spline is to be evaluated
c    x,y = the arrays of data abscissas and ordinates
c    b,c,d = arrays of spline coefficients computed by spline
c
c  if  u  is not in the same interval as the previous call, then a
c  binary search is performed to determine the proper interval.
c========================================================================
      integer n
      double precision  u, x(n), y(n), b(n), c(n), d(n)
      integer i, j, k
      double precision dx
      logical debug
      
      i=1
      
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 4510
      if ( u .le. x(i+1) ) go to 4530
c
c  binary search
c
4510  i = 1
      j = n+1
4520  k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 4520
c
c  evaluate spline
c
4530  dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      if (debug) then
          write(*,*) 'i,THETA,seval:', i,u,dx,seval
          write(*,8879) 'b     :', b
          write(*,8879) 'c     :', c
          write(*,8879) 'd     :', d
          write(*,8879) 'szas  :', x
          write(*,8879) 'svj_tj:', y
      endif
8879     FORMAT(1A6,100000(E26.17))
      return
      end

      subroutine polint(f,a,n,x,r)
c----------------------------------------------------------
c service program for fintr
c----------------------------------------------------------
        real*8 a, al, f, x, r
        integer*4 i, j, n
        dimension f(n),a(n)
        r=0.0
        do 1 j=1,n
        al=1.0
        do 2 i=1,n
          if (i-j) 3,2,3
 3        al=al*(x-a(i))/(a(j)-a(i))
 2      continue
 1      r=r+al*f(j)
        return
      end

      subroutine set_up_photol(iyear, imonth, iday, O3col, ralt, rlat,
     &           rlon, box_temp, raird, bs, cs, ds, sza, svj_tj)
        incLude 'params'
        real*8 b(19),c(19),d(19)
        real*8 bs(19,kj), cs(19,kj), ds(19,kj)
        REAL*8 O3col, ralt, box_temp, rlat, rlon, raird
        REAL*4 sO3col, sralt, sbox_temp, srlat, srlon, sraird
        REAL*8 y,dy,x
        integer i, n, j, iyear, imonth, iday
        real svj_tj(kt,kj), sza(kt) 
        real*8 temp2(19), temp(19)
        character*50 jlabels(kj)

               
        so3col = real(o3col)
        sralt = real(ralt)
        srlat = real(rlat)
        srlon = real(rlon)
        sbox_temp = real(box_temp)
        sraird = real(raird)
        call tuv(iyear, imonth, iday, so3col, sralt, 
     &    srlat, srlon, sbox_temp, sraird,
     &    sza, svj_tj, jlabels)
        do j=1,kj
          do i=1,19 
            temp(i)=sza(i)
            temp2(i)=svj_tj(i,j)
          enddo    
                                
          n=19
          call spline(n,temp,temp2,b,c,d)
          do i=1,19
            bs(i,j)=b(i)
            cs(i,j)=c(i)
            ds(i,j)=d(i)    
          enddo
        enddo
     
        return 
      end
