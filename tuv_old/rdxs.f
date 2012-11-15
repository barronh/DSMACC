* This file contains the following subroutines, related to reading the
* absorption cross sections of atmospheric gases:
*     rdno2xs
*     rdo2xs
*     rdo3xs
*       o3xs_mm
*       o3xs_mal
*       o3xs_bass
*     rdso2xs
*=============================================================================*

      SUBROUTINE rdno2xs(nw,wl,no2xs)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read NO2 molecular absorption cross section.  Re-grid data to match      =*
*=  specified wavelength working grid.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  NO2XS  - REAL, molecular absoprtion cross section (cm^2) of NO2 at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input: (altitude working grid)
      INTEGER nw
      REAL wl(kw)

* output:

      REAL no2xs(kw)

* local:
      REAL x1(kdata)
      REAL y1(kdata)
      REAL yg(kw)
      REAL dum
      INTEGER ierr
      INTEGER i, l, n, idum
      CHARACTER*40 fil
*_______________________________________________________________________

************* absorption cross sections:
*     measurements by:
* Davidson, J. A., C. A. Cantrell, A. H. McDaniel, R. E. Shetter,
* S. Madronich, and J. G. Calvert, Visible-ultraviolet absorption
* cross sections for NO2 as a function of temperature, J. Geophys.
* Res., 93, 7105-7112, 1988.
*  Values at 273K from 263.8 to 648.8 nm in approximately 0.5 nm intervals

      fil = 'DATAE1/NO2/NO2_ncar_00.abs'
      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      n = 750
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), dum, dum, idum
      ENDDO
      CLOSE(kin)

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
         no2xs(l) = yg(l)
   13 CONTINUE

*_______________________________________________________________________

      RETURN
      END

*=============================================================================*

      SUBROUTINE rdo2xs(nw,wl,o2xs1)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Compute equivalent O2 cross section, except                              =*
*=  the SR bands and the Lyman-alpha line.                                   =*
*-----------------------------------------------------------------------------* 
*=  PARAMETERS:                                   
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
*=            working wavelength grid           
*=            vertical layer at each specified wavelength                    =*
*=  O2XS1   - REAL, O2 molecular absorption cross section                    =*
*=
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* Input

      INTEGER nw
      REAL wl(kw)

* Output O2 xsect, temporary, will be over-written in Lyman-alpha and 
*   Schumann-Runge wavelength bands.

      REAL o2xs1(kw)

* Internal

      INTEGER i, n, kdata
      PARAMETER (kdata = 200)
      REAL x1(kdata), y1(kdata)
      REAL x, y
      INTEGER ierr

*-----------------------------------------------------

* Read O2 absorption cross section data:
*  116.65 to 203.05 nm = from Brasseur and Solomon 1986
*  205 to 240 nm = Yoshino et al. 1988

* Note that subroutine la_srb.f will over-write values in the spectral regions
*   corresponding to:
* - Lyman-alpha (LA: 121.4-121.9 nm, Chabrillat and Kockarts parameterization) 
* - Schumann-Runge bands (SRB: 174.4-205.8 nm, Koppers parameteriaztion)

      n = 0

      OPEN(UNIT=kin,FILE='DATAE1/O2/O2_brasseur.abs')
      DO i = 1, 7
         READ(kin,*)
      ENDDO
      DO i = 1, 78
         READ(kin,*) x, y
         IF (x .LE. 204.) THEN
            n = n + 1
            x1(n) = x
            y1(n) = y
         ENDIF
      ENDDO
      CLOSE(kin)

      OPEN(UNIT=kin,FILE='DATAE1/O2/O2_yoshino.abs',STATUS='old')
      DO i = 1, 8
         READ(kin,*)
      ENDDO
      DO i = 1, 36
         n = n + 1
         READ(kin,*) x, y
         y1(n) = y*1.E-24
         x1(n) = x
      END DO
      CLOSE (kin)

* Add termination points and interpolate onto the 
*  user grid (set in subroutine gridw):

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,0.               ,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,              1.E+38,0.)
      CALL inter2(nw,wl,o2xs1, n,x1,y1, ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O2 -> O + O'
         STOP
      ENDIF

*------------------------------------------------------

      RETURN
      END

*=============================================================================*

      SUBROUTINE rdo3xs(nw,wl,nz,tlev,o3xs)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read ozone molecular absorption cross section.  Re-grid data to match    =*
*=  specified wavelength working grid.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  O3XS   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
*=           each specified wavelength (WMO value at 273)                    =*
*=  S226   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
*=           each specified wavelength (value from Molina and Molina at 226K)=*
*=  S263   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
*=           each specified wavelength (value from Molina and Molina at 263K)=*
*=  S298   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
*=           each specified wavelength (value from Molina and Molina at 298K)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input: (altitude working grid)

      INTEGER nw
      REAL wl(kw)

      INTEGER nz
      REAL tlev(kz)

* internal

      INTEGER mopt

* output:
* ozone absorption cross section at three different 
* temperatures: 226, 263, 298 Kelvin.  Can interpolate
* to different temperatures. Units are cm2 molecule-1

      REAL o3xs(kz,kw)

*_______________________________________________________________________
* cross sections from WMO 1985 Ozone Assessment
* from 175.439 to 847.500 nm. Using value at 273 K.
* Values are over-written in Hartly and Huggins bands, using different
* options depending on value of mopt:
*     mopt = 1 = Molina and Molina
*     mopt = 2 = Malicet et al.
*     mopt = 3 = Bass et al.

      mopt = 1

      IF(mopt .eq. 1) CALL o3xs_mm(nw,wl,nz,tlev, o3xs)
      IF(mopt .eq. 2) CALL o3xs_mal(nw,wl,nz,tlev, o3xs)
      IF(mopt .eq. 3) CALL o3xs_bass(nw,wl,nz,tlev, o3xs)

      RETURN
      END

*=============================================================================*

      SUBROUTINE o3xs_mm(nw,wl,nz,tlev, xs)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read and interpolate the O3 cross section                                =*
*=  Combined data from WMO 85 Ozone Assessment (use 273K value from          =*
*=  175.439-847.5 nm) and:                                                   =*
*=  For Hartley and Huggins bands, use temperature-dependent values from     =*
*=  Molina, L. T., and M. J. Molina, Absolute absorption cross sections      =*
*=  of ozone in the 185- to 350-nm wavelength range, J. Geophys. Res.,       =*
*=  vol. 91, 14501-14508, 1986.                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  XS     - REAL, cross section (cm^2) for O3                            (O)=*
*=           at each defined wavelength and each defined altitude level      =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

*  input

      INTEGER nw, iw
      REAL wl(kw)

      INTEGER nz, iz
      REAL tlev(kz)

* output

      REAL xs(kz,kw)

* internal

      INTEGER kdata
      PARAMETER (kdata = 250)

* for wmo values

      INTEGER i, n, idum
      REAL a1, a2, dum
      REAL x1(kdata), y1(kdata)
      
      INTEGER ierr
      REAL yg(kw), o3xs(kw)

* for Molina & Molina values

      INTEGER n1, n2, n3
      REAL x2(kdata), y2(kdata), x3(kdata), y3(kdata)
      REAL s226(kw), s263(kw), s298(kw)

*----------------------------------------------------------
* cross sections from WMO 1985 Ozone Assessment
* from 175.439 to 847.500 nm
* use value at 273 K

      OPEN(UNIT=kin,FILE='DATAE1/wmo85',STATUS='old')
      DO i = 1, 3
         read(kin,*)
      ENDDO
      n = 158
      DO i = 1, n
         READ(kin,*) idum, a1, a2, dum, dum, dum, dum, y1(i)
         x1(i) = (a1+a2)/2.
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 cross section - WMO'
         STOP
      ENDIF

      DO iw = 1, nw-1
         o3xs(iw) = yg(iw)
      ENDDO

* For Hartley and Huggins bands, use temperature-dependent values from
* Molina, L. T., and M. J. Molina, Absolute absorption cross sections
* of ozone in the 185- to 350-nm wavelength range,
* J. Geophys. Res., vol. 91, 14501-14508, 1986.

      OPEN(UNIT=kin,FILE='DATAE1/O3/O3.molina.abs',STATUS='old')
      DO i = 1, 5
         READ(kin,*)
      ENDDO
      n1 = 220
      n2 = 220
      n3 = 220
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i), y2(i), y3(i)
         x2(i) = x1(i)
         x3(i) = x1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 xsect - 226K Molina'
         STOP
      ENDIF
      DO iw = 1, nw-1
         s226(iw) = yg(iw)*1.E-20
      ENDDO

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 xsect - 263K Molina'
         STOP
      ENDIF
      DO iw = 1, nw-1
         s263(iw) = yg(iw)*1.E-20
      ENDDO

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,               0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 xsect - 298K Molina'
         STOP
      ENDIF
      DO iw = 1, nw-1
         s298(iw) = yg(iw)*1.E-20
      ENDDO

* assign:

      DO 10 iw = 1, nw-1
         DO 20 iz = 1, nz

            xs(iz,iw) = o3xs(iw)
            IF ( wl(iw) .GT. 240.5  .AND. wl(iw+1) .LT. 350. ) THEN
               IF (tlev(iz) .LT. 263.) THEN
                  xs(iz,iw) = s226(iw) + 
     $                 (s263(iw)-s226(iw))*(tlev(iz)-226.)/(263.-226.) 
               ELSE
                  xs(iz,iw) = s263(iw) + 
     $                 (s298(iw)-s263(iw))*(tlev(iz)-263.)/(298.-263.)
               ENDIF
            ENDIF

 20      CONTINUE

 10   CONTINUE

      RETURN
      END

*=============================================================================*

      SUBROUTINE o3xs_mal(nw,wl,nz,tlev, xs)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read and interpolate the O3 cross section                                =*
*=  Combined data from WMO 85 Ozone Assessment (use 273K value from          =*
*=  175.439-847.5 nm) and:                                                   =*
*=  For Hartley and Huggins bands, use temperature-dependent values from     =*
*=  Malicet et al., J. Atmos. Chem.  v.21, pp.263-273, 1995.                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  XS     - REAL, cross section (cm^2) for O3                            (O)=*
*=           at each defined wavelength and each defined altitude level      =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

*  input

      INTEGER nw, iw
      REAL wl(kw)

      INTEGER nz, iz
      REAL tlev(kz)

* output

      REAL xs(kz,kw)

* internal

      INTEGER kdata
      PARAMETER (kdata = 16000)

* for wmo values

      INTEGER i, n, idum
      REAL a1, a2, dum
      REAL x1(kdata), y1(kdata)
      
      INTEGER ierr
      REAL yg(kw), o3xs(kw)

* for Malicet values

      INTEGER n1, n2, n3, n4
      REAL x2(kdata), x3(kdata), x4(kdata)
      REAL y2(kdata), y3(kdata), y4(kdata)
      REAL s218(kw), s228(kw), s243(kw), s295(kw)

*----------------------------------------------------------
* cross sections from WMO 1985 Ozone Assessment
* from 175.439 to 847.500 nm
* use value at 273 K

      OPEN(UNIT=kin,FILE='DATAE1/wmo85',STATUS='old')
      DO i = 1, 3
         read(kin,*)
      ENDDO
      n = 158
      DO i = 1, n
         READ(kin,*) idum, a1, a2, dum, dum, dum, dum, y1(i)
         x1(i) = (a1+a2)/2.
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 cross section - WMO'
         STOP
      ENDIF

      DO iw = 1, nw-1
         o3xs(iw) = yg(iw)
      ENDDO

*=  For Hartley and Huggins bands, use temperature-dependent values from     =*
*=  Malicet et al., J. Atmos. Chem.  v.21, pp.263-273, 1995.                 =*

      OPEN(UNIT=kin,FILE='DATAE1/O3/o3absqs.dat',STATUS='old')
      DO i = 1, 1
         READ(kin,*)
      ENDDO
      n1 = 15001
      n2 = 15001
      n3 = 15001
      n4 = 15001

      DO i = 1, n1
         READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i)
         x2(i) = x1(i)
         x3(i) = x1(i)
         x4(i) = x1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 xsect - 295K Malicet'
         STOP
      ENDIF
      DO iw = 1, nw-1
         s295(iw) = yg(iw)
      ENDDO

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 xsect - 243K Malicet'
         STOP
      ENDIF
      DO iw = 1, nw-1
         s243(iw) = yg(iw)
      ENDDO

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,               0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 xsect - 228K Malicet'
         STOP
      ENDIF
      DO iw = 1, nw-1
         s228(iw) = yg(iw)
      ENDDO

      CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,               0.,0.)
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n4,x4,y4,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 xsect - 218K Malicet'
         STOP
      ENDIF
      DO iw = 1, nw-1
         s218(iw) = yg(iw)
      ENDDO

* assign:

      DO 10 iw = 1, nw-1
         DO 20 iz = 1, nz

            xs(iz,iw) = o3xs(iw)
            IF ( wl(iw) .GT. 195.  .AND. wl(iw+1) .LT. 345. ) THEN
               IF(tlev(iz) .GE. 243.) THEN
                  xs(iz,iw) = s243(iw) + 
     $                 (s295(iw)-s243(iw))*(tlev(iz)-243.)/(295.-243.)
               ENDIF
               IF(tlev(iz) .LT. 254. .AND. tlev(iz) .GE. 228.) THEN
                  xs(iz,iw) = s228(iw) + 
     $                 (s243(iw)-s228(iw))*(tlev(iz)-228.)/(243.-228.)
               ENDIF
               IF(tlev(iz) .LT. 228.) THEN
                  xs(iz,iw) = s218(iw) + 
     $                 (s228(iw)-s218(iw))*(tlev(iz)-218.)/(228.-218.)
               ENDIF
            ENDIF

 20      CONTINUE

 10   CONTINUE

      RETURN
      END

*=============================================================================*

      SUBROUTINE o3xs_bass(nw,wl,nz,tlev, xs)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read and interpolate the O3 cross section                                =*
*=  Combined data from WMO 85 Ozone Assessment (use 273K value from          =*
*=  175.439-847.5 nm) and:                                                   =*
*=  For Hartley and Huggins bands, use temperature-dependent values from     =*
*=  Bass
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  XS     - REAL, cross section (cm^2) for O3                            (O)=*
*=           at each defined wavelength and each defined altitude level      =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

*  input

      INTEGER nw, iw
      REAL wl(kw)

      INTEGER nz, iz
      REAL tlev(kz)

* output

      REAL xs(kz,kw)

* internal

      INTEGER kdata
      PARAMETER (kdata = 2000)

* for wmo values

      INTEGER i, n, idum
      REAL a1, a2, dum
      REAL x1(kdata), y1(kdata)
      
      INTEGER ierr
      REAL yg(kw), o3xs(kw)

* for Bass values

      INTEGER n1, n2, n3
      REAL x2(kdata), y2(kdata), x3(kdata), y3(kdata)
      REAL c0(kw), c1(kw), c2(kw)
      REAL tc

*----------------------------------------------------------
* cross sections from WMO 1985 Ozone Assessment
* from 175.439 to 847.500 nm
* use value at 273 K

      OPEN(UNIT=kin,FILE='DATAE1/wmo85',STATUS='old')
      DO i = 1, 3
         read(kin,*)
      ENDDO
      n = 158
      DO i = 1, n
         READ(kin,*) idum, a1, a2, dum, dum, dum, dum, y1(i)
         x1(i) = (a1+a2)/2.
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 cross section - WMO'
         STOP
      ENDIF

      DO iw = 1, nw-1
         o3xs(iw) = yg(iw)
      ENDDO

* For Hartley and Huggins bands, use temperature-dependent values from
*  Bass et al.

      OPEN(UNIT=kin,FILE='DATAE1/O3/O3_bass.abs',STATUS='old')
      DO i = 1, 8
         READ(kin,*)
      ENDDO
      n1 = 1915
      n2 = 1915
      n3 = 1915
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i), y2(i), y3(i)
         x2(i) = x1(i)
         x3(i) = x1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 xsect - c0 Bass'
         STOP
      ENDIF
      DO iw = 1, nw-1
         c0(iw) = yg(iw)
      ENDDO

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 xsect - c1 Bass'
         STOP
      ENDIF
      DO iw = 1, nw-1
         c1(iw) = yg(iw)
      ENDDO

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,               0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O3 xsect - c2 Bass'
         STOP
      ENDIF
      DO iw = 1, nw-1
         c2(iw) = yg(iw)
      ENDDO

* assign:

      DO 10 iw = 1, nw-1
         DO 20 iz = 1, nz

            tc = tlev(iz) - 273.15

            xs(iz,iw) = o3xs(iw)
            IF ( wl(iw) .GT. 245. .AND. wl(iw+1) .LT. 341. ) THEN

               xs(iz,iw) = 1.e-20 * (c0(iw) + c1(iw)*tc + 
     $              c2(iw)*tc*tc)

            ENDIF

 20      CONTINUE

 10   CONTINUE

      RETURN
      END

*=============================================================================*

      SUBROUTINE rdso2xs(nw,wl,so2xs)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read SO2 molecular absorption cross section.  Re-grid data to match      =*
*=  specified wavelength working grid.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  SO2XS  - REAL, molecular absoprtion cross section (cm^2) of SO2 at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Changed offset for grid-end interpolation to relative number      =*
*=         (x * (1 +- deltax)                                                =*
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

      REAL so2xs(kw)

* local:
      REAL x1(kdata)
      REAL y1(kdata)
      REAL yg(kw)
      REAL dum
      INTEGER ierr
      INTEGER i, l, n, idum
      CHARACTER*40 fil
*_______________________________________________________________________

************* absorption cross sections:
* SO2 absorption cross sections from J. Quant. Spectrosc. Radiat. Transfer
* 37, 165-182, 1987, T. J. McGee and J. Burris Jr.
* Angstrom vs. cm2/molecule, value at 221 K

      fil = 'DATA/McGee87'
      OPEN(UNIT=kin,FILE='DATAE1/SO2/SO2xs.all',STATUS='old')
      DO 11, i = 1,3 
         read(kin,*)
   11 CONTINUE
c      n = 681 
      n = 704 
      DO 12, i = 1, n
         READ(kin,*) x1(i), y1(i)
         x1(i) = x1(i)/10.
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
         so2xs(l) = yg(l)
   13 CONTINUE

*_______________________________________________________________________

      RETURN
      END
