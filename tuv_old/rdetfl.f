* This file contains the following subroutines, related to reading the
* extraterrestrial spectral irradiances:
*     rdetfl
*     read1
*     read2
*=============================================================================*

      SUBROUTINE rdetfl(nw,wl,f)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read and re-grid extra-terrestrial flux data.                            =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      integer kdata
      parameter(kdata=20000)

* input: (wavelength grid)
      INTEGER nw
      REAL wl(kw)
      INTEGER iw

* output: (extra terrestrial solar flux)
      REAL f(kw)

* INTERNAL:

* work arrays for input data files:

      CHARACTER*40 fil
      REAL x1(kdata)
      REAL y1(kdata)
      INTEGER nhead, n, i, ierr
      REAL dum

* data gridded onto wl(kw) grid:

      REAL yg1(kw)
      REAL yg2(kw)
      REAL yg3(kw)

      REAL hc
      PARAMETER(hc = 6.62E-34 * 2.998E8)

      INTEGER msun

*_______________________________________________________________________
* select desired extra-terrestrial solar irradiance, using msun:

*  1 =   extsol.flx:  De Luisi, JGR 80, 345-354, 1975
*                     280-400 nm, 1 nm steps.
*  2 =   lowsun3.flx:  Lowtran (John Bahr, priv. comm.)
*                      173.974-500000 nm, ca. 0.1 nm steps in UV-B
*  3 =   modtran1.flx:  Modtran (Gail Anderson, priv. comm.)
*                       200.55-949.40, 0.05 nm steps
*  4 =   nicolarv.flx:  wvl<300 nm from Nicolet, Plan. Sp. Sci., 29,  951-974, 1981.
*                       wvl>300 nm supplied by Thekaekera, Arvesen Applied Optics 8, 
*                       11, 2215-2232, 1969 (also see Thekaekera, Applied Optics, 13,
*                       3, 518, 1974) but with corrections recommended by:
*                       Nicolet, Plan. Sp. Sci., 37, 1249-1289, 1989.
*                       270.0-299.0 nm in 0.5 nm steps
*                       299.6-340.0 nm in ca. 0.4 nm steps
*                       340.0-380.0 nm in ca. 0.2 nm steps
*                       380.0-470.0 nm in ca. 0.1 nm steps   
*  5 =  solstice.flx:  From:   MX%"ROTTMAN@virgo.hao.ucar.edu" 12-OCT-1994 13:03:01.62
*                      Original data gave Wavelength in vacuum
*                      (Converted to wavelength in air using Pendorf, 1967, J. Opt. Soc. Am.)
*                      279.5 to 420 nm, 0.24 nm spectral resolution, approx 0.07 nm steps
*  6 =  suntoms.flx: (from TOMS CD-ROM).  280-340 nm, 0.05 nm steps.
*  7 =  neckel.flx:  H.Neckel and D.Labs, "The Solar Radiation Between 3300 and 12500 A",
*                    Solar Physics v.90, pp.205-258 (1984).
*                    1 nm between 330.5 and 529.5 nm
*                    2 nm between 631.0 and 709.0 nm
*                    5 nm between 872.5 and 1247.4 nm
*                    Units: must convert to W m-2 nm-1 from photons cm-2 s-1 nm-1
*  8 =  atlas3.flx:  ATLAS3-SUSIM 13 Nov 94 high resolution (0.15 nm FWHM)
*                    available by ftp from susim.nrl.navy.mil
*                    atlas3_1994_317_a.dat, downloaded 30 Sept 98.
*                    150-407.95 nm, in 0.05 nm steps
*                    (old version from Dianne Prinz through Jim Slusser)
*                    orig wavelengths in vac, correct here to air.
*  9 =  solstice.flx:  solstice 1991-1996, average
*                    119.5-420.5 nm in 1 nm steps

* 10 =  susim_hi.flx:  SUSIM SL2 high resolution
*                      120.5-400.0 in 0.05 nm intervals (0.15 nm resolution)
* 11 =  wmo85.flx: from WMO 1995 Ozone Atmospheric Ozone (report no. 16)
*                  on variable-size bins.  Original values are per bin, not
*                  per nm.
* 12 = combine susim_hi.flx for .lt. 350 nm, neckel.flx for .gt. 350 nm.
*
* 13 = combine 
*     for wl(iw) .lt. 150.01                                susim_hi.flx
*     for wl(iw) .ge. 150.01 and wl(iw) .le. 400            atlas3.flx 
*     for wl(iw) .gt. 400                                   Neckel & Labs 

      msun = 13

* simple files are read and interpolated here in-line. Reading of 
* more complex files may be done with longer code in a read#.f subroutine.

 1    IF (msun .EQ. 1) THEN
         fil = 'DATAE1/SUN/extsol.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 3
         n =121
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

 2    ELSEIF (msun .EQ. 2) THEN
         fil = 'DATAE1/SUN/lowsun3.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 3
         n = 4327
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

 3    ELSEIF (msun .EQ. 3) THEN
         fil = 'DATAE1/SUN/modtran1.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 6
         n = 14980
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

 4    ELSEIF (msun .EQ. 4) THEN
         fil = 'DATAE1/SUN/nicolarv.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 8
         n = 1260
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

 5    ELSEIF (msun .EQ. 5) THEN
* unofficial - do not use
         fil = 'DATAE2/SUN/solstice.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 11
         n = 2047
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

 6    ELSEIF (msun .EQ. 6) THEN
* unofficial - do not use
         fil = 'DATAE2/SUN/suntoms.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 3
         n = 1200
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i)* 1.e-3
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

 7    ELSEIF (msun .EQ. 7) THEN
         fil = 'DATAE1/SUN/neckel.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 11
         n = 496
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) dum, y1(i)
            if (dum .lt. 630.0) x1(i) = dum - 0.5
            if (dum .gt. 630.0 .and. dum .lt. 870.0) x1(i) = dum - 1.0
            if (dum .gt. 870.0) x1(i) = dum - 2.5
            y1(i) = y1(i) * 1.E4 * hc / (dum * 1.E-9)
         ENDDO
         CLOSE (kin)
         x1(n+1) = x1(n) + 2.5
         do i = 1, n
            y1(i) = y1(i) * (x1(i+1)-x1(i))
         enddo
         call inter3(nw,wl,yg2,n+1,x1,y1,0)
         do iw = 1, nw-1
            yg1(iw) = yg1(iw) / (wl(iw+1)-wl(iw))
         enddo
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

 8    ELSEIF (msun .EQ. 8) THEN
         nhead = 5
         fil = 'DATAE1/SUN/atlas3_1994_317_a.dat'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 13
         n = 5160
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-3
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

 9    ELSEIF (msun .EQ. 9) THEN
         fil = 'DATAE1/SUN/solstice.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 2
         n = 302
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

 10   ELSEIF (msun .EQ. 10) THEN
         WRITE(kout,*) 'DATAE1/SUN/susim_hi.flx'
         CALL read1(nw,wl,yg1)
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO


 11   ELSEIF (msun .EQ. 11) THEN
         WRITE(kout,*) 'DATAE1/SUN/wmo85.flx'
         CALL read2(nw,wl,yg1)
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

 12   ELSEIF (msun .EQ. 12) THEN
         WRITE(kout,*) 'DATAE1/SUN/susim_hi.flx'
         CALL read1(nw,wl,yg1)
         fil = 'DATAE1/SUN/neckel.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 11
         n = 496
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) dum, y1(i)
            if (dum .lt. 630.0) x1(i) = dum - 0.5
            if (dum .gt. 630.0 .and. dum .lt. 870.0) x1(i) = dum - 1.0
            if (dum .gt. 870.0) x1(i) = dum - 2.5
            y1(i) = y1(i) * 1.E4 * hc / (dum * 1.E-9)
         ENDDO
         CLOSE (kin)
         x1(n+1) = x1(n) + 2.5
         do i = 1, n
            y1(i) = y1(i) * (x1(i+1)-x1(i))
         enddo
         call inter3(nw,wl,yg2,n+1,x1,y1,0)
         do iw = 1, nw-1
            yg2(iw) = yg2(iw) / (wl(iw+1)-wl(iw))
         enddo

         DO iw = 1, nw-1
            IF (wl(iw) .GT. 350.) THEN
               f(iw) = yg2(iw)
            ELSE
               f(iw) = yg1(iw)
            ENDIF
         ENDDO

 13   ELSEIF (msun .EQ. 13) THEN

         WRITE(kout,*) 'DATAE1/SUN/susim_hi.flx'
         CALL read1(nw,wl,yg1)

         nhead = 5
         fil = 'DATAE1/SUN/atlas3_1994_317_a.dat'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         n = 5160
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-3
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         

         fil = 'DATAE1/SUN/neckel.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 11
         n = 496
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) dum, y1(i)
            if (dum .lt. 630.0) x1(i) = dum - 0.5
            if (dum .gt. 630.0 .and. dum .lt. 870.0) x1(i) = dum - 1.0
            if (dum .gt. 870.0) x1(i) = dum - 2.5
            y1(i) = y1(i) * 1.E4 * hc / (dum * 1.E-9)
         ENDDO
         CLOSE (kin)

         x1(n+1) = x1(n) + 2.5
         call inter4(nw,wl,yg3,n+1,x1,y1,0)

         DO iw = 1, nw-1

            IF (wl(iw) .LT. 150.01) THEN
               f(iw) = yg1(iw)
            ELSE IF ((wl(iw) .GE. 150.01) .AND. wl(iw) .LE. 400.) THEN
               f(iw) = yg2(iw)
            ELSE IF (wl(iw) .GT. 400.) THEN
               f(iw) = yg3(iw)
            ENDIF

         ENDDO


      ENDIF

      RETURN
      END

*=============================================================================*

      SUBROUTINE read1(nw,wl,f)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read extra-terrestrial flux data.  Re-grid data to match specified       =*
*=  working wavelength grid.                                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input: (wavelength grid)
      INTEGER nw
      REAL wl(kw)

* output: (extra terrestrial solar flux)
      REAL f(kw)

* local:

      REAL lambda_hi(10000),irrad_hi(10000)
      REAL lambda
      INTEGER ierr
      INTEGER i, j, n
      CHARACTER*40 FIL

*_______________________________________________________________________

******* SUSIM irradiance 
*_______________________________________________________________________
* VanHoosier, M. E., J.-D. F. Bartoe, G. E. Brueckner, and
* D. K. Prinz, Absolute solar spectral irradiance 120 nm -
* 400 nm (Results from the Solar Ultraviolet Spectral Irradiance
* Monitor - SUSIM- Experiment on board Spacelab 2), 
* Astro. Lett. and Communications, 1988, vol. 27, pp. 163-168.
*     SUSIM SL2 high resolution (0.15nm) Solar Irridance data.
*     Irradiance values are given in milliwatts/m^2/nanomenters
*     and are listed at 0.05nm intervals.  The wavelength given is
*     the center wavelength of the 0.15nm triangular bandpass.
*     Normalized to 1 astronomical unit.
*  DATA for wavelengths > 350 nm are unreliable
* (Van Hoosier, personal communication, 1994).
*_______________________________________________________________________

** high resolution

      fil = 'DATAE1/SUN/susim_hi.flx'
      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO 11, i = 1, 7
         READ(kin,*)
   11 CONTINUE
      DO 12, i = 1, 559
         READ(kin,*)lambda,(irrad_hi(10*(i-1)+j), j=1, 10)
   12 CONTINUE
      CLOSE (kin)

* compute wavelengths, convert from mW to W

      n = 559*10
      DO 13, i = 1, n
         lambda_hi(i)=120.5 + FLOAT(i-1)*.05
         irrad_hi(i) = irrad_hi(i)  /  1000.
   13 CONTINUE
*_______________________________________________________________________

      CALL addpnt(lambda_hi,irrad_hi,10000,n,
     >            lambda_hi(1)*(1.-deltax),0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,                 0.,0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,
     >            lambda_hi(n)*(1.+deltax),0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,              1.e38,0.)
      CALL inter2(nw,wl,f,n,lambda_hi,irrad_hi,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      RETURN
      END

*=============================================================================*

      SUBROUTINE read2(nw,wl,f)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read extra-terrestrial flux data.  Re-grid data to match specified       =*
*=  working wavelength grid.                                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input: (wavelength grid)
      INTEGER nw
      REAL wl(kw)
      REAL yg(kw)

*
      INTEGER iw

* output: (extra terrestrial solar flux)
      REAL f(kw)

* local:

      REAL x1(1000), y1(1000) 
      REAL x2(1000)
      REAL x3(1000)
      INTEGER i, n
      REAL DUM
      INTEGER IDUM

*_______________________________________________________________________

*********WMO 85 irradiance

      OPEN(UNIT=kin,FILE='DATAE1/SUN/wmo85.flx',STATUS='old')
      DO 11, i = 1, 3
         READ(kin,*)
   11 CONTINUE
      n = 158
      DO 12, i = 1, n
         READ(kin,*) idum, x1(i),x2(i),y1(i), dum, dum, dum
         x3(i) = 0.5 * (x1(i) + x2(i))

C average value needs to be calculated only if inter2 is
C used to interpolate onto wavelength grid (see below)
C        y1(i) =  y1(i) / (x2(i) - x1(i)) 

   12 CONTINUE
      CLOSE (kin)

      x1(n+1) = x2(n)

C inter2: INPUT : average value in each bin 
C         OUTPUT: average value in each bin
C inter3: INPUT : total area in each bin
C         OUTPUT: total area in each bin

      CALL inter3(nw,wl,yg, n+1,x1,y1,0)
C      CALL inter2(nw,wl,yg,n,x3,y1,ierr)

      DO 10,  iw = 1, nw-1
* from quanta s-1 cm-2 bin-1 to  watts m-2 nm-1
* 1.e4 * ([hc =] 6.62E-34 * 2.998E8)/(wc*1e-9) 
         
C the scaling by bin width needs to be done only if
C inter3 is used for interpolation

         yg(iw) = yg(iw) / (wl(iw+1)-wl(iw))
         f(iw) = yg(iw) * 1.e4 * (6.62E-34 * 2.998E8) / 
     $        ( 0.5 * (wl(iw+1)+wl(iw)) * 1.e-9)

   10 CONTINUE

      RETURN
      END
