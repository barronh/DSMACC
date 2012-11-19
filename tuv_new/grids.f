* This file contains the following subroutine, related to setting up
* grids for numerical calculations:
*     gridw
*     gridz
*     gridt
*     gridck
*=============================================================================*

      SUBROUTINE gridw(wstart, wstop, nwint,
     $     nw,wl,wc,wu)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create the wavelength grid for all interpolations and radiative transfer =*
*=  calculations.  Grid may be irregularly spaced.  Wavelengths are in nm.   =*
*=  No gaps are allowed within the wavelength grid.                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW  - INTEGER, number of wavelength grid _points_                     (O)=*
*=  WL  - REAL, vector carrying the lower limit of each wavel. interval   (O)=*
*=  WC  - REAL, vector carrying the center wavel of each wavel. interval  (O)=*
*=              (wc(i) = 0.5*(wl(i)+wu(i), i = 1..NW-1)                      =*
*=  WU  - REAL, vector carrying the upper limit of each wavel. interval   (O)=*
*=
*=  MOPT- INTEGER OPTION for wave-length IF 3 good for JO2                (O)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'params'

* input:

      REAL wstart, wstop
      INTEGER nwint

* output:

      REAL wl(kw), wc(kw), wu(kw)
      INTEGER nw

* local:

      INTEGER mopt
      REAL wincr
      INTEGER iw

      CHARACTER*40 fi
      CHARACTER*20 wlabel

      REAL dum

      LOGICAL ok
*_______________________________________________________________________

**** chose wavelength grid

* some pre-set options
*     mopt = 1    equal spacing
*     mopt = 2    grid defined in data table
*     mopt = 3    user-defined
*     mopt = 4    fast-TUV, troposheric wavelengths only

      mopt = 1
      IF(nwint .EQ. -156) mopt = 2
      IF(nwint .EQ. -7) mopt = 4

      IF (mopt .EQ. 1) GO TO 1
      IF (mopt .EQ. 2) GO TO 2
      IF (mopt .EQ. 3) GO TO 3
      IF (mopt .EQ. 4) GO TO 4

*_______________________________________________________________________

 1    CONTINUE

      wlabel = 'equal spacing'
      nw = nwint + 1
      wincr = (wstop - wstart) / FLOAT (nwint)
      DO 10, iw = 1, nw-1
         wl(iw) = wstart + wincr*FLOAT(iw-1)
         wu(iw) = wl(iw) + wincr
         wc(iw) = ( wl(iw) + wu(iw) )/2.
 10   CONTINUE
      wl(nw) = wu(nw-1)
      GO TO 9

*_______________________________________________________________________

 2    CONTINUE

* Input from table.  In this example:
* Wavelength grid will be read from a file.
* First line of table is:  nw = number of wavelengths (no. of intervals + 1)
* Then, nw wavelengths are read in, and assigned to wl(iw)
* Finally, wu(iw) and wc(iw) are computed from wl(iw)

c      wlabel = 'isaksen.grid'
      wlabel = 'combined.grid'

      fi = 'DATAE1/GRIDS/'//wlabel
      OPEN(unit=kin,file=fi,status='old')
      READ(kin,*) nw
      DO iw = 1, nw
         READ(kin,*) wl(iw)
      ENDDO
      CLOSE(kin)
      DO iw = 1, nw-1
         wu(iw) = wl(iw+1)
         wc(iw) = 0.5*(wl(iw) + wu(iw))
      ENDDO
      GO TO 9

*_______________________________________________________________________

 3    CONTINUE

* user-defined grid.  In this example, a single calculation is used to 
* obtain results for two 1 nm wide intervals centered at 310 and 400 nm:
* interval 1 : 1 nm wide, centered at 310 nm
* interval 3 : 2 nm wide, centered at 400 nm
* (inteval 2 : 310.5 - 399.5 nm, required to connect intervals 1 & 3)

      nw = 4
      wl(1) = 309.5
      wl(2) = 310.5
      wl(3) = 399.5
      wl(4) = 400.5
      DO iw = 1, nw-1
         wu(iw) = wl(iw+1)
         wc(iw) = 0.5*(wl(iw) + wu(iw))
      ENDDO
      GO TO 9

*_______________________________________________________________________

 4    CONTINUE
      wlabel = 'fast-TUV tropospheric grid'
      
      fi = 'DATAE1/GRIDS/fast_tuv.grid'
      OPEN(UNIT=kin,FILE=fi,STATUS='old')
      DO iw = 1, 4
         READ(kin,*)
      ENDDO

* skip wavelength shorter than 289.9 nm

      DO iw = 1, 10
         READ(kin,*)
      ENDDO
      nw = ABS(nwint) + 1
      DO iw = 1, nw-1
         READ(kin,*) dum, wl(iw), dum, dum
      ENDDO
      wl(nw) = dum
      DO iw = 1, nw-1
         wu(iw) = wl(iw+1)
         wc(iw) = 0.5*(wl(iw) + wu(iw))
      ENDDO



      GO TO 9

*_______________________________________________________________________

 9    CONTINUE

* check grid for assorted improprieties:

      CALL gridck(kw,nw,wl,ok)

      IF (.NOT. ok) THEN
         WRITE(*,*)'STOP in GRIDW:  The w-grid does not make sense'
         STOP
      ENDIF

*_______________________________________________________________________

      RETURN
      END

*=============================================================================*

      SUBROUTINE gridz(zstart, zstop, nz, z, zout, izout)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create the altitude grid for all interpolations and radiative transfer   =*
*=  calculations.  Grid may be irregularly spaced.  All altitudes are in     =*
*=  kilometers (km).  The altitude at index 1 specifies the surface elevation=*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  nz  - INTEGER, number of altitude points (levels)                     (O)=*
*=  z   - REAL, vector of altitude levels (in km)                         (O)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'params'

* input surface elevation, altitude for output

      REAL zstart, zstop, zout

* output: altitude working grid, index for output

      REAL z(kz)
      INTEGER nz
      INTEGER izout

* local:

      REAL zincr, tol
      INTEGER i, n, nlev
      LOGICAL ok
*_______________________________________________________________________

* Set vertical grid of the atmosphere: atmospheric level altitudes 
* (in real km), including top-most level.
* User specifies grid (surface at lowest km value), increasing
* upwards:
*     -  nz = total number of user levels
*     -  z(I) = altitude in km for each level.
* z(1) is the elevation of the surface (km asl), and can be specified either
* here or in the main program.
* Non-uniform spacing is possible:
*     - zincr = altitude increment between current and previous level (km)
*     - nlev = number of levels in current equally-spaced section
*     - n = index of top level of equally-spaced section
* Note "levels" are vertical points
*      "layers" are vertical distances between levels


* Grid selection options:
* 1 = standard equally spaced grid, manual
* 2 = standard equally spaced grid, auto generated
* 3 = variable spacing grid, example for snow
* 4 = mirage z-grid for Mexico City
* 5 = arbitrary user-defined grid

      GO TO 1

*-----grid option 2: manual -----------------
* entire grid (nz levels) in increments zincr 

 1    CONTINUE
      WRITE(*,*) 'equally spaced z-grid'
      zincr = (zstop - zstart) / FLOAT(nz - 1)
      z(1) = zstart
      DO i = 2, nz
         z(i) = z(1) + zincr*FLOAT(i-1)
      ENDDO
      GOTO 10

*-----grid option 3: automatic -----------------
* entire grid (nz levels) in increments zincr 

 2    CONTINUE
      WRITE(*,*) 'equally spaced z-grid'
      zincr = (zstop - zstart) / FLOAT(nz - 1)
      nlev = nz-1
      n = 1
      CALL buildz(zincr, nlev, n, z)
      GOTO 10

*-----grid option 4: variable grid example----------------------
*-----copy & edit this section for non-uniform grid----
* the example provided below is high vertical resolution in 
*   snow, with atmosphere above it.

 3    CONTINUE
      WRITE(*,*) 'snow-atmosphere grid'
* 0.-10. cm from ground, in 1 cm increments ( 1 cm = 1e-5 km):
      zincr = 1.e-5
      nlev = 10
      n = 1
      CALL buildz(zincr,nlev,n,z)

* 10-90 cm from ground, in 10 cm increments ( 1 cm = 1e-5 km):
      zincr = 1.e-4
      nlev = 8
      CALL buildz(zincr,nlev,n,z)

* 90-95 cm from ground, in 1x 5 cm increment ( 1 cm = 1e-5 km):
      zincr = 5.e-5
      nlev = 1
      CALL buildz(zincr,nlev,n,z)

* 95-99 cm from ground, in 4x 1 cm increments ( 1 cm = 1e-5 km):
      zincr = 1.e-5
      nlev = 4
      CALL buildz(zincr,nlev,n,z)

* 99-99.5 cm from ground, in 1x 0.5 cm increment ( 1 cm = 1e-5 km):
      zincr = 5.e-6
      nlev = 1
      CALL buildz(zincr,nlev,n,z)

* 99.5 centimeters - 1m, in 0.1 cm increments (1 cm = 1e-5 km):
      zincr = 1.e-6
      nlev = 5
      CALL buildz(zincr,nlev,n,z)

*atmosphere
* 1.-10. m in 1 m increments
      zincr = 1.e-3
      nlev = 9
      CALL buildz(zincr,nlev,n,z)

* 10.-100 m in 10 m increments
      zincr = 1.e-2
      nlev = 9
      CALL buildz(zincr,nlev,n,z)

* 100.- 1000. meters, in 100 m increments
      zincr = 1.e-1
      nlev = 9
      CALL buildz(zincr,nlev,n,z)

* 1.-2. km in 1 km increments
      zincr = 1.
      nlev =  1
      CALL buildz(zincr,nlev,n,z)

* 2.-80. km in 2 km increments
      zincr = 2.
      nlev =  39
      CALL buildz(zincr,nlev,n,z)

      GOTO 10

*-----grid option 4:  grid for Mexico City

 4    CONTINUE
      WRITE(*,*) 'mirage z-grid'

* grid for mirage km: incr(range)i 
* 0.1(0-4)   2-41
* 0.2(4-8)   42-61
*   1(8-30)  62-83
*   2(30-50) 84-93 
*   5(50-80) 94-99

      nz = 99
      z(1) = zstart
      DO i = 2, 41
         z(i) = z(1) + 0.1*FLOAT(i-1)
      ENDDO
      DO i = 42, 61
         z(i) = z(41) + 0.2*FLOAT(i-41)
      ENDDO
      DO i = 62, 83
         z(i) = z(61) + 1.*FLOAT(i-61)
      ENDDO
      DO i = 84, 93
         z(i) = z(83) + 2.*FLOAT(i-83)
      ENDDO
      DO i = 94, 99
         z(i) = z(93) + 5.*FLOAT(i-93)
      ENDDO
      GOTO 10
		
*-----grid option 5: user defined

 5    CONTINUE

* insert your grid values here:
* specify:
*  nz = total number of altutudes
* Table:  z(iz), where iz goes from 1 to nz
* trivial example of 2-layer (3-altitudes) shown below, user should modify
      WRITE(*,*) 'user-defined grid, named...'
      nz = 3
      z(1) = 0.
      z(2) = 10.
      z(3) = 80.

      GOTO 10

*-----end of user options.
*------------------------------------------------

 10   CONTINUE      

* Insert additional altitude for selected outputs.

      DO i = 1, nz
         IF(ABS(z(i) - zout) .LT. 0.001) THEN
            izout = i
            GO TO 24
         ENDIF
      ENDDO

* locate index for new altitude

      izout = 0
      DO i = 1, nz
         IF(z(i) .GT. zout) THEN
            izout = i
            GO TO 22
         ENDIF
      ENDDO
 22   CONTINUE
      IF(izout .LE. 1) STOP 'zout not in range - '

* shift overlying levels and insert new point

      nz = nz + 1
      DO i = nz, izout + 1, -1
         z(i) = z(i-1)
      ENDDO
      z(izout) = zout

 24   CONTINUE

* check grid for assorted improprieties:

 99	CONTINUE
      CALL gridck(kz,nz,z,ok)

      IF (.NOT. ok) THEN
         WRITE(*,*)'STOP in GRIDZ:  The z-grid does not make sense'
         STOP
      ENDIF

*_______________________________________________________________________

      RETURN
      END

*=============================================================================*
      SUBROUTINE  buildz(zincr,nlev,n,z)
*-----------------------------------------------------------------------------*
*= Purpose: to construct the altitude grid from parameters in gridz          =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      REAL zincr                   ! i
      INTEGER nlev                 ! i
      INTEGER n                    ! i/o
      REAL z(kz)                   ! i/o
      INTEGER i, j                    ! internal

      j = 0
      DO i = n + 1, n + nlev
        j = j + 1
        z(i) = z(n) + FLOAT(j)*zincr
      ENDDO
      n = n + nlev

      RETURN
      END

*=============================================================================*

      SUBROUTINE gridt(lat, lon, tmzone,
     $     iyear, imonth, iday,
     $     lzenit, tstart, tstop,
     $     nt, t, sza, esrm2)

*-----------------------------------------------------------------------------*
*=  Subroutine to create time (or solar zenith angle) grid                   =*
*=  Also computes earth-sun distance (1/R**2) correction.                    =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'params'

* INPUTS

      REAL lat, lon, tmzone
      INTEGER iyear, imonth, iday
      LOGICAL lzenit
      INTEGER nt
      REAL tstart, tstop


* OUTPUTS

      REAL t(kt), sza(kt), esrm2(kt)

* INTERNAL

      INTEGER it
      REAL ut, dt

      INTEGER jday, nday
      LOGICAL oky, okm, okd

      REAL az, el, soldia, soldst

*  switch for refraction correction to solar zenith angle. Because
* this is only for the observed sza at the surface, do not use.

      LOGICAL lrefr
      DATA lrefr /.FALSE./

***************

      IF(nt .EQ. 1) THEN
         dt = 0.
      ELSE
         dt = (tstop - tstart) / FLOAT(nt - 1)
      ENDIF

      DO 10 it = 1, nt
         t(it) = tstart + dt * FLOAT(it - 1)

* solar zenith angle calculation:
*  If lzenit = .TRUE., use selected solar zenith angles, also
*  set Earth-Sun distance to 1 AU.

         IF (lzenit) THEN
            sza(it) = t(it)
            esrm2(it) = 1.

*  If lzenit = .FALSE., compute solar zenith angle for specified
* location, date, time of day.  Assume no refraction (lrefr = .FALSE.)
*  Also calculate corresponding
* Earth-Sun correcton factor. 

         ELSE
            CALL calend(iyear, imonth, iday,
     $           jday, nday, oky, okm, okd)
            IF( oky .AND. okm .AND. okd) THEN

               ut = t(it) - tmzone
               CALL sunae(iyear, jday, ut, lat, lon, lrefr,
     &              az, el, soldia, soldst )
               sza(it) = 90. - el
               esrm2(it) = 1./(soldst*soldst)
            ELSE
               WRITE(*,*) '**** incorrect date specification'
               STOP ' in gridt '
            ENDIF
         ENDIF
            
 10   CONTINUE
      RETURN
      END

*=============================================================================*

      SUBROUTINE gridck(k,n,x,ok)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Check a grid X for various improperties.  The values in X have to comply =*
*=  with the following rules:                                                =*
*=  1) Number of actual points cannot exceed declared length of X            =*
*=  2) Number of actual points has to be greater than or equal to 2          =*
*=  3) X-values must be non-negative                                         =*
*=  4) X-values must be unique                                               =*
*=  5) X-values must be in ascending order                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  K  - INTEGER, length of X as declared in the calling program          (I)=*
*=  N  - INTEGER, number of actual points in X                            (I)=*
*=  X  - REAL, vector (grid) to be checked                                (I)=*
*=  OK - LOGICAL, .TRUE. -> X agrees with rules 1)-5)                     (O)=*
*=                .FALSE.-> X violates at least one of 1)-5)                 =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'params'

* input:
      INTEGER k, n
      REAL x(k)

* output:
      LOGICAL ok

* local:
      INTEGER i
*_______________________________________________________________________

      ok = .TRUE.

* check if dimension meaningful and within bounds

      IF (n .GT. k) THEN
         ok = .false.
         WRITE(kout,100)
         RETURN
      ENDIF         
  100 FORMAT('Number of data exceeds dimension')

      IF (n .LT. 2) THEN
         ok = .FALSE.
         WRITE(kout,101)
         RETURN
      ENDIF
  101 FORMAT('Too few data, number of data points must be >= 2')

* disallow negative grid values

      IF(x(1) .LT. 0.) THEN
         ok = .FALSE.
         WRITE(kout,105)
         RETURN
      ENDIF
  105 FORMAT('Grid cannot start below zero')

* check sorting

      DO 10, i = 2, n
         IF( x(i) .LE. x(i-1)) THEN
            ok = .FALSE.
            WRITE(kout,110)
            RETURN
         ENDIF
   10 CONTINUE
  110 FORMAT('Grid is not sorted or contains multiple values')
*_______________________________________________________________________

      RETURN
      END
