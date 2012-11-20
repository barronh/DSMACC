#INLINE F90_RATES
  FUNCTION GREGDATETIME(JDAY)
      IMPLICIT NONE
      REAL(dp), DIMENSION(4) :: GREGDATETIME
      REAL(dp), INTENT(IN) :: JDAY
      INTEGER IYEAR, IMTH, IDAY
      REAL(dp) GMT
      INTEGER      IMN( 12 )          ! number of days in each month
      DATA IMN /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      SAVE IMN
      IYEAR = INT(JDAY / 1000.)
      IDAY = INT(MOD(JDAY, 1000.)) ! Temporary day
      GMT = MOD(JDAY, 1.) * 24.
      IF ( MOD( IYEAR, 4 ) .EQ. 0 ) THEN
         IMN( 2 ) = 29
      ELSE
         IMN( 2 ) = 28
      END IF
      print *, IYEAR, IDAY, GMT
      IMTH = 1
      DO WHILE (SUM(IMN(1:IMTH)) .LT. IDAY)
        IMTH = IMTH + 1
      ENDDO
      IDAY = IDAY - SUM(IMN(1:IMTH-1))
      print *, IYEAR, IMTH, IDAY, GMT
      GREGDATETIME(1) = DBLE(IYEAR)      
      GREGDATETIME(2) = DBLE(IMTH)      
      GREGDATETIME(3) = DBLE(IDAY)      
      GREGDATETIME(4) = GMT
      
  END FUNCTION GREGDATETIME

  REAL(kind=dp) FUNCTION ZENITH ( LAT, LONG, JDAY )

!***********************************************************************
!   Portions of Models-3/CMAQ software were developed or based on      *
!   information from various groups: Federal Government employees,     *
!   contractors working on a United States Government contract, and    *
!   non-Federal sources (including research institutions).  These      *
!   research institutions have given the Government permission to      *
!   use, prepare derivative works, and distribute copies of their      *
!   work in Models-3/CMAQ to the public and to permit others to do     *
!   so.  EPA therefore grants similar permissions for use of the       *
!   Models-3/CMAQ software, but users are requested to provide copies  *
!   of derivative works to the Government without restrictions as to   *
!   use by others.  Users are responsible for acquiring their own      *
!   copies of commercial software associated with Models-3/CMAQ and    *
!   for complying with vendor requirements.  Software copyrights by    *
!   the MCNC Environmental Modeling Center are used with their         *
!   permissions subject to the above restrictions.                     *
!***********************************************************************


! RCS file, release, date & time of last delta, author, state, [and locker]
! $Header: /project/work/rep//JPROC/src/driver/jproc_table/calczen.F,v 1.4 2002/04/15 18:00:43 yoj Exp $ 

! what(1) key, module and SID; SCCS file; date and time of last delta:
! @(#)calczen.F	1.2 /project/mod3/JPROC/src/driver/jproc_table/SCCS/s.calczen.F 04 Jun 1997 10:48:01

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        
!*********************************************************************
!
!  this subroutine calculates solar zenith angle for a
!     particulat time and location.  Must specify:
!  INPUT:
!       LAT - latitude in decimal degrees
!       LONG - longitude in decimal degrees (West is positive; East is negative)
!       JDAY - Fractional days at Greenwich - specify year (yyyy), 
!              day (jjj), time (decimal) day format is six-digit integer:  yyyyjjj
!  OUTPUT
!       Zenith
!
!*********************************************************************

      IMPLICIT NONE

!      INCLUDE SUBST_CONST        ! commonly used constants
      REAL*8      PI ! pi (single precision 3.141593)
      PARAMETER ( PI = 3.14159265358979324 )
 
      REAL        PI180 ! pi/180 [ rad/deg ]
      PARAMETER ( PI180  = PI / 180.0 )


!...........ARGUMENTS and their descriptions:
      REAL*8 JDAY

      REAL         LAT                ! latitude (decimal degrees)
      REAL         LONG               ! longitude (decimal degrees)
      REAL         GMT                ! Greenwich mean time (dec.milt)
!      REAL         ZENITH             ! zenith angle (degrees)
            
!...........LOCAL VARIABLES and their descriptions:

      INTEGER      IMN( 12 )          ! number of days in each month
      DATA IMN /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      SAVE         IMN

      INTEGER      IIYEAR             ! year (yyyy)
      INTEGER      LEAP               ! # leap years since 1974
      INTEGER      NYEARS             ! years since 1974
      INTEGER      NOLEAP             ! # nonleap years since 1974
      INTEGER      IJD                ! julian days (ddd)
      INTEGER      IN                 ! month pointer (mm-1)
      INTEGER      IMTH               ! month (mm)
      INTEGER      I                  ! month index
      INTEGER      IDAY               ! day (dd)
      INTEGER      JD                 ! julian days since 1974 ref
        
      REAL         LBGMT              ! 
      REAL         LZGMT              ! 
      REAL         ML                 ! geometric mean longitude (deg)
      REAL         RLT                ! latitude (radians)
      REAL         YREF               ! number of days to 1974
      REAL         YR                 ! number of days to yyyy
      REAL         D                  ! jul. days since 1974 + hh frac
      REAL         RML                ! geometric mean long (rad)
      REAL         W                  ! mean longitude of perigee (deg)
      REAL         EC                 ! eccentricity 
      REAL         EPSI               ! mean obliquity of ecliptic (deg)
      REAL         YT                 ! tan**2 (mean obl. of eclipt.)
      REAL         PEPSI              ! mean obliquity of ecliptic (rad)
      REAL         CW                 ! cosine mean long. perigee
      REAL         WR                 ! mean longitude of perigee (rad)
      REAL         SW                 ! sine mean long. perigee
      REAL         SSW                ! sine 2*mean long. perigee
      REAL         EYT                !
      REAL         FEQT               !
      REAL         FEQT1              !
      REAL         FEQT2              !
      REAL         FEQT3              !
      REAL         FEQT4              !
      REAL         FEQT5              !
      REAL         FEQT6              !
      REAL         FEQT7              !
      REAL         REQT               !
      REAL         EQT                !
      REAL         RA                 ! right ascension (degrees)
      REAL         TAB                !
      REAL         RRA                ! right ascension (radians)
      REAL         RDECL              ! declination angle (rad)
      REAL         CSZ                ! cosine (zenith angle)
      REAL         ZPT                ! 
      REAL         ZR                 ! zenith angle (radians)
      
!*********************************************************************
!     begin body of subroutine CALZEN2        

!...convert to radians

      RLT = LAT * PI180

!...parse year and time
      IIYEAR = int(JDAY / 1000)
      IDAY = INT(MOD(JDAY, DBLE(IIYEAR))) ! Temporary Julian day (1-365)
      GMT = MOD(JDAY, 1.) * 24.

!...identify and correct leap years

      IF ( MOD( IIYEAR, 4 ) .EQ. 0 ) THEN
         IMN( 2 ) = 29
      ELSE
         IMN( 2 ) = 28
      END IF

!...parse month and gregorian day
      IMTH = 1
      DO WHILE (SUM(IMN(1:IMTH+1)) .LT. IDAY)
        IMTH = IMTH + 1
      ENDDO

      IDAY = IDAY - SUM(IMN(1:IMTH)) ! Remove month days      

!...count days from Dec.31,1973 to Jan 1, YEAR, then add to 2,442,047.5

      YREF =  2442047.5
      NYEARS = IIYEAR - 1974
      LEAP = ( NYEARS + 1 ) / 4
      IF ( NYEARS .LE. -1 ) LEAP = ( NYEARS - 2 ) / 4
      NOLEAP = NYEARS - LEAP
      YR = YREF + 365.0 * NOLEAP + 366.0 * LEAP

      IJD = 0
      IN = IMTH - 1
      
      IF ( IN .EQ. 0 ) THEN
        IJD = IDAY
      ELSE IF ( IN .GT. 0 ) THEN      
        DO 101 I = 1, IN
          IJD = IJD + IMN( I )
101     CONTINUE
        IJD = IJD + IDAY
      END IF

!...print julian days current "ijd"

      JD = IJD + ( YR - YREF )
      D = JD + GMT / 24.0

!...calc geom mean longitude

      ML = 279.2801988 + 0.9856473354 * D + 2.267E-13 * D * D
      RML = ML * PI180

!...calc equation of time in sec
!... w = mean long of perigee
!... e = eccentricity
!... epsi = mean obliquity of ecliptic

      W = 282.4932328 + 4.70684E-5 * D + 3.39E-13 * D * D
      WR = W * PI180
      EC = 1.6720041E-2 - 1.1444E-9 * D - 9.4E-17 * D * D
      EPSI = 23.44266511 - 3.5626E-7 * D - 1.23E-15 * D * D
      PEPSI = EPSI * PI180
      YT = ( TAN( PEPSI / 2.0 ) )**2
      CW = COS( WR )
      SW = SIN( WR )
      SSW = SIN( 2.0 * WR )
      EYT = 2.0 * EC * YT
      FEQT1 = SIN( RML ) * ( -EYT * CW - 2.0 * EC * CW )
      FEQT2 = COS( RML ) * ( 2.0 * EC * SW - EYT * SW )
      FEQT3 = SIN( 2.0 * RML ) * ( YT - ( 5.0 * EC**2 / 4.0 ) &
     &      * ( CW**2 - SW**2 ) )
      FEQT4 = COS( 2.0 * RML ) * ( 5.0 * EC**2 * SSW / 4.0 )
      FEQT5 = SIN( 3.0 * RML ) * ( EYT * CW )
      FEQT6 = COS( 3.0 * RML ) * ( -EYT * SW )
      FEQT7 = -SIN( 4.0 * RML ) * ( 0.5 * YT**2 )
      FEQT = FEQT1 + FEQT2 + FEQT3 + FEQT4 + FEQT5 + FEQT6 + FEQT7
      EQT = FEQT * 13751.0

!...convert eq of time from sec to deg

      REQT = EQT / 240.0

!...calc right ascension in rads

      RA = ML - REQT
      RRA = RA * PI180

!...calc declination in rads, deg

      TAB = 0.43360 * SIN( RRA )
      RDECL = ATAN( TAB )

!...calc local hour angle

      LBGMT = 12.0 - EQT / 3600.0 + LONG * 24.0 / 360.0
      LZGMT = 15.0 * ( GMT - LBGMT )
      ZPT = LZGMT * PI180
      CSZ = SIN( RLT ) * SIN( RDECL ) + COS( RLT ) * COS( RDECL ) &
     &    * COS( ZPT )
      ZR = ACOS( CSZ )
      ZENITH = DBLE(ZR / PI180)
     
      RETURN
  END FUNCTION ZENITH

  SUBROUTINE UPDATE_TUV() 
     REAL(kind=dp) alta
     INTEGER :: IJDAY, OLD_JDAY = -9999 ! JDAY for use in photolysis
     REAL(dp) :: OLD_o3col = -9999 ! Ozone Column (Dobson Units)
     REAL(dp) :: OLD_press = -9999 ! Pressure (hPa)
     REAL(dp) :: OLD_TEMP = -9999 ! Temperature (K)
     REAL(dp) :: OLD_lat = -9999, OLD_lon = -9999 ! Latitude and Longitude (Decimal degrees)
     REAL(dp) :: OLD_M = -9999! Constants for concentrations that are not in the mechanism (molecule cm**-3); provided by input file
     REAL(dp), DIMENSION(4) :: tdate
     SAVE OLD_M, OLD_LAT, OLD_LON, OLD_PRESS, OLD_O3COL, OLD_JDAY, OLD_TEMP
     !    Set up the photolysis rates
     !    First calculate pressure altitude from altitude
     IJDAY = INT(JDAY)
     IF ((OLD_TEMP .eq. -9999._dp) .or. &
        &(ABS(TEMP - OLD_TEMP) .gt. 1.)) THEN
     IF (OLD_M .ne. M .or. &
  &        OLD_LAT .ne. LAT .or. &
  &        OLD_LON .ne. LON .or. &
  &        OLD_PRESS .ne. PRESS .or. &
  &        OLD_TEMP .ne. TEMP .or. &
  &        OLD_O3COL .ne. O3COL .or. &
  &        OLD_JDAY .ne. IJDAY) THEN
         WRITE (6,*) 'hvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhv'
         WRITE (6,*) 'Using TUV to calculate photolysis rates as a function of O3 column, SZA, ALT, and TEMP'
         alta=max(0., (1-(press/1013.25)**0.190263)*288.15/0.00198122*0.304800/1000.)
     
         if (o3col .eq. 0) then 
           o3col=260.
           write (6,*) 'Ozone column not specified using 260 Dobsons'
         else
           write (6,*) 'Ozone column =', o3col
         endif
           
    !    Calculate the photolysis rates for the run
         tdate = GREGDATETIME(JDAY)
         write(6,*) int(tdate(1)), int(tdate(2)), int(tdate(3))
         call set_up_photol(int(tdate(1)), int(tdate(2)), int(tdate(3)),&
                            O3col, alta, lat, lon, temp, M, bs,cs,ds,szas,svj_tj)
         WRITE (6,*) 'Photolysis rates calculated'
         WRITE (6,*) 'hvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhv'
        OLD_M = M
        OLD_LAT = LAT
        OLD_LON = LON
        OLD_PRESS = PRESS
        OLD_TEMP = TEMP
        OLD_O3COL = O3COL
        OLD_JDAY = IJDAY
     endif
     ENDIF
  END SUBROUTINE UPDATE_TUV
  REAL(kind=dp) FUNCTION J(IJ)
    USE dsmacc_Global, ONLY : THETA
    TYPE JMAP
       CHARACTER(50) :: label 
       INTEGER       :: tuvj
       INTEGER       :: mcmj
       REAL(dp)          :: factor
    END TYPE JMAP
    INTEGER i, IJ, mcmi
    INTEGER, PARAMETER :: NMCMJ = 42
!      JMAP(label                      , tuvj, mcmj, fac)
    TYPE( JMAP ), SAVE :: JMAPPING( NMCMJ ) = (/&
     & JMAP('O3->O1D                                           ',    2,    1, 1.00),&
     & JMAP('O3->O3P                                           ',    3,    2, 1.00),&
     & JMAP('H2O2->2*OH                                        ',    5,    3, 1.00),&
     & JMAP('NO2->NO+O3P                                       ',    6,    4, 1.00),&
     & JMAP('NO3->NO+O2                                        ',    7,    5, 1.00),&
     & JMAP('NO3->NO2+O3P                                      ',    8,    6, 1.00),&
     & JMAP('HNO2->OH+NO                                       ',   12,    7, 1.00),&
     & JMAP('HNO3->NO2+OH                                      ',   13,    8, 1.00),&
     & JMAP('CH2O -> H + HCO                                   ',   17,   11, 1.00),&
     & JMAP('CH2O -> H2 + CO                                   ',   18,   12, 1.00),&
     & JMAP('CH3CHO -> CH3 + HCO                               ',   19,   13, 1.00),&
     & JMAP('C2H5CHO -> C2H5 + HCO                             ',   22,   14, 1.00),&
     & JMAP('nC3H7CHO -> nC3H7 + HCO                           ',   84,   15, 1.00),&
     & JMAP('nC3H7CHO -> C2H4 + CH3CHO                         ',   85,   16, 1.00),&
     & JMAP('iC3H7CHO -> iC3H7 + HCO                           ',   86,   17, 1.00),&
     & JMAP('CH2=C(CH3)CHO -> Products                         ',   25,   18, 0.50),&
     & JMAP('CH2=C(CH3)CHO -> Products                         ',   25,   19, 0.50),&
     & JMAP('CH3COCH3 -> CH3CO + CH3                           ',   26,   21, 1.00),&
     & JMAP('CH3COC2H5 -> C2H5 + CH3CO                         ',   96,   22, 1.00),&
     & JMAP('CH3COCHCH2 -> Products                            ',   27,   23, 0.50),&
     & JMAP('CH3COCHCH2 -> Products                            ',   27,   24, 0.50),&
     & JMAP('CHOCHO -> 2CO + H2                                ',   44,   31, 1.00),&
     & JMAP('CHOCHO -> CH2O + CO                               ',   46,   32, 1.00),&
     & JMAP('CHOCHO -> HCO + HCO                               ',   45,   33, 1.00),&
     & JMAP('CH3COCHO -> CH3CO + HCO                           ',   47,   34, 1.00),&
     & JMAP('CH3COCOCH3 -> Products                            ',   48,   35, 1.00),&
     & JMAP('ClONO2 -> ClO + NO2                               ',   54,   41, 1.00),&
     & JMAP('CH3ONO2 -> CH3O + NO2                             ',   34,   51, 1.00),&
     & JMAP('CH3CH2ONO2 -> CH3CH2O + NO2                       ',   35,   52, 1.00),&
     & JMAP('nC3H7ONO2 -> nC3H7O + NO2                         ',  100,   53, 1.00),&
     & JMAP('CH3CH2ONO2 -> CH3CH2O + NO2                       ',   35,   54, 1.00),&
     & JMAP('2-C4H9ONO2 -> 2-C4H9O + NO2                       ',  103,   55, 1.00),&
     & JMAP('CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2              ',   38,   56, 0.75),&
     & JMAP('CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2              ',   38,   57, 0.25),&
     & JMAP('HOBr                                              ',   57, 1001, 1.00),&
     & JMAP('BrO                                               ',   56, 1002, 1.00),&
     & JMAP('Br2                                               ',   55, 1003, 1.00),&
     & JMAP('BrONO2 -> Br + NO3                                ',   59, 1004, 1.00),&
     & JMAP('BrONO2 -> BrO + NO2                               ',   58, 1005, 1.00),&
     & JMAP('ClONO2 -> Cl + NO3                                ',   53, 1006, 1.00),&
     & JMAP('ClONO2 -> ClO + NO2                               ',   54, 1007, 1.00),&
     & JMAP('Cl2 -> Cl + Cl                                    ',   50, 1008, 1.00)&
    /)
    J = 0._dp    
    DO i=1,NMCMJ
      mcmi = JMAPPING(i)%%%%mcmj
      IF (mcmi .eq. IJ) THEN
        J = TUV_J(JMAPPING(i)%%%%tuvj, THETA, JMAPPING(i)%%%%factor)
        RETURN
      ENDIF
    ENDDO
  END FUNCTION J

  REAL(kind=dp) FUNCTION TUV_J(IJ, THETA, JFACT)
    USE dsmacc_Global,  only: BS, CS, DS, SZAS, SVJ_TJ, &
  &               PRESS, TEMP, kj
    INTEGER, parameter :: NSZAS = 19
    REAL(kind=dp) B(NSZAS), C(NSZAS), D(NSZAS), TMP_SVJ_TJ(NSZAS), &
  &               TMP_SZAS(NSZAS), THETA, JFACT
  ! IJ is the integer index of the TUV photolysis calculation
  ! THETA is the current solar zenith angle
    INTEGER IJ, THIS_CSZA
    
    REAL*8 SEVAL ! function from TUV
    EXTERNAL SEVAL ! function from TUV
    INTEGER THIS_SZA ! temporary variables

    call update_tuv()
    IF (THETA .LT. 90.D0) THEN
      DO THIS_CSZA = 1, NSZAS
        B(THIS_CSZA) = BS(THIS_CSZA, IJ)
        C(THIS_CSZA) = CS(THIS_CSZA, IJ)
        D(THIS_CSZA) = DS(THIS_CSZA, IJ)
        TMP_SZAS(THIS_CSZA) = SZAS(THIS_CSZA)
        TMP_SVJ_TJ(THIS_CSZA) = SVJ_TJ(THIS_CSZA, IJ)
      ENDDO
      
      TUV_J = SEVAL(NSZAS, THETA, TMP_SZAS, TMP_SVJ_TJ, B, C, D, .false.) * JFACT
      if (.false.) then
          write(*,*) 'MP'
          write(*,*) 'I,THETA,J:', IJ, THETA, TUV_J
          write(*,8879) 'B     :', B
          write(*,8879) 'C     :', C
          write(*,8879) 'D     :', D
          write(*,8879) 'SZAS  :', TMP_SZAS
          write(*,8879) 'SVJ_TJ:', TMP_SVJ_TJ
          
8879     FORMAT(1A6,100000(E26.17))
          TUV_J = SEVAL(NSZAS, THETA, TMP_SZAS, TMP_SVJ_TJ, B, C, D, .true.)

      endif
    
      IF (TUV_J .LT. 0.d0) TUV_J = 0.d0
    ELSE
      TUV_J = 0.d0
    ENDIF

  END FUNCTION TUV_J

#ENDINLINE
