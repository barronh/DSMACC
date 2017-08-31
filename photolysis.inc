#INLINE F90_GLOBAL
    CHARACTER*50 JLABELS(kj)
    TYPE JMAP
       CHARACTER(50) :: label 
       INTEGER       :: tuvj
       INTEGER       :: mcmj
       REAL(dp)          :: factor
    END TYPE JMAP
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
     & JMAP('nC3H7CHO -> nC3H7 + HCO                           ',   85,   15, 1.00),&
     & JMAP('nC3H7CHO -> C2H4 + CH3CHO                         ',   86,   16, 1.00),&
     & JMAP('iC3H7CHO -> iC3H7 + HCO                           ',   87,   17, 1.00),&
     & JMAP('CH2=C(CH3)CHO -> Products                         ',   25,   18, 0.50),&
     & JMAP('CH2=C(CH3)CHO -> Products                         ',   25,   19, 0.50),&
     & JMAP('CH3COCH3 -> CH3CO + CH3                           ',   26,   21, 1.00),&
     & JMAP('CH3COC2H5 -> C2H5 + CH3CO                         ',   29,   22, 1.00),&
     & JMAP('CH3COCHCH2 -> Products                            ',   28,   23, 0.50),&
     & JMAP('CH3COCHCH2 -> Products                            ',   28,   24, 0.50),&
     & JMAP('CHOCHO -> 2CO + H2                                ',   45,   31, 1.00),&
     & JMAP('CHOCHO -> CH2O + CO                               ',   47,   32, 1.00),&
     & JMAP('CHOCHO -> HCO + HCO                               ',   46,   33, 1.00),&
     & JMAP('CH3COCHO -> CH3CO + HCO                           ',   48,   34, 1.00),&
     & JMAP('CH3COCOCH3 -> Products                            ',   49,   35, 1.00),&
     & JMAP('ClONO2 -> ClO + NO2                               ',   55,   41, 1.00),&
     & JMAP('CH3ONO2 -> CH3O + NO2                             ',   35,   51, 1.00),&
     & JMAP('CH3CH2ONO2 -> CH3CH2O + NO2                       ',   36,   52, 1.00),&
     & JMAP('nC3H7ONO2 -> nC3H7O + NO2                         ',  101,   53, 1.00),&
     & JMAP('CH3CH2ONO2 -> CH3CH2O + NO2                       ',   36,   54, 1.00),&
     & JMAP('2-C4H9ONO2 -> 2-C4H9O + NO2                       ',  104,   55, 1.00),&
     & JMAP('CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2              ',   39,   56, 0.75),&
     & JMAP('CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2              ',   39,   57, 0.25),&
     & JMAP('HOBr                                              ',   58, 1001, 1.00),&
     & JMAP('BrO                                               ',   57, 1002, 1.00),&
     & JMAP('Br2                                               ',   56, 1003, 1.00),&
     & JMAP('BrONO2 -> Br + NO3                                ',   60, 1004, 1.00),&
     & JMAP('BrONO2 -> BrO + NO2                               ',   59, 1005, 1.00),&
     & JMAP('ClONO2 -> Cl + NO3                                ',   54, 1006, 1.00),&
     & JMAP('ClONO2 -> ClO + NO2                               ',   55, 1007, 1.00),&
     & JMAP('Cl2 -> Cl + Cl                                    ',   51, 1008, 1.00)&
    /)
#ENDINLINE
#INLINE F90_RCONST
  THETA = ZENITH ()
#ENDINLINE F90_RCONST
#INLINE F90_RATES
  FUNCTION GREGDATETIME(JDAY_GMT)
      ! ARGUMENTS
      ! - JDAY_GMT (YYYYJJJ.F where F is fractional day and time is in GMT)
      ! RETURNS
      ! - YEAR, MONTH, DAY and TIME in seconds from 0 UTC
      IMPLICIT NONE
      LOGICAL :: DEBUG = .false.
      INTEGER, DIMENSION(4) :: GREGDATETIME
      REAL(dp), INTENT(IN) :: JDAY_GMT
      INTEGER IYEAR, IMTH, IDAY
      REAL(dp) GMT
      INTEGER      IMN( 12 )          ! number of days in each month
      DATA IMN /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      SAVE IMN
      IYEAR = INT(JDAY_GMT / 1000.)
      IDAY = INT(MOD(JDAY_GMT, 1000.)) ! Temporary day
      GMT = MOD(JDAY_GMT, 1.) * 24.
      IF ( MOD( IYEAR, 4 ) .EQ. 0 ) THEN
         IMN( 2 ) = 29
      ELSE
         IMN( 2 ) = 28
      END IF
      IF (DEBUG) THEN
      WRITE(*, *) IYEAR, IDAY, GMT
      ENDIF
      IMTH = 1
      DO WHILE (SUM(IMN(1:IMTH)) .LT. IDAY)
        IMTH = IMTH + 1
      ENDDO
      IDAY = IDAY - SUM(IMN(1:IMTH-1))
      IF (DEBUG) THEN
      write(*,*) IYEAR, IMTH, IDAY, GMT
      ENDIF
      GREGDATETIME(1) = IYEAR
      GREGDATETIME(2) = IMTH
      GREGDATETIME(3) = IDAY
      GREGDATETIME(4) = INT(GMT * 3600)
      
  END FUNCTION GREGDATETIME
  FUNCTION ZENITH()
      ! ZENITH uses libtuv to calculate
      ! Solar Zenith Angle
      ! Relies on:
      ! - LAT (degrees N)
      ! - LON (degrees E)
      ! - JDAY_GMT (YYYYJJJ.F) where F is a fractional day
      ! Returns:
      ! ZENITH - degrees from zenith
      USE dsmacc_Global, ONLY : LAT, LON, JDAY_GMT
      LOGICAL :: DEBUG = .false.
      REAL rlat, rlon, rut
      REAL az, el, soldia, soldst
      REAL*8 ZENITH
      INTEGER iyear, imonth, iday, myjday, nday
      INTEGER, DIMENSION(4) :: YYYYMMDD_TS
      LOGICAL oky, okm, okd
      LOGICAL :: lrefr = .false.
      rlat = REAL(LAT)
      rlon = REAL(LON)
      YYYYMMDD_TS = GREGDATETIME(JDAY_GMT)
      IYEAR = INT(YYYYMMDD_TS(1))
      IMONTH = INT(YYYYMMDD_TS(2))
      IDAY = INT(YYYYMMDD_TS(3))
      RUT = REAL(YYYYMMDD_TS(4) / 3600.)
      CALL calend(iyear, imonth, iday, myjday, nday, oky, okm, okd)
      IF (DEBUG) THEN
        print *,IYEAR,IMONTH, IDAY,MYJDAY,RUT,RLAT,RLON
      ENDIF
      IF( oky .AND. okm .AND. okd) THEN
         CALL sunae(iyear, myjday, rut, rlat, rlon, lrefr, az, el, soldia, soldst )
         ZENITH = DBLE(90. - el)
         RETURN
      ELSE
         WRITE(*,*) 'Bad date'
         WRITE(*,*) IYEAR, OKY
         WRITE(*,*) IMONTH, OKM
         WRITE(*,*) IDAY, OKD
      ENDIF
  END FUNCTION ZENITH

  SUBROUTINE UPDATE_TUV() 
     USE dsmacc_Global, ONLY : JLABELS
     LOGICAL :: DEBUG = .false.
     REAL(kind=dp) alta
     INTEGER :: IJDAY, OLD_JDAY = -9999 ! JDAY for use in photolysis
     REAL(dp) :: OLD_o3col = -9999 ! Ozone Column (Dobson Units)
     REAL(dp) :: OLD_press = -9999 ! Pressure (hPa)
     REAL(dp) :: OLD_TEMP = -9999 ! Temperature (K)
     REAL(dp) :: OLD_lat = -9999, OLD_lon = -9999 ! Latitude and Longitude (Decimal degrees)
     REAL(dp) :: OLD_M = -9999! Constants for concentrations that are not in the mechanism (molecule cm**-3); provided by input file
     INTEGER, DIMENSION(4) :: tdate
     SAVE OLD_M, OLD_LAT, OLD_LON, OLD_PRESS, OLD_O3COL, OLD_JDAY, OLD_TEMP
     !    Set up the photolysis rates
     !    First calculate pressure altitude from altitude
     IJDAY = INT(JDAY_GMT)
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
         tdate = GREGDATETIME(JDAY_GMT)
         call set_up_photol(tdate(1), tdate(2), tdate(3),&
                            O3col, alta, lat, lon, temp, &
                            M, bs,cs,ds,szas,svj_tj, jlabels,DEBUG)
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
    USE dsmacc_Global, ONLY : THETA, JMAPPING, NMCMJ, JLABELS
    LOGICAL, SAVE :: FIRST = .true.
    INTEGER i, IJ, mcmi, JL
    J = 0._dp    
    DO i=1,NMCMJ
      mcmi = JMAPPING(i)%%%%mcmj
      IF (mcmi .eq. IJ) THEN
        J = TUV_J(JMAPPING(i)%%%%tuvj, THETA, JMAPPING(i)%%%%factor)
        IF (FIRST) THEN
          WRITE(*,*) '* START TUV -> MCM Mapping Confirmation *'
          DO JL = 1,NMCMJ
            mcmi = JMAPPING(JL)%%%%mcmj
            write(*,*) mcmi,' MCM ', JMAPPING(JL)%%%%LABEL
            write(*,*) mcmi,' TUV ', JLABELS(JMAPPING(JL)%%%%tuvj)
          ENDDO
          WRITE(*,*) '* END   TUV -> MCM Mapping Confirmation *'
          FIRST = .FALSE.
        ENDIF
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
    LOGICAL :: DEBUG = .false.
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
      if (DEBUG) then
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
