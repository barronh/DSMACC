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
#ENDINLINE F90_RATES
