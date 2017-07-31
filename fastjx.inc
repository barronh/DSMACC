#INLINE F90_RATES
    REAL*8 FUNCTION FASTJX(JIDX, THETA)
       USE dsmacc_Global, ONLY : JDAY, PRESS, TEMP
       USE CMN_FJX_MOD
       IMPLICIT NONE
       REAL*8 XLNG,YLAT,THETA
       REAL*8 OLDLON, OLDLAT, OLDPRESS,OLDTEMP
       LOGICAL RECALCULATE
       integer, parameter :: CLDFLAG=1
       integer :: J,JIDX
       integer, parameter :: DEBUG = 1
       real*8, dimension(JVN_) :: JVALUES
       SAVE OLDPRESS, OLDTEMP, OLDLON, OLDLAT, RECALCULATE
       RECALCULATE = (ABS(PRESS/OLDPRESS - 1) > 0.05) .OR. &
                     (ABS(TEMP /OLDTEMP - 1) > 0.05) .OR.
                     (ABS(LON - OLDLON) > 5.) .OR.
                     (ABS(LAT - OLDLAT - 1) > 5.)
       IF (RECALCULATE) THEN
           XLNG = LON
           YLAT = LAT
           CALL DSMACC_FAST_JX(XLNG,YLAT,JDAY_GMT,PRESS,TEMP,CLDFLAG,JVALUES, DEBUG)
       ENDIF
       DO J=JIDX,JIDX !1,NRATJ
            WRITE(*,'(I2,x,A,x,A,x,e8.2,x,e8.2)') J,JLABEL(J),JVMAP(J),JFACTA(J),JVALUES(JIND(J))
       ENDDO
       FASTJX = JFACTA(JIDX)*JVALUES(JIDX)
  END FUNCTION FASTJX
#ENDINLINE
