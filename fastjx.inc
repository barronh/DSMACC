#INLINE F90_RATES
    REAL*8 FUNCTION FASTJX(JIDX, THETA)
       USE dsmacc_Global, ONLY : JDAY_GMT, PRESS, TEMP
       IMPLICIT NONE
       REAL*8 XLNG,YLAT,THETA
       REAL*8 OLDLON, OLDLAT, OLDPRESS,OLDTEMP
       LOGICAL RECALCULATE
       integer, parameter :: CLDFLAG=1
       integer :: J,JIDX,NJ
       integer, parameter :: DEBUG = 1
       integer, parameter :: JVN_ = 101 ! see cmn_fjx_mod.F90
       real*8, dimension(JVN_) :: JFACTORS,JVALUES
       integer, dimension(JVN_) :: JINDS
       character*6, dimension(JVN_) :: JVMAPS
       character*50, dimension(JVN_) :: JLABELS
       
       SAVE OLDPRESS, OLDTEMP, OLDLON, OLDLAT, RECALCULATE
       RECALCULATE = (ABS(PRESS/OLDPRESS - 1) > 0.05) .OR. &
                     (ABS(TEMP /OLDTEMP - 1) > 0.05) .OR. &
                     (ABS(LON - OLDLON) > 5.) .OR. &
                     (ABS(LAT - OLDLAT - 1) > 5.)
       IF (RECALCULATE) THEN
           XLNG = LON
           YLAT = LAT
           CALL DSMACC_FAST_JX(XLNG,YLAT,JDAY_GMT,PRESS,TEMP,CLDFLAG,JLABELS,JVMAPS,JFACTORS,JVALUES, JINDS, NJ,DEBUG)
       ENDIF
       IF (JIDX > NJ) THEN
        write(*,*) 'ERROR', JIDX, ' greater than available ', NJ
        STOP
       ENDIF
       DO J=JIDX,JIDX !1,NRATJ
            WRITE(*,'(I2,x,A,x,A,x,e8.2,x,e8.2)') J,JLABELS(J),JVMAPS(J),JFACTORS(J),JVALUES(JINDS(J))
       ENDDO
       FASTJX = JFACTORS(JIDX)*JVALUES(JINDS(JIDX))
  END FUNCTION FASTJX
#ENDINLINE
