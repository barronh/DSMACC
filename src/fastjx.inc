#INLINE F90_GLOBAL
    REAL(kind=dp) :: JSCALE = 1.0d0
#ENDINLINE F90_GLOBAL

#INLINE F90_RATES
    REAL(kind=dp) FUNCTION J(JIDX)
       IMPLICIT NONE
       REAL*8 XLNG,YLAT
       REAL*8 OLDLON, OLDLAT, OLDPRESS,OLDTEMP,OLDJDAY
       LOGICAL RECALCULATE
       integer, parameter :: CLDFLAG=1
       integer :: JIDX,NJ
       integer, parameter :: DEBUG = 0
       integer, parameter :: JVN_ = 2000 ! see cmn_fjx_mod.F90
       real*8, dimension(JVN_) :: JFACTORS,JVALUES
       integer, dimension(JVN_) :: JINDS
       character*6, dimension(JVN_) :: JVMAPS
       character*50, dimension(JVN_) :: JLABELS
       
       SAVE OLDPRESS, OLDTEMP, OLDLON, OLDLAT, OLDJDAY, RECALCULATE
       RECALCULATE = (ABS(PRESS/OLDPRESS - 1) > 0.05) .OR. &
                     (ABS(TEMP /OLDTEMP - 1) > 0.05) .OR. &
                     (ABS(LON - OLDLON) > 5.) .OR. &
                     (ABS(JDAY_GMT - OLDJDAY) > .01)
                     
       IF (RECALCULATE) THEN
           XLNG = LON
           YLAT = LAT
           CALL DSMACC_FAST_JX(XLNG,YLAT,JDAY_GMT,PRESS,TEMP,CLDFLAG,JLABELS,JVMAPS,JFACTORS,JVALUES, JINDS, NJ,DEBUG)
           OLDPRESS = PRESS
           OLDTEMP = TEMP
           OLDLON = LON
           OLDLAT = LAT
           OLDJDAY = JDAY_GMT
       ENDIF
       IF ((JIDX .gt. NJ) .or. (JIDX .eq. 0) .or. (JINDS(JIDX).eq.0)) THEN
        write(*,*) 'ERROR', JIDX, ' greater than available ', NJ
        STOP
       ENDIF
       ! For debugging FAST-JX
       J = JFACTORS(JIDX)*JVALUES(JINDS(JIDX))*JSCALE
       ! WRITE(*,'(I2,x,A30,x,A6,x,es8.2,x,es8.2,x,es8.2)') JIDX,JLABELS(JIDX),JVMAPS(JIDX),JFACTORS(JIDX),JVALUES(JINDS(JIDX)),J
       RETURN 
  END FUNCTION J
#ENDINLINE
