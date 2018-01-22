#INLINE F90_GLOBAL
    REAL(kind=dp) :: JSCALE = 1.0d0
    INTEGER :: NJREACT
    integer, parameter :: JVN_ = 2000 ! see cmn_fjx_mod.F90
    character*50, dimension(JVN_) :: JLABELS
#ENDINLINE F90_GLOBAL

#INLINE F90_RATES
    REAL(kind=dp) FUNCTION J(JIDX)
       IMPLICIT NONE
       REAL*8 :: PRESS_hPA, XLNG,YLAT, DPRESSFRAC, DTEMPFRAC, DLON,DJDAY,DLAT
       REAL*8 :: OLDLON = -999, OLDLAT = -999, OLDPRESS = -999
       REAL*8 :: OLDTEMP = -999,OLDJDAY = -999
       LOGICAL :: RECALCULATE
       integer, parameter :: CLDFLAG=1
       integer :: JIDX
       integer, parameter :: DEBUG = 0
       real*8, dimension(JVN_) :: JFACTORS,JVALUES
       integer, dimension(JVN_) :: JINDS
       character*6, dimension(JVN_) :: JVMAPS
       
       SAVE OLDPRESS, OLDTEMP, OLDLON, OLDLAT, OLDJDAY, RECALCULATE, JFACTORS, JVMAPS, JVALUES, JINDS
       PRESS_hPA = PRESS / 100.
       DPRESSFRAC = PRESS_hPA/OLDPRESS - 1
       DTEMPFRAC = TEMP /OLDTEMP - 1
       DLON = LON - OLDLON
       DLAT = LAT - OLDLAT
       DJDAY = JDAY_GMT - OLDJDAY
       RECALCULATE = (ABS(DPRESSFRAC) > 0.05) .OR. &
                     (ABS(DTEMPFRAC) > 0.05) .OR. &
                     (ABS(DLON) > 5.) .OR. &
                     (ABS(DLAT) > 5.) .OR. &
                     (ABS(DJDAY) > .01)
                     
       IF (RECALCULATE) THEN
           XLNG = LON
           YLAT = LAT
           IF (DEBUG>0) THEN
             print*,'DPRESSFRAC','DTEMPFRAC','DLON','DLAT','DJDAY'
             print*,DPRESSFRAC,DTEMPFRAC,DLON,DLAT,DJDAY
           ENDIF
           CALL DSMACC_FAST_JX(XLNG,YLAT,JDAY_GMT,PRESS_hPA,TEMP,&
                               CLDFLAG,JLABELS,JVMAPS,JFACTORS,&
                               JVALUES, JINDS, NJREACT,THETA,DEBUG)
           OLDPRESS = PRESS_hPA
           OLDTEMP = TEMP
           OLDLON = LON
           OLDLAT = LAT
           OLDJDAY = JDAY_GMT
       ENDIF
       IF ((JIDX .gt. NJREACT) .or. (JIDX .eq. 0)) THEN
        write(*,*) 'ERROR', JIDX, ' greater than available ', NJREACT
        STOP
       ELSEIF (JINDS(JIDX).eq.0) THEN
        write(*,*) 'ERROR', JIDX, ' not mapped; returning 0. for RXN'
        write(*,*) ' - ', JLABELS(JIDX)
        J = 0.
        RETURN
       ENDIF
       J = JFACTORS(JIDX)*JVALUES(JINDS(JIDX))*JSCALE
       IF (TRIM(ADJUSTL(JVMAPS(JIDX))) .eq. 'O3(1D)') THEN
          J = J * JFACTO1D
       ELSEIF (TRIM(ADJUSTL(JVMAPS(JIDX))) .eq. 'NO2') THEN
          J = J * JFACTNO2
       ENDIF
       ! For debugging FAST-JX
       IF (DEBUG>0) THEN
         WRITE(*,'(I2,x,A30,x,A6,x,es8.2,x,es8.2,x,es8.2)') JIDX,&
                 JLABELS(JIDX),JVMAPS(JIDX),JFACTORS(JIDX),&
                 JVALUES(JINDS(JIDX)),J
       ENDIF
       RETURN 
  END FUNCTION J
#ENDINLINE
