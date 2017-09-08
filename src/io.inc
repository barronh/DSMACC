#INLINE F90_UTIL

    SUBROUTINE NewInitVal (counter )
    
      USE dsmacc_Global
      use dsmacc_Monitor
    
      INTEGER :: i,j,found, counter, STATUS, ARGLEN
      REAL(kind=dp) :: x
      CHARACTER(LEN=10000) :: line, INPATH
      CHARACTER(LEN=15) :: test
      REAL(kind=dp) :: concs(10000)
      REAL(kind=dp) :: oldvar(10000)
      LOGICAL :: SPECIAL
      INTEGER :: COUNT_NOX_CONSTRAINTS, IOSTAT
      SAVE SPECIAL, INPATH
    
    ! Open the file with the info
    
      write (6,*) 'Initializing model point', counter
      call flush(6)
    
    
      DO I=1,NVAR
         OLDVAR(I)=C(I)
         VAR(I)=0.
         C(I)=0.
      ENDDO
      SAREA=0.
      ALBEDO=0.
      RP1=0.
    
      IF (counter .eq. 0) then 
        IF (COMMAND_ARGUMENT_COUNT() == 0) THEN
           WRITE(6,*) "There were no command arguments; using Init_cons.dat"
           INPATH = 'Init_cons.dat'
        ELSE
           CALL GET_COMMAND_ARGUMENT(1, VALUE = INPATH, LENGTH = ARGLEN, STATUS = STATUS)	
           IF (STATUS .ne. 0) THEN
               STOP 'INPATH was truncated'
           ENDIF
        ENDIF
            open(UNIT=21,FILE=INPATH, IOSTAT = IOSTAT)
    
    ! set everything to zero for the first iteration
    ! after that IF doing a constrained run don't
            LINECOUNT=0
            DO WHILE (.NOT.IS_IOSTAT_END(IOSTAT))
            READ (21,'(A1000000)', IOSTAT = IOSTAT) line
            IF ((TRIM(ADJUSTL(line)) .ne. '') .and. .NOT.IS_IOSTAT_END(IOSTAT)) THEN
                LINECOUNT=LINECOUNT+1
            ENDIF
            ENDDO
            
         CLOSE(21)
         WRITE (6,*) 'Input file has ',LINECOUNT,' lines'
      ENDIF
    
      IF (COUNTER .EQ. 1) THEN 
        open(UNIT=21,FILE=INPATH)
      ENDIF    
      time=tstart
    ! Set everything to zero
      
      IF (counter .eq. 1) then 
         SPEC_CH4=.FALSE.
         SPEC_H2=.FALSE.
         READ(21,'(i10)') IntTime
         CONSTRAIN_RUN=.FALSE.
    
         IF (INTTIME .LE. 0) THEN 
            SPECIAL=.TRUE.
         ENDIF
    
         IF (INTTIME .EQ. -1) THEN 
            WRITE (6,*) 'Integration to convergence'
            CONSTRAIN_RUN=.TRUE.
            INTTIME=50*24.*60.*60.
            OUTPUT_LAST=.FALSE.
            SPECIAL=.FALSE.
         ENDIF
    
         IF (INTTIME .EQ. -2) THEN 
            WRITE (6,*) 'Integration to convergence'
            CONSTRAIN_RUN=.TRUE.
            INTTIME=50*24.*60.*60.
            OUTPUT_LAST=.TRUE.
            SPECIAL=.FALSE.
         ENDIF
    
         IF (SPECIAL .EQV. .TRUE.) THEN
            WRITE (6,*) 'Negative Integration Time', INTTIME
            WRITE (6,*) 'But not a special case', SPECIAL
            STOP
         ENDIF
    
         
         READ(21,'(10000(a15,x))') spec_name
         READ(21,'(10000(i15,x))') const_method     
      ENDIF
    
      IF (COUNTER .NE. 0) THEN 
    
      READ (21,'(10000(e15.4,x))', IOSTAT = IOSTAT) concs
      IF (IS_IOSTAT_END(IOSTAT)) THEN
         LAST_POINT=.TRUE.
      ENDIF
      DO I=1,10000
         FOUND=0
         IF (SPEC_NAME(I) .NE. '') THEN
            
            DO J=1,NVAR 
               
               SPEC_NAME(I)=ADJUSTL(SPEC_NAME(I))
               
               IF (TRIM(SPEC_NAME(I)) .EQ. TRIM(SPC_NAMES(J))) THEN 
                  FOUND=1
                  VAR(J)=CONCS(I)
                  IF (const_method(i) .EQ. 1) CONSTRAIN(J)=CONCS(I)
                  IF (SPEC_NAME(I) .EQ. 'CH4') SPEC_CH4=.TRUE.
                  IF (SPEC_NAME(I) .EQ. 'H2') SPEC_H2=.TRUE.
               ENDIF
            ENDDO
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'H2O') THEN
               FOUND=1
               H2O=CONCS(I)
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'PRESS_hPa') THEN 
               FOUND=1
               PRESS=CONCS(I)
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'NOx') THEN
               FOUND=1
               CONSTRAIN_NOX=.TRUE.
               WRITE (6,*) 'Constraining total NOx concentration'
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'LAT_degN') THEN
               FOUND=1
               LAT=CONCS(I)
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'LON_degE') THEN 
               FOUND=1
               LON=CONCS(I)
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'TEMP_K') THEN
               FOUND=1
               TEMP=CONCS(I)
            ENDIF
            
            IF (TRIM(SPEC_NAME(I)) .EQ. 'JDAY_GMT') THEN 
               FOUND=1
               JDAY_GMT=CONCS(I)
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'TIME(GMTs)') THEN
               FOUND=1
               
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'O3COL_DU') THEN
               FOUND=1
               O3COL=CONCS(I)
            ENDIF
    
        IF (TRIM(SPEC_NAME(I)) .EQ. 'ALBEDO') THEN 
           FOUND=1
               ALBEDO=CONCS(I)
        ENDIF
    
        IF (TRIM(SPEC_NAME(I)) .EQ. 'SAREA') THEN 
           FOUND=1
           SAREA=CONCS(I)
        ENDIF
            
        IF (TRIM(SPEC_NAME(I)) .EQ. 'RP1') THEN 
           FOUND=1
           RP1=CONCS(I)
        ENDIF
    
        IF (TRIM(SPEC_NAME(I)) .EQ. 'JNO2_per_s') THEN
           FOUND=1
               IF (CONST_METHOD(I) .GE. 1) THEN 
                   JNO2=CONCS(I)
           ENDIF 
               JREPEAT=0
           IF (CONST_METHOD(I) .EQ. 2) THEN 
            JREPEAT=1
           ENDIF        
        ENDIF
            
        IF (TRIM(SPEC_NAME(I)) .EQ. 'JO1D_per_s') THEN
           FOUND=1
           IF (CONST_METHOD(I) .GE. 1) THEN 
               JO1D=CONCS(I)
           ENDIF
           JREPEAT=0
           IF (CONST_METHOD(I) .EQ. 2) THEN 
               JREPEAT=1
               ENDIF    
            ENDIF
         ENDIF
          IF (TRIM(SPEC_NAME(I)) .NE. '' .AND. FOUND .EQ. 0) THEN
             WRITE (6,*) SPEC_NAME(I),' NOT FOUND'
             IF (SPEC_NAME(I)(1:1) .NE. 'X') STOP
             WRITE (6,*) 'Starts with an X so ignored and continued'   
          ENDIF
      ENDDO
         
      CFACTOR=PRESS*1e2*1e-6/(8.314*TEMP)*6.022E23
      
      H2O=H2O*CFACTOR 
      DO I=1,NVAR
         VAR(I)=VAR(I)*CFACTOR
    !     IF (OLDVAR(I) .NE. CONSTRAIN(I)) THEN 
          CONSTRAIN(I)=CONSTRAIN(I)*CFACTOR
    !  ENDIF
      ENDDO
    
     DO I=1,NVAR
    !    C(I)=0.
    !     IF (VAR(I) .EQ. 0) VAR(I)=OLDVAR(I)  
      ENDDO
    
    ! FIND NOX species
      IF (CONSTRAIN_NOX) THEN 
      COUNT_NOX_CONSTRAINTS=0
      DO I=1,NVAR 
         IF (TRIM(SPC_NAMES(I)) .EQ. 'NO2') THEN
            NOX(I)=1
            CONSTRAIN(I)=-1.*CONSTRAIN(I)
            IF (CONSTRAIN(I) .NE. 0) THEN 
               COUNT_NOX_CONSTRAINTS=COUNT_NOX_CONSTRAINTS+1
            ENDIF
         ENDIF 
    
         IF (TRIM(SPC_NAMES(I)) .EQ. 'NO') THEN
            NOX(I)=1
            CONSTRAIN(I)=-1.*CONSTRAIN(I)
            IF (CONSTRAIN(I) .NE. 0) THEN 
               COUNT_NOX_CONSTRAINTS=COUNT_NOX_CONSTRAINTS+1
            ENDIF
         ENDIF
    
         IF (TRIM(SPC_NAMES(I)) .EQ. 'NO3') THEN 
            NOX(I)=1
            CONSTRAIN(I)=-1.*CONSTRAIN(I)
              IF (CONSTRAIN(I) .NE. 0) THEN 
               COUNT_NOX_CONSTRAINTS=COUNT_NOX_CONSTRAINTS+1
            ENDIF
         ENDIF
    
         IF (TRIM(SPC_NAMES(I)) .EQ. 'N2O5') THEN
            NOX(I)=2
            CONSTRAIN(I)=-1.*CONSTRAIN(I)
              IF (CONSTRAIN(I) .NE. 0) THEN 
               COUNT_NOX_CONSTRAINTS=COUNT_NOX_CONSTRAINTS+1
            ENDIF
         ENDIF
    
         IF (TRIM(SPC_NAMES(I)) .EQ. 'HONO') THEN
            NOX(I)=1
            CONSTRAIN(I)=-1.*CONSTRAIN(I)
              IF (CONSTRAIN(I) .NE. 0) THEN 
               COUNT_NOX_CONSTRAINTS=COUNT_NOX_CONSTRAINTS+1
            ENDIF
         ENDIF
    
         IF (TRIM(SPC_NAMES(I)) .EQ. 'HNO2') THEN
            NOX(I)=1
            CONSTRAIN(I)=-1.*CONSTRAIN(I)
              IF (CONSTRAIN(I) .NE. 0) THEN 
               COUNT_NOX_CONSTRAINTS=COUNT_NOX_CONSTRAINTS+1
            ENDIF
         ENDIF
    
         IF (TRIM(SPC_NAMES(I)) .EQ. 'HO2NO2') THEN
            NOX(I)=1
            CONSTRAIN(I)=-1.*CONSTRAIN(I)
              IF (CONSTRAIN(I) .NE. 0) THEN 
               COUNT_NOX_CONSTRAINTS=COUNT_NOX_CONSTRAINTS+1
            ENDIF
         ENDIF
    
         IF (TRIM(SPC_NAMES(I)) .EQ. 'PNA') THEN 
            NOX(I)=1
            CONSTRAIN(I)=-1*CONSTRAIN(I)
         ENDIF
    
         IF (TRIM(SPC_NAMES(I)) .EQ. 'HNO4') THEN 
            NOX(I)=1
            CONSTRAIN(I)=-1.*CONSTRAIN(I)
            IF (CONSTRAIN(I) .NE. 0) THEN 
               COUNT_NOX_CONSTRAINTS=COUNT_NOX_CONSTRAINTS+1
            ENDIF
         ENDIF
    
         IF (NOX(I) .NE. 0) WRITE (6,*) SPC_NAMES(I),' IN NOX FAMILY'
      ENDDO
      
      IF (COUNT_NOX_CONSTRAINTS .GT. 1) THEN 
         WRITE (6,*) 'You can not contrains multiple NOX species'
         STOP
      ENDIF
      ENDIF
    
    ! FIND CH4 and H2 species
      DO I=1,NVAR 
         IF (TRIM(SPC_NAMES(I)) .EQ. 'CH4' .AND. (SPEC_CH4 .EQV. .FALSE.)) THEN
            WRITE (6,*) 'No CH4 specified assuming 1770 ppbv'
            VAR(I)=1770e-9*CFACTOR
            CONSTRAIN(I)=VAR(I)
         ENDIF
    
         IF (TRIM(SPC_NAMES(I)) .EQ. 'H2' .AND. (SPEC_H2 .EQV. .FALSE.)) THEN 
            WRITE (6,*) 'No H2 specified assuming 550 ppbv'
            VAR(I)= 550e-9*CFACTOR
            CONSTRAIN(I)=VAR(I)
         ENDIF
      ENDDO
    ! INLINED initialisations
    
    ! END INLINED initialisations
    
      ENDIF
          
    END SUBROUTINE NewInitval

! ****************************************************************
!                            
! InitSaveData - Opens the data file for writing
!   Parameters :                                                  
!
! ****************************************************************

    SUBROUTINE NewInitSaveData(counter)
        USE dsmacc_Parameters
        USE dsmacc_monitor
        USE dsmacc_GLOBAL
        INTEGER counter, length, ierr, STATUS, ARGLEN
        CHARACTER*20 filename1
        CHARACTER(LEN=10000) :: SPECROOT, RATEROOT, filename
        INTEGER i
        IF (COMMAND_ARGUMENT_COUNT() < 2) THEN
           WRITE(6,*) "There were fewer than 2 command arguments; using Spec_"
           SPECROOT = 'Spec_'
        ELSE
           CALL GET_COMMAND_ARGUMENT(2, VALUE = SPECROOT, LENGTH = ARGLEN, STATUS = STATUS)	
           IF (STATUS .ne. 0) THEN
               STOP 'SPECROOT was truncated'
           ENDIF
        ENDIF
        IF (COMMAND_ARGUMENT_COUNT() < 3) THEN
           WRITE(6,*) "There were fewer than 3 command arguments; using Rate_"
           RATEROOT = 'Rate_'
        ELSE
           CALL GET_COMMAND_ARGUMENT(3, VALUE = RATEROOT, LENGTH = ARGLEN, STATUS = STATUS)	
           IF (STATUS .ne. 0) THEN
               STOP 'RATEROOT was truncated'
           ENDIF
        ENDIF

        WRITE(filename1,'(I20)') counter
        filename=TRIM(ADJUSTL(SPECROOT))//TRIM(ADJUSTL(filename1))//'.dat'
        open(SPEC_UNIT, file=filename)
 
        filename=TRIM(ADJUSTL(RATEROOT))//TRIM(ADJUSTL(filename1))//'.dat'
        open(RATE_UNIT, file=filename)

        WRITE(SPEC_UNIT,'(100000(a25,"!"))') 'TIME','JDAY_GMT','LAT_degN','LON_degE','PRESS_hPa','TEMP_K','THETA',&
                'H2O', 'CFACTOR', 'M', 'N2', 'O2', 'JNO2FACT','JO1DFACT','RO2',&
                (SPC_NAMES(LOOKAT(i)), i=1,NLOOKAT)

        WRITE(RATE_UNIT,'(100000(a50,"!"))') 'TIME','JDAY_GMT','LAT_degN','LON_degE','PRESS_hPa','TEMP_K','THETA','H2O', &
                'CFACTOR', 'M', 'N2', 'O2', 'JNO2FACT','JO1DFACT','RO2',&
                                              EQN_NAMES

    END SUBROUTINE NewInitSaveData

! END of InitSaveData FUNCTION
! ****************************************************************
! ****************************************************************
!                            
! SaveData - Write LOOKAT species in the data file 
!   Parameters :                                                  
!
! ****************************************************************

    SUBROUTINE NewSaveData ()

        USE dsmacc_Global
        USE dsmacc_Monitor
        USE dsmacc_Integrator
        
        INTEGER i
        IF (CONSTRAIN_RUN) THEN
            WRITE(SPEC_UNIT,'(100000(E25.16E3,"!"))') JDAY_GMT,LAT,LON,PRESS,TEMP,THETA,H2O,&
            CFACTOR,M,N2,O2,JFACTNO2, JFACTO1D, RO2,C
            WRITE(RATE_UNIT,'(100000(E50.16E3,"!"))') JDAY_GMT,LAT,LON,PRESS,TEMP,THETA,H2O,&
            CFACTOR,M,N2,O2,JFACTNO2, JFACTO1D, RO2,RCONST
        ELSE
            WRITE(SPEC_UNIT,'(100000(E25.16E3,"!"))') TIME,JDAY_GMT,LAT,LON,PRESS,TEMP,THETA,H2O,&
            CFACTOR,M,N2,O2,JFACTNO2, JFACTO1D,RO2,(C(LOOKAT(i)), i=1,NLOOKAT)
            WRITE(RATE_UNIT,'(100000(E50.16E3,"!"))') TIME,JDAY_GMT,LAT,LON,PRESS,TEMP,THETA,H2O,&
            CFACTOR,M,N2,O2,JFACTNO2, JFACTO1D,RO2,RCONST
        ENDIF

    END SUBROUTINE NewSaveData

! END of SaveData FUNCTION
! ****************************************************************

! ****************************************************************
!                            
! CloseSaveData - Close the data file 
!   Parameters :                                                  
!
! ****************************************************************

    SUBROUTINE NewCloseData ()
        
        USE dsmacc_Parameters
        USE dsmacc_Global, ONLY : SPEC_UNIT, RATE_UNIT
        
        CLOSE(SPEC_UNIT)
        CLOSE(RATE_UNIT)
        
    END SUBROUTINE NewCloseData

! END of CloseSaveData FUNCTION
! ****************************************************************

         
#ENDINLINE
