#INLINE F90_UTIL
    SUBROUTINE ShowHelp()
        WRITE(6,*) "There are 3 optional arguments:"
        WRITE(6,*) " - initpath: path to initialization file with 3 header lines"
        WRITE(6,*) "   and N data lines."
        WRITE(6,*) "      - First line (i10) indicates "
        WRITE(6,*) "        integration time (>0), constrained (-1) or diurnal"
        WRITE(6,*) "        steady-state constrained (-2) run. "
        WRITE(6,*) "      - The second line (10000(a15,x)) has the names of "
        WRITE(6,*) "        species or global properties (e.g., JDAY_GMT, TEMP_K,"
        WRITE(6,*) "        PRESS_Pa, LON_degE, LAT_degN). "
        WRITE(6,*) "      - The third line (10000(i15,x)) has 1 (yes) or 0 (no)"
        WRITE(6,*) "        indicating if the species will be constrained"
        WRITE(6,*) "      - The data lines (10000(e15.4,x)) specify data values"
        WRITE(6,*) " - concpath: prefix path to write out concentrations in volume "
        WRITE(6,*) "   mixing ratios. Run number and .dat will be appended)"
        WRITE(6,*) " - ratepath: prefix path to write out reaction rates (1/s, "
        WRITE(6,*) "   cm3/molecule/s). Run number and .dat will be appended)"
        WRITE(6,*) ""
        WRITE(6,*) "This warning is shown if any are not provided, even though"
        WRITE(6,*) "default values will be tested."
    END SUBROUTINE ShowHelp
    SUBROUTINE NewInitVal (counter )
    
      INTEGER :: i,JJI,found, counter, STATUS, ARGLEN
      REAL(kind=dp) :: x
      CHARACTER(LEN=10000) :: line, INPATH
      CHARACTER(LEN=15) :: test
      CHARACTER(LEN=24) :: THIS_SPEC_NAME
      REAL(kind=dp) :: concs(10000)
      REAL(kind=dp) :: oldvar(10000)
      LOGICAL :: SPECIAL
      INTEGER :: COUNT_NOX_CONSTRAINTS, IOSTAT
      SAVE SPECIAL, INPATH
    
    ! Open the file with the info
    
      write (6,*) 'Initializing model point', counter
      call flush(6)
    
    
      DO I=1,NSPEC
         OLDVAR(I)=C(I)
         C(I)=0.
      ENDDO
      SAREA=0.
      ALBEDO=0.
      RP1=0.
    
      IF (counter .eq. 0) then 
        IF (COMMAND_ARGUMENT_COUNT() == 0) THEN
           WRITE(6,*) "There were no command arguments; using Init_cons.dat."
           WRITE(6,*) "Use -h for more details"
           INPATH = 'Init_cons.dat'
        ELSE
           CALL GET_COMMAND_ARGUMENT(1, VALUE = INPATH, LENGTH = ARGLEN, STATUS = STATUS)
           SELECT CASE (TRIM(ADJUSTL(INPATH)))
                  CASE ('-h', '--help')
                    CALL ShowHelp
                    STOP 0
           END SELECT
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
            STOP 1
         ENDIF
    
         
         READ(21,'(10000(a15,x))') spec_name
         READ(21,'(10000(i15,x))') const_method     
      ENDIF
      DO JJI=1,NSPEC 
         SPEC_NAME(I)=ADJUSTL(SPEC_NAME(I))
      ENDDO
    
      IF (COUNTER .NE. 0) THEN 
    
      READ (21,'(10000(e15.4,x))', IOSTAT = IOSTAT) concs
      IF (IS_IOSTAT_END(IOSTAT)) THEN
         LAST_POINT=.TRUE.
      ENDIF
      DO I=1,10000
         IF (TRIM(SPEC_NAME(I)).eq.'') CYCLE
         FOUND=0
         DO JJI=1,NSPEC 
            IF (TRIM(SPEC_NAME(I)) .EQ. TRIM(SPC_NAMES(JJI))) THEN 
               FOUND=1
               C(JJI)=CONCS(I)
               IF (const_method(i) .EQ. 1) CONSTRAIN(JJI)=CONCS(I)
               IF (SPEC_NAME(I) .EQ. 'CH4') SPEC_CH4=.TRUE.
               IF (SPEC_NAME(I) .EQ. 'H2') SPEC_H2=.TRUE.
            ENDIF
         ENDDO
         THIS_SPEC_NAME = TRIM(SPEC_NAME(I))
         ! Check for Global Variables
         ! by name
         SELECT CASE (THIS_SPEC_NAME)
            CASE('H2O')
               H2O=CONCS(I)
            CASE('PRESS_hPa')
               PRESS=CONCS(I)*100.
            CASE('PRESS_Pa')
               PRESS=CONCS(I)
            CASE('PRESS')
               PRESS=CONCS(I)
            CASE('NOx')
               CONSTRAIN_NOX=.TRUE.
               WRITE (6,*) 'Constraining total NOx concentration'
            CASE('LAT_degN')
               LAT=CONCS(I)
            CASE('LON_degE') 
               LON=CONCS(I)
            CASE('TEMP_C')
               TEMP=CONCS(I) + 273.15
            CASE('TEMP_K')
               TEMP=CONCS(I)
            CASE('TEMP')
               TEMP=CONCS(I)
            CASE('JDAY_GMT') 
               JDAY_GMT=CONCS(I)
            CASE('TIME(GMTs)')
            CASE('O3COL_DU')
               O3COL=CONCS(I)
            CASE('ALBEDO') 
               ALBEDO=CONCS(I)
            CASE('SAREA') 
               SAREA=CONCS(I)
            CASE('RP1') 
               RP1=CONCS(I)
            CASE('JNO2_per_s')
               IF (CONST_METHOD(I) .GE. 1) THEN
                   JNO2=CONCS(I)
               ENDIF 
               JREPEAT=0
               IF (CONST_METHOD(I) .EQ. 2) THEN
                   JREPEAT=1
               ENDIF
            CASE('JO1D_per_s')
               IF (CONST_METHOD(I) .GE. 1) THEN
                   JO1D=CONCS(I)
               ENDIF
               JREPEAT=0
               IF (CONST_METHOD(I) .EQ. 2) THEN
                   JREPEAT=1
               ENDIF    
            CASE DEFAULT
                IF (FOUND .eq. 0) THEN
                    WRITE (6,*) SPEC_NAME(I),' NOT FOUND'
                    IF (SPEC_NAME(I)(1:1) .NE. 'X') STOP 1
                    WRITE (6,*) 'Starts with an X so ignored and continued'   
                ENDIF
        END SELECT
      ENDDO
         
      CFACTOR=PRESS*1e-6/(8.314*TEMP)*6.022E23
      
      H2O=H2O*CFACTOR 
      DO I=1,NSPEC
         C(I)=C(I)*CFACTOR
    !     IF (OLDVAR(I) .NE. CONSTRAIN(I)) THEN 
          CONSTRAIN(I)=CONSTRAIN(I)*CFACTOR
    !  ENDIF
      ENDDO
    
     DO I=1,NSPEC
    !    C(I)=0.
    !     IF (C(I) .EQ. 0) C(I)=OLDVAR(I)  
     ENDDO
    
    ! FIND NOX species
     IF (CONSTRAIN_NOX) THEN 
       COUNT_NOX_CONSTRAINTS=0
       DO I=1,NSPEC 
          SELECT CASE (TRIM(SPC_NAMES(I)))
              CASE ('NO2', 'NO', 'NO3', 'N2O5', 'HONO', 'HNO2', 'HO2NO2', 'PNA', 'HNO4')
                IF (TRIM(SPC_NAMES(I)).EQ. 'N2O5') THEN
                    NOX(I)=2
                ELSE
                    NOX(I)=1
                ENDIF
                CONSTRAIN(I)=-1.*CONSTRAIN(I)
                IF (CONSTRAIN(I) .NE. 0) THEN 
                   COUNT_NOX_CONSTRAINTS=COUNT_NOX_CONSTRAINTS+1
                ENDIF
                WRITE (6,*) SPC_NAMES(I),' IN NOX FAMILY'
              CASE DEFAULT
          END SELECT
       ENDDO
      
       IF (COUNT_NOX_CONSTRAINTS .GT. 1) THEN 
         WRITE (6,*) 'You can not contrains multiple NOX species'
         STOP 1
       ENDIF
     ENDIF
    
    ! FIND CH4 and H2 species
      DO I=1,NSPEC 
         IF (TRIM(SPC_NAMES(I)) .EQ. 'CH4' .AND. (SPEC_CH4 .EQV. .FALSE.)) THEN
            WRITE (6,*) 'No CH4 specified assuming 1770 ppbv'
            C(I)=1770e-9*CFACTOR
            CONSTRAIN(I)=C(I)
         ENDIF
    
         IF (TRIM(SPC_NAMES(I)) .EQ. 'H2' .AND. (SPEC_H2 .EQV. .FALSE.)) THEN 
            WRITE (6,*) 'No H2 specified assuming 550 ppbv'
            C(I)= 550e-9*CFACTOR
            CONSTRAIN(I)=C(I)
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
        INTEGER counter, length, ierr, STATUS, ARGLEN
        CHARACTER*20 filename1
        CHARACTER(LEN=10000) :: SPECROOT, RATEROOT, filename
        INTEGER i
        IF (COMMAND_ARGUMENT_COUNT() .lt. 2) THEN
           WRITE(6,*) "There were fewer than 2 command arguments; using Spec_"
           WRITE(6,*) "Use -h for more details"
           SPECROOT = 'Spec_'
        ELSE
           CALL GET_COMMAND_ARGUMENT(2, VALUE = SPECROOT, LENGTH = ARGLEN, STATUS = STATUS)
           IF (STATUS .ne. 0) THEN
               WRITE(*,*)'SPECROOT was truncated'
               STOP 1
           ENDIF
        ENDIF
        IF (COMMAND_ARGUMENT_COUNT() .lt. 3) THEN
           WRITE(6,*) "There were fewer than 3 command arguments; using Rate_"
           WRITE(6,*) "Use -h for more details"
           RATEROOT = 'Rate_'
        ELSE
           CALL GET_COMMAND_ARGUMENT(3, VALUE = RATEROOT, LENGTH = ARGLEN, STATUS = STATUS)
           IF (STATUS .ne. 0) THEN
               WRITE(*,*)'RATEROOT was truncated'
               STOP 1
           ENDIF
        ENDIF

        WRITE(filename1,'(I20)') counter
        filename=TRIM(ADJUSTL(SPECROOT))//TRIM(ADJUSTL(filename1))//'.dat'
        open(SPEC_UNIT, file=filename)
 
        filename=TRIM(ADJUSTL(RATEROOT))//TRIM(ADJUSTL(filename1))//'.dat'
        open(RATE_UNIT, file=filename)

        WRITE(SPEC_UNIT,'(100000(a25,"!"))') 'TIME','JDAY_GMT','LAT_degN','LON_degE','PRESS_Pa','TEMP_K','THETA',&
                'H2O', 'CFACTOR', 'M', 'N2', 'O2', 'JNO2FACT','JO1DFACT','RO2',&
                (SPC_NAMES(LOOKAT(i)), i=1,NLOOKAT)

        WRITE(RATE_UNIT,'(100000(a50,"!"))') 'TIME','JDAY_GMT','LAT_degN','LON_degE','PRESS_Pa','TEMP_K','THETA','H2O', &
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
        
        CLOSE(SPEC_UNIT)
        CLOSE(RATE_UNIT)
        
    END SUBROUTINE NewCloseData

! END of CloseSaveData FUNCTION
! ****************************************************************

         
#ENDINLINE
