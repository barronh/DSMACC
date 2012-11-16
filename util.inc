#INLINE F90_UTIL

    SUBROUTINE NewInitVal (counter )
    
      USE dsmacc_Global
      use dsmacc_Monitor
    
      INTEGER :: i,j,found, counter
      REAL(kind=dp) :: x
      CHARACTER(LEN=10000) :: line
      CHARACTER(LEN=15) :: test
      REAL(kind=dp) :: concs(10000)
      REAL(kind=dp) :: oldvar(10000)
      LOGICAL :: SPECIAL
      INTEGER :: COUNT_NOX_CONSTRAINTS, IOSTAT
    
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
            open(UNIT=21,FILE='Init_cons.dat', IOSTAT = IOSTAT)
    
    ! set everything to zero for the first iteration
    ! after that IF doing a constrained run don't
            LINECOUNT=0
            DO WHILE (.NOT.IS_IOSTAT_END(IOSTAT))
            READ (21,*, IOSTAT = IOSTAT)
            LINECOUNT=LINECOUNT+1
            ENDDO
         CLOSE(21)
         WRITE (6,*) 'Input file has ',LINECOUNT,' lines'
      ENDIF
    
      IF (COUNTER .EQ. 1) THEN 
        open(UNIT=21,FILE='Init_cons.dat')
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
            WRITE (6,*) 'Negative Integration Time'
            WRITE (6,*) 'But not a special case'
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
                  IF (const_method(i) .EQ. 1) THEN 
                             CONSTRAIN(J)=CONCS(I)
                  IF (SPEC_NAME(I) .EQ. 'CH4') SPEC_CH4=.TRUE.
                  IF (SPEC_NAME(I) .EQ. 'H2') SPEC_H2=.TRUE.
                  ENDIF
               ENDIF
            ENDDO
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'H2O') THEN
               FOUND=1
               H2O=CONCS(I)
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'PRESS') THEN 
               FOUND=1
               PRESS=CONCS(I)
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'NOx') THEN
               FOUND=1
               CONSTRAIN_NOX=.TRUE.
               WRITE (6,*) 'Constraining total NOx concentation'
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'LAT') THEN
               FOUND=1
               LAT=CONCS(I)
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'LON') THEN 
               FOUND=1
               LON=CONCS(I)
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'TEMP') THEN
               FOUND=1
               TEMP=CONCS(I)
            ENDIF
            
            IF (TRIM(SPEC_NAME(I)) .EQ. 'JDAY') THEN 
               FOUND=1
               JDAY=CONCS(I)
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'TIME(GMTs)') THEN
               FOUND=1
               
            ENDIF
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'O3COL') THEN
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
    
            IF (TRIM(SPEC_NAME(I)) .EQ. 'JNO2') THEN
           FOUND=1
               IF (CONST_METHOD(I) .GE. 1) THEN 
                   JNO2=CONCS(I)
           ENDIF 
               JREPEAT=0
           IF (CONST_METHOD(I) .EQ. 2) THEN 
            JREPEAT=1
           ENDIF        
        ENDIF
            
        IF (TRIM(SPEC_NAME(I)) .EQ. 'JO1D') THEN
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
        INTEGER counter, length, ierr
        CHARACTER*20 filename1, filename
        INTEGER i
      
            
        CALL  ICNVRT(0, counter, filename1, length, ierr)

        filename='Spec_'//TRIM(filename1)//'.dat'
        open(10, file=filename)
 
        filename='Rate_'//TRIM(filename1)//'.dat'
        open(12, file=filename)

        WRITE(SPEC_UNIT,'(100000(a25,"!"))') 'TIME','LAT','LON','PRESS','TEMP','H2O', &
                'M', 'JNO2FACT','JO1DFACT','RO2', (SPC_NAMES(LOOKAT(i)), i=1,NLOOKAT)

        WRITE(RATE_UNIT,'(100000(a50,"!"))') 'TIME','LAT','LON','PRESS','TEMP','M', &
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
            WRITE(SPEC_UNIT,'(100000(E25.16,"!"))') JDAY,LAT,LON,PRESS,TEMP,H2O,CFACTOR,&
            JFACTNO2, JFACTO1D, RO2,C
            WRITE(RATE_UNIT,'(100000(E50.16,"!"))') JDAY,LAT,LON,PRESS,TEMP,CFACTOR, RCONST
        ELSE
            WRITE(SPEC_UNIT,'(100000(E25.16,"!"))') TIME,LAT,LON,PRESS,TEMP,H2O,CFACTOR,&
            JFACTNO2, JFACTO1D,RO2,(C(LOOKAT(i)), i=1,NLOOKAT)
            WRITE(RATE_UNIT,'(100000(E50.16,"!"))') TIME,LAT,LON,PRESS,TEMP,CFACTOR, RCONST
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
        
        CLOSE(10)
        CLOSE(12)
        
    END SUBROUTINE NewCloseData

! END of CloseSaveData FUNCTION
! ****************************************************************

    SUBROUTINE ICNVRT(WAY,NUM,STRING,LENGTH,IERR)
!
!       FUNCTION:
!F
!F        This subroutine does an INTEGER-to-CHARACTER conversion
!F        or a characater-to-INTEGER conversion depending on the
!F        INTEGER WAY:
!F                IF WAY = 0 then an INTEGER-to-CHARACTER conversion
!F                is done. IF WAY .NE. 0 then a CHARACTER-to-INTEGER
!F                conversion is done.
!F
!       USAGE:
!U
!U        !ALL I!NVRT(WAY,NUM,STRING)
!U             where WAY, NUM, STRING, and LENGTH are defined below.
!U
!U        Example: !ALL I!NVRT(0,1000,STRING,LENGTH)
!U                 on return STRING = '1000' and
!U                 LENGTH = 4.
!U         
!       INPUTS:
!I
!I        WAY - INTEGER; Determines which way the conversion goes:
!I              IF WAY = 0 then an INTEGER-to-CHARACTER conversion
!I                         is performed;
!I              IF WAY.NE.0 then a CHARACTER-to-INTEGER conversion
!I                         is performed.
!I
!I         NUM - INTEGER; an input only IF WAY = 0. NUM is the INTEGER
!I               number to be converted to a CHARACTER expression.
!I
!I         STRING - !HARACTER; an input only IF WAY .NE. 0. STRING
!I                is the CHARACTER expression to be converted to an
!CI                INTEGER value. It contain no decimal points or 
!I                non-numeric characters other than possibly a
!I                sign. IF STRING contains  a '+' sign, it will be
!I                stripped of it on return.
!I
!       OUTPUTS:
!O
!O         NUM - INTEGER; contains the INTEGER representation of 
!O                STRING.
!O
!O         STRING - !HARA!TER; contains the CHARACTER representation of
!O                  NUM.
!CO
!CO         LENGTH - INTEGER; The length of STRING to the first blank.
!CO                  The significant part of STRING can be accessed with
!CO                  the declaration STRING(1:LENGTH).
!CO
!CO         IERR - INTEGER variable giving return condition:
!CO                IERR = 0 for normal return;
!CO                IERR = 1 IF NUM cannot be converted to STRING because
!CO                       STRING is too short or STRING cannot be
!CO                       converted to NUM because STRING is too long.
!CO                IERR = 2 IF STRING contained a non-numeric CHARACTER
!CO                       other than a leading sign or something went
!CO                       wrong with an INTEGER-to-CHARACTER conversion.
!CO
!C       ALGORITHM:
!CA
!CA         Nothing noteworthy, except that this subroutine will work
!CA          for strange CHARACTER sets where the CHARACTER '1' doesn't
!CA          follow '0', etc.
!CA
!C       MACHINE DEPENDENCIES:
!CM
!CM          The parameter MAXINT (below) should be set to the
!CM          number of digits that an INTEGER data type can have
!CM          not including leading signs. For VAX FORTRAN V4.4-177
!CM          MAXINT = 10.
!CM
!CM          NOTE: Under VAX FORTRAN V4.4-177, the
!CM          error condition IERR = 1 will never occur for an
!CM          INTEGER-to-CHARACTER conversion IF STRING
!CM          is allocated at least 11 bytes (CHARACTER*11).
!CM
!C       HISTORY:
!CH
!CH      written by:             bobby bodenheimer
!CH      date:                   september 1986
!CH      current version:        1.0
!CH      modifications:          none
!CH
!C       ROUTINES CALLED:
!CC
!CC          NONE.
!CC
!C----------------------------------------------------------------------
!C       written for:    The CASCADE Project
!C                       Oak Ridge National Laboratory
!C                       U.S. Department of Energy
!C                       contract number DE-AC05-840R21400
!C                       subcontract number 37B-7685 S13
!C                       organization:  The University of Tennessee
!C----------------------------------------------------------------------
!C       THIS SOFTWARE IS IN THE PUBLIC DOMAIN
!C       NO RESTRICTIONS ON ITS USE ARE IMPLIED
!C----------------------------------------------------------------------
!C
!C Global Variables.
!C
        INTEGER       WAY
        INTEGER       LENGTH
        INTEGER       NUM
        INTEGER       IERR
        CHARACTER*(*) STRING
!C
!C Local Variables
!C
        INTEGER       I
        INTEGER       MAXINT
        INTEGER       MNUM
        INTEGER       M
        LOGICAL       NEG
!C
        PARAMETER(MAXINT=10)
!C
        NEG = .FALSE.
        IERR = 0
!C
!C  INTEGER-to-CHARACTER conversion.
!C
        IF (WAY.EQ.0) THEN
            STRING = ' '
            IF (NUM.LT.0) THEN
                NEG = .TRUE.
                MNUM = -NUM
                LENGTH = INT(LOG10(REAL(MNUM))) + 1
            ELSE IF (NUM.EQ.0) THEN
                MNUM = NUM
                LENGTH = 1
            ELSE
                MNUM = NUM
                LENGTH = INT(LOG10(REAL(MNUM))) + 1
            END IF
            IF (LENGTH.GT.LEN(STRING)) THEN
                IERR = 1
                RETURN
            END IF
            DO 10, I=LENGTH,1,-1
                M=INT(REAL(MNUM)/10**(I-1))
                IF (M.EQ.0) THEN
                   STRING(LENGTH-I+1:LENGTH-I+1) = '0'
                ELSE IF (M.EQ.1) THEN
                   STRING(LENGTH-I+1:LENGTH-I+1) = '1'
                ELSE IF (M.EQ.2) THEN
                   STRING(LENGTH-I+1:LENGTH-I+1) = '2'
                ELSE IF (M.EQ.3) THEN
                   STRING(LENGTH-I+1:LENGTH-I+1) = '3'
                ELSE IF (M.EQ.4) THEN
                   STRING(LENGTH-I+1:LENGTH-I+1) = '4'
                ELSE IF (M.EQ.5) THEN
                   STRING(LENGTH-I+1:LENGTH-I+1) = '5'
                ELSE IF (M.EQ.6) THEN
                   STRING(LENGTH-I+1:LENGTH-I+1) = '6'
                ELSE IF (M.EQ.7) THEN
                   STRING(LENGTH-I+1:LENGTH-I+1) = '7'
                ELSE IF (M.EQ.8) THEN
                   STRING(LENGTH-I+1:LENGTH-I+1) = '8'
                ELSE IF (M.EQ.9) THEN
                   STRING(LENGTH-I+1:LENGTH-I+1) = '9'
                ELSE
                   IERR = 2
                   RETURN
                END IF
                MNUM = MNUM - M*10**(I-1)
            10       CONTINUE
                IF (NEG) THEN
                    STRING = '-'//STRING
                    LENGTH = LENGTH + 1
                END IF
!C
!C  CHARACTER-to-INTEGER conversion.
!!C
        ELSE
            IF (STRING(1:1).EQ.'-') THEN
                NEG = .TRUE.
                STRING = STRING(2:LEN(STRING))
            END IF

            IF (STRING(1:1).EQ.'+') STRING = STRING(2:LEN(STRING))
                NUM = 0
                LENGTH = INDEX(STRING,' ') - 1
                IF (LENGTH.GT.MAXINT) THEN
                    IERR = 1
                    RETURN
                END IF
                DO 20, I=LENGTH,1,-1
                    IF (STRING(LENGTH-I+1:LENGTH-I+1).EQ.'0') THEN
                        M = 0
                    ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1).EQ.'1') THEN
                        M = 1
                    ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1).EQ.'2') THEN
                        M = 2
                    ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1).EQ.'3') THEN
                        M = 3
                    ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1).EQ.'4') THEN
                        M = 4
                    ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1).EQ.'5') THEN
                        M = 5
                    ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1).EQ.'6') THEN
                        M = 6
                    ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1).EQ.'7') THEN
                        M = 7
                    ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1).EQ.'8') THEN
                        M = 8
                    ELSE IF (STRING(LENGTH-I+1:LENGTH-I+1).EQ.'9') THEN
                        M = 9
                    ELSE
                        IERR = 2
                        RETURN
                    END IF
                    NUM = NUM + INT(10**(I-1))*M
     20         CONTINUE
                IF (NEG) THEN
                    NUM = -NUM
                    STRING = '-'//STRING
                    LENGTH = LENGTH + 1
                END IF
            END IF
!C
!C  Last lines of ICNVRT
!C
        RETURN
    END SUBROUTINE
         
#ENDINLINE
