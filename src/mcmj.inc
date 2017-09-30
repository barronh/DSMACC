// mcmj is designed to read a photolysis.txt file
// following the format 1 header line followed by N data lines
// A data line has the format JIDX<tab>L<tab>M<tab>N
// where
// - JIDX is the integer index for this j value, used for J function
// - L, M, and N are floating point parameters to calculate J using
//   the function J = 
#INLINE F90_GLOBAL
  INTEGER, parameter :: KJ=61
  REAL(dp) :: JSCALE = 1.d0
  ! variables for zenith routine which calculates zenith angle
  REAL(dp) cosx
#ENDINLINE 
#INCLUDE theta.inc
#INLINE F90_RATES
  RECURSIVE FUNCTION J(IJ) RESULT(TMPJ)
    INTEGER, intent(in) :: IJ
    REAL(dp), SAVE :: L, M, N
    REAL(dp), SAVE :: MCM_L(KJ)
    REAL(dp), SAVE :: MCM_M(KJ)
    REAL(dp), SAVE :: MCM_N(KJ)
    INTEGER, SAVE :: MCM_J(KJ)
    REAL(dp) COSX, TMPJ
    INTEGER, parameter :: PHOT_UNIT=100
    REAL(DP), parameter :: DEG2RAD = 3.141592653589793/180.d0
    INTEGER :: IOSTAT, tj, mcmj
    LOGICAL, SAVE :: FIRST = .true.
    IF (THETA .ge. 90) THEN
       TMPJ = 0.
       RETURN
    ENDIF
    COSX = COS(THETA*DEG2RAD)
    IF (FIRST) THEN
        ! open and read photolysis.txt
        open(UNIT=PHOT_UNIT,FILE='photolysis.txt', IOSTAT = IOSTAT, STATUS = 'OLD')
        IF (IOSTAT.ne.0) THEN
          WRITE(*, *) 'photolysis.txt file not found. Using default MCM v3.3 photolysis J,l,m,n'
          MCM_J(1:KJ) = &
            (/ 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
               21, 22, 23, 24, 0, 0, 0, 0, 0, 0, 31, 32, 33, 34, 35, 0, 0, 0, 0, 0,&
               41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 51, 52, 53, 54, 55, 56, 0, 0, 0, 0, 61 /)
          MCM_L(1:KJ) = &
            (/ 6.073D-05, 4.775D-04, 1.041D-05, 1.165D-02, 2.485D-02, 1.747D-01,&
               2.644D-03, 9.312D-07, 0.000D+00, 0.000D+00, 4.642D-05, 6.853D-05,&
               7.344D-06, 2.879D-05, 2.792D-05, 1.675D-05, 7.914D-05, 1.482D-06,&
               1.482D-06, 7.600D-04, 7.992D-07, 5.804D-06, 2.4246D-06, 2.424D-06,&
               0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00,&
               6.845D-05, 1.032D-05, 3.802D-05, 1.537D-04, 3.326D-04, 0.000D+00,&
               0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, 7.649D-06, 0.000D+00,&
               0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00,&
               0.000D+00, 0.000D+00, 1.588D-06, 1.907D-06, 2.485D-06, 4.095D-06,&
               1.135D-05, 4.365D-05, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00,&
               7.537D-04 /)
          MCM_M(1:KJ) = &
            (/ 1.743, 0.298, 0.723, 0.244, 0.168, 0.155, 0.261, 1.230, 0.000, 0.000,&
               0.762, 0.477, 1.202, 1.067, 0.805, 0.805, 0.764, 0.396, 0.396, 0.396,&
               1.578, 1.092, 0.395, 0.395, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,&
               0.130, 0.130, 0.644, 0.170, 0.148, 0.000, 0.000, 0.000, 0.000, 0.000,&
               0.682, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,&
               1.154, 1.244, 1.196, 1.111, 0.974, 1.089, 0.000, 0.000, 0.000, 0.000,&
               0.499 /)
          MCM_N(1:KJ) = &
            (/ 0.474, 0.080, 0.279, 0.267, 0.108, 0.125, 0.288, 0.307, 0.000, 0.000,&
               0.353, 0.323, 0.417, 0.358, 0.338, 0.338, 0.364, 0.298, 0.298, 0.298,&
               0.271, 0.377, 0.296, 0.296, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,&
               0.201, 0.201, 0.312, 0.208, 0.215, 0.000, 0.000, 0.000, 0.000, 0.000,&
               0.279, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,&
               0.318, 0.335, 0.328, 0.316, 0.309, 0.323, 0.000, 0.000, 0.000, 0.000,&
               0.266 /)
        ELSE
          READ(PHOT_UNIT, *) 
          !J	l	m	n
          DO WHILE (.NOT.IS_IOSTAT_END(IOSTAT))
            READ(PHOT_UNIT, *, IOSTAT = IOSTAT) mcmj, L, M, N
            IF (.NOT.IS_IOSTAT_END(IOSTAT)) THEN
            MCM_J(mcmj) = mcmj
            MCM_L(mcmj) = L
            MCM_M(mcmj) = M
            MCM_N(mcmj) = N
            END IF
          ENDDO
        END IF
        FIRST = .false.
        DO TJ=1,KJ
          MCMJ=MCM_J(TJ)
          IF (MCMJ.ne.0) THEN
          TMPJ = J(IJ)
          WRITE(*,'(I2,x,ES9.3,x,F7.3,x,F7.3x,F7.3)') MCMJ, MCM_L(MCMJ), MCM_M(MCMJ), MCM_N(MCMJ), TMPJ
          ENDIF
        ENDDO
    ENDIF
    TMPJ = MCM_L(IJ) * COSX**MCM_M(IJ) * EXP(-MCM_N(KJ) / COSX) *  JSCALE
    IF ((IJ .eq. 1) .and. (JFACTO1D.ne.1.)) THEN
       TMPJ = TMPJ * JFACTO1D
    ELSEIF ((IJ .eq. 4) .and. (JFACTO1D.ne.1.)) THEN
       TMPJ = TMPJ * JFACTNO2
    ENDIF

  END FUNCTION J
#ENDINLINE
