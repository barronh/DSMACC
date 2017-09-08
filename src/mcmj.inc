// mcmj is designed to read a photolysis.txt file
// following the format 1 header line followed by N data lines
// A data line has the format JIDX<tab>L<tab>M<tab>N
// where
// - JIDX is the integer index for this j value, used for J function
// - L, M, and N are floating point parameters to calculate J using
//   the function J = 
#INLINE F90_GLOBAL
  INTEGER, parameter :: KJ=61, NTHETA=19
  ! variables for zenith routine which calculates zenith angle
  REAL(dp) secx, cosx
  ! generic reaction rate variables
  ! variables for calculation of photolysis reaction rates
  real*8 bs(NTHETA,kj), cs(NTHETA,kj), ds(NTHETA,kj)
#ENDINLINE 
#INCLUDE theta.inc
#INLINE F90_RATES
  REAL(kind=dp) FUNCTION J(IJ)
    USE dsmacc_Global, ONLY : THETA, KJ
    INTEGER, intent(in) :: IJ
    REAL(dp), SAVE :: L, M, N
    REAL(dp), SAVE :: MCM_L(KJ), MCM_M(KJ), MCM_N(KJ)
    INTEGER, SAVE :: MCM_J(KJ)
    REAL(dp) COSX
    INTEGER, parameter :: PHOT_UNIT=100
    REAL(DP), parameter :: DEG2RAD = 3.141592653589793/180.d0
    INTEGER :: IOSTAT, tj, mcmj
    LOGICAL, SAVE :: FIRST = .true.
    IF (THETA .ge. 90) THEN
       J = 0.
       RETURN
    ENDIF
    COSX = COS(THETA*DEG2RAD)
    IF (FIRST) THEN
        ! open and read photolysis.txt
        open(UNIT=PHOT_UNIT,FILE='photolysis.txt', IOSTAT = IOSTAT)
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
        FIRST = .false.
        DO TJ=1,KJ
          MCMJ=MCM_J(TJ)
          IF (MCMJ.ne.0) THEN
          J = MCM_L(IJ) * COSX**MCM_M(IJ) * EXP(-MCM_N(KJ) / COSX)
          WRITE(*,*) MCMJ, MCM_L(MCMJ), MCM_M(MCMJ), MCM_N(MCMJ), J
          ENDIF
        ENDDO
    ENDIF
    J = MCM_L(IJ) * COSX**MCM_M(IJ) * EXP(-MCM_N(KJ) / COSX)
  END FUNCTION J
#ENDINLINE
