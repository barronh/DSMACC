#INLINE F90_GLOBAL
    REAL(dp) :: PRESS, LAT, LON, O3COL, JO1D, JNO2
    REAL(dp) :: JDAY, JREPEAT, ALBEDO, SAREA, RP1
    INTEGER :: INIT_TIME, NOX(NVAR)
    REAL(dp):: CONSTRAIN(NVAR)
    LOGICAL :: SPEC_CH4, SPEC_H2
    INTEGER :: IntTime
    LOGICAL :: CONSTRAIN_NOX
    LOGICAL :: CONSTRAIN_RUN = .false., LAST_POINT = .false., OUTPUT_LAST = .false.
    CHARACTER(LEN=15) :: spec_name(10000)
    INTEGER :: const_method(10000), LINECOUNT
    ! C - Concentration of all species
    REAL(kind=dp) :: COLD(NSPEC)
    INTEGER, PARAMETER :: OUTPUT_UNIT = 6, ERROR_UNIT = 0,&
     SPEC_UNIT = 10, RATE_UNIT = 12
    !Photolysis variables
    include '../tuv_old/params'
    REAL*8::bs(19,kj), cs(19,kj),ds(19,kj)
    REAL(dp) :: THETA
    REAL::svj_tj(kt,kj), szas(kt), jfactno2, jfacto1d
    !END photolysis variables
#ENDINLINE {above lines go into MODULE dsmacc_Global}
