#INLINE F90_GLOBAL
    REAL(dp) :: PRESS, LAT, LON, O3COL, JO1D, JNO2
    REAL(dp) :: JDAY_GMT, JDAY_LOCAL, JREPEAT, ALBEDO, SAREA, RP1
    REAL(dp) :: TSTART_GMT, TSTART_LOCAL
    INTEGER :: INIT_TIME, NOX(NSPEC)
    REAL(dp):: CONSTRAIN(NSPEC)
    LOGICAL :: SPEC_CH4, SPEC_H2
    INTEGER :: IntTime
    LOGICAL :: JDAYISGMT = .true.
    LOGICAL :: DEBUG = .true.
    LOGICAL :: CONSTRAIN_NOX
    LOGICAL :: CONSTRAIN_RUN = .false., LAST_POINT = .false., OUTPUT_LAST = .false.
    CHARACTER(LEN=15) :: spec_name(10000)
    INTEGER :: const_method(10000), LINECOUNT
    ! C - Concentration of all species
    REAL(kind=dp) :: COLD(NSPEC)
    INTEGER, PARAMETER :: OUTPUT_UNIT = 6, ERROR_UNIT = 0,&
    SPEC_UNIT = 10, RATE_UNIT = 12
    REAL:: jfactno2 = 1., jfacto1d = 1.
    REAL(dp) :: THETA
    !END photolysis variables
#ENDINLINE {above lines go into MODULE KPP_ROOT_Global}
