#INLINE F90_RATES
  REAL(kind=dp) FUNCTION GEOS_STD(A0, B0, C0)
    REAL A0, B0, C0

!   GEOS Chem reaction rates have the form K = A * (300 / T)**B * EXP(C / T)
!   KPP ARR reaction rates have the form   K = A * (T/300.0)**C * EXP(-B/T) 
!   
!   Translation reorders B and C and changes both their signs
    GEOS_STD = A0 * (300. / TEMP)**B0 * EXP(C0 / TEMP)
  END FUNCTION
  
  !---------------------------------------------------------------------------
  
  REAL(kind=dp) FUNCTION GEOS_P(A0, B0, C0, A1, B1, C1, &
                                FCV, FCT1, FCT2)
    REAL A0, B0, C0, A1, B1, C1 ,CF
    REAL FCV, FCT1, FCT2
    REAL(kind=dp) K0M, K1

    IF (FCT2 .NE. 0.000000e+00) THEN 
      CF = EXP(-TEMP / FCT1) + EXP(-FCT2 / TEMP) 
    ELSEIF (FCT1 .NE. 0.000000e+00) THEN 
      CF = EXP(-TEMP / FCT1)
    ELSE
      CF = FCV
    ENDIF
    
!   GEOS Chem reaction rates have the form K = A * (300 / T)**B * EXP(C / T)
!   KPP ARR reaction rates have the form   K = A * (T/300.0)**C * EXP(-B/T) 
!   
!   Translation reorders B and C and changes both their signs

    K0M = GEOS_STD(A0, B0, C0) * M

    K1 = GEOS_STD(A1, B1, C1)
    K1 = K0M / K1

    GEOS_P = (K0M / (1.0_dp + K1))*   &
           DBLE(CF)**(1.0_dp / (1.0_dp + (LOG10(K1))**2))
  END FUNCTION

  !---------------------------------------------------------------------------
  
  REAL(kind=dp) FUNCTION GEOS_Z(A0, B0, C0, A1, B1, C1, A2, B2, C2)
    REAL A0, B0, C0, A1, B1, C1, A2, B2, C2
    REAL(kind=dp) K0, K1, K2

    K0 = GEOS_STD(A0, B0, C0)
    K1 = GEOS_STD(A1, B1, C1)*M
    K2 = GEOS_STD(A2, B2, C2)

    GEOS_Z = (K0 + K1) * (1 + H2O * K2)
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_Y(A0, B0, C0)
    REAL A0, B0, C0
    REAL(kind=dp) K0
    REAL(kind=dp) KHI1,KLO1,XYRAT1,BLOG1,FEXP1,KHI2,KLO2,XYRAT2,BLOG2
    REAL(kind=dp) FEXP2,KCO1,KCO2,KCO

    !IGNORES INPUTS per v08-02-04 update
    !K0 = GEOS_STD(A0, B0, C0)
    !GEOS_Y = K0 * (1 + .6 * (PRESS) / 101325.)
    KLO1 = 5.9D-33 * (300 / TEMP)**(1.4D0) 
    KHI1 = 1.1D-12 * (300 / TEMP)**(-1.3D0)
    XYRAT1 = KLO1 * M /KHI1
    BLOG1 = LOG10(XYRAT1)
    FEXP1 = 1.D0 / (1.D0 + BLOG1 * BLOG1)
    KCO1 = KLO1 * M * 0.6**FEXP1 /(1.d0 + XYRAT1)
    KLO2 = 1.5D-13 * (300 / TEMP)**(-0.6D0)
    KHI2 = 2.1D09 * (300 / TEMP)**(-6.1D0)
    XYRAT2 = KLO2 * M /KHI2
    BLOG2 = LOG10(XYRAT2)
    FEXP2 = 1.D0 / (1.D0 + BLOG2 * BLOG2)
    KCO2 = KLO2 * 0.6**FEXP2 / (1.d0 + XYRAT2)
    GEOS_Y = KCO1 + KCO2
  END FUNCTION
  
  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_X(A0, B0, C0, A1, B1, C1, A2, B2, C2)
    REAL A0, B0, C0, A1, B1, C1, A2, B2, C2
    REAL(kind=dp) K0, K2, K3            
    K0 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    K3 = GEOS_STD(A2, B2, C2)
    K3 = K3 * M
    GEOS_X = K0 + K3 / (1.0_dp + K3 / K2 )
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_C(A0, B0, C0)
    REAL A0, B0, C0, A1, B1, C1, A2, B2, C2
    REAL(kind=dp) K1
    K1 = GEOS_STD(A0, B0, C0)
    GEOS_C = K1 * (O2 + 3.5D18) / (2.0_dp * O2 + 3.5D18)
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_K(A0, B0, C0)
    REAL A0, B0, C0
    GEOS_K = 0
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_V(A0, B0, C0, A1, B1, C1)
    REAL A0, B0, C0, A1, B1, C1
    REAL(kind=dp) K1, K2
    K1 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    GEOS_V = K1 / (1 + K2)
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_E(A0, B0, C0, Kf)
    REAL A0, B0, C0
    REAL(kind=dp) K1, Kf
    K1 = GEOS_STD(A0, B0, C0)
    GEOS_E = Kf / K1
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL*4 FUNCTION FYRNO3(CN)
    REAL*4, PARAMETER :: Y300 = .826, ALPHA = 1.94E-22
    REAL*4, PARAMETER :: BETA = .97, XM0 = 0., XMINF = 8.1
    REAL*4, PARAMETER :: XF = .411
    
    REAL*4 CN
    REAL*4 XCARBN, ZDNUM, TT, XXYN, YYYN, AAA, ZZYN, RARB
    XCARBN = CN
    ZDNUM = M
    TT = TEMP
    
    XXYN = ALPHA * exp(BETA * XCARBN) * ZDNUM * ((300. / TT)**XM0)
    YYYN = Y300 * ((300. / TT)**XMINF)
    AAA = log10(XXYN / YYYN)
    ZZYN = 1. / (1. + AAA / AAA)
    RARB = (XXYN / (1. + (XXYN / YYYN))) * (XF**ZZYN)
    FYRNO3 = RARB / (1. + RARB)
  END FUNCTION
  
  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_A(A0, B0, C0, A1, B1, C1 )
    REAL A0, B0, C0, A1, B1, C1
    REAL TMP_A0
    TMP_A0 = A0 * FYRNO3(A1)
    GEOS_A = GEOS_STD(TMP_A0, B0, C0)
  END FUNCTION
  
  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_B(A0, B0, C0, A1, B1, C1 )
    REAL A0, B0, C0, A1, B1, C1
    REAL TMP_A0
    TMP_A0 = A0 * ( 1. - FYRNO3(A1) )
    GEOS_B = GEOS_STD(TMP_A0, B0, C0)
  END FUNCTION
  
  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_JO3(O3J)
    REAL(kind=dp) O3J, T3I, RO1DplH2O, RO1DplH2, RO1DplN2, RO1DplO2, RO1D
    T3I = 1.0_dp/TEMP
    RO1DplH2O = 1.63d-10 * EXP( 60.d0 * T3I) * H2O
    RO1DplH2 = 1.2e-10 * H2
    RO1DplN2 = 2.15d-11 * EXP(110.d0 * T3I) * N2
    RO1DplO2 = 3.30d-11 * EXP( 55.d0 * T3I) * O2
    RO1D = RO1DplH2O + RO1DplH2 + RO1DplN2 + RO1DplO2
    GEOS_JO3 = O3J * RO1DplH2O / RO1D
  END FUNCTION
  
  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_JO3_2(O3J)
    REAL(kind=dp) O3J, T3I, RO1DplH2O, RO1DplH2, RO1DplN2, RO1DplO2, RO1D
    T3I = 1.0_dp/TEMP
    RO1DplH2O = 1.63d-10 * EXP( 60.d0 * T3I) * H2O
    RO1DplH2 = 1.2e-10 * H2
    RO1DplN2 = 2.15d-11 * EXP(110.d0 * T3I) * N2
    RO1DplO2 = 3.30d-11 * EXP( 55.d0 * T3I) * O2
    RO1D = RO1DplH2O + RO1DplH2 + RO1DplN2 + RO1DplO2
    GEOS_JO3_2 = O3J * RO1DplH2 / RO1D
  END FUNCTION
  
  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_G(A0, B0, C0, A1, B1, C1)
    REAL A0, B0, C0, A1, B1, C1
    REAL(kind=dp) K1, K2
    K1 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    GEOS_G = K1 / ( 1.0_dp + K1 * O2 )
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION JHNO4_NEAR_IR(HNO4J)
    REAL(kind=dp) HNO4J
    IF (HNO4J .GT. 0.d0) THEN
      JHNO4_NEAR_IR = HNO4J + 1D-5
    ELSE
      JHNO4_NEAR_IR = HNO4J
    ENDIF
    
  END FUNCTION

! I have to add a bunch of new reactions here(jmao, 07/21/2014):

  !---------------------------------------------------------------------------

   REAL(kind=dp) FUNCTION GEOS_F(A0, B0, C0)
     REAL A0, B0, C0
     REAL(kind=dp) FRAC
!     HAC + OH
      FRAC=1-23.7*exp(-1./(60*TEMP))
      IF (FRAC<0d0) FRAC=0.01d0
      GEOS_F=GEOS_STD(A0, B0, C0)*FRAC

   END FUNCTION


  !---------------------------------------------------------------------------

   REAL(kind=dp) FUNCTION GEOS_L(A0, B0, C0)
     REAL A0, B0, C0
     REAL(kind=dp) FRAC
!     HAC + OH
      FRAC=1-23.7*exp(-1./(60*TEMP))
      IF (FRAC<0d0) FRAC=0.01d0
      GEOS_L=GEOS_STD(A0, B0, C0)*(1.0_dp-FRAC)      
   END FUNCTION
   
   REAL(kind=dp) FUNCTION GEOS_N(A0, B0, C0)
     REAL A0, B0, C0
     REAL(kind=dp) FRAC
!        GLYC+OH
!FRAC=1!>11.0729*exp(!>1/73T),Butkovskaya 2006
      FRAC=1-11.0729*exp(-1./(73.*TEMP))
      IF (FRAC<0d0) FRAC=0.01d0
      GEOS_N=GEOS_STD(A0, B0, C0)*FRAC

   END FUNCTION

  !---------------------------------------------------------------------------


   REAL(kind=dp) FUNCTION GEOS_O(A0, B0, C0)
     REAL A0, B0, C0
     REAL(kind=dp) FRAC
!  GLYC+OH
!  FRAC=1!>11.0729*exp(!>1/73T),Butkovskaya 2006
      FRAC=1-11.0729*exp(-1./(73.*TEMP))
      IF (FRAC<0d0) FRAC=0.01d0
      GEOS_O=GEOS_STD(A0, B0, C0)*(1.0_dp-FRAC)      

   END FUNCTION

  !---------------------------------------------------------------------------

   REAL(kind=dp) FUNCTION GEOS_HR(A0, B0, C0, A1, B1, C1 )
      REAL A0, B0, C0, A1, B1, C1
      REAL(kind=dp) FRAC
      FRAC = 1.-EXP(-0.245*A1)
      GEOS_HR = GEOS_STD(A0, B0, C0)*FRAC
   END FUNCTION
#ENDINLINE
