#INLINE F90_RATES
  SUBROUTINE UPDATE_TUV() 
     REAL(kind=dp) alta
     !    Set up the photolysis rates
     !    First calculate pressure altitude from altitude
     WRITE (6,*) 'hvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhv'
     WRITE (6,*) 'Using TUV to calculate photolysis rates as a function of O3 column, SZA, ALT, and TEMP'
     alta=(1-(press/1013.25)**0.190263)*288.15/0.00198122*0.304800/1000.
 
     if (o3col .eq. 0) then 
       o3col=260.
       write (6,*) 'Ozone column not specified using 260 Dobsons'
     else
       write (6,*) 'Ozone column =', o3col
     endif
       
!    Calculate the photolysis rates for the run
     call set_up_photol(O3col, alta, temp, bs,cs,ds,szas,svj_tj)
     WRITE (6,*) 'Photolysis rates calculated'
     WRITE (6,*) 'hvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhvhv'
  END SUBROUTINE

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION RACM_TROE(A0, B0, A1, B1)
    REAL A0, B0, A1, B1
    RACM_TROE = CMAQ_10 ( A0, -B0, 0.0, A1, -B1, 0.0, 0.6, 1.)
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION RACM_TROE_EQUIL(A0, B0, A1, B1, A2, C2)
    REAL A0, B0, A1, B1, A2, C2
    RACM_TROE_EQUIL = CMAQ_1to4(A2, 0.0, C2) * CMAQ_10 ( A0, -B0, 0.0, A1, -B1, 0.0, 0.6, 1.)
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION RACM_THERMAL(A0, B0)
    REAL A0, B0
    REAL, PARAMETER :: C0 = 0.
!   RACM2 reaction rates have the form K = A * EXP(-B / T)
!   KPP ARR reaction rates have the form   K = A * (T/300.0)**C * EXP(-B/T) 
!   
!   Translation adds a 0 C
    RACM_THERMAL = ARR(A0, B0, C0)
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION RACM_THERMAL_T2(A0, B0)
    REAL A0, B0
    REAL, PARAMETER :: C0 = 0.
!   RACM2 reaction rates have the form K = A * EXP(-B / T)
!   KPP ARR reaction rates have the form   K = A * (T/300.0)**C * EXP(-B/T) 
!   
!   Translation adds a 0 C
    RACM_THERMAL_T2 = DBLE(A0)*temp**2*exp(-DBLE(B0)/temp)
  END FUNCTION

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_STD(A0, B0, C0)
    REAL A0, B0, C0

!   GEOS Chem reaction rates have the form K = A * (300 / T)**B * EXP(C / T)
!   KPP ARR reaction rates have the form   K = A * (T/300.0)**C * EXP(-B/T) 
!   
!   Translation reorders B and C and changes both their signs
    GEOS_STD = ARR(A0, -C0, -B0)
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
    K0 = GEOS_STD(A0, B0, C0)
    GEOS_Y = K0 * (1 + .6 * (PRESS * 100.) / 101325.)
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
    REAL(kind=dp) O3J, T3I
    T3I = 1.0_dp/TEMP
    GEOS_JO3 = O3J * &
               1.45d-10 * EXP( 89.0_dp * T3I) * H2O / &
               ( 1.45d-10 * EXP( 89.0_dp * T3I) * H2O + &
                 2.14d-11 * EXP(110.0_dp * T3I) * N2 + &
                 3.20d-11 * EXP( 70.0_dp * T3I) * O2 &
               )
  END FUNCTION
  
  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION GEOS_JHNO4(HNO4J)
    REAL(kind=dp) HNO4J
    IF (HNO4J .GT. 0.d0) THEN
      GEOS_JHNO4 = HNO4J + 1D-5
    ELSE
      GEOS_JHNO4 = HNO4J
    ENDIF
    
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

  REAL(kind=dp) FUNCTION CMAQ_1to4(A0, B0, C0)
    REAL A0, B0, C0

!   CMAQ reaction rates have the form K = A * (T/300.0)**B * EXP(-C/T) 
!   KPP ARR reaction rates have the form K = A * (T/300.0)**C * EXP(-B/T) 
!   
!   Translation reorders B and C
    CMAQ_1to4 = ARR(A0, C0, B0)
  END FUNCTION
  
  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION CMAQ_5(A0, B0, C0, Kf)
    REAL A0, B0, C0
    REAL(kind=dp) K1, Kf
    K1 = CMAQ_1to4(A0, B0, C0)
    CMAQ_5 = Kf / K1
  END FUNCTION

  !---------------------------------------------------------------------------
  
  REAL(kind=dp) FUNCTION CMAQ_6(A0, B0, C0, Kf)
    REAL A0, B0, C0
    REAL(kind=dp) K1, Kf
    K1 = CMAQ_1to4(A0, B0, C0)
    CMAQ_6 = Kf * K1
  END FUNCTION

  !---------------------------------------------------------------------------
  
  REAL(kind=dp) FUNCTION CMAQ_7(A0, B0, C0)
    REAL A0, B0, C0
    REAL(kind=dp) K0
    K0 = CMAQ_1to4(A0, B0, C0)
    CMAQ_7 = K0 * (1 + .6 * PRESS * 100.)
  END FUNCTION

  !---------------------------------------------------------------------------
  
   REAL(kind=dp) FUNCTION CMAQ_8(A0, C0, A2, C2, A3, C3)
      REAL A0, C0, A2, C2, A3, C3
      REAL(kind=dp) K0, K2, K3            
      K0 = DBLE(A0) * EXP(-DBLE(C0) / TEMP)
      K2 = DBLE(A2) * EXP(-DBLE(C2) / TEMP)
      K3 = DBLE(A3) * EXP(-DBLE(C3) / TEMP)
      K3 = K3 * M
      CMAQ_8 = K0 + K3 / (1.0_dp + K3 / K2 )
   END FUNCTION CMAQ_8

  !---------------------------------------------------------------------------

   REAL(kind=dp) FUNCTION CMAQ_9(A1, C1, A2, C2) 
      REAL*8 A1, C1, A2, C2
      REAL(kind=dp) K1, K2      
      K1 = DBLE(A1) * EXP(-DBLE(C1) / TEMP)
      K2 = DBLE(A2) * EXP(-DBLE(C2) / TEMP)
      CMAQ_9 = K1 + K2 * M
   END FUNCTION CMAQ_9 

  !---------------------------------------------------------------------------

   REAL(kind=dp) FUNCTION CMAQ_10 ( A0, B0, C0, A1, B1, C1, CF, N)
      REAL A0, B0, C0, A1, B1, C1, CF, N
      REAL(kind=dp) K0, K1     
      K0 = CMAQ_1to4(A0, B0, C0)
      K1 = CMAQ_1to4(A1, B1, C1)
      K0 = K0 * M
      K1 = K0 / K1
      CMAQ_10 = (K0 / (1.0_dp + K1))*   &
           DBLE(CF)**(1.0_dp / (1.0_dp / DBLE(N) + (LOG10(K1))**2))
   END FUNCTION CMAQ_10

  !---------------------------------------------------------------------------

   REAL(kind=dp) FUNCTION OH_CO ( A0, B0, C0, A1, B1, C1, CF, N)
      REAL A0, B0, C0, A1, B1, C1, CF, N
      REAL(kind=dp) K0, K1     
      K0 = CMAQ_1to4(A0, B0, C0)
      K1 = CMAQ_1to4(A1, B1, C1)
      K0 = K0
      K1 = K0 / (K1 / M)
      OH_CO = (K0 / (1.0_dp + K1))*   &
           DBLE(CF)**(1.0_dp / (1.0_dp / DBLE(N) + (LOG10(K1))**2))
   END FUNCTION OH_CO

  !---------------------------------------------------------------------------

  REAL(kind=dp) FUNCTION TUV_J(IJ, THETA)
    USE model_Global,  only: BS, CS, DS, SZAS, SVJ_TJ, NSZAS
    REAL(kind=dp) B(NSZAS), C(NSZAS), D(NSZAS), TMP_SVJ_TJ(NSZAS), &
  &               TMP_SZAS(NSZAS), THETA
  ! IJ is the integer index of the TUV photolysis calculation
  ! THETA is the current solar zenith angle
    INTEGER IJ, THIS_CSZA

    REAL SEVAL ! function from TUV
    EXTERNAL SEVAL ! function from TUV
    INTEGER THIS_SZA ! temporary variables
    IF (THETA .LT. 90.D0) THEN
      DO THIS_CSZA = 1, NSZAS
        B(THIS_CSZA) = BS(THIS_CSZA, IJ)
        C(THIS_CSZA) = CS(THIS_CSZA, IJ)
        D(THIS_CSZA) = DS(THIS_CSZA, IJ)
        TMP_SZAS(THIS_CSZA) = SZAS(THIS_CSZA)
        TMP_SVJ_TJ(THIS_CSZA) = SVJ_TJ(THIS_CSZA, IJ)
      ENDDO
      
      TUV_J = SEVAL(NSZAS, THETA, TMP_SZAS, TMP_SVJ_TJ, B, C, D, .false.)
      if (.false.) then
          write(*,*) 'MP'
          write(*,*) 'I,THETA,J:', IJ, THETA, TUV_J
          write(*,8879) 'B     :', B
          write(*,8879) 'C     :', C
          write(*,8879) 'D     :', D
          write(*,8879) 'SZAS  :', TMP_SZAS
          write(*,8879) 'SVJ_TJ:', TMP_SVJ_TJ
          
8879     FORMAT(1A6,100000(E26.17))
          TUV_J = SEVAL(NSZAS, THETA, TMP_SZAS, TMP_SVJ_TJ, B, C, D, .true.)

      endif
    
      IF (TUV_J .LT. 0.d0) TUV_J = 0.d0
    ELSE
      TUV_J = 0.d0
    ENDIF

  END FUNCTION

  REAL(kind=dp) FUNCTION MZ4_TROE ( A0, B0, A1, B1, factor)
      REAL(kind=dp) A0, B0, factor, A1, B1
      REAL(kind=dp) ko, kinf, xpo
      ko = A0 * (300.d0 / TEMP)**B0
      kinf = A1 * (300.d0 / TEMP)**B1
      xpo  = ko * m / kinf
      MZ4_TROE = ko / (1. + xpo)
      xpo  = log10( xpo )
      xpo  = 1. / (1. + xpo*xpo)
      MZ4_TROE = MZ4_TROE * factor**xpo
  END FUNCTION MZ4_TROE

  REAL(kind=dp) FUNCTION MZ4_USR1()
      MZ4_USR1 = 6.d-34 * (300.d0/TEMP)**2.4
  END FUNCTION MZ4_USR1

  REAL(kind=dp) FUNCTION MZ4_USR2()
      MZ4_USR2 = MZ4_TROE(8.5d-29, 6.5d0, 1.1d-11, 1.d0, .6d0)
  END FUNCTION MZ4_USR2

  REAL(kind=dp) FUNCTION MZ4_USR3()
    MZ4_USR3 = MZ4_USR2() * 3.333d26 * exp( -10990.d0/TEMP )
  END FUNCTION MZ4_USR3

  REAL(kind=dp) FUNCTION MZ4_USR4()
    MZ4_USR4 = MZ4_TROE(2.0d-30, 3.0d0, 2.5d-11, 0.d0, .6d0)
  END FUNCTION MZ4_USR4

  REAL(kind=dp) FUNCTION MZ4_USR5()
    REAL(kind=dp) KO, TINV
    
    TINV = 1/TEMP
    ko = M * 6.5e-34 * exp( 1335.*tinv )
    ko = ko / (1. + ko/(2.7e-17*exp( 2199.*tinv )))
    MZ4_USR5 = ko + 2.4e-14*exp( 460.*tinv )
  END FUNCTION MZ4_USR5

  REAL(kind=dp) FUNCTION MZ4_USR6()
    MZ4_USR6 = MZ4_TROE(1.8d-31, 3.2d0, 4.7d-12, 1.4d0, .6d0)
  END FUNCTION MZ4_USR6

  REAL(kind=dp) FUNCTION MZ4_USR7()
    MZ4_USR7 = MZ4_USR6() * exp( -10900./TEMP )/ 2.1d-27
  END FUNCTION MZ4_USR7

  REAL(kind=dp) FUNCTION MZ4_USR8()
    real, parameter ::  boltz    = 1.38044e-16      ! erg/k
    MZ4_USR8 = 1.5d-13 * (1. + 6.d-7 * BOLTZ * M * TEMP)
  END FUNCTION MZ4_USR8

  REAL(kind=dp) FUNCTION MZ4_USR9()
    real(kind = dp) ko, kinf, fc, tinv
    tinv = 1.d0/TEMP
    ko   = 2.3d-13 * exp( 600.d0*tinv )
    kinf = 1.7d-33 * m * exp( 1000.d0*tinv )
    fc   = 1.d0 + 1.4d-21 * H2O * exp( 2200.d0*tinv )
    
    MZ4_USR9 = (ko + kinf) * fc
  END FUNCTION MZ4_USR9

  REAL(kind=dp) FUNCTION MZ4_USR10()
    MZ4_USR10 = MZ4_TROE(8.d-27, 3.5d0, 3.d-11, 0d0, .5d0)
  END FUNCTION MZ4_USR10

  REAL(kind=dp) FUNCTION MZ4_USR11()
    MZ4_USR11 = MZ4_TROE(8.5d-29, 6.5d0, 1.1d-11, 1.d0, .6d0)
  END FUNCTION MZ4_USR11

  REAL(kind=dp) FUNCTION MZ4_USR12()
    MZ4_USR12 =  MZ4_USR11() * 1.111d28 * exp( -14000.d0 / TEMP )
  END FUNCTION MZ4_USR12

  REAL(kind=dp) FUNCTION MZ4_USR14()
    MZ4_USR14 = 1.1d-11 * 300.d0/ TEMP / M
  END FUNCTION MZ4_USR14

  REAL(kind=dp) FUNCTION MZ4_USR15()
    MZ4_USR15 = MZ4_USR14() * 1.111d28 *  exp( -14000.d0 / TEMP )
  END FUNCTION MZ4_USR15

  REAL(kind=dp) FUNCTION MZ4_USR16()
    MZ4_USR16 = 0d0 ! aerosol
  END FUNCTION MZ4_USR16

  REAL(kind=dp) FUNCTION MZ4_USR17()
    MZ4_USR17 = 0d0 ! aerosol
  END FUNCTION MZ4_USR17

  REAL(kind=dp) FUNCTION MZ4_USR17a()
    MZ4_USR17a = 0d0 ! aerosol
  END FUNCTION MZ4_USR17a

  REAL(kind=dp) FUNCTION MZ4_USR21()
    MZ4_USR21 = temp**2 * 7.69d-17 * exp( 253.d0/temp )
  END FUNCTION MZ4_USR21

  REAL(kind=dp) FUNCTION MZ4_USR22()
    MZ4_USR22 = 3.82d-11 * exp( -2000.d0/TEMP ) + 1.33d-13
  END FUNCTION MZ4_USR22

  REAL(kind=dp) FUNCTION MZ4_USR23()
    REAL(kind=dp) ko, fc
    
    fc = 3.0d-31 *(300.d0/TEMP)**3.3d0
    ko = fc * m / (1.d0 + fc * m / 1.5d-12) 
    MZ4_USR23 = ko * .6d0**(1. + (log10(fc * m / 1.5d-12))**2.d0)**(-1.d0)
  END FUNCTION MZ4_USR23

  REAL(kind=dp) FUNCTION MZ4_USR24()
    REAL(kind=dp) ko
    ko = 1.d0 + 5.5d-31 * exp( 7460.d0/TEMP ) * M * 0.21d0
    MZ4_USR24 = 1.7d-42 * exp( 7810.d0/TEMP ) * M * 0.21d0 / ko
  END FUNCTION MZ4_USR24

  REAL(kind=dp) FUNCTION MZ4_USR25()
    MZ4_USR25 = 0d0 ! unknown
  END FUNCTION MZ4_USR25

  REAL(kind=dp) FUNCTION MZ4_USR26()
    MZ4_USR26 = 0d0 ! aerosol
  END FUNCTION MZ4_USR26

  REAL(kind=dp) FUNCTION MZ4_USR27()
    MZ4_USR27 = 0d0 ! unknown
  END FUNCTION MZ4_USR27

  REAL(kind=dp) FUNCTION MZ4_USR28()
    MZ4_USR28 = 0d0 ! unknown
  END FUNCTION MZ4_USR28

  !---------------------------------------------------------------------------

  SUBROUTINE SHOW_RCONST()
    INTEGER THIS_RK
    do THIS_RK = 1, 288
      write(*,*) THIS_RK, RCONST(THIS_RK)
    enddo
  END SUBROUTINE

#ENDINLINE
