#INLINE F90_GLOBAL 
 REAL(dp):: PRESS_PA, PRESS_hPA, M, N2, O2, RO2, SO4GASFRAC = 1.0
#ENDINLINE {above lines go into MODULE KPP_ROOT_Global}

#INLINE F90_RATES
 SUBROUTINE UPDATEJ()
     !
     ! UPDATEJ calls FAST-JX for all GEOS-Chem species
     ! and applies corrections in PHOTRATES_ADJ following
     ! v11-01-public-release code
     !
     INTEGER :: JIDX
     INTEGER, SAVE :: RXN_O2 = 0, RXN_O3_1 = 0, RXN_O3_2a = 0,&
                     RXN_O3_2b = 0, RXN_H2SO4 = 0, RXN_NO2 = 0,&
                     RXN_NO = 0, RXN_NO3 = 0, RXN_N2O = 0
     REAL*8 :: ITEMPK, C_O2, C_N2, C_H2, C_H2O, RO1DplH2O, RO1DplH2, RO1D
     LOGICAL, SAVE :: FIRST = .true.
     IF (FIRST) THEN
         PHOTOL(1) = J(1)
         DO JIDX=1,NJREACT
             SELECT CASE (TRIM(ADJUSTL(JLABELS(JIDX))))
             CASE ('O2        PHOTON    O         O')
               RXN_O2 = JIDX
               print *,'Found RXN_O2'
               print *,JLABELS(JIDX)
             CASE ('SO4       PHOTON    SO2       OH        OH')
               print *,'Found RXN_H2SO4'
               print *,JLABELS(JIDX)
               RXN_H2SO4 = JIDX
             CASE ('O3        PHOTON    O2        O')
               print *,'Found O3 + hv = O2 + O'
               print *,JLABELS(JIDX)
               RXN_O3_1 = JIDX
             CASE ('O3        PHOTON    O2        OH       OH')
               print *,'Found O3 + hv = 2 OH via O(1D) + H2O'
               print *,JLABELS(JIDX)
               RXN_O3_2a = JIDX
             CASE ('O3        PHOTON    O2        HO2       OH')
               print *,'Found O3 + hv = HO2 + OH via O(1D) + H2'
               print *,JLABELS(JIDX)
               RXN_O3_2b = JIDX
             END SELECT
         END DO
         FIRST = .false.
     ENDIF
     DO JIDX=1,NJREACT
         PHOTOL(JIDX) = J(JIDX)       
     END DO
     IF (UCX) THEN
         IF (RXN_H2SO4.gt.0) THEN
             !PHOTOL(RXN_H2SO4) = PHOTOL(RXN_H2SO4) * FRAC
             PHOTOL(RXN_H2SO4) = J(RXN_H2SO4) * SO4GASFRAC
         ENDIF
         IF (RXN_O3_1.gt.0) THEN
           !ZPJ(L,RXN_O3_1,I,J) = ZPJ(L,RXN_O3_1,I,J) - ZPJ(L,RXN_O3_2a,I,J)
           PHOTOL(RXN_O3_1) = J(RXN_O3_1) - J(RXN_O3_2a)
         ENDIF
     ELSE
         ! Inverse temperature [K-1]
         ITEMPK    = 1.0d0 / TEMP

         ! Set species concentrations [molec/m3] ???
         C_O2      = 0.2095d+0 * NUMDEN
         C_N2      = 0.7808d+0 * NUMDEN

         ! Added H2 concentration (bhh, jmao, eam, 7/18/11)
         ! Seasonal variability of H2 may be important, 
         ! but not included in this update (bhh, jmao, eam, 7/18/11)
         C_H2      = 0.5000d-6 * NUMDEN
         C_H2O = H2O
         RO1DplH2O = 1.63d-10 * EXP(  60.0d0 * ITEMPK ) * C_H2O

         RO1DplH2  = 1.2e-10                                * C_H2

         RO1D      = RO1DplH2O + &
                     RO1DplH2 + &
                     2.15d-11 * EXP( 110.0d0 * ITEMPK ) * C_N2 + &
                     3.30d-11 * EXP(  55.0d0 * ITEMPK ) * C_O2

         ! Prevent div-by-zero
         IF ( RO1D > 0.0d0 ) THEN
             ! RXN_O3_2a specifies: O3 + hv -> O2 + O(1D) #1
             !ZPJ(L,RXN_O3_2a,I,J) = ZPJ(L,RXN_O3_2a,I,J) * RO1DplH2O / RO1D
             PHOTOL(RXN_O3_2a) = PHOTOL(RXN_O3_2a) * RO1DplH2O / RO1D
             
             ! RXN_O3_2b specifies: O3 + hv -> O2 + O(1D) #2
             !ZPJ(L,RXN_O3_2b,I,J) = ZPJ(L,RXN_O3_2b,I,J) * RO1DplH2  / RO1D
             PHOTOL(RXN_O3_2b) = PHOTOL(RXN_O3_2b) * RO1DplH2  / RO1D
         ENDIF
     ENDIF
 END SUBROUTINE UPDATEJ
#ENDINLINE

#INLINE F90_RCONST
 M = PRESS/(8.314*TEMP)*6.022E23*1e-6 ! molecule/cm3
 NUMDEN = M
 PRESS_PA = PRESS
 PRESS_hPA = PRESS / 100.
 CALL UPDATEJ()
#ENDINLINE

#LOOKATALL;
