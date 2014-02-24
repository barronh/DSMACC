#INLINE F90_RATES
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

#ENDINLINE
