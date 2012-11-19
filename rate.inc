#INLINE F90_GLOBAL

  ! generic reaction rate variables
  REAL(dp) kro2no, kro2ho2, kapho2, kapno, kro2no3, kno3al, kdec, &
    krosec, kalkoxy, kalkpxy, kroprim
  REAL(dp) k298ch3o2
  REAL(dp) kch3o2
  ! variables for calculation of kfpan and kbpan
  REAL(dp) kfpan, kbpan
  REAL(dp) kc0, kci, krc, fcc, fc
  REAL(dp) kd0, kdi, krd, fcd, fd
  ! variables for calculation of kmt01
  REAL(dp) kmt01
  REAL(dp) k10, k1i, kr1, fc1, f1
  ! variables for calculation of kmt02
  REAL(dp) kmt02
  REAL(dp) k20, k2i, kr2, fc2, f2
  ! variables for calculation of kmt03
  REAL(dp) kmt03
  REAL(dp) k30, k3i, kr3, fc3, f3
  ! variables for calculation of kmt04
  REAL(dp) kmt04
  REAL(dp) k40, k4i, kr4, fc4, f4
  ! variables for calculation of kmt05
  REAL(dp) kmt05
  ! variables for calculation of kmt06
  REAL(dp) kmt06
  ! variables for calculation of kmt07
  REAL(dp) kmt07
  REAL(dp) k70, k7i, kr7, fc7, f7
  ! variables for calculation of kmt08
  REAL(dp) kmt08
  REAL(dp) k80, k8i, kr8, fc8, f8
  ! variables for calculation of kmt09
  REAL(dp) kmt09
  REAL(dp) k90, k9i, kr9, fc9, f9
  ! variables for calculation of kmt10
  REAL(dp) kmt10
  REAL(dp) k100, k10i, kr10, fc10, f10
  ! variables for calculation of kmt11
  REAL(dp) kmt11
  REAL(dp) k1,k2,k3,k4
  ! variables for calculation of kmt12
  REAL(dp) kmt12
  REAL(dp) k0, ki, ssign,f
  ! variables for calculation of kmt13
  REAL(dp) kmt13
  REAL(dp) k130, k13i, kr13, fc13, f13
  ! variables for calculation of kmt14
  REAL(dp) kmt14
  REAL(dp) k140, k14i, kr14, fc14, f14
  ! variables for calculation of kmt15
  REAL(dp) kmt15
  ! variables for calculation of kmt16
  REAL(dp) kmt16
  REAL(dp) k160, k16i, kr16, fc16, f16
  ! variables for calculation of kmt17
  REAL(dp) kmt17
  REAL(dp) kmt18
  REAL(dp) f12, f15, f17, fc12, fc15, fc17, k120
  REAL(dp) k12i, k150, k15i, k170, k17i, kr12, kr15
  REAL(dp) kr17, nc, nc1, nc10, nc12, nc13, nc14
  REAL(dp) nc15, nc16, nc17, nc2
  REAL(dp) nc3, nc4, nc7, nc8, nc9, ncd
#ENDINLINE
#INLINE F90_RCONST
  THETA = DBLE( ZENITH ( REAL(LAT), -REAL(LON), JDAY))
#ENDINLINE F90_RCONST
