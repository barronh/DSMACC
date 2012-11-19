MODULE constants

  USE dsmacc_Precision, ONLY: dp
  IMPLICIT NONE
  PRIVATE dp

  INTEGER, PARAMETER :: mnsp=250, mre=2000
  INTEGER i
  ! variables for zenith routine which calculates zenith angle
  REAL(dp) theta, secx, cosx
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
  REAL(dp) k0, ki, x, ssign,f
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
  ! variables for calculation of photolysis reaction rates
  ! J increased to 1100. Upto 1000 for inorganics/organics, 1000 onwards halogens
  ! GEOS_CHEM paramterization
  REAL(dp)  l(1500), mm(1500), nn(1500), j(1500)
  INTEGER k

CONTAINS

  !***************************************************************************

  SUBROUTINE mcm_constants(time, temp, M, N2, O2, RO2, H2O)
    ! calculates rate constants from arrhenius informtion
    USE dsmacc_Global,  ONLY:LAT, LON, JDAY, SZAS, SVJ_TJ, bs,cs,ds, jfactno2, jfacto1d
    REAL(dp) time, temp, M, N2, O2, RO2, H2O, THETA, TIME2, LAT2
    REAL*8 y,dy,x,tmp(19), tmp2(19),b(19),c(19),d(19)
    integer i,n,jl
    INTEGER LK
    include '../tuv_old/params'
  !  OPEN(13,file='photolysis.txt', status='old')

!    PI=3.1415926
    ! calculate zenith angle for time of day
!    CALL zenith(theta, secx, cosx, time)
! Time2 is local time in hours
    Time2=mod(Time/(60.*60.), 24.)
    IF (TIME2 .LT. 0) TIME2=TIME2+24.
    LAT2=LAT
    THETA=ZENANG(int(jday),Time2,LAT2)*180./PI
 
    ! ************************************************************************
    ! define generic reaction rates.
    ! ************************************************************************

    ! constants used in calculation of reaction rates
    !M  = 2.55E+19
    N2 = 0.79*M
    O2 = 0.2095*M

    ! kro2no : ro2      + no      = ro      + no2
    !        : ro2      + no      = rono2
    ! mcm 2009
    kro2no    = 2.54d-12*EXP(360.0/temp)

    ! kro2ho2: ro2      + ho2     = rooh    + o2
    ! mcm protocol [1997]
    kro2ho2   = 2.91d-13*EXP(1300.0/temp)

    ! kapho2 : rcoo2    + ho2     = products
    ! mcm protocol [1997]
    kapho2    = 4.30d-13*EXP(1040.0/temp)

    ! kapno  : rcoo2    + no      = products
    ! mej [1998]
    kapno = 8.10d-12*EXP(270.0/temp)

    ! kro2no3: ro2      + no3     = products
    ! mcm protocol [1997]
    kro2no3   = 2.50d-12

    ! kno3al : no3      + rcho    = rcoo2   + hno3
    ! mcm protocol [1997]
    kno3al    = 1.44d-12*EXP(-1862.0/temp)

    ! kdec   : ro                 = products
    ! mcm protocol [1997]
    kdec      = 1.00d+06
    krosec = 1.50e-14*EXP(-200.0/temp)

    kalkoxy=6.00d-14*EXP(-550.0/temp)*o2
    kalkpxy=1.50d-14*EXP(-200.0/temp)*o2
 
    ! -------------------------------------------------------------------
    ! complex reactions
    ! -------------------------------------------------------------------

    ! kfpan kbpan
    ! formation and decomposition of pan
    ! iupac 1997
    kc0     = 2.70d-28*m*(temp/298.0)**(-7.1)
    kci     = 1.21d-11*(temp/298.0)**(-0.9)
    krc     = kc0/kci
    fcc     = 0.30
    fc      = 10**(dlog10(fcc)/(1+(dlog10(krc))**2))
    kfpan   = (kc0*kci)*fc/(kc0+kci)

    kd0     = 4.90d-03*m*EXP(-12100.0/temp)
    kdi     = 5.40d+16*EXP(-13830.0/temp)
    krd     = kd0/kdi
    fcd     = 0.30
    fd      = 10**(dlog10(fcd)/(1+(dlog10(krd))**2))
    kbpan   = (kd0*kdi)*fd/(kd0+kdi)

    ! kmt01  : o        + no      = no2
    ! iupac 1997, 2006 gives Fc=0.85, otherwise the same
    k10     = 1.00d-31*m*(temp/300.0)**(-1.6)

    k1i     = 3.00d-11*(temp/300.0)**(0.3)
    kr1     = k10/k1i
    fc1     = EXP(-temp/1850.0)
    f1      = 10**(dlog10(fc1)/(1+(dlog10(kr1))**2))
    kmt01   = (k10*k1i)*f1/(k10+k1i)

    ! kmt02  : o        + no2     = no3
    ! iupac 2006
    k20     = 1.30d-31*m*(temp/300.0)**(-1.5)
    k2i     = 2.30d-11*(temp/300.0)**(0.24)
    kr2     = k20/k2i
    fc2     = 0.6
    f2      = 10**(dlog10(fc2)/(1+(dlog10(kr2))**2))
    kmt02   = (k20*k2i)*f2/(k20+k2i)

    ! kmt03  : no2      + no3     = n2o5
    ! iupac 2006
    k30     = 3.60d-30*m*(temp/300.0)**(-4.1)
    k3i     = 1.90d-12*(temp/300.0)**(0.2)
    kr3     = k30/k3i
    fc3     = 0.35
    f3      = 10**(dlog10(fc3)/(1+(dlog10(kr3))**2))
    kmt03   = (k30*k3i)*f3/(k30+k3i)

    ! kmt04  : n2o5               = no2     + no3
    ! iupac 2006
    k40     = 1.30d-03*m*(temp/300.0)**(-3.5)*EXP(-11000.0/temp)
    k4i     = 9.70d+14*(temp/300.0)**(0.1)*EXP(-11080.0/temp)
    kr4     = k40/k4i
    fc4     = 0.35
    f4      = 10**(dlog10(fc4)/(1+(dlog10(kr4))**2))
    kmt04   = (k40*k4i)*f4/(k40+k4i)

    ! kmt05  : oh       + co(+o2) = ho2     + co2
    ! iupac 2006
    kmt05  = (1 + m/4.2d19)

    ! kmt06  : ho2      + ho2     = h2o2    + o2
    ! water enhancement factor
    ! iupac 1992

    kmt06  = 1 + (1.40d-21*EXP(2200.0/temp)*h2o)

    ! kmt06  = 1 + (2.00d-25*EXP(4670.0/temp)*h2o)
    ! S+R 2005 values

    ! kmt07  : oh       + no      = hono

    ! iupac 2006
    k70     = 7.40d-31*m*(temp/300.0)**(-2.4)
    k7i     = 3.30d-11*(temp/300.0)**(-0.3)
    kr7     = k70/k7i
    fc7     = 0.81
    f7      = 10**(dlog10(fc7)/(1+(dlog10(kr7))**2))
    kmt07   = (k70*k7i)*f7/(k70+k7i)

    ! kmt08  : oh       + no2     = hno3

    ! iupac 2006
    k80     = 3.30d-30*m*(temp/300.0)**(-3.0)
    k8i     = 4.10d-11
    kr8     = k80/k8i
    fc8     = 0.4
    f8      = 10**(dlog10(fc8)/(1+(dlog10(kr8))**2))
    kmt08   = (k80*k8i)*f8/(k80+k8i)

    ! kmt09  : ho2      + no2     = ho2no2
    ! iupac 1997

    k90     = 1.80d-31*m*(temp/300.0)**(-3.2)
    ! k90     = 2.0d-31*m*(temp/300.0)**(-3.4)
    k9i     = 4.70d-12
    ! k9i     = 2.90d-12*(temp/300.0)**(-1.1)
    kr9     = k90/k9i
    fc9     = 0.6
    f9      = 10**(dlog10(fc9)/(1+(dlog10(kr9))**2))
    kmt09   = (k90*k9i)*f9/(k90+k9i)

    ! kmt10  : ho2no2             = ho2     + no2
    ! iupac 1997

    k100     = 4.10d-05*m*EXP(-10650.0/temp)
    k10i     = 5.70d+15*EXP(-11170.0/temp)
    kr10     = k100/k10i
    fc10     = 0.5
    f10      = 10**(dlog10(fc10)/(1+(dlog10(kr10))**2))
    kmt10    = (k100*k10i)*f10/(k100+k10i)

    ! kmt10    = kmt09/(3.9d-27*EXP(10125/temp))

    ! kmt11  : oh       + hno3    = h2o     + no3
    ! iupac 2006

    k1     = 2.40d-14*EXP(460.0/temp)
    k3     = 6.50d-34*EXP(1335.0/temp)
    k4     = 2.70d-17*EXP(2199.0/temp)
    k2     = (k3*m)/(1+(k3*m/k4))
    kmt11  = k1 + k2

    ! kmt12 iupac 2006

    k0 = 4.50d-31*((temp/300.0)**(-3.9))*m
    ki = 1.30d-12*((temp/300.0)**(-0.7))
    fc = 0.525
    x = 1.0d+0
    ssign = dsign(x,(k0-ki))
    f=10**(dlog10(fc)/(1.0+(ssign*(ABS(dlog10(k0/ki)))**(2.0))))
    kmt12=(k0*ki*f)/(k0+ki)

    ! kmt13  : ch3o2    + no2     = ch3o2no2
    ! iupac 2006

    k130     = 2.50d-30*((temp/300.0)**(-5.5))*m
    k13i     = 1.80d-11
    kr13     = k130/k13i
    fc13     = 0.36
    f13      = 10**(dlog10(fc13)/(1+(dlog10(kr13))**2))
    kmt13    = (k130*k13i)*f13/(k130+k13i)

    ! kmt14  : ch3o2no2           = ch3o2   + no2
    ! iupac 2006

    k140     = 9.00d-05*EXP(-9690.0/temp)*m
    k14i     = 1.10d+16*EXP(-10560.0/temp)
    kr14     = k140/k14i
    fc14     = 0.6
    f14      = 10**(dlog10(fc14)/(1+(dlog10(kr14))**2))
    kmt14    = (k140*k14i)*f14/(k140+k14i)

    ! kmt15 iupac 2006

    k0=8.60d-29*((temp/300.0)**(-3.1))*m
    ki=9.00d-12*((temp/300.0)**(-0.85))
    fc=0.48
    x = 1.0d+0
    ssign = dsign(x,(k0-ki))
    f=10**(dlog10(fc)/(1.0+(ssign*(ABS(dlog10(k0/ki)))**(2.0))))
    kmt15=(k0*ki*f)/(k0+ki)

    ! kmt16  :  oh  +  c3h6
    ! iupac 2006

    k160     = 8.00d-27*((temp/300.0)**(-3.5))*m
    k16i     = 3.00d-11*((temp/300.0)**(-1.0))
    kr16     = k160/k16i
    fc16     = 0.5
    f16      = 10**(dlog10(fc16)/(1+(dlog10(kr16))**2))
    kmt16    = (k160*k16i)*f16/(k160+k16i)

    ! kmt17 iupac 2006

    k0 = 5.00d-30*((temp/300.0)**(-1.5))*m
    ki = 1.00d-12
    fc = 0.37
    x = 1.0d+0
    ssign = dsign(x,(k0-ki))
    f=10**(dlog10(fc)/(1.0+(ssign*(ABS(dlog10(k0/ki)))**(2.0))))
    kmt17=(k0*ki*f)/(k0+ki)

    !       mcm 2001

    kroprim  = 3.70d-14*EXP(-460.0/temp)
    krosec   = 1.80d-14*EXP(-260.0/temp)

    !    added by Mr.Wei Zhou 2012/10/09
    !    ref. mcm.leeds.ac.uk/MCM/parameters/simple.htt
    k298ch3o2 = 3.5D-13

    !    added by Mr.Wei Zhou
    !    ref. mcm.leeds.ac.uk/MCM/parameters/simple.htt
    kch3o2 = 1.03D-13*EXP(365/temp)
    ! ************************************************************************
    ! define photolysis reaction rates from cubic splines of the TUV output
    ! ************************************************************************

	if (theta .le. 90) then 
	 n=19
	 do jl=1,kj
                       do i=1,19 
                               tmp(i)=szas(i)
                               tmp2(i)=svj_tj(i,jl)

                               b(i)=bs(i,jl)
                               c(i)=cs(i,jl)
                               d(i)=ds(i,jl)
                       enddo

                       if (jl .eq.  2)  j(1)=seval(n,theta,tmp, tmp2, b,c,d) ! O3->O1D
                       if (jl .eq.  3)  j(2)=seval(n,theta,tmp, tmp2, b,c,d) ! O3->O3P
                       if (jl .eq. 11)  j(3)=seval(n,theta,tmp, tmp2, b,c,d) ! H2O2->2*OH
                       if (jl .eq.  4)  j(4)=seval(n,theta,tmp, tmp2, b,c,d) ! NO2->NO+O3P
	                   if (jl .eq.  5)  j(5)=seval(n,theta,tmp, tmp2, b,c,d) ! NO3->NO+O2
                       if (jl .eq.  6)  j(6)=seval(n,theta,tmp, tmp2, b,c,d) ! NO3->NO2+O3P
                       if (jl .eq. 12)  j(7)=seval(n,theta,tmp, tmp2, b,c,d) ! HNO2->OH+NO
                       if (jl .eq. 13)  j(8)=seval(n,theta,tmp, tmp2, b,c,d) ! HNO3->NO2+OH
                       if (jl .eq. 14)  j(1300)=seval(n,theta,tmp, tmp2, b,c,d) 
                       if (jl .eq. 15)  j(11)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 16)  j(12)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 17)  j(13)=seval(n,theta,tmp, tmp2, b,c,d)
 		               if (jl .eq. 20)  j(14)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 75)  j(15)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 76)  j(16)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 77)  j(17)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 62)  j(18)=seval(n,theta,tmp, tmp2, b,c,d)*0.5
                       if (jl .eq. 62)  j(19)=seval(n,theta,tmp, tmp2, b,c,d)*0.5
                       if (jl .eq. 25)  j(21)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 87)  j(22)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 61)  j(23)=seval(n,theta,tmp, tmp2, b,c,d)*0.5
                       if (jl .eq. 61)  j(24)=seval(n,theta,tmp, tmp2, b,c,d)*0.5
                       if (jl .eq. 21)  j(31)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 23)  j(32)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 22)  j(33)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 24)  j(34)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 60)  j(35)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 26)  j(41)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 27)  j(51)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 64)  j(52)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 91)  j(53)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 64)  j(54)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 94)  j(55)=seval(n,theta,tmp, tmp2, b,c,d)
                       if (jl .eq. 67)  j(56)=seval(n,theta,tmp, tmp2, b,c,d)*0.75
                       if (jl .eq. 67)  j(57)=seval(n,theta,tmp, tmp2, b,c,d)*0.25
!Halogens
                       if (jl .eq. 72)  j(1001)=seval(n,theta,tmp, tmp2, b,c,d) ! HOBr
                       if (jl .eq. 73)  j(1002)=seval(n,theta,tmp, tmp2, b,c,d) ! BrO
                       if (jl .eq. 74)  j(1003)=seval(n,theta,tmp, tmp2, b,c,d) ! Br2
                       if (jl .eq. 50)  j(1004)=seval(n,theta,tmp, tmp2, b,c,d) ! BrNO3->Br+NO3
                       if (jl .eq. 51)  j(1005)=seval(n,theta,tmp, tmp2, b,c,d) ! BrNO3->BrO+NO2
                       if (jl .eq. 30)  j(1006)=seval(n,theta,tmp, tmp2, b,c,d) ! ClNO3->Cl+NO3
                       if (jl .eq. 31)  j(1007)=seval(n,theta,tmp, tmp2, b,c,d) ! ClNO3->ClO+NO2
                       if (jl .eq. 58)  j(1008)=seval(n,theta,tmp, tmp2, b,c,d) ! Cl2->2Cl
 

                       
          enddo
         
	else
           do i=1,1500
              j(i)=0.
           enddo
	endif

        
        if (jfactno2 .eq. 0) jfactno2=1.
        if (jfacto1d .eq. 0) jfacto1d=1. 
        do i=1,1500
           if (i .ne. 1) then 
              j(i)=j(i)*jfactno2
           else
              j(i)=j(i)*jfacto1d
           endif
        enddo
      
	do i=1,1500
		if (j(i) .lt. 0.e0) j(i)=0e0
	enddo
	
!    REWIND(13)
!    READ(13,*)
!    DO lk = 1, 35
!      READ(13,*) k, l(k), mm(k), nn(k)
!     
!      IF (cosx<1.e-10) THEN
!        j(k) = 1.e-30
!      ELSE
!        j(k) = (l(k)*cosx**( mm(k))*EXP(-nn(k)*secx))
!      ENDIF
!  
!    END DO
  
  END SUBROUTINE mcm_constants

  !***************************************************************************

!  SUBROUTINE zenith(theta, secx, cosx, ttime)
!
!    REAL(dp), INTENT(IN)  :: ttime
!    REAL(dp), INTENT(OUT) :: theta, secx, cosx
!    REAL(dp) lat, pi, radian, dec, lha, sinld, cosld
!
!    ! solar declination angle from jday 172
!    dec = 23.44
!    ! latitude
!    lat = 14.0
!    pi = 4.0*ATAN(1.0)
!    ! local hour angle - representing time of day
!    lha = (1.0+ttime/4.32d+4)*pi
!    radian = 180.0/pi
!    lat = lat/radian
!    dec = dec/radian
!    theta = ACOS(COS(lha)*COS(dec)*COS(lat)+SIN(dec)*SIN(lat))
!    sinld = SIN(lat)*SIN(dec)
!    cosld = COS(lat)*COS(dec)
!    cosx = (COS(lha)*cosld)+sinld
!    cosx = COS(theta)
!    secx = 1.0d+0/(cosx+1.0d-30)
!
!  END SUBROUTINE zenith

  !***************************************************************************

      REAL FUNCTION UPTAKE(GAMMA,TEMP,AREA,RP,MASS)
  ! INPUT GAMMA=UPTAKE COEFFICIENT (0..1)
  !       TEMP=TEMPERTURE(K)
  !       AREA=SURFACE AREA (M^2/M^3)
  !       MASS=MOLECULAR MASS OF MOLECULE
  !       R=BOLTZMAN'S CONSTANT
  !       RP=RADIUS OF PARTICLE (M)
  !       DG=GAS-PHASE DIFFUSION COEFFICIENT (M^2 S-1)

        REAL*8 GAMMA, TEMP, AREA, RP, MASS
        REAL*8 V, V2, R, DG



        R=8.314 

        DG=2.47D-5

  !     DG TAKEN FROM MOZUREVICH ET EL 1987
       
        V2=3*R*TEMP/MASS
        V=(8*V2/(3*3.1415))**0.5

        UPTAKE=(((RP/DG)+(4/(GAMMA*V)))**-1)*AREA
        
  !      UPTAKE=v*AREA*GAMMA/4

     END FUNCTION UPTAKE

       FUNCTION ZENANG(Nday,Time,Geolat)
 
        REAL(dp) :: Geolat , Time
        INTEGER :: Nday
        REAL(dp) :: ZENANG
        INTENT (IN) Geolat , Time
  !
  ! Local variables
  !
       REAL :: arg , dec , georad , w,torad
!       REAL :: SOLDEC
  !
  !-----------------------------------------------------------------------
  !
  !     a function to calculate the solar zenith angle for a given day,
  !     nday, local time, ntime & geographical latitude, geolat.
  !
  !     the solar zenith angle, also called the zenith distance, is the
  !     angle between the local zenith and the line joining the observer
  !     and the sun.
  !
  !     it is returned in radians.
  !
  !-----------------------------------------------------------------------
  !
  !     the hour angle, w, changes by 15 degrees for every hour, it is 0
  !     at local noon. it is computed here in radians.
  !
  !     the declination, dec, is the angular position of the sun at
  !     solar noon with respect to the plane of the equator, north is
  !     positive, it is a function of the time of year.
  !     it is required in radians.
  !
  !-----------------------------------------------------------------------
  !
  !     calculate the hour angle in radians.
	torad=4.0*ATAN(1.0)/180.
        w = torad*15.0*(Time-12.)
  !
  !     calculate declination in radians.
         dec = SOLDEC(Nday)
  !
  !     geographic latitude in degrees converted to radians.
         georad = Geolat*torad
  !
  !     we know the geographical latitude, so from a simple geometrical
  !     relationship we can now calculate the solar zenith angle.
         arg = SIN(dec)*SIN(georad) + COS(dec)*COS(georad)*COS(w)
  !
         ZENANG = ACOS(arg)
  !
  !-----------------------------------------------------------------------
  !
       END FUNCTION ZENANG
 
       FUNCTION SOLDEC(Nday)
 
        IMPLICIT NONE
  !
  ! Dummy arguments
  !
         INTEGER :: Nday
         REAL :: SOLDEC, pi
         INTENT (IN) Nday
  !
  ! Local variables
  !
         REAL :: dayang
         REAL :: REAL

         PI=3.14159265
  !
  !-----------------------------------------------------------------------
  !
  !     a function to calculate the solar declination using the
  !     relation given by j. w. spencer, fourier series representation
  !     of the position of the sun, search vol. 2 (5), pp. 172, 1971.
  !
  !     the declination, dec, is the angular position of the sun at
  !     solar noon with respect to the plane of the equator, north is
 !     positive, it is a function of the time of year.
  !
  !     it is returned in radians.
  !
  !-----------------------------------------------------------------------
  !
  !     nday is the day number, 1 on january 1st, 365 on december 31st.
  !
  !     dayang is the day angle, dayang=2*pi*(nday-1)/365, in radians,
  !     from spencer.
  !
  !-----------------------------------------------------------------------
  !
  !     calculate day angle in radians.
         dayang = 2.0*pi*REAL(Nday-1)/365.0
  !
  !     calculate declination in radians.
         SOLDEC = (0.006918-0.399912*COS(dayang)+0.070257*SIN(dayang)      &
                & -0.006758*COS(2.0*dayang)+0.000907*SIN(2.0*dayang)       &
                & -0.002697*COS(3.0*dayang)+0.00148*SIN(3.0*dayang))
  !
  !-----------------------------------------------------------------------
  !
         END FUNCTION SOLDEC

                   double precision function seval(n, u, x, y, b, c, d)
!======================================================================
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!========================================================================
      integer n
      double precision  u, x(n), y(n), b(n), c(n), d(n)
      integer i, j, k
      double precision dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
!
!  binary search
!
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
!
!  evaluate spline
!
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end function seval

      subroutine polint(f,a,n,x,r)
!----------------------------------------------------------
! service program for fintr
!----------------------------------------------------------
       real*8 r,al,x,a,f 
       integer  i,j,k,l,n	
       dimension f(n),a(n)
       r=0.0
       do 1 j=1,n
       al=1.0
       do 2 i=1,n
          if (i-j) 3,2,3
 3        al=al*(x-a(i))/(a(j)-a(i))
 2     continue
 1     r=r+al*f(j)
       return
       end subroutine polint


END MODULE constants
