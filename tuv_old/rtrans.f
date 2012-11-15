* This file contains the following subroutines, related to the solution of
* the equation of radiative transfer in multiple homogeneous layers.
*     rtlink
*     ps2str
*        tridag
*     psndo
*        asymtx
*        chekin
*        fluxes
*        lepoly
*        pravin
*        prtinp
*        prtint
*        qgausn
*        setdis
*        setmtx
*        soleig
*        solve0
*        surfac
*        solvec
*        upbeam
*        zeroal
*        zeroit
*        errmsg
*        sgbco
*        sgbfa
*        sgbsl
*        sgeco
*        sgefa
*        sgesl
*        saxpy
*        sscal
*        sswap
*        t665d
*        t665r
* and the functions
*        dref
*        ratio
*        wrtbad
*        wrtdim
*        tstbad
*        sasum
*        sdot
*        isamax
*        d1mach
*        r1mach
*=============================================================================*

      SUBROUTINE rtlink(nstr, nz, 
     $     iw, ag, zen,
     $     dsdh, nid,
     $     dtrl, 
     $     dto3, 
     $     dto2,
     $     dtso2,
     $     dtno2, 
     $     dtcld, omcld, gcld,
     $     dtaer, omaer, gaer,
     $     edir, edn, eup, fdir, fdn, fup)

      IMPLICIT NONE

      INCLUDE 'params'

* input

      INTEGER nstr
      INTEGER nz, iw
      REAL ag
      REAL zen
      REAL dsdh(0:kz,kz)
      INTEGER nid(0:kz)

      REAL dtrl(kz,kw)
      REAL dto3(kz,kw), dto2(kz,kw), dtso2(kz,kw), dtno2(kz,kw)
      REAL dtcld(kz,kw), omcld(kz,kw), gcld(kz,kw)
      REAL dtaer(kz,kw), omaer(kz,kw), gaer(kz,kw)

* output

      REAL edir(kz), edn(kz), eup(kz)
      REAL fdir(kz), fdn(kz), fup(kz)

* program constants:

      REAL dr
      PARAMETER (dr = pi/180.)

* local:

      REAL dt(kz), om(kz), g(kz)
      REAL dtabs, dtsct, dscld, dacld, dsaer, daaer
      INTEGER i, ii


* specific two ps2str

      REAL ediri(kz), edni(kz), eupi(kz)
      REAL fdiri(kz), fdni(kz), fupi(kz)
      LOGICAL delta
      DATA delta /.true./

*  specific to psndo:

      REAL pmcld, pmray, pmaer
      REAL om1
      INTEGER istr, iu

      INTEGER MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXPHI
      PARAMETER(MAXCLY=151,MAXULV=151)
      PARAMETER(MAXUMU=32,MAXCMU=32)
      PARAMETER(MAXPHI=3)

      INTEGER   NLYR, NUMU
      REAL     ALBEDO, DTAUC( MAXCLY ),
     $     PMOM( 0:MAXCMU, MAXCLY ),
     $     SSALB( MAXCLY ),  
     $     UMU( MAXUMU ), CWT( MAXUMU ), UMU0

      REAL     RFLDIR( MAXULV ), RFLDN( MAXULV ), FLUP( MAXULV ),
     $     U0U( MAXUMU, MAXULV ),
     $     uavgso( maxulv ), uavgup( maxulv ), uavgdn( maxulv ),
     $     sindir( maxulv ), sinup( maxulv ), sindn( maxulv )

*bm  added array LDIF for convenience (sky radiance)
*bm  sine weighted intensity

      REAL irrad
      REAL ldif(MAXUMU, kz)
      REAL sdir(kz), sdn(kz), sup(kz)

*_______________________________________________________________________

* initialize:

      DO 5 i = 1, nz
         fdir(i) = 0.
         fup(i) = 0.
         fdn(i) = 0.
         edir(i) = 0.
         eup(i) = 0.
         edn(i) = 0.
         sdir(i) = 0.
         sup(i) = 0.
         sdn(i) = 0.
 5    CONTINUE

      UMU0 = cos(zen*dr)
      NLYR = nz - 1
      ALBEDO = ag

      DO 10 i = 1, nz - 1

         dscld = dtcld(i,iw)*omcld(i,iw)
         dacld = dtcld(i,iw)*(1.-omcld(i,iw))

         dsaer = dtaer(i,iw)*omaer(i,iw)
         daaer = dtaer(i,iw)*(1.-omaer(i,iw))

         dtsct = dtrl(i,iw) + dscld + dsaer
         dtabs = dtso2(i,iw) + dto2(i,iw) + dto3 (i,iw) +
     >           dtno2(i,iw) + dacld + daaer

 	 dtabs = amax1(dtabs,1./largest)
 	 dtsct = amax1(dtsct,1./largest)

* invert z-coordinate:

         ii = nz - i
         dt(ii) = dtsct + dtabs
         om(ii) = dtsct/(dtsct + dtabs)
           IF(dtsct .EQ. 1./largest) om(ii) = 1./largest
         g(ii) = (gcld(i,iw)*dscld + gaer(i,iw)*dsaer)/dtsct

         IF(nstr .LT. 2) GO TO 10

* DISORD parameters

         OM1 = AMIN1(OM(ii),1.-PRECIS)
         SSALB( II ) = AMAX1(OM1,PRECIS)
         DTAUC( II ) = AMAX1(DT(ii),PRECIS)

*  phase function - assume Henyey-Greenstein for cloud and aerosol
*  and Rayleigh for molecular scattering

         PMOM(0,II) = 1.0
         DO 15 ISTR = 1, NSTR
            PMCLD = GCLD(i,iw)**(ISTR)
            PMAER = GAER(i,iw)**(ISTR)
            IF(ISTR .EQ. 2) THEN
               PMRAY = 0.1
            ELSE
               PMRAY = 0.
            ENDIF
            PMOM(ISTR,II) = (PMCLD*DSCLD + PMAER*DSAER + 
     $           PMRAY*DTRL(i,iw)) / DTSCT
 15      CONTINUE

 10   CONTINUE

* call rt routine:

      IF( nstr .LT. 2 ) THEN

         CALL ps2str(nz,zen,ag,dt,om,g,
     $        dsdh, nid, delta,
     $        fdiri, fupi, fdni, ediri, eupi, edni)

      ELSE

         CALL  PSNDO( dsdh, nid,
     $        NLYR, DTAUC, SSALB, PMOM, 
     $        ALBEDO, NSTR, 
     $        NUMU, UMU, CWT, UMU0,
     $        MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXPHI, 
     $        RFLDIR,RFLDN, FLUP, U0U,
     $        uavgso, uavgup, uavgdn,
     $        sindir, sinup, sindn)

      ENDIF

* output (invert z-coordinate)

      DO 20 i = 1, nz
         ii = nz - i + 1

         IF( nstr .LT. 2 ) THEN
            fdir(i) = fdiri(ii)
            fup(i) = fupi(ii)
            fdn(i) = fdni(ii)
            edir(i) = ediri(ii)
            eup(i) = eupi(ii)
            edn(i) = edni(ii)
         ELSE
            edir(i) = RFLDIR(II)
            edn(i)  = RFLDN(II)
            eup(i)  = FLUP(II)
            fdir(i) = 4.* pi * uavgso(ii)
            fdn(i)  = 4.* pi * uavgdn(ii)
            fup(i)  = 4.* pi * uavgup(ii)
            sdir(i) = sindir(ii)
            sdn(i)  = sindn(ii)
            sup(i)  = sinup(ii)

*bm  azimutally averaged radiances at computational angles:
*bm  ldif(iu,i) is the radiance at level i and cosine of polar angle UMU(iu);
*bm  the polar angle is measured from the upward direction, implying that 
*bm  positive mu is upwelling and negative mu down-welling radiation.

            DO iu = 1, numu
               ldif(iu,i) = u0u (iu, ii)
            ENDDO

         ENDIF

 20   CONTINUE

*bm  example output:
c         DO iu = 1, numu
c            WRITE (*,*) 'rad',iu,umu(iu),ldif(iu,1)
c         ENDDO

*bm  example for an integral, irradiance

c         irrad = 0.
c         DO iu = 1, NUMU/2
c            irrad = irrad + ldif(iu,1)*cwt(iu)*umu(iu)
c         ENDDO
c         irrad = irrad * 2. * pi

c         WRITE (*,*) edn(1),' = ',irrad,' ?'

      RETURN
      END

*=============================================================================*

      SUBROUTINE ps2str(nlevel,zen,rsfc,tauu,omu,gu,
     $     dsdh, nid, delta,
     $     fdr, fup, fdn, edr, eup, edn)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Solve two-stream equations for multiple layers.  The subroutine is based =*
*=  on equations from:  Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.=*
*=  It contains 9 two-stream methods to choose from.  A pseudo-spherical     =*
*=  correction has also been added.                                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NLEVEL  - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
*=  RSFC    - REAL, surface albedo at current wavelength                  (I)=*
*=  TAUU    - REAL, unscaled optical depth of each layer                  (I)=*
*=  OMU     - REAL, unscaled single scattering albedo of each layer       (I)=*
*=  GU      - REAL, unscaled asymmetry parameter of each layer            (I)=*
*=  DSDH    - REAL, slant path of direct beam through each layer crossed  (I)=*
*=            when travelling from the top of the atmosphere to layer i;     =*
*=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
*=  NID     - INTEGER, number of layers crossed by the direct beam when   (I)=*
*=            travelling from the top of the atmosphere to layer i;          =*
*=            NID(i), i = 0..NZ-1                                            =*
*=  DELTA   - LOGICAL, switch to use delta-scaling                        (I)=*
*=            .TRUE. -> apply delta-scaling                                  =*
*=            .FALSE.-> do not apply delta-scaling                           =*
*=  FDR     - REAL, contribution of the direct component to the total     (O)=*
*=            actinic flux at each altitude level                            =*
*=  FUP     - REAL, contribution of the diffuse upwelling component to    (O)=*
*=            the total actinic flux at each altitude level                  =*
*=  FDN     - REAL, contribution of the diffuse downwelling component to  (O)=*
*=            the total actinic flux at each altitude level                  =*
*=  EDR     - REAL, contribution of the direct component to the total     (O)=*
*=            spectral irradiance at each altitude level                     =*
*=  EUP     - REAL, contribution of the diffuse upwelling component to    (O)=*
*=            the total spectral irradiance at each altitude level           =*
*=  EDN     - REAL, contribution of the diffuse downwelling component to  (O)=*
*=            the total spectral irradiance at each altitude level           =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER nrows
      PARAMETER(nrows=2*kz)

*******
* input:
*******
      INTEGER nlevel
      REAL zen, rsfc
      REAL tauu(kz), omu(kz), gu(kz)
      REAL dsdh(0:kz,kz)
      INTEGER nid(0:kz)
      LOGICAL delta

*******
* output:
*******
      REAL fup(kz),fdn(kz),fdr(kz)
      REAL eup(kz),edn(kz),edr(kz)

*******
* local:
*******
      REAL tausla(0:kz), tauc(0:kz)
      REAL mu2(0:kz), mu, sum

* internal coefficients and matrix
      REAL lam(kz),taun(kz),bgam(kz)
      REAL e1(kz),e2(kz),e3(kz),e4(kz)
      REAL cup(kz),cdn(kz),cuptn(kz),cdntn(kz)
      REAL mu1(kz)
      INTEGER row
      REAL a(nrows),b(nrows),d(nrows),e(nrows),y(nrows)

*******
* other:
*******
      REAL pifs, fdn0
      REAL gi(kz), omi(kz), tempg
      REAL f, g, om
      REAL gam1, gam2, gam3, gam4

* For calculations of Associated Legendre Polynomials for GAMA1,2,3,4
* in delta-function, modified quadrature, hemispheric constant,
* Hybrid modified Eddington-delta function metods, p633,Table1.
* W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
* W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440, 
* uncomment the following two lines and the appropriate statements further
* down.
C     REAL YLM0, YLM2, YLM4, YLM6, YLM8, YLM10, YLM12, YLMS, BETA0,
C    >     BETA1, BETAn, amu1, subd

      REAL expon, expon0, expon1, divisr, temp, up, dn
      REAL ssfc
      INTEGER nlayer, mrows, lev

      INTEGER i, j

* Some additional program constants:

      REAL eps
      PARAMETER (eps = 1.E-3)
*_______________________________________________________________________

* MU = cosine of solar zenith angle
* RSFC = surface albedo
* TAUU =  unscaled optical depth of each layer
* OMU  =  unscaled single scattering albedo
* GU   =  unscaled asymmetry factor
* KLEV = max dimension of number of layers in atmosphere
* NLAYER = number of layers in the atmosphere
* NLEVEL = nlayer + 1 = number of levels

* initial conditions:  pi*solar flux = 1;  diffuse incidence = 0

      pifs = 1.      
      fdn0 = 0.

      nlayer = nlevel - 1

       mu = COS(zen*pi/180.)

************** compute coefficients for each layer:
* GAM1 - GAM4 = 2-stream coefficients, different for different approximations
* EXPON0 = calculation of e when TAU is zero
* EXPON1 = calculation of e when TAU is TAUN
* CUP and CDN = calculation when TAU is zero
* CUPTN and CDNTN = calc. when TAU is TAUN
* DIVISR = prevents division by zero

        do j = 0, kz
           tauc(j) = 0.
           tausla(j) = 0.
           mu2(j) = 1./SQRT(largest)

        end do

       IF( .NOT. delta ) THEN
         DO i = 1, nlayer
           gi(i) = gu(i)
           omi(i) = omu(i)
           taun(i) = tauu(i)
         ENDDO
       ELSE 

* delta-scaling. Have to be done for delta-Eddington approximation, 
* delta discrete ordinate, Practical Improved Flux Method, delta function,
* and Hybrid modified Eddington-delta function methods approximations

         DO i = 1, nlayer
           f = gu(i)*gu(i)
           gi(i) = (gu(i) - f)/(1 - f)
           omi(i) = (1 - f)*omu(i)/(1 - omu(i)*f)       
           taun(i) = (1 - omu(i)*f)*tauu(i)
         ENDDO
        END IF

*
* calculate slant optical depth at the top of the atmosphere when zen>90.
* in this case, higher altitude of the top layer is recommended which can 
* be easily changed in gridz.f.
*
         IF(zen .GT. 90.0) THEN
           IF(nid(0) .LT. 0) THEN
             tausla(0) = largest
           ELSE
             sum = 0.0
             DO j = 1, nid(0)
              sum = sum + 2.*taun(j)*dsdh(0,j)
             END DO
             tausla(0) = sum 
           END IF
         END IF
  
*
        DO 11, i = 1, nlayer

         g = gi(i)
         om = omi(i)
         tauc(i) = tauc(i-1) + taun(i)

* stay away from 1 by precision.  For g, also stay away from -1

         tempg = AMIN1(abs(g),1. - precis)
         g = SIGN(tempg,g)
         om = AMIN1(om,1.-precis)


* calculate slant optical depth
*              
          IF(nid(i) .LT. 0) THEN
            tausla(i) = largest
          ELSE
            sum = 0.0
            DO j = 1, MIN(nid(i),i)
               sum = sum + taun(j)*dsdh(i,j)
            ENDDO
            DO j = MIN(nid(i),i)+1,nid(i)
               sum = sum + 2.*taun(j)*dsdh(i,j)
            ENDDO
            tausla(i) = sum 
            IF(tausla(i) .EQ. tausla(i-1)) THEN
              mu2(i) = SQRT(largest)
            ELSE
              mu2(i) = (tauc(i)-tauc(i-1))/(tausla(i)-tausla(i-1))
              mu2(i) = SIGN( AMAX1(ABS(mu2(i)),1./SQRT(largest)),
     $                     mu2(i) )
            END IF
          END IF
*
*** the following gamma equations are from pg 16,289, Table 1
*** save mu1 for each approx. for use in converting irradiance to actinic flux

* Eddington approximation(Joseph et al., 1976, JAS, 33, 2452):

        gam1 =  (7. - om*(4. + 3.*g))/4.
        gam2 = -(1. - om*(4. - 3.*g))/4.
        gam3 = (2. - 3.*g*mu)/4.
        gam4 = 1. - gam3
        mu1(i) = 0.5

* quadrature (Liou, 1973, JAS, 30, 1303-1326; 1974, JAS, 31, 1473-1475):

c          gam1 = 1.7320508*(2. - om*(1. + g))/2.
c          gam2 = 1.7320508*om*(1. - g)/2.
c          gam3 = (1. - 1.7320508*g*mu)/2.
c          gam4 = 1. - gam3
c          mu1(i) = 1./sqrt(3.)
         
* hemispheric mean (Toon et al., 1089, JGR, 94, 16287):

c          gam1 = 2. - om*(1. + g)
c          gam2 = om*(1. - g)
c          gam3 = (2. - g*mu)/4.
c          gam4 = 1. - gam3
c          mu1(i) = 0.5

* PIFM  (Zdunkovski et al.,1980, Conrib.Atmos.Phys., 53, 147-166):
c         GAM1 = 0.25*(8. - OM*(5. + 3.*G))
c         GAM2 = 0.75*OM*(1.-G)
c         GAM3 = 0.25*(2.-3.*G*MU)
c         GAM4 = 1. - GAM3
c         mu1(i) = 0.5

* delta discrete ordinates  (Schaller, 1979, Contrib.Atmos.Phys, 52, 17-26):
c         GAM1 = 0.5*1.7320508*(2. - OM*(1. + G))
c         GAM2 = 0.5*1.7320508*OM*(1.-G)
c         GAM3 = 0.5*(1.-1.7320508*G*MU)
c         GAM4 = 1. - GAM3
c         mu1(i) = 1./sqrt(3.)

* Calculations of Associated Legendre Polynomials for GAMA1,2,3,4
* in delta-function, modified quadrature, hemispheric constant,
* Hybrid modified Eddington-delta function metods, p633,Table1.
* W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
* W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440
c      YLM0 = 2.
c      YLM2 = -3.*G*MU
c      YLM4 = 0.875*G**3*MU*(5.*MU**2-3.)
c      YLM6=-0.171875*G**5*MU*(15.-70.*MU**2+63.*MU**4)
c     YLM8=+0.073242*G**7*MU*(-35.+315.*MU**2-693.*MU**4
c    *+429.*MU**6)
c     YLM10=-0.008118*G**9*MU*(315.-4620.*MU**2+18018.*MU**4
c    *-25740.*MU**6+12155.*MU**8)
c     YLM12=0.003685*G**11*MU*(-693.+15015.*MU**2-90090.*MU**4
c    *+218790.*MU**6-230945.*MU**8+88179.*MU**10)
c      YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
c      YLMS=0.25*YLMS
c      BETA0 = YLMS
c
c         amu1=1./1.7320508
c      YLM0 = 2.
c      YLM2 = -3.*G*amu1
c      YLM4 = 0.875*G**3*amu1*(5.*amu1**2-3.)
c      YLM6=-0.171875*G**5*amu1*(15.-70.*amu1**2+63.*amu1**4)
c     YLM8=+0.073242*G**7*amu1*(-35.+315.*amu1**2-693.*amu1**4
c    *+429.*amu1**6)
c     YLM10=-0.008118*G**9*amu1*(315.-4620.*amu1**2+18018.*amu1**4
c    *-25740.*amu1**6+12155.*amu1**8)
c     YLM12=0.003685*G**11*amu1*(-693.+15015.*amu1**2-90090.*amu1**4
c    *+218790.*amu1**6-230945.*amu1**8+88179.*amu1**10)
c      YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
c      YLMS=0.25*YLMS
c      BETA1 = YLMS
c
c         BETAn = 0.25*(2. - 1.5*G-0.21875*G**3-0.085938*G**5
c    *-0.045776*G**7)


* Hybrid modified Eddington-delta function(Meador and Weaver,1980,JAS,37,630):
c         subd=4.*(1.-G*G*(1.-MU))
c         GAM1 = (7.-3.*G*G-OM*(4.+3.*G)+OM*G*G*(4.*BETA0+3.*G))/subd
c         GAM2 =-(1.-G*G-OM*(4.-3.*G)-OM*G*G*(4.*BETA0+3.*G-4.))/subd
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = (1. - g*g*(1.- mu) )/(2. - g*g)

*****
* delta function  (Meador, and Weaver, 1980, JAS, 37, 630):
c         GAM1 = (1. - OM*(1. - beta0))/MU
c         GAM2 = OM*BETA0/MU
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = mu
*****
* modified quadrature (Meador, and Weaver, 1980, JAS, 37, 630):
c         GAM1 = 1.7320508*(1. - OM*(1. - beta1))
c         GAM2 = 1.7320508*OM*beta1
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = 1./sqrt(3.)

* hemispheric constant (Toon et al., 1989, JGR, 94, 16287):
c         GAM1 = 2.*(1. - OM*(1. - betan))
c         GAM2 = 2.*OM*BETAn
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = 0.5

*****

* lambda = pg 16,290 equation 21
* big gamma = pg 16,290 equation 22
* if gam2 = 0., then bgam = 0. 

         lam(i) = sqrt(gam1*gam1 - gam2*gam2)

         IF( gam2 .NE. 0.) THEN
            bgam(i) = (gam1 - lam(i))/gam2
         ELSE
            bgam(i) = 0.
         ENDIF

         expon = EXP(-lam(i)*taun(i))

* e1 - e4 = pg 16,292 equation 44
         
         e1(i) = 1. + bgam(i)*expon
         e2(i) = 1. - bgam(i)*expon
         e3(i) = bgam(i) + expon
         e4(i) = bgam(i) - expon

* the following sets up for the C equations 23, and 24
* found on page 16,290
* prevent division by zero (if LAMBDA=1/MU, shift 1/MU^2 by EPS = 1.E-3
* which is approx equiv to shifting MU by 0.5*EPS* (MU)**3

         expon0 = EXP(-tausla(i-1))
         expon1 = EXP(-tausla(i))
          
         divisr = lam(i)*lam(i) - 1./(mu2(i)*mu2(i))
         temp = AMAX1(eps,abs(divisr))
         divisr = SIGN(temp,divisr)

         up = om*pifs*((gam1 - 1./mu2(i))*gam3 + gam4*gam2)/divisr
         dn = om*pifs*((gam1 + 1./mu2(i))*gam4 + gam2*gam3)/divisr
         
* cup and cdn are when tau is equal to zero
* cuptn and cdntn are when tau is equal to taun

         cup(i) = up*expon0
         cdn(i) = dn*expon0
         cuptn(i) = up*expon1
         cdntn(i) = dn*expon1
 
   11 CONTINUE

***************** set up matrix ******
* ssfc = pg 16,292 equation 37  where pi Fs is one (unity).

      ssfc = rsfc*mu*EXP(-tausla(nlayer))*pifs

* MROWS = the number of rows in the matrix

      mrows = 2*nlayer     
      
* the following are from pg 16,292  equations 39 - 43.
* set up first row of matrix:

      i = 1
      a(1) = 0.
      b(1) = e1(i)
      d(1) = -e2(i)
      e(1) = fdn0 - cdn(i)

      row=1

* set up odd rows 3 thru (MROWS - 1):

      i = 0
      DO 20, row = 3, mrows - 1, 2
         i = i + 1
         a(row) = e2(i)*e3(i) - e4(i)*e1(i)
         b(row) = e1(i)*e1(i + 1) - e3(i)*e3(i + 1)
         d(row) = e3(i)*e4(i + 1) - e1(i)*e2(i + 1)
         e(row) = e3(i)*(cup(i + 1) - cuptn(i)) + 
     $        e1(i)*(cdntn(i) - cdn(i + 1))
   20 CONTINUE

* set up even rows 2 thru (MROWS - 2): 

      i = 0
      DO 30, row = 2, mrows - 2, 2
         i = i + 1
         a(row) = e2(i + 1)*e1(i) - e3(i)*e4(i + 1)
         b(row) = e2(i)*e2(i + 1) - e4(i)*e4(i + 1)
         d(row) = e1(i + 1)*e4(i + 1) - e2(i + 1)*e3(i + 1)
         e(row) = (cup(i + 1) - cuptn(i))*e2(i + 1) - 
     $        (cdn(i + 1) - cdntn(i))*e4(i + 1)
   30 CONTINUE

* set up last row of matrix at MROWS:

      row = mrows
      i = nlayer
      
      a(row) = e1(i) - rsfc*e3(i)
      b(row) = e2(i) - rsfc*e4(i)
      d(row) = 0.
      e(row) = ssfc - cuptn(i) + rsfc*cdntn(i)

* solve tri-diagonal matrix:

      CALL tridag(a, b, d, e, y, mrows)

**** unfold solution of matrix, compute output fluxes:

      row = 1 
      lev = 1
      j = 1
      
* the following equations are from pg 16,291  equations 31 & 32

      fdr(lev) = EXP( -tausla(0) )
      edr(lev) = mu * fdr(lev)
      edn(lev) = fdn0
      eup(lev) =  y(row)*e3(j) - y(row + 1)*e4(j) + cup(j)
      fdn(lev) = edn(lev)/mu1(lev)
      fup(lev) = eup(lev)/mu1(lev)

      DO 60, lev = 2, nlayer + 1
         fdr(lev) = EXP(-tausla(lev-1))
         edr(lev) =  mu *fdr(lev)
         edn(lev) =  y(row)*e3(j) + y(row + 1)*e4(j) + cdntn(j)
         eup(lev) =  y(row)*e1(j) + y(row + 1)*e2(j) + cuptn(j)
         fdn(lev) = edn(lev)/mu1(j)
         fup(lev) = eup(lev)/mu1(j)

         row = row + 2
         j = j + 1
   60 CONTINUE
*_______________________________________________________________________

      RETURN
      END

*=============================================================================*

      SUBROUTINE tridag(a,b,c,r,u,n)

*_______________________________________________________________________
* solves tridiagonal system.  From Numerical Recipies, p. 40
*_______________________________________________________________________

      IMPLICIT NONE

* input:
      INTEGER n
      REAL a, b, c, r
      DIMENSION a(n),b(n),c(n),r(n)

* output:
      REAL u
      DIMENSION u(n)

* local:
      INTEGER j

      INCLUDE 'params'
      REAL bet, gam
      DIMENSION gam(2*kz)
*_______________________________________________________________________

      IF (b(1) .EQ. 0.) STOP 1001
      bet   = b(1)
      u(1) = r(1)/bet
      DO 11, j = 2, n   
         gam(j) = c(j - 1)/bet
         bet = b(j) - a(j)*gam(j)
         IF (bet .EQ. 0.) STOP 2002 
         u(j) = (r(j) - a(j)*u(j - 1))/bet
   11 CONTINUE
      DO 12, j = n - 1, 1, -1  
         u(j) = u(j) - gam(j + 1)*u(j + 1)
   12 CONTINUE
*_______________________________________________________________________

      RETURN
      END

*=============================================================================*

c  Note: CDIR$ and CFPP$ comment lines are relevant only when running
c        on Cray computers.  They cause better optimization of loops
c        immediately following.


C      SUBROUTINE DISORT( dsdh, nid,
      SUBROUTINE PSNDO( dsdh, nid,
     &                   NLYR, DTAUC, SSALB, PMOM,
     &                   ALBEDO, NSTR, 
     &                   NUMU, UMU, CWT, UMU0,
     &                   MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXPHI, 
     $                   RFLDIR, RFLDN, FLUP, U0U,
     &                   uavgso, uavgup, uavgdn, 
     &                   sindir, sinup, sindn )

c Improved handling of numerical instabilities. Bernhard Mayer on 5/3/99.
c  disort seems to produce unstable results for certain combinations
c  of single scattering albedo and phase function. A temporary fix has been 
c  introduced to avoid this problem: The original instability check in 
c  UPBEAM fails on certain compiler/machine combinations (e.g., gcc/LINUX, 
c  or xlf/IBM RS6000). This check has therefore been replaced by a new one.
c  Whenever UPBEAM reports an instability, the single scattering albedo 
c  of the respective layer is changed by a small amount, and the 
c  calculation is repeated until numerically stable conditions are reached 
c  (all the necessary changes are confined to the new subroutine SOLVEC 
c  and the slighly changed subroutine UPBEAM). To check for potential 
c  instabilities, the variable 'RCOND' returned by SGECO is compared to
c  a machine-dependent constant, 'MINRCOND'. The value of this constant 
c  determines (a) if really all instabilities are caught; and (b) the 
c  amount by which the single scattering albedo has to be changed. The 
c  value of 'MINRCOND' is therefore a compromise between numerical 
c  stability on the one hand and uncertainties introduced by changing 
c  the atmospheric conditions and increased computational time on the 
c  other hand (an increase of MINRCOND will lead to the detection of
c  more potential numerical instabilities, and thus to an increase in 
c  computational time; by changing the atmospheric conditions, that is,
c  the single scattering albedo, the result might however be changed 
c  unfavourably, if the change is too large). From a limited number 
c  of experiments we found that 'MINRCOND = 5000. * R1MACH(4)' seems 
c  to be a good choice if high accuracy is required (more tests are 
c  definitely neccessary!). If an instability is encountered, a message 
c  is printed telling about neccessary changes to the single scattering 
c  albedo. This message may be switched off by setting 'DEBUG = .FALSE.' 
c  in subroutine SOLVEC. 
c
c
c modified to calculate sine-weighted intensities. Bernhard Mayer on 2/12/99.
c modified to handle some numerical instabilities. Chris Fischer on 1/22/99.
c modified by adding pseudo-spherical correction. Jun Zeng on 3/11/97.
c dsdh: slant path of direct beam through each layer crossed  
c       when travelling from the top of the atmosphere to layer i;    
c       dsdh(i,j), i = 0..nlyr, j = 1..nlyr;
c nid:  number of layers crossed by the direct beam when   
c       travelling from the top of the atmosphere to layer i; 
c       NID(i), i = 0..nlyr.
c uavgso, uvagup, and uvagdn are direct, downward diffuse, and upward
c diffuse actinic flux (mean intensity).
c u0u is the azimuthally averaged intensity, check DISORT.doc for details.
c *******************************************************************
c       Plane-parallel discrete ordinates radiative transfer program
c                      V E R S I O N    1.1
c             ( see DISORT.DOC for complete documentation )
c *******************************************************************


c +------------------------------------------------------------------+
c  Calling Tree (omitting calls to ERRMSG):
c  (routines in parentheses are not in this file)

c  DISORT-+-(R1MACH)
c         +-ZEROIT
c         +-CHEKIN-+-(WRTBAD)
c         |        +-(WRTDIM)
c         |        +-DREF
c         +-ZEROAL
c         +-SETDIS-+-QGAUSN (1)-+-(D1MACH)
c         +-PRTINP
c         +-LEPOLY see 2
c         +-SURFAC-+-QGAUSN see 1
c         |        +-LEPOLY see 2
c         |        +-ZEROIT
c         +-SOLEIG see 3
c         +-UPBEAM-+-(SGECO)
c         |        +-(SGESL)
c         +-TERPEV
c         +-TERPSO
c         +-SETMTX see 4
c         +-SOLVE0-+-ZEROIT
c         |        +-(SGBCO)
c         |        +-(SGBSL)
c         +-FLUXES--ZEROIT
c         +-PRAVIN
c         +-RATIO--(R1MACH)
c         +-PRTINT

c *** Intrinsic Functions used in DISORT package which take
c     non-negligible amount of time:

c    EXP :  Called by- ALBTRN, ALTRIN, CMPINT, FLUXES, SETDIS,
c                      SETMTX, SPALTR, USRINT, PLKAVG

c    SQRT : Called by- ASYMTX, LEPOLY, SOLEIG

c +-------------------------------------------------------------------+

c  Index conventions (for all DO-loops and all variable descriptions):

c     IU     :  for user polar angles

c  IQ,JQ,KQ  :  for computational polar angles ('quadrature angles')

c   IQ/2     :  for half the computational polar angles (just the ones
c               in either 0-90 degrees, or 90-180 degrees)

c     J      :  for user azimuthal angles

c     K,L    :  for Legendre expansion coefficients or, alternatively,
c               subscripts of associated Legendre polynomials

c     LU     :  for user levels

c     LC     :  for computational layers (each having a different
c               single-scatter albedo and/or phase function)

c    LEV     :  for computational levels

c    MAZIM   :  for azimuthal components in Fourier cosine expansion
c               of intensity and phase function

c +------------------------------------------------------------------+

c               I N T E R N A L    V A R I A B L E S

c   AMB(IQ/2,IQ/2)    First matrix factor in reduced eigenvalue problem
c                     of Eqs. SS(12), STWJ(8E)  (used only in SOLEIG)

c   APB(IQ/2,IQ/2)    Second matrix factor in reduced eigenvalue problem
c                     of Eqs. SS(12), STWJ(8E)  (used only in SOLEIG)

c   ARRAY(IQ,IQ)      Scratch matrix for SOLEIG, UPBEAM and UPISOT
c                     (see each subroutine for definition)

c   B()               Right-hand side vector of Eq. SC(5) going into
c                     SOLVE0,1;  returns as solution vector
c                     vector  L, the constants of integration

c   BDR(IQ/2,0:IQ/2)  Bottom-boundary bidirectional reflectivity for a
c                     given azimuthal component.  First index always
c                     refers to a computational angle.  Second index:
c                     if zero, refers to incident beam angle UMU0;
c                     if non-zero, refers to a computational angle.

c   BEM(IQ/2)         Bottom-boundary directional emissivity at compu-
c                     tational angles.

c   BPLANK            Intensity emitted from bottom boundary

c   CBAND()           Matrix of left-hand side of the linear system
c                     Eq. SC(5), scaled by Eq. SC(12);  in banded
c                     form required by LINPACK solution routines

c   CC(IQ,IQ)         C-sub-IJ in Eq. SS(5)

c   CMU(IQ)           Computational polar angles (Gaussian)

c   CWT(IQ)           Quadrature weights corresponding to CMU

c   DELM0             Kronecker delta, delta-sub-M0, where M = MAZIM
c                     is the number of the Fourier component in the
c                     azimuth cosine expansion

c   DITHER            Small quantity subtracted from single-scattering
c                     albedos of unity, in order to avoid using special
c                     case formulas;  prevents an eigenvalue of exactly
c                     zero from occurring, which would cause an
c                     immediate overflow

c   DTAUCP(LC)        Computational-layer optical depths (delta-M-scaled
c                     if DELTAM = TRUE, otherwise equal to DTAUC)

c   EMU(IU)           Bottom-boundary directional emissivity at user
c                     angles.

c   EVAL(IQ)          Temporary storage for eigenvalues of Eq. SS(12)

c   EVECC(IQ,IQ)      Complete eigenvectors of SS(7) on return from
c                     SOLEIG; stored permanently in  GC

c   EXPBEA(LC)        Transmission of direct beam in delta-M optical
c                     depth coordinates

c   FLYR(LC)          Truncated fraction in delta-M method

c   GL(K,LC)          Phase function Legendre polynomial expansion
c                     coefficients, calculated from PMOM by
c                     including single-scattering albedo, factor
c                     2K+1, and (if DELTAM=TRUE) the delta-M
c                     scaling

c   GC(IQ,IQ,LC)      Eigenvectors at polar quadrature angles,
c                     g  in Eq. SC(1)

c   GU(IU,IQ,LC)      Eigenvectors interpolated to user polar angles
c                     ( g  in Eqs. SC(3) and S1(8-9), i.e.
c                       G without the L factor )

c   HLPR()            Legendre coefficients of bottom bidirectional
c                     reflectivity (after inclusion of 2K+1 factor)

c   IPVT(LC*IQ)       Integer vector of pivot indices for LINPACK
c                     routines

c   KK(IQ,LC)         Eigenvalues of coeff. matrix in Eq. SS(7)

c   KCONV             Counter in azimuth convergence test

c   LAYRU(LU)         Computational layer in which user output level
c                     UTAU(LU) is located

c   LL(IQ,LC)         Constants of integration L in Eq. SC(1),
c                     obtained by solving scaled version of Eq. SC(5)

c   LYRCUT            TRUE, radiation is assumed zero below layer
c                     NCUT because of almost complete absorption

c   NAZ               Number of azimuthal components considered

c   NCUT              Computational layer number in which absorption
c                     optical depth first exceeds ABSCUT

c   OPRIM(LC)         Single scattering albedo after delta-M scaling

c   PASS1             TRUE on first entry, FALSE thereafter

c   PKAG(0:LC)        Integrated Planck function for internal emission

c   PSI(IQ)           Sum just after square bracket in  Eq. SD(9)

c   RMU(IU,0:IQ)      Bottom-boundary bidirectional reflectivity for a
c                     given azimuthal component.  First index always
c                     refers to a user angle.  Second index:
c                     if zero, refers to incident beam angle UMU0;
c                     if non-zero, refers to a computational angle.

c   TAUC(0:LC)        Cumulative optical depth (un-delta-M-scaled)

c   TAUCPR(0:LC)      Cumulative optical depth (delta-M-scaled if
c                     DELTAM = TRUE, otherwise equal to TAUC)

c   TPLANK            Intensity emitted from top boundary

c   UUM(IU,LU)        Expansion coefficients when the intensity
c                     (u-super-M) is expanded in Fourier cosine series
c                     in azimuth angle

c   U0C(IQ,LU)        Azimuthally-averaged intensity

c   UTAUPR(LU)        Optical depths of user output levels in delta-M
c                     coordinates;  equal to  UTAU(LU) if no delta-M

c   WK()              scratch array

c   XR0(LC)           X-sub-zero in expansion of thermal source func-
c                     tion preceding Eq. SS(14) (has no mu-dependence)

c   XR1(LC)           X-sub-one in expansion of thermal source func-
c                     tion;  see  Eqs. SS(14-16)

c   YLM0(L)           Normalized associated Legendre polynomial
c                     of subscript L at the beam angle (not saved
c                     as function of superscipt M)

c   YLMC(L,IQ)        Normalized associated Legendre polynomial
c                     of subscript L at the computational angles
c                     (not saved as function of superscipt M)

c   YLMU(L,IU)        Normalized associated Legendre polynomial
c                     of subscript L at the user angles
c                     (not saved as function of superscipt M)

c   Z()               scratch array used in  SOLVE0,1  to solve a
c                     linear system for the constants of integration

c   Z0(IQ)            Solution vectors Z-sub-zero of Eq. SS(16)

c   Z0U(IU,LC)        Z-sub-zero in Eq. SS(16) interpolated to user
c                     angles from an equation derived from SS(16)

c   Z1(IQ)            Solution vectors Z-sub-one  of Eq. SS(16)

c   Z1U(IU,LC)        Z-sub-one in Eq. SS(16) interpolated to user
c                     angles from an equation derived from SS(16)

c   ZBEAM(IU,LC)      Particular solution for beam source

c   ZJ(IQ)            Right-hand side vector  X-sub-zero in
c                     Eq. SS(19), also the solution vector
c                     Z-sub-zero after solving that system

c   ZZ(IQ,LC)         Permanent storage for the beam source vectors ZJ

c   ZPLK0(IQ,LC)      Permanent storage for the thermal source
c                     vectors  Z0  obtained by solving  Eq. SS(16)

c   ZPLK1(IQ,LC)      Permanent storage for the thermal source
c                     vectors  Z1  obtained by solving  Eq. SS(16)

c +-------------------------------------------------------------------+

c  LOCAL SYMBOLIC DIMENSIONS (have big effect on storage requirements):

c       MXCLY  = Max no. of computational layers
c       MXULV  = Max no. of output levels
c       MXCMU  = Max no. of computation polar angles
c       MXUMU  = Max no. of output polar angles
c       MXPHI  = Max no. of output azimuthal angles

c +-------------------------------------------------------------------+

      INCLUDE 'params'

c     .. Parameters ..

      INTEGER   MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI
      PARAMETER ( MXCLY = 151, MXULV = 151, MXCMU = 32, MXUMU = 32,
     &          MXPHI = 3, MI = MXCMU / 2, MI9M2 = 9*MI - 2,
     &          NNLYRI = MXCMU*MXCLY )
c     ..
c     .. Scalar Arguments ..

      CHARACTER HEADER*127
      LOGICAL   DELTAM, LAMBER, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND, MAXCLY, MAXCMU, MAXPHI, MAXULV, MAXUMU, NLYR,
     &          NPHI, NSTR, NTAU, NUMU
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0, WVNMHI, WVNMLO

c     sherical geometry
      REAL dsdh(0:kz,kz)
      INTEGER nid(0:kz)
      REAL tausla(0:kz), tauslau(0:kz), mu2(0:kz)
c     ..
c     .. Array Arguments ..

      LOGICAL   PRNT( 7 )
      REAL      ALBMED( MAXUMU ), DFDT( MAXULV ), DTAUC( MAXCLY ),
     &          FLUP( MAXULV ), HL( 0:MAXCMU ), PHI( MAXPHI ),
     &          PMOM( 0:MAXCMU, MAXCLY ), RFLDIR( MAXULV ),
     &          RFLDN( MAXULV ), SSALB( MAXCLY ), TEMPER( 0:MAXCLY ),
     &          TRNMED( MAXUMU ), U0U( MAXUMU, MAXULV ), UAVG( MAXULV ),
     &          UMU( MAXUMU ), CWT( MAXCMU ), UTAU( MAXULV ),
     &          UU( MAXUMU, MAXULV, MAXPHI ), 
     &          uavgso( maxulv ), uavgup( maxulv ), uavgdn( maxulv ),
     &          sindir( maxulv ), sinup( maxulv ),  sindn ( maxulv )
c     ..
c     .. Local Scalars ..

      LOGICAL   COMPAR, LYRCUT, PASS1
      INTEGER   IQ, IU, J, KCONV, L, LC, LEV, LU, MAZIM, NAZ, NCOL,
     &          NCOS, NCUT, NN
      REAL      ANGCOS, AZERR, AZTERM, BPLANK, COSPHI, DELM0, DITHER,
     &          DUM, RPD, SGN, TPLANK
c     ..
c     .. Local Arrays ..

      INTEGER   IPVT( NNLYRI ), LAYRU( MXULV )

      REAL      AMB( MI, MI ), APB( MI, MI ), ARRAY( MXCMU, MXCMU ),
     &          B( NNLYRI ), BDR( MI, 0:MI ), BEM( MI ),
     &          CBAND( MI9M2, NNLYRI ), CC( MXCMU, MXCMU ),
     &          CMU( MXCMU ), DTAUCP( MXCLY ),
     &          EMU( MXUMU ), EVAL( MI ), EVECC( MXCMU, MXCMU ),
     &          EXPBEA( 0:MXCLY ), FLDIR( MXULV ), FLDN( MXULV ),
     &          FLYR( MXCLY ), GC( MXCMU, MXCMU, MXCLY ),
     &          GL( 0:MXCMU, MXCLY ), GU( MXUMU, MXCMU, MXCLY ),
     &          HLPR( 0:MXCMU ), KK( MXCMU, MXCLY ), LL( MXCMU, MXCLY ),
     &          OPRIM( MXCLY ), PHIRAD( MXPHI ), PKAG( 0:MXCLY ),
     &          PSI( MXCMU ), RMU( MXUMU, 0:MI ), TAUC( 0:MXCLY ),
     &          TAUCPR( 0:MXCLY ), U0C( MXCMU, MXULV ), UTAUPR( MXULV ),
     &          UUM( MXUMU, MXULV ), WK( MXCMU ), XR0( MXCLY ),
     &          XR1( MXCLY ), YLM0( 0:MXCMU ), YLMC( 0:MXCMU, MXCMU ),
     &          YLMU( 0:MXCMU, MXUMU ), Z( NNLYRI ), Z0( MXCMU ),
     &          Z0U( MXUMU, MXCLY ), Z1( MXCMU ), Z1U( MXUMU, MXCLY ),
     &          ZBEAM( MXUMU, MXCLY ), ZJ( MXCMU ),
     &          ZPLK0( MXCMU, MXCLY ), ZPLK1( MXCMU, MXCLY ),
     &          ZZ( MXCMU, MXCLY )

cgy added glsave and dgl to allow adjustable dimensioning in SOLVEC
      REAL GLSAVE( 0:MXCMU ), DGL( 0:MXCMU )

      DOUBLE PRECISION AAD( MI, MI ), EVALD( MI ), EVECCD( MI, MI ),
     &                 WKD( MXCMU )
c     ..
c     .. External Functions ..

      REAL      PLKAVG, R1MACH, RATIO
      EXTERNAL  PLKAVG, R1MACH, RATIO
c     ..
c     .. External Subroutines ..

      EXTERNAL  CHEKIN, FLUXES, LEPOLY, PRAVIN, PRTINP,
     &          PRTINT, SETDIS, SETMTX, SOLEIG, SOLVE0, SURFAC,
     &          UPBEAM, ZEROAL, ZEROIT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, ASIN, COS, LEN, MAX
c     ..
      SAVE      PASS1, DITHER, RPD
      DATA      PASS1 / .TRUE. /

* Discrete ordinate constants:
* For pseudo-spherical DISORT, PLANK, USRTAU and USRANG must be .FALSE.;
* ONLYFL must be .TRUE.; FBEAM = 1.; FISOT = 0.; IBCND = 0

      data LAMBER /.TRUE./
      data USRTAU /.FALSE./
      data PLANK /.FALSE./
      data USRANG /.FALSE./
      data ONLYFL /.TRUE./
      data PRNT /.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,
     $		 .FALSE.,.FALSE./
      data ACCUR /0.0001/
      data HEADER /' '/
      data NPHI /0/
      data IBCND /0/
      data FBEAM /1./
      data FISOT /0.0/
      data PHI0 /0.0/

* delat-M scaling option

      data DELTAM /.true./





      IF( PASS1 ) THEN

         DITHER = 10.*R1MACH( 4 )

c                            ** Must dither more on Cray (14-digit prec)

         IF( DITHER.LT.1.E-10 ) DITHER = 10.*DITHER

         RPD  = PI / 180.0
         PASS1 = .FALSE.
      END IF
 
   10 CONTINUE

c                                  ** Calculate cumulative optical depth
c                                     and dither single-scatter albedo
c                                     to improve numerical behavior of
c                                     eigenvalue/vector computation
      CALL ZEROIT( TAUC, MXCLY + 1 )

      DO 20 LC = 1, NLYR

         IF( SSALB( LC ).EQ.1.0 ) SSALB( LC ) = 1.0 - DITHER
         TAUC( LC ) = TAUC( LC - 1 ) + DTAUC( LC )

   20 CONTINUE
c                                ** Check input dimensions and variables

      CALL CHEKIN( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,
     &             USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI,
     &             PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO,
     &             HL, BTEMP, TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, TAUC,
     &             MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXPHI, MXCLY, MXULV,
     &             MXUMU, MXCMU, MXPHI )

c                                 ** Zero internal and output arrays

      CALL  ZEROAL( MXCLY, EXPBEA(1), FLYR, OPRIM, TAUCPR(1), XR0, XR1,
     $              MXCMU, CMU, CWT, PSI, WK, Z0, Z1, ZJ,
     $              MXCMU+1, HLPR, YLM0,
     $              MXCMU**2, ARRAY, CC, EVECC,
     $              (MXCMU+1)*MXCLY, GL,
     $              (MXCMU+1)*MXCMU, YLMC,
     $              (MXCMU+1)*MXUMU, YLMU,
     $              MXCMU*MXCLY, KK, LL, ZZ, ZPLK0, ZPLK1,
     $              MXCMU**2*MXCLY, GC,
     $              MXULV, LAYRU, UTAUPR,
     $              MXUMU*MXCMU*MXCLY, GU,
     $              MXUMU*MXCLY, Z0U, Z1U, ZBEAM,
     $              MI, EVAL,
     $              MI**2, AMB, APB,
     $              NNLYRI, IPVT, Z,
     $              MAXULV, RFLDIR, RFLDN, FLUP, UAVG, DFDT,
     $              MAXUMU, ALBMED, TRNMED,
     $              MAXUMU*MAXULV, U0U,
     $              MAXUMU*MAXULV*MAXPHI, UU )

c                                 ** Perform various setup operations

      CALL SETDIS( dsdh, nid, tausla, tauslau, mu2,
     &             CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM, FLYR,
     &             GL, HL, HLPR, IBCND, LAMBER, LAYRU, LYRCUT, MAXUMU,
     &             MAXCMU, MXCMU, NCUT, NLYR, NTAU, NN, NSTR, PLANK,
     &             NUMU, ONLYFL, OPRIM, PMOM, SSALB, TAUC, TAUCPR, UTAU,
     &             UTAUPR, UMU, UMU0, USRTAU, USRANG )

c                                 ** Print input information
      IF ( PRNT(1) )
     $     CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, PMOM, TEMPER,
     $                  WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,
     $                  NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT,
     $                  LAMBER, ALBEDO, HL, BTEMP, TTEMP, TEMIS,
     $                  DELTAM, PLANK, ONLYFL, ACCUR, FLYR, LYRCUT,
     $                  OPRIM, TAUC, TAUCPR, MAXCMU, PRNT(7) )

c                              ** Handle special case for getting albedo
c                                 and transmissivity of medium for many
c                                 beam angles at once
c                                   ** Calculate Planck functions

         BPLANK = 0.0
         TPLANK = 0.0
         CALL ZEROIT( PKAG, MXCLY + 1 )

c ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
c           (EQ STWJ 5)

      KCONV  = 0
      NAZ  = NSTR - 1
c                                    ** Azimuth-independent case

      IF( FBEAM.EQ.0.0 .OR. ( 1.- UMU0 ).LT.1.E-5 .OR. ONLYFL .OR.
     &      ( NUMU.EQ.1 .AND. ( 1.- UMU(1) ).LT.1.E-5 ) )
     &   NAZ = 0

      DO 160 MAZIM = 0, NAZ

         IF( MAZIM.EQ.0 ) DELM0  = 1.0
         IF( MAZIM.GT.0 ) DELM0  = 0.0

c                             ** Get normalized associated Legendre
c                                polynomials for
c                                (a) incident beam angle cosine
c                                (b) computational and user polar angle
c                                    cosines
         IF( FBEAM.GT.0.0 ) THEN

            NCOS   = 1
            ANGCOS = -UMU0

            CALL LEPOLY( NCOS, MAZIM, MXCMU, NSTR - 1, ANGCOS, YLM0 )

         END IF


         IF( .NOT.ONLYFL .AND. USRANG )
     &       CALL LEPOLY( NUMU, MAZIM, MXCMU, NSTR-1, UMU, YLMU )

         CALL LEPOLY( NN, MAZIM, MXCMU, NSTR-1, CMU, YLMC )

c                       ** Get normalized associated Legendre polys.
c                          with negative arguments from those with
c                          positive arguments; Dave/Armstrong Eq. (15)
         SGN  = - 1.0

         DO 50 L = MAZIM, NSTR - 1

            SGN  = - SGN

            DO 40 IQ = NN + 1, NSTR
               YLMC( L, IQ ) = SGN*YLMC( L, IQ - NN )
   40       CONTINUE

   50    CONTINUE
c                                 ** Specify users bottom reflectivity
c                                    and emissivity properties
      IF ( .NOT.LYRCUT )
     $   CALL  SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER,
     $                 MI, MAZIM, MXCMU, MXUMU, NN, NUMU, NSTR, ONLYFL,
     $                 UMU, USRANG, YLM0, YLMC, YLMU, BDR, EMU, BEM,
     $                 RMU )


c ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============

         DO 60 LC = 1, NCUT

            CALL SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL( 0,LC ), MI,
     &           MAZIM, MXCMU, NN, NSTR, YLM0, YLMC, CC, 
     &           EVECC, EVAL, KK( 1,LC ), GC( 1,1,LC ), AAD, EVECCD, 
     &           EVALD, WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0,
     &           ZJ, ZZ(1,LC), OPRIM(LC), LC, DITHER, mu2(lc),
     &           glsave, dgl)
cgy added glsave and dgl to call to allow adjustable dimensioning

 60      CONTINUE


c ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============


c                      ** Set coefficient matrix of equations combining
c                         boundary and layer interface conditions

         CALL SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,
     &                LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT,
     &                NNLYRI, NN, NSTR, TAUCPR, WK )

c                      ** Solve for constants of integration in homo-
c                         geneous solution (general boundary conditions)

         CALL SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,
     &                FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM, MI,
     &                MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, NNLYRI, PI,
     &                TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )

c                                  ** Compute upward and downward fluxes

      IF ( MAZIM.EQ.0 )
     $     CALL FLUXES( tausla, tauslau,
     $                  CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,
     $                  MAXULV, MXCMU, MXULV, NCUT, NN, NSTR, NTAU,
     $                  PI, PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,
     $                  XR0, XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP,
     $                  FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C,
     $                  uavgso, uavgup, uavgdn,
     $                  sindir, sinup, sindn)

         IF( ONLYFL ) THEN

            IF( MAXUMU.GE.NSTR ) THEN
c                                     ** Save azimuthal-avg intensities
c                                        at quadrature angles
               DO 80 LU = 1, NTAU

                  DO 70 IQ = 1, NSTR
                     U0U( IQ, LU ) = U0C( IQ, LU )
   70             CONTINUE

   80          CONTINUE

            END IF

            GO TO  170

         END IF


         CALL ZEROIT( UUM, MXUMU*MXULV )

         IF( MAZIM.EQ.0 ) THEN
c                               ** Save azimuthally averaged intensities

            DO 110 LU = 1, NTAU

               DO 100 IU = 1, NUMU
                  U0U( IU, LU ) = UUM( IU, LU )

                  DO 90 J = 1, NPHI
                     UU( IU, LU, J ) = UUM( IU, LU )
   90             CONTINUE

  100          CONTINUE

  110       CONTINUE
c                              ** Print azimuthally averaged intensities
c                                 at user angles

            IF( PRNT( 4 ) ) CALL PRAVIN( UMU, NUMU, MAXUMU, UTAU, NTAU,
     &                                   U0U )
            IF( NAZ.GT.0 ) THEN

               CALL ZEROIT( PHIRAD, MXPHI )
               DO 120 J = 1, NPHI
                  PHIRAD( J ) = RPD*( PHI( J ) - PHI0 )
  120          CONTINUE

            END IF


         ELSE
c                                ** Increment intensity by current
c                                   azimuthal component (Fourier
c                                   cosine series);  Eq SD(2)
            AZERR  = 0.0

            DO 150 J = 1, NPHI

               COSPHI = COS( MAZIM*PHIRAD( J ) )

               DO 140 LU = 1, NTAU

                  DO 130 IU = 1, NUMU
                     AZTERM = UUM( IU, LU )*COSPHI
                     UU( IU, LU, J ) = UU( IU, LU, J ) + AZTERM
                     AZERR = MAX( AZERR,
     &                       RATIO( ABS(AZTERM), ABS(UU(IU,LU,J)) ) )
  130             CONTINUE

  140          CONTINUE

  150       CONTINUE

            IF( AZERR.LE.ACCUR ) KCONV  = KCONV + 1

            IF( KCONV.GE.2 ) GO TO  170

         END IF

  160 CONTINUE

c ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============


c                                          ** Print intensities
  170 CONTINUE
      IF( PRNT( 5 ) .AND. .NOT.ONLYFL ) CALL PRTINT( UU, UTAU, NTAU,
     &    UMU, NUMU, PHI, NPHI, MAXULV, MAXUMU )

      END

      SUBROUTINE ASYMTX( AA, EVEC, EVAL, M, IA, IEVEC, IER, WKD, AAD,
     &                   EVECD, EVALD )

c    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

c       Solves eigenfunction problem for real asymmetric matrix
c       for which it is known a priori that the eigenvalues are real.

c       This is an adaptation of a subroutine EIGRF in the IMSL
c       library to use real instead of complex arithmetic, accounting
c       for the known fact that the eigenvalues and eigenvectors in
c       the discrete ordinate solution are real.  Other changes include
c       putting all the called subroutines in-line, deleting the
c       performance index calculation, updating many DO-loops
c       to Fortran77, and in calculating the machine precision
c       TOL instead of specifying it in a data statement.

c       EIGRF is based primarily on EISPACK routines.  The matrix is
c       first balanced using the Parlett-Reinsch algorithm.  Then
c       the Martin-Wilkinson algorithm is applied.

c       References:
c          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
c             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
c             Sources and Development of Mathematical Software,
c             Prentice-Hall, Englewood Cliffs, NJ
c         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
c             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
c         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
c             Clarendon Press, Oxford

c   I N P U T    V A R I A B L E S:

c       AA    :  input asymmetric matrix, destroyed after solved
c        M    :  order of  AA
c       IA    :  first dimension of  AA
c    IEVEC    :  first dimension of  EVEC

c   O U T P U T    V A R I A B L E S:

c       EVEC  :  (unnormalized) eigenvectors of  AA
c                   ( column J corresponds to EVAL(J) )

c       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M )

c       IER   :  if .NE. 0, signals that EVAL(IER) failed to converge;
c                   in that case eigenvalues IER+1,IER+2,...,M  are
c                   correct but eigenvalues 1,...,IER are set to zero.

c   S C R A T C H   V A R I A B L E S:

c       WKD   :  work area ( dimension at least 2*M )
c       AAD   :  double precision stand-in for AA
c       EVECD :  double precision stand-in for EVEC
c       EVALD :  double precision stand-in for EVAL

c   Called by- SOLEIG
c   Calls- D1MACH, ERRMSG
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..

      INTEGER   IA, IER, IEVEC, M
c     ..
c     .. Array Arguments ..

      REAL      AA( IA, M ), EVAL( M ), EVEC( IEVEC, M )
      DOUBLE PRECISION AAD( IA, M ), EVALD( M ), EVECD( IA, M ),
     &                 WKD( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   NOCONV, NOTLAS
      INTEGER   I, II, IN, J, K, KA, KKK, L, LB, LLL, N, N1, N2
      DOUBLE PRECISION C1, C2, C3, C4, C5, C6, COL, DISCRI, F, G, H,
     &                 ONE, P, Q, R, REPL, RNORM, ROW, S, SCALE, SGN, T,
     &                 TOL, UU, VV, W, X, Y, Z, ZERO
c     ..
c     .. External Functions ..

      DOUBLE PRECISION D1MACH
      EXTERNAL  D1MACH
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, DBLE, MIN, SIGN, SQRT
c     ..
      DATA      C1 / 0.4375D0 / , C2 / 0.5D0 / , C3 / 0.75D0 / ,
     &          C4 / 0.95D0 / , C5 / 16.D0 / , C6 / 256.D0 / ,
     &          ZERO / 0.D0 / , ONE / 1.D0 /


      IER  = 0
      TOL  = D1MACH( 4 )

      IF( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M )
     &    CALL ERRMSG( 'ASYMTX--bad input variable(s)', .TRUE. )

c                           ** Handle 1x1 and 2x2 special cases

      IF( M.EQ.1 ) THEN

         EVAL( 1 )    = AA( 1, 1 )
         EVEC( 1, 1 ) = 1.0
         RETURN

      ELSE IF( M.EQ.2 ) THEN

         DISCRI = ( AA( 1,1 ) - AA( 2,2 ) )**2 +
     &              4.*AA( 1, 2 )*AA( 2, 1 )

         IF( DISCRI.LT.0.0 )
     &       CALL ERRMSG( 'ASYMTX--complex evals in 2x2 case',.TRUE. )

         SGN  = 1.0

         IF( AA( 1,1 ).LT.AA( 2,2 ) ) SGN  = - 1.0

         EVAL( 1 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) + SGN*SQRT( DISCRI ) )
         EVAL( 2 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) - SGN*SQRT( DISCRI ) )
         EVEC( 1, 1 ) = 1.0
         EVEC( 2, 2 ) = 1.0

         IF( AA( 1,1 ).EQ.AA( 2,2 ) .AND.
     &       ( AA( 2,1 ).EQ.0.0 .OR. AA( 1,2 ).EQ.0.0 ) ) THEN

            RNORM  = ABS( AA( 1,1 ) ) + ABS( AA( 1,2 ) ) +
     &               ABS( AA( 2,1 ) ) + ABS( AA( 2,2 ) )
            W  = TOL*RNORM
            EVEC( 2, 1 ) =   AA( 2, 1 ) / W
            EVEC( 1, 2 ) = - AA( 1, 2 ) / W

         ELSE

            EVEC( 2, 1 ) = AA( 2, 1 ) / ( EVAL( 1 ) - AA( 2,2 ) )
            EVEC( 1, 2 ) = AA( 1, 2 ) / ( EVAL( 2 ) - AA( 1,1 ) )

         END IF

         RETURN

      END IF
c                               ** Put s.p. matrix into d.p. matrix
      DO 20 J = 1, M

         DO 10 K = 1, M
            AAD( J, K ) = DBLE( AA( J,K ) )
   10    CONTINUE

   20 CONTINUE

c                                ** Initialize output variables
      IER  = 0

      DO 40 I = 1, M
         EVALD( I ) = ZERO

         DO 30 J = 1, M
            EVECD( I, J ) = ZERO
   30    CONTINUE

         EVECD( I, I ) = ONE
   40 CONTINUE

c                  ** Balance the input matrix and reduce its norm by
c                     diagonal similarity transformation stored in WK;
c                     then search for rows isolating an eigenvalue
c                     and push them down
      RNORM  = ZERO
      L  = 1
      K  = M

   50 CONTINUE
      KKK  = K

      DO 90 J = KKK, 1, -1

         ROW  = ZERO

         DO 60 I = 1, K

            IF( I.NE.J ) ROW  = ROW + ABS( AAD( J,I ) )

   60    CONTINUE

         IF( ROW.EQ.ZERO ) THEN

            WKD( K ) = J

            IF( J.NE.K ) THEN

               DO 70 I = 1, K
                  REPL        = AAD( I, J )
                  AAD( I, J ) = AAD( I, K )
                  AAD( I, K ) = REPL
   70          CONTINUE

               DO 80 I = L, M
                  REPL        = AAD( J, I )
                  AAD( J, I ) = AAD( K, I )
                  AAD( K, I ) = REPL
   80          CONTINUE

            END IF

            K  = K - 1
            GO TO  50

         END IF

   90 CONTINUE
c                                ** Search for columns isolating an
c                                   eigenvalue and push them left
  100 CONTINUE
      LLL  = L

      DO 140 J = LLL, K

         COL  = ZERO

         DO 110 I = L, K

            IF( I.NE.J ) COL  = COL + ABS( AAD( I,J ) )

  110    CONTINUE

         IF( COL.EQ.ZERO ) THEN

            WKD( L ) = J

            IF( J.NE.L ) THEN

               DO 120 I = 1, K
                  REPL        = AAD( I, J )
                  AAD( I, J ) = AAD( I, L )
                  AAD( I, L ) = REPL
  120          CONTINUE

               DO 130 I = L, M
                  REPL        = AAD( J, I )
                  AAD( J, I ) = AAD( L, I )
                  AAD( L, I ) = REPL
  130          CONTINUE

            END IF

            L  = L + 1
            GO TO  100

         END IF

  140 CONTINUE

c                           ** Balance the submatrix in rows L through K
      DO 150 I = L, K
         WKD( I ) = ONE
  150 CONTINUE

  160 CONTINUE
      NOCONV = .FALSE.

      DO 220 I = L, K

         COL  = ZERO
         ROW  = ZERO

         DO 170 J = L, K

            IF( J.NE.I ) THEN

               COL  = COL + ABS( AAD( J,I ) )
               ROW  = ROW + ABS( AAD( I,J ) )

            END IF

  170    CONTINUE

         F  = ONE
         G  = ROW / C5
         H  = COL + ROW

  180    CONTINUE
         IF( COL.LT.G ) THEN

            F    = F*C5
            COL  = COL*C6
            GO TO  180

         END IF

         G  = ROW*C5

  190    CONTINUE
         IF( COL.GE.G ) THEN

            F    = F / C5
            COL  = COL / C6
            GO TO  190

         END IF
c                                                ** Now balance
         IF( ( COL + ROW ) / F.LT.C4*H ) THEN

            WKD( I ) = WKD( I )*F
            NOCONV = .TRUE.

            DO 200 J = L, M
               AAD( I, J ) = AAD( I, J ) / F
  200       CONTINUE

            DO 210 J = 1, K
               AAD( J, I ) = AAD( J, I )*F
  210       CONTINUE

         END IF

  220 CONTINUE


      IF( NOCONV ) GO TO  160
c                                   ** Is A already in Hessenberg form?
      IF( K-1 .LT. L+1 ) GO TO  370

c                                   ** Transfer A to a Hessenberg form
      DO 310 N = L + 1, K - 1

         H  = ZERO
         WKD( N + M ) = ZERO
         SCALE  = ZERO
c                                                 ** Scale column
         DO 230 I = N, K
            SCALE  = SCALE + ABS( AAD( I,N - 1 ) )
  230    CONTINUE

         IF( SCALE.NE.ZERO ) THEN

            DO 240 I = K, N, -1
               WKD( I + M ) = AAD( I, N - 1 ) / SCALE
               H  = H + WKD( I + M )**2
  240       CONTINUE

            G    = - SIGN( SQRT( H ), WKD( N + M ) )
            H    = H - WKD( N + M )*G
            WKD( N + M ) = WKD( N + M ) - G
c                                            ** Form (I-(U*UT)/H)*A
            DO 270 J = N, M

               F  = ZERO

               DO 250 I = K, N, -1
                  F  = F + WKD( I + M )*AAD( I, J )
  250          CONTINUE

               DO 260 I = N, K
                  AAD( I, J ) = AAD( I, J ) - WKD( I + M )*F / H
  260          CONTINUE

  270       CONTINUE
c                                    ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 300 I = 1, K

               F  = ZERO

               DO 280 J = K, N, -1
                  F  = F + WKD( J + M )*AAD( I, J )
  280          CONTINUE

               DO 290 J = N, K
                  AAD( I, J ) = AAD( I, J ) - WKD( J + M )*F / H
  290          CONTINUE

  300       CONTINUE

            WKD( N + M ) = SCALE*WKD( N + M )
            AAD( N, N - 1 ) = SCALE*G

         END IF

  310 CONTINUE


      DO 360 N = K - 2, L, -1

         N1   = N + 1
         N2   = N + 2
         F  = AAD( N + 1, N )

         IF( F.NE.ZERO ) THEN

            F  = F*WKD( N + 1 + M )

            DO 320 I = N + 2, K
               WKD( I + M ) = AAD( I, N )
  320       CONTINUE

            IF( N + 1.LE.K ) THEN

               DO 350 J = 1, M

                  G  = ZERO

                  DO 330 I = N + 1, K
                     G  = G + WKD( I + M )*EVECD( I, J )
  330             CONTINUE

                  G  = G / F

                  DO 340 I = N + 1, K
                     EVECD( I, J ) = EVECD( I, J ) + G*WKD( I + M )
  340             CONTINUE

  350          CONTINUE

            END IF

         END IF

  360 CONTINUE


  370 CONTINUE

      N  = 1

      DO 390 I = 1, M

         DO 380 J = N, M
            RNORM  = RNORM + ABS( AAD( I,J ) )
  380    CONTINUE

         N  = I

         IF( I.LT.L .OR. I.GT.K ) EVALD( I ) = AAD( I, I )

  390 CONTINUE

      N  = K
      T  = ZERO

c                                      ** Search for next eigenvalues
  400 CONTINUE
      IF( N.LT.L ) GO TO  550

      IN  = 0
      N1  = N - 1
      N2  = N - 2
c                          ** Look for single small sub-diagonal element
  410 CONTINUE

      DO 420 I = L, N
         LB  = N + L - I

         IF( LB.EQ.L ) GO TO  430

         S  = ABS( AAD( LB - 1,LB - 1 ) ) + ABS( AAD( LB,LB ) )

         IF( S.EQ.ZERO ) S  = RNORM

         IF( ABS( AAD( LB, LB-1 ) ).LE. TOL*S ) GO TO  430

  420 CONTINUE


  430 CONTINUE
      X  = AAD( N, N )

      IF( LB.EQ.N ) THEN
c                                        ** One eigenvalue found
         AAD( N, N ) = X + T
         EVALD( N ) = AAD( N, N )
         N  = N1
         GO TO  400

      END IF

C next line has been included to avoid run time error caused by xlf

      IF ( ( N1.LE.0 ).OR.( N.LE.0 ) ) THEN
        WRITE(0,*) 'Subscript out of bounds in ASYMTX'
        STOP 9999
      ENDIF

      Y  = AAD( N1, N1 )
      W  = AAD( N, N1 )*AAD( N1, N )

      IF( LB.EQ.N1 ) THEN
c                                        ** Two eigenvalues found
         P  = ( Y - X )*C2
         Q  = P**2 + W
         Z  = SQRT( ABS( Q ) )
         AAD( N, N ) = X + T
         X  = AAD( N, N )
         AAD( N1, N1 ) = Y + T
c                                        ** Real pair
         Z  = P + SIGN( Z, P )
         EVALD( N1 ) = X + Z
         EVALD( N ) = EVALD( N1 )

         IF( Z.NE.ZERO ) EVALD( N ) = X - W / Z

         X  = AAD( N, N1 )
c                                  ** Employ scale factor in case
c                                     X and Z are very small
         R  = SQRT( X*X + Z*Z )
         P  = X / R
         Q  = Z / R
c                                             ** Row modification
         DO 440 J = N1, M
            Z  = AAD( N1, J )
            AAD( N1, J ) = Q*Z + P*AAD( N, J )
            AAD( N, J ) = Q*AAD( N, J ) - P*Z
  440    CONTINUE
c                                             ** Column modification
         DO 450 I = 1, N
            Z  = AAD( I, N1 )
            AAD( I, N1 ) = Q*Z + P*AAD( I, N )
            AAD( I, N ) = Q*AAD( I, N ) - P*Z
  450    CONTINUE
c                                          ** Accumulate transformations
         DO 460 I = L, K
            Z  = EVECD( I, N1 )
            EVECD( I, N1 ) = Q*Z + P*EVECD( I, N )
            EVECD( I, N ) = Q*EVECD( I, N ) - P*Z
  460    CONTINUE

         N  = N2
         GO TO  400

      END IF


      IF( IN.EQ.30 ) THEN

c                    ** No convergence after 30 iterations; set error
c                       indicator to the index of the current eigenvalue
         IER  = N
         GO TO  700

      END IF
c                                                  ** Form shift
      IF( IN.EQ.10 .OR. IN.EQ.20 ) THEN

         T  = T + X

         DO 470 I = L, N
            AAD( I, I ) = AAD( I, I ) - X
  470    CONTINUE

         S  = ABS( AAD( N,N1 ) ) + ABS( AAD( N1,N2 ) )
         X  = C3*S
         Y  = X
         W  = -C1*S**2

      END IF


      IN  = IN + 1

c                ** Look for two consecutive small sub-diagonal elements

C inhibit vectorization by CF77, as this will cause a run time error

CDIR$ NEXTSCALAR
      DO 480 J = LB, N2
         I  = N2 + LB - J
         Z  = AAD( I, I )
         R  = X - Z
         S  = Y - Z
         P  = ( R*S - W ) / AAD( I + 1, I ) + AAD( I, I + 1 )
         Q  = AAD( I + 1, I + 1 ) - Z - R - S
         R  = AAD( I + 2, I + 1 )
         S  = ABS( P ) + ABS( Q ) + ABS( R )
         P  = P / S
         Q  = Q / S
         R  = R / S

         IF( I.EQ.LB ) GO TO  490

         UU   = ABS( AAD( I, I-1 ) )*( ABS( Q ) + ABS( R ) )
         VV   = ABS( P ) * ( ABS( AAD( I-1, I-1 ) ) + ABS( Z ) +
     &                       ABS( AAD( I+1, I+1 ) ) )

         IF( UU .LE. TOL*VV ) GO TO  490

  480 CONTINUE

  490 CONTINUE
      AAD( I+2, I ) = ZERO

c                      ** fpp vectorization of this loop triggers
c                         array bounds errors, so inhibit
CFPP$ NOVECTOR L
      DO 500 J = I + 3, N
         AAD( J, J - 2 ) = ZERO
         AAD( J, J - 3 ) = ZERO
  500 CONTINUE

c             ** Double QR step involving rows K to N and columns M to N

      DO 540 KA = I, N1

         NOTLAS = KA.NE.N1

         IF( KA.EQ.I ) THEN

            S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )

            IF( LB.NE.I ) AAD( KA, KA - 1 ) = -AAD( KA, KA - 1 )

         ELSE

            P  = AAD( KA, KA - 1 )
            Q  = AAD( KA + 1, KA - 1 )
            R  = ZERO

            IF( NOTLAS ) R  = AAD( KA + 2, KA - 1 )

            X  = ABS( P ) + ABS( Q ) + ABS( R )

            IF( X.EQ.ZERO ) GO TO  540

            P  = P / X
            Q  = Q / X
            R  = R / X
            S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
            AAD( KA, KA - 1 ) = -S*X

         END IF

         P  = P + S
         X  = P / S
         Y  = Q / S
         Z  = R / S
         Q  = Q / P
         R  = R / P
c                                              ** Row modification
         DO 510 J = KA, M

            P  = AAD( KA, J ) + Q*AAD( KA + 1, J )

            IF( NOTLAS ) THEN

               P  = P + R*AAD( KA + 2, J )
               AAD( KA + 2, J ) = AAD( KA + 2, J ) - P*Z

            END IF

            AAD( KA + 1, J ) = AAD( KA + 1, J ) - P*Y
            AAD( KA, J ) = AAD( KA, J ) - P*X
  510    CONTINUE
c                                                 ** Column modification
         DO 520 II = 1, MIN( N, KA + 3 )

            P  = X*AAD( II, KA ) + Y*AAD( II, KA + 1 )

            IF( NOTLAS ) THEN

               P  = P + Z*AAD( II, KA + 2 )
               AAD( II, KA + 2 ) = AAD( II, KA + 2 ) - P*R

            END IF

            AAD( II, KA + 1 ) = AAD( II, KA + 1 ) - P*Q
            AAD( II, KA ) = AAD( II, KA ) - P
  520    CONTINUE
c                                          ** Accumulate transformations
         DO 530 II = L, K

            P  = X*EVECD( II, KA ) + Y*EVECD( II, KA + 1 )

            IF( NOTLAS ) THEN

               P  = P + Z*EVECD( II, KA + 2 )
               EVECD( II, KA + 2 ) = EVECD( II, KA + 2 ) - P*R

            END IF

            EVECD( II, KA + 1 ) = EVECD( II, KA + 1 ) - P*Q
            EVECD( II, KA ) = EVECD( II, KA ) - P
  530    CONTINUE

  540 CONTINUE

      GO TO  410
c                     ** All evals found, now backsubstitute real vector
  550 CONTINUE

      IF( RNORM.NE.ZERO ) THEN

         DO 580 N = M, 1, -1
            N2   = N
            AAD( N, N ) = ONE

            DO 570 I = N - 1, 1, -1
               W  = AAD( I, I ) - EVALD( N )

               IF( W.EQ.ZERO ) W  = TOL*RNORM

               R  = AAD( I, N )

               DO 560 J = N2, N - 1
                  R  = R + AAD( I, J )*AAD( J, N )
  560          CONTINUE

               AAD( I, N ) = -R / W
               N2   = I
  570       CONTINUE

  580    CONTINUE
c                      ** End backsubstitution vectors of isolated evals
         DO 600 I = 1, M

            IF( I.LT.L .OR. I.GT.K ) THEN

               DO 590 J = I, M
                  EVECD( I, J ) = AAD( I, J )
  590          CONTINUE

            END IF

  600    CONTINUE
c                                   ** Multiply by transformation matrix
         IF( K.NE.0 ) THEN

            DO 630 J = M, L, -1

               DO 620 I = L, K
                  Z  = ZERO

                  DO 610 N = L, MIN( J, K )
                     Z  = Z + EVECD( I, N )*AAD( N, J )
  610             CONTINUE

                  EVECD( I, J ) = Z
  620          CONTINUE

  630       CONTINUE

         END IF

      END IF


      DO 650 I = L, K

         DO 640 J = 1, M
            EVECD( I, J ) = EVECD( I, J )*WKD( I )
  640    CONTINUE
  650 CONTINUE

c                           ** Interchange rows if permutations occurred
      DO 670 I = L-1, 1, -1

         J  = WKD( I )

         IF( I.NE.J ) THEN

            DO 660 N = 1, M
               REPL   = EVECD( I, N )
               EVECD( I, N ) = EVECD( J, N )
               EVECD( J, N ) = REPL
  660       CONTINUE

         END IF

  670 CONTINUE


      DO 690 I = K + 1, M

         J  = WKD( I )

         IF( I.NE.J ) THEN

            DO 680 N = 1, M
               REPL   = EVECD( I, N )
               EVECD( I, N ) = EVECD( J, N )
               EVECD( J, N ) = REPL
  680       CONTINUE

         END IF

  690 CONTINUE

c                         ** Put results into output arrays
  700 CONTINUE

      DO 720 J = 1, M

         EVAL( J ) = EVALD( J )

         DO 710 K = 1, M
            EVEC( J, K ) = EVECD( J, K )
  710    CONTINUE

  720 CONTINUE


      END

      SUBROUTINE CHEKIN( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
     &                   UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
     &                   FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP, TEMIS,
     &                   PLANK, ONLYFL, ACCUR, TAUC, MAXCLY, MAXULV,
     &                   MAXUMU, MAXCMU, MAXPHI, MXCLY, MXULV, MXUMU,
     &                   MXCMU, MXPHI )

c           Checks the input dimensions and variables

c   Calls- WRTBAD, WRTDIM, DREF, ERRMSG
c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      LOGICAL   LAMBER, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND, MAXCLY, MAXCMU, MAXPHI, MAXULV, MAXUMU, MXCLY,
     &          MXCMU, MXPHI, MXULV, MXUMU, NLYR, NPHI, NSTR, NTAU, NUMU
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0, WVNMHI, WVNMLO
c     ..
c     .. Array Arguments ..

      REAL      DTAUC( MAXCLY ), HL( 0:MAXCMU ), PHI( MAXPHI ),
     &          PMOM( 0:MAXCMU, MAXCLY ), SSALB( MAXCLY ),
     &          TAUC( 0:MXCLY ), TEMPER( 0:MAXCLY ), UMU( MAXUMU ),
     &          UTAU( MAXULV )
c     ..
c     .. Local Scalars ..

      LOGICAL   INPERR
      INTEGER   IRMU, IU, J, K, LC, LU
      REAL      FLXALB, RMU
c     ..
c     .. External Functions ..

      LOGICAL   WRTBAD, WRTDIM
      REAL      DREF
      EXTERNAL  WRTBAD, WRTDIM, DREF
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MOD
c     ..


      INPERR = .FALSE.

      IF( NLYR.LT.1 ) INPERR = WRTBAD( 'NLYR' )

      IF( NLYR.GT.MAXCLY ) INPERR = WRTBAD( 'MAXCLY' )

      DO 20 LC = 1, NLYR

         IF( DTAUC( LC ).LT.0.0 ) INPERR = WRTBAD( 'DTAUC' )

         IF( SSALB( LC ).LT.0.0 .OR. SSALB( LC ).GT.1.0 )
     &       INPERR = WRTBAD( 'SSALB' )

         IF( PLANK .AND. IBCND.NE.1 ) THEN

            IF( LC.EQ.1 .AND. TEMPER( 0 ).LT.0.0 )
     &          INPERR = WRTBAD( 'TEMPER' )

            IF( TEMPER( LC ).LT.0.0 ) INPERR = WRTBAD( 'TEMPER' )

         END IF

         DO 10 K = 0, NSTR

            IF( PMOM( K,LC ).LT.-1.0 .OR. PMOM( K,LC ).GT.1.0 )
     &          INPERR = WRTBAD( 'PMOM' )

   10    CONTINUE

   20 CONTINUE


      IF( IBCND.EQ.1 ) THEN

         IF( MAXULV.LT.2 ) INPERR = WRTBAD( 'MAXULV' )

      ELSE IF( USRTAU ) THEN

         IF( NTAU.LT.1 ) INPERR = WRTBAD( 'NTAU' )

         IF( MAXULV.LT.NTAU ) INPERR = WRTBAD( 'MAXULV' )

         DO 30 LU = 1, NTAU

            IF( ABS( UTAU( LU ) - TAUC( NLYR ) ).LE. 1.E-4 )
     &          UTAU( LU ) = TAUC( NLYR )

            IF( UTAU( LU ).LT.0.0 .OR. UTAU( LU ).GT. TAUC( NLYR ) )
     &          INPERR = WRTBAD( 'UTAU' )

   30    CONTINUE

      ELSE

         IF( MAXULV.LT.NLYR + 1 ) INPERR = WRTBAD( 'MAXULV' )

      END IF


      IF( NSTR.LT.2 .OR. MOD( NSTR,2 ).NE.0 ) INPERR = WRTBAD( 'NSTR' )

c     IF( NSTR.EQ.2 )
c    &    CALL ERRMSG( 'CHEKIN--2 streams not recommended;'//
c    &                 ' use specialized 2-stream code instead',.False.)

      IF( NSTR.GT.MAXCMU ) INPERR = WRTBAD( 'MAXCMU' )

      IF( USRANG ) THEN

         IF( NUMU.LT.0 ) INPERR = WRTBAD( 'NUMU' )

         IF( .NOT.ONLYFL .AND. NUMU.EQ.0 ) INPERR = WRTBAD( 'NUMU' )

         IF( NUMU.GT.MAXUMU ) INPERR = WRTBAD( 'MAXUMU' )

         IF( IBCND.EQ.1 .AND. 2*NUMU.GT.MAXUMU )
     &       INPERR = WRTBAD( 'MAXUMU' )

         DO 40 IU = 1, NUMU

            IF( UMU( IU ).LT.-1.0 .OR. UMU( IU ).GT.1.0 .OR.
     &          UMU( IU ).EQ.0.0 ) INPERR = WRTBAD( 'UMU' )

            IF( IBCND.EQ.1 .AND. UMU( IU ).LT.0.0 )
     &          INPERR = WRTBAD( 'UMU' )

            IF( IU.GT.1 ) THEN

               IF( UMU( IU ).LT.UMU( IU-1 ) ) INPERR = WRTBAD( 'UMU' )

            END IF

   40    CONTINUE

      ELSE

         IF( MAXUMU.LT.NSTR ) INPERR = WRTBAD( 'MAXUMU' )

      END IF


      IF( .NOT.ONLYFL .AND. IBCND.NE.1 ) THEN

         IF( NPHI.LE.0 ) INPERR = WRTBAD( 'NPHI' )

         IF( NPHI.GT.MAXPHI ) INPERR = WRTBAD( 'MAXPHI' )

         DO 50 J = 1, NPHI

            IF( PHI( J ).LT.0.0 .OR. PHI( J ).GT.360.0 )
     &          INPERR = WRTBAD( 'PHI' )

   50    CONTINUE

      END IF


      IF( IBCND.LT.0 .OR. IBCND.GT.1 ) INPERR = WRTBAD( 'IBCND' )

      IF( IBCND.EQ.0 ) THEN

         IF( FBEAM.LT.0.0 ) INPERR = WRTBAD( 'FBEAM' )

         IF( FBEAM.GT.0.0 .AND. abs(UMU0).GT.1.0 )
     &       INPERR = WRTBAD( 'UMU0' )

         IF( FBEAM.GT.0.0 .AND. ( PHI0.LT.0.0 .OR.PHI0.GT.360.0 ) )
     &       INPERR = WRTBAD( 'PHI0' )

         IF( FISOT.LT.0.0 ) INPERR = WRTBAD( 'FISOT' )

         IF( LAMBER ) THEN

            IF( ALBEDO.LT.0.0 .OR. ALBEDO.GT.1.0 )
     &          INPERR = WRTBAD( 'ALBEDO' )

         ELSE
c                    ** Make sure flux albedo at dense mesh of incident
c                       angles does not assume unphysical values

            DO 60 IRMU = 0, 100
               RMU  = IRMU*0.01
               FLXALB = DREF( RMU, HL, NSTR )

               IF( FLXALB.LT.0.0 .OR. FLXALB.GT.1.0 )
     &             INPERR = WRTBAD( 'HL' )

   60       CONTINUE

         END IF


      ELSE IF( IBCND.EQ.1 ) THEN

         IF( ALBEDO.LT.0.0 .OR. ALBEDO.GT.1.0 )
     &       INPERR = WRTBAD( 'ALBEDO' )

      END IF


      IF( PLANK .AND. IBCND.NE.1 ) THEN

         IF( WVNMLO.LT.0.0 .OR. WVNMHI.LE.WVNMLO )
     &       INPERR = WRTBAD( 'WVNMLO,HI' )

         IF( TEMIS.LT.0.0 .OR. TEMIS.GT.1.0 ) INPERR = WRTBAD( 'TEMIS' )

         IF( BTEMP.LT.0.0 ) INPERR = WRTBAD( 'BTEMP' )

         IF( TTEMP.LT.0.0 ) INPERR = WRTBAD( 'TTEMP' )

      END IF


      IF( ACCUR.LT.0.0 .OR. ACCUR.GT.1.E-2 ) INPERR = WRTBAD( 'ACCUR' )

      IF( MXCLY.LT.NLYR ) INPERR = WRTDIM( 'MXCLY', NLYR )

      IF( IBCND.NE.1 ) THEN

         IF( USRTAU .AND. MXULV.LT.NTAU )
     &       INPERR = WRTDIM( 'MXULV',NTAU )

         IF( .NOT.USRTAU .AND. MXULV .LT. NLYR + 1 )
     &       INPERR = WRTDIM( 'MXULV', NLYR + 1 )

      ELSE

         IF( MXULV.LT.2 ) INPERR = WRTDIM( 'MXULV', 2 )

      END IF

      IF( MXCMU.LT.NSTR ) INPERR = WRTDIM( 'MXCMU', NSTR )

      IF( USRANG .AND. MXUMU.LT.NUMU ) INPERR = WRTDIM( 'MXUMU', NUMU )

      IF( USRANG .AND. IBCND.EQ.1 .AND.MXUMU.LT.2*NUMU )
     &    INPERR = WRTDIM( 'MXUMU', NUMU )

      IF( .NOT.USRANG .AND. MXUMU.LT.NSTR )
     &    INPERR = WRTDIM( 'MXUMU', NSTR )

      IF( .NOT.ONLYFL .AND. IBCND.NE.1 .AND. MXPHI.LT.NPHI )
     &    INPERR = WRTDIM( 'MXPHI', NPHI )

      IF( INPERR )
     &    CALL ERRMSG( 'DISORT--input and/or dimension errors',.True.)

      IF( PLANK ) THEN

         DO 70 LC = 1, NLYR

            IF( ABS( TEMPER( LC ) - TEMPER( LC-1 ) ).GT. 20.0 )
     &          CALL ERRMSG('CHEKIN--vertical temperature step may'
     &                      // ' be too large for good accuracy',
     &                      .False.)
   70    CONTINUE

      END IF

      END

      SUBROUTINE FLUXES( tausla, tauslau,
     &                   CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,
     &                   MAXULV, MXCMU, MXULV, NCUT, NN, NSTR, NTAU, PI,
     &                   PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR, XR0,
     &                   XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP, FLDN, FLDIR,
     &                   RFLDIR, RFLDN, UAVG, U0C,
     &                   uavgso, uavgup, uavgdn,
     $                   sindir, sinup, sindn)

c       Calculates the radiative fluxes, mean intensity, and flux
c       derivative with respect to optical depth from the m=0 intensity
c       components (the azimuthally-averaged intensity)

c    I N P U T     V A R I A B L E S:

c       CMU      :  Abscissae for Gauss quadrature over angle cosine
c       CWT      :  Weights for Gauss quadrature over angle cosine
c       GC       :  Eigenvectors at polar quadrature angles, SC(1)
c       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
c       LAYRU    :  Layer number of user level UTAU
c       LL       :  Constants of integration in Eq. SC(1), obtained
c                     by solving scaled version of Eq. SC(5);
c                     exponential term of Eq. SC(12) not included
c       LYRCUT   :  Logical flag for truncation of comput. layer
c       NN       :  Order of double-Gauss quadrature (NSTR/2)
c       NCUT     :  Number of computational layer where absorption
c                     optical depth exceeds ABSCUT
c       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
c       UTAUPR   :  Optical depths of user output levels in delta-M
c                     coordinates;  equal to UTAU if no delta-M
c       XR0      :  Expansion of thermal source function in Eq. SS(14)
c       XR1      :  Expansion of thermal source function Eqs. SS(16)
c       ZZ       :  Beam source vectors in Eq. SS(19)
c       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
c       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
c       (remainder are DISORT input variables)


c                   O U T P U T     V A R I A B L E S:

c       U0C      :  Azimuthally averaged intensities
c                   ( at polar quadrature angles )
c       (RFLDIR, RFLDN, FLUP, DFDT, UAVG are DISORT output variables)


c                   I N T E R N A L       V A R I A B L E S:

c       DIRINT   :  Direct intensity attenuated
c       FDNTOT   :  Total downward flux (direct + diffuse)
c       FLDIR    :  Direct-beam flux (delta-M scaled)
c       FLDN     :  Diffuse down-flux (delta-M scaled)
c       FNET     :  Net flux (total-down - diffuse-up)
c       FACT     :  EXP( - UTAUPR / UMU0 )
c       PLSORC   :  Planck source function (thermal)
c       ZINT     :  Intensity of m = 0 case, in Eq. SC(1)

c   Called by- DISORT
c   Calls- ZEROIT
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..

      LOGICAL   LYRCUT
      INTEGER   MAXULV, MXCMU, MXULV, NCUT, NN, NSTR, NTAU
      REAL      FBEAM, PI, UMU0
c     ..
c     .. Array Arguments ..

      LOGICAL   PRNT( * )
      INTEGER   LAYRU( MXULV )
      REAL      CMU( MXCMU ), CWT( MXCMU ), DFDT( MAXULV ),
     &          FLDIR( MXULV ), FLDN( MXULV ), FLUP( MAXULV ),
     &          GC( MXCMU, MXCMU, * ), KK( MXCMU, * ), LL( MXCMU, * ),
     &          RFLDIR( MAXULV ), RFLDN( MAXULV ), SSALB( * ),
     &          TAUCPR( 0:* ), U0C( MXCMU, MXULV ), UAVG( MAXULV ),
     &          UTAU( MAXULV ), UTAUPR( MXULV ), XR0( * ), XR1( * ),
     &          ZPLK0( MXCMU, * ), ZPLK1( MXCMU, * ), ZZ( MXCMU, * ),
     &          uavgso(*),uavgup(*), uavgdn(*),
     &          sindir(*),sinup(*), sindn(*)
      REAL tausla(0:*), tauslau(0:*)
c     ..
c     .. Local Scalars ..

      INTEGER   IQ, JQ, LU, LYU
      REAL      ANG1, ANG2, DIRINT, FACT, FDNTOT, FNET, PLSORC, ZINT
c     ..
c     .. External Subroutines ..

      EXTERNAL  ZEROIT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ACOS, EXP
c     ..


      IF( PRNT( 2 ) ) WRITE( *, 9000 )
c                                          ** Zero DISORT output arrays
      CALL ZEROIT( U0C, MXULV*MXCMU )
      CALL ZEROIT( FLDIR, MXULV )
      CALL ZEROIT( FLDN, MXULV )
      call  zeroit( uavgso,   maxulv )
      call  zeroit( uavgup,   maxulv )
      call  zeroit( uavgdn,   maxulv )
      call  zeroit( sindir,   maxulv )
      call  zeroit( sinup,    maxulv )
      call  zeroit( sindn,    maxulv )

c                                        ** Loop over user levels
      DO 80 LU = 1, NTAU

         LYU  = LAYRU( LU )

         IF( LYRCUT .AND. LYU.GT.NCUT ) THEN
c                                                ** No radiation reaches
c                                                ** this level
            FDNTOT = 0.0
            FNET   = 0.0
            PLSORC = 0.0
            GO TO  70

         END IF

         IF( FBEAM.GT.0.0 ) THEN
 
            FACT  = EXP( - tausla(LU-1) )
            DIRINT       = FBEAM*FACT
            FLDIR( LU )  = UMU0*( FBEAM*FACT )
            RFLDIR( LU ) = UMU0*FBEAM * EXP( -tauslau(lu-1) )
            sindir( LU ) = SQRT(1.-UMU0*UMU0)*FBEAM * 
     $                     EXP( -tauslau(lu-1) )

         ELSE

            DIRINT       = 0.0
            FLDIR( LU )  = 0.0
            RFLDIR( LU ) = 0.0
            sindir( LU ) = 0.0

         END IF


         DO 30 IQ = 1, NN

            ZINT   = 0.0

            DO 10 JQ = 1, NN
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU ) ) )
   10       CONTINUE

            DO 20 JQ = NN + 1, NSTR
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU - 1 ) ) )
   20       CONTINUE

            U0C( IQ, LU ) = ZINT

            IF( FBEAM.GT.0.0 ) U0C( IQ, LU ) = ZINT + ZZ( IQ, LYU )*FACT

            U0C( IQ, LU ) = U0C( IQ, LU ) + ZPLK0( IQ, LYU ) +
     &                      ZPLK1( IQ, LYU )*UTAUPR( LU )
            UAVG( LU ) = UAVG( LU ) + CWT( NN + 1 - IQ )*U0C( IQ, LU )
            uavgdn(lu) = uavgdn(lu) + cwt(nn+1-iq) * u0c( iq,lu )
            sindn(lu)  = sindn(lu)  + cwt(nn+1-iq) * 
     &                   SQRT(1.-CMU(NN+1-IQ)*CMU(NN+1-IQ))*
     &                   U0C( IQ, LU )
            FLDN( LU ) = FLDN( LU ) + CWT( NN + 1 - IQ )*
     &                   CMU( NN + 1 - IQ )*U0C( IQ, LU )
   30    CONTINUE


         DO 60 IQ = NN + 1, NSTR

            ZINT   = 0.0

            DO 40 JQ = 1, NN
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU ) ) )
   40       CONTINUE

            DO 50 JQ = NN + 1, NSTR
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU - 1 ) ) )
   50       CONTINUE

            U0C( IQ, LU ) = ZINT

            IF( FBEAM.GT.0.0 ) U0C( IQ, LU ) = ZINT + ZZ( IQ, LYU )*FACT

            U0C( IQ, LU ) = U0C( IQ, LU ) + ZPLK0( IQ, LYU ) +
     &                      ZPLK1( IQ, LYU )*UTAUPR( LU )
            UAVG( LU ) = UAVG( LU ) + CWT( IQ - NN )*U0C( IQ, LU )
            uavgup(lu) = uavgup(lu) + cwt(iq-nn) * u0c( iq,lu )
            sinup (lu) = sinup(lu)  + cwt(iq-nn) * 
     &                   SQRT(1.-CMU(IQ-NN)*CMU(IQ-NN))*
     &                   U0C( IQ, LU )
            FLUP( LU ) = FLUP( LU ) + CWT( IQ - NN )*CMU( IQ - NN )*
     &                   U0C( IQ, LU )
   60    CONTINUE


         FLUP( LU )  = 2.*PI*FLUP( LU )
         FLDN( LU )  = 2.*PI*FLDN( LU )
         FDNTOT      = FLDN( LU ) + FLDIR( LU )
         FNET        = FDNTOT - FLUP( LU )
         RFLDN( LU ) = FDNTOT - RFLDIR( LU )
         UAVG( LU )  = ( 2.*PI*UAVG( LU ) + DIRINT ) / ( 4.*PI )
         uavgso( lu ) = dirint / (4.*pi)
         uavgup( lu ) = (2.0 * pi * uavgup(lu) )/ (4.*pi)
         uavgdn( lu)  = (2.0 * pi * uavgdn(lu) )/ (4.*pi)
         sindn ( lu ) = 2.*PI*sindn ( LU )
         sinup ( lu ) = 2.*PI*sinup ( LU )

         PLSORC      = XR0( LYU ) + XR1( LYU )*UTAUPR( LU )
         DFDT( LU )  = ( 1.- SSALB( LYU ) ) * 4.*PI *
     &                 ( UAVG( LU ) - PLSORC )

   70    CONTINUE
         IF( PRNT( 2 ) ) WRITE( *, FMT = 9010 ) UTAU( LU ), LYU,
     &       RFLDIR( LU ), RFLDN( LU ), FDNTOT, FLUP( LU ), FNET,
     &       UAVG( LU ), PLSORC, DFDT( LU )

   80 CONTINUE


      IF( PRNT( 3 ) ) THEN

         WRITE( *, FMT = 9020 )

         DO 100 LU = 1, NTAU

            WRITE( *, FMT = 9030 ) UTAU( LU )

            DO 90 IQ = 1, NN
               ANG1   = 180./ PI* ACOS( CMU( 2*NN - IQ + 1 ) )
               ANG2   = 180./ PI* ACOS( CMU( IQ ) )
               WRITE( *, 9040 ) ANG1, CMU(2*NN-IQ+1), U0C(IQ,LU),
     $                          ANG2, CMU(IQ),        U0C(IQ+NN,LU)
   90       CONTINUE

  100    CONTINUE

      END IF


 9000 FORMAT( //, 21X,
     $ '<----------------------- FLUXES ----------------------->', /,
     $ '   Optical  Compu    Downward    Downward    Downward     ',
     $ ' Upward                    Mean      Planck   d(Net Flux)', /,
     $ '     Depth  Layer      Direct     Diffuse       Total     ',
     $ 'Diffuse         Net   Intensity      Source   / d(Op Dep)', / )
 9010 FORMAT( F10.4, I7, 1P, 7E12.3, E14.3 )
 9020 FORMAT( / , / , ' ******** AZIMUTHALLY AVERAGED INTENSITIES',
     &      ' ( at polar quadrature angles ) *******' )
 9030 FORMAT( /, ' Optical depth =', F10.4, //,
     $  '     Angle (deg)   cos(Angle)     Intensity',
     $  '     Angle (deg)   cos(Angle)     Intensity' )
 9040 FORMAT( 2( 0P,F16.4,F13.5,1P,E14.3 ) )

      END

      SUBROUTINE LEPOLY( NMU, M, MAXMU, TWONM1, MU, YLM )

c       Computes the normalized associated Legendre polynomial,
c       defined in terms of the associated Legendre polynomial
c       Plm = P-sub-l-super-m as

c             Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)

c       for fixed order m and all degrees from l = m to TWONM1.
c       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available
c       from a prior call to the routine.

c       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
c                  High-Order Associated Legendre Polynomials,
c                  J. Quant. Spectrosc. Radiat. Transfer 10,
c                  557-562, 1970.  (hereafter D/A)

c       METHOD: Varying degree recurrence relationship.

c       NOTE 1: The D/A formulas are transformed by
c               setting  M = n-1; L = k-1.
c       NOTE 2: Assumes that routine is called first with  M = 0,
c               then with  M = 1, etc. up to  M = TWONM1.
c       NOTE 3: Loops are written in such a way as to vectorize.

c  I N P U T     V A R I A B L E S:

c       NMU    :  Number of arguments of YLM
c       M      :  Order of YLM
c       MAXMU  :  First dimension of YLM
c       TWONM1 :  Max degree of YLM
c       MU(i)  :  Arguments of YLM (i = 1 to NMU)

c       If M.GT.0, YLM(M-1,i) for i = 1 to NMU is assumed to exist
c       from a prior call.

c  O U T P U T     V A R I A B L E:

c       YLM(l,i) :  l = M to TWONM1, normalized associated Legendre
c                   polynomials evaluated at argument MU(i)

c   Called by- DISORT, ALBTRN, SURFAC
c   Calls- ERRMSG
c +-------------------------------------------------------------------+

c     .. Parameters ..

      INTEGER   MAXSQT
      PARAMETER ( MAXSQT = 1000 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   M, MAXMU, NMU, TWONM1
c     ..
c     .. Array Arguments ..

      REAL      MU( * ), YLM( 0:MAXMU, * )
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      INTEGER   I, L, NS
      REAL      TMP1, TMP2
c     ..
c     .. Local Arrays ..

      REAL      SQT( MAXSQT )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC FLOAT, SQRT
c     ..
      SAVE      SQT, PASS1
      DATA      PASS1 / .TRUE. /


      IF( PASS1 ) THEN

         PASS1  = .FALSE.

         DO 10 NS = 1, MAXSQT
            SQT( NS ) = SQRT( FLOAT( NS ) )
   10    CONTINUE

      END IF

      IF( 2*TWONM1.GT.MAXSQT )
     &    CALL ERRMSG('LEPOLY--need to increase param MAXSQT',.True.)


      IF( M.EQ.0 ) THEN
c                             ** Upward recurrence for ordinary
c                                Legendre polynomials
         DO 20 I = 1, NMU
            YLM( 0, I ) = 1.0
            YLM( 1, I ) = MU( I )
   20    CONTINUE


         DO 40 L = 2, TWONM1

            DO 30 I = 1, NMU
               YLM( L, I ) = ( ( 2*L - 1 )*MU( I )*YLM( L - 1, I ) -
     &                         ( L - 1 )*YLM( L - 2, I ) ) / L
   30       CONTINUE

   40    CONTINUE


      ELSE

         DO 50 I = 1, NMU
c                               ** Y-sub-m-super-m; derived from
c                               ** D/A Eqs. (11,12)

            YLM( M, I ) = - SQT( 2*M - 1 ) / SQT( 2*M )*
     &                      SQRT( 1.- MU(I)**2 )*YLM( M - 1, I )

c                              ** Y-sub-(m+1)-super-m; derived from
c                              ** D/A Eqs.(13,14) using Eqs.(11,12)

            YLM( M + 1, I ) = SQT( 2*M + 1 )*MU( I )*YLM( M, I )

   50    CONTINUE

c                                   ** Upward recurrence; D/A EQ.(10)
         DO 70 L = M + 2, TWONM1

            TMP1  = SQT( L - M )*SQT( L + M )
            TMP2  = SQT( L - M - 1 )*SQT( L + M - 1 )

            DO 60 I = 1, NMU
               YLM( L, I ) = ( ( 2*L - 1 )*MU( I )*YLM( L-1, I ) -
     &                         TMP2*YLM( L-2, I ) ) / TMP1
   60       CONTINUE

   70    CONTINUE

      END IF


      END

      SUBROUTINE PRAVIN( UMU, NUMU, MAXUMU, UTAU, NTAU, U0U )

c        Print azimuthally averaged intensities at user angles

c   Called by- DISORT

c     LENFMT   Max number of polar angle cosines UMU that can be
c                printed on one line, as set in FORMAT statement
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   MAXUMU, NTAU, NUMU
c     ..
c     .. Array Arguments ..

      REAL      U0U( MAXUMU, NTAU ), UMU( NUMU ), UTAU( NTAU )
c     ..
c     .. Local Scalars ..

      INTEGER   IU, IUMAX, IUMIN, LENFMT, LU, NP, NPASS
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..


      IF( NUMU.LT.1 )  RETURN

      WRITE( *, '(//,A)' )
     &   ' *******  AZIMUTHALLY AVERAGED INTENSITIES ' //
     &   '(at user polar angles)  ********'

      LENFMT = 8
      NPASS  = 1 + (NUMU-1) / LENFMT

      WRITE( *,'(/,A,/,A)') '   Optical   Polar Angle Cosines',
     &                      '     Depth'

      DO 20 NP = 1, NPASS

         IUMIN  = 1 + LENFMT * ( NP - 1 )
         IUMAX  = MIN( LENFMT*NP, NUMU )
         WRITE( *,'(/,10X,8F14.5)') ( UMU(IU), IU = IUMIN, IUMAX )

         DO 10 LU = 1, NTAU
            WRITE( *, '(0P,F10.4,1P,8E14.4)' ) UTAU( LU ),
     &           ( U0U( IU,LU ), IU = IUMIN, IUMAX )
   10    CONTINUE

   20 CONTINUE


      END

      SUBROUTINE PRTINP( NLYR, DTAUC, DTAUCP, SSALB, PMOM, TEMPER,
     &                   WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,
     &                   NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT,
     &                   LAMBER, ALBEDO, HL, BTEMP, TTEMP, TEMIS,
     &                   DELTAM, PLANK, ONLYFL, ACCUR, FLYR, LYRCUT,
     &                   OPRIM, TAUC, TAUCPR, MAXCMU, PRTMOM )

c        Print values of input variables

c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      LOGICAL   DELTAM, LAMBER, LYRCUT, ONLYFL, PLANK, PRTMOM
      INTEGER   IBCND, MAXCMU, NLYR, NPHI, NSTR, NTAU, NUMU
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0, WVNMHI, WVNMLO
c     ..
c     .. Array Arguments ..

      REAL      DTAUC( * ), DTAUCP( * ), FLYR( * ), HL( 0:MAXCMU ),
     &          OPRIM( * ), PHI( * ), PMOM( 0:MAXCMU, * ), SSALB( * ),
     &          TAUC( 0:* ), TAUCPR( 0:* ), TEMPER( 0:* ), UMU( * ),
     &          UTAU( * )
c     ..
c     .. Local Scalars ..

      INTEGER   IU, J, K, LC, LU
      REAL      YESSCT
c     ..


      WRITE( *, '(/,A,I4,A,I4)' ) ' No. streams =', NSTR,
     &       '     No. computational layers =', NLYR

      IF( IBCND.NE.1 ) WRITE( *, '(I4,A,10F10.4,/,(26X,10F10.4))' )
     &    NTAU,' User optical depths :', ( UTAU(LU), LU = 1, NTAU )

      IF( .NOT.ONLYFL ) WRITE( *, '(I4,A,10F9.5,/,(31X,10F9.5))' )
     &    NUMU,' User polar angle cosines :',( UMU(IU), IU = 1, NUMU )

      IF( .NOT.ONLYFL .AND. IBCND.NE.1 )
     &    WRITE( *, '(I4,A,10F9.2,/,(28X,10F9.2))' )
     &           NPHI,' User azimuthal angles :',( PHI(J), J = 1, NPHI )

      IF( .NOT.PLANK .OR. IBCND.EQ.1 )
     &    WRITE( *, '(A)' ) ' No thermal emission'


      WRITE( *, '(A,I2)' ) ' Boundary condition flag: IBCND =', IBCND

      IF( IBCND.EQ.0 ) THEN

         WRITE( *, '(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P,E11.3)' )
     &          '    Incident beam with intensity =', FBEAM,
     &          ' and polar angle cosine = ', UMU0,
     &          '  and azimuth angle =', PHI0,
     &          '    plus isotropic incident intensity =', FISOT

         IF( LAMBER ) WRITE( *, '(A,0P,F8.4)' )
     &                '    Bottom albedo (Lambertian) =', ALBEDO

         IF( .NOT.LAMBER ) WRITE( *, '(A,/,(10X,10F9.5))' )
     &     '    Legendre coeffs of bottom bidirectional reflectivity :',
     &         ( HL( K ), K = 0, NSTR )

         IF( PLANK ) WRITE( *, '(A,2F14.4,/,A,F10.2,A,F10.2,A,F8.4)' )
     &       '    Thermal emission in wavenumber interval :', WVNMLO,
     &       WVNMHI,
     &       '    Bottom temperature =', BTEMP,
     &       '    Top temperature =', TTEMP,
     &       '    Top emissivity =',TEMIS

      ELSE IF( IBCND.EQ.1 ) THEN

         WRITE(*,'(A)') '    Isotropic illumination from top and bottom'
         WRITE( *, '(A,0P,F8.4)' )
     &          '    Bottom albedo (Lambertian) =', ALBEDO
      END IF


      IF( DELTAM ) WRITE( *, '(A)' ) ' Uses delta-M method'
      IF( .NOT.DELTAM ) WRITE( *, '(A)' ) ' Does not use delta-M method'


      IF( IBCND.EQ.1 ) THEN

         WRITE( *, '(A)' ) ' Calculate albedo and transmissivity of'//
     &                     ' medium vs. incident beam angle'

      ELSE IF( ONLYFL ) THEN

         WRITE( *, '(A)' )
     &          ' Calculate fluxes and azim-averaged intensities only'

      ELSE

         WRITE( *, '(A)' ) ' Calculate fluxes and intensities'

      END IF


      WRITE( *, '(A,1P,E11.2)' )
     &       ' Relative convergence criterion for azimuth series =',
     &       ACCUR

      IF( LYRCUT ) WRITE( *, '(A)' )
     &    ' Sets radiation = 0 below absorption optical depth 10'


c                                        ** Print layer variables
      IF( PLANK ) WRITE( *, FMT = 9180 )
      IF( .NOT.PLANK ) WRITE( *, FMT = 9190 )

      YESSCT = 0.0

      DO 10 LC = 1, NLYR

         YESSCT = YESSCT + SSALB( LC )

         IF( PLANK )
     &       WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4,F14.3)')
     &             LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),
     &             DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM(1,LC),
     &             TEMPER( LC-1 )

         IF( .NOT.PLANK )
     &       WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4)')
     &             LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),
     &             DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM( 1,LC )
   10 CONTINUE

      IF( PLANK ) WRITE( *, '(85X,F14.3)' ) TEMPER( NLYR )


      IF( PRTMOM .AND. YESSCT.GT.0.0 ) THEN

         WRITE( *, '(/,A)' ) ' Layer   Phase Function Moments'

         DO 20 LC = 1, NLYR

            IF( SSALB( LC ).GT.0.0 )
     &          WRITE( *, '(I6,10F11.6,/,(6X,10F11.6))' )
     &                 LC, ( PMOM( K, LC ), K = 0, NSTR )
   20    CONTINUE

      END IF

c                ** (Read every other line in these formats)

 9180 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,
     &'                   Total    Single                           ',
     &               'Total    Single', /,
     &'       Optical   Optical   Scatter   Truncated   ',
     &   'Optical   Optical   Scatter    Asymm', /,
     &'         Depth     Depth    Albedo    Fraction     ',
     &     'Depth     Depth    Albedo   Factor   Temperature' )
 9190 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,
     &'                   Total    Single                           ',
     &               'Total    Single', /,
     &'       Optical   Optical   Scatter   Truncated   ',
     &   'Optical   Optical   Scatter    Asymm', /,
     &'         Depth     Depth    Albedo    Fraction     ',
     &     'Depth     Depth    Albedo   Factor' )

      END

      SUBROUTINE PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                   MAXUMU )

c         Prints the intensity at user polar and azimuthal angles

c     All arguments are DISORT input or output variables

c   Called by- DISORT

c     LENFMT   Max number of azimuth angles PHI that can be printed
c                on one line, as set in FORMAT statement
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      INTEGER   MAXULV, MAXUMU, NPHI, NTAU, NUMU
c     ..
c     .. Array Arguments ..

      REAL      PHI( * ), UMU( * ), UTAU( * ), UU( MAXUMU, MAXULV, * )
c     ..
c     .. Local Scalars ..

      INTEGER   IU, J, JMAX, JMIN, LENFMT, LU, NP, NPASS
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..


      IF( NPHI.LT.1 )  RETURN

      WRITE( *, '(//,A)' )
     &   ' *********  I N T E N S I T I E S  *********'

      LENFMT = 10
      NPASS  = 1 + (NPHI-1) / LENFMT

      WRITE( *, '(/,A,/,A,/,A)' )
     &   '             Polar   Azimuth angles (degrees)',
     &   '   Optical   Angle',
     &   '    Depth   Cosine'

      DO 30 LU = 1, NTAU

         DO 20 NP = 1, NPASS

            JMIN   = 1 + LENFMT * ( NP - 1 )
            JMAX   = MIN( LENFMT*NP, NPHI )

            WRITE( *, '(/,18X,10F11.2)' ) ( PHI(J), J = JMIN, JMAX )

            IF( NP.EQ.1 ) WRITE( *, '(F10.4,F8.4,1P,10E11.3)' )
     &             UTAU(LU), UMU(1), (UU(1, LU, J), J = JMIN, JMAX)
            IF( NP.GT.1 ) WRITE( *, '(10X,F8.4,1P,10E11.3)' )
     &                       UMU(1), (UU(1, LU, J), J = JMIN, JMAX)

            DO 10 IU = 2, NUMU
               WRITE( *, '(10X,F8.4,1P,10E11.3)' ) 
     &                 UMU( IU ), ( UU( IU, LU, J ), J = JMIN, JMAX )
   10       CONTINUE

   20    CONTINUE

   30 CONTINUE


      END

      SUBROUTINE QGAUSN( M, GMU, GWT )

c       Compute weights and abscissae for ordinary Gaussian quadrature
c       on the interval (0,1);  that is, such that

c           sum(i=1 to M) ( GWT(i) f(GMU(i)) )

c       is a good approximation to

c           integral(0 to 1) ( f(x) dx )

c   INPUT :    M       order of quadrature rule

c   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M)
c             GWT(I)   array of weights (I = 1 TO M)

c   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
c                   Integration, Academic Press, New York, pp. 87, 1975

c   METHOD:  Compute the abscissae as roots of the Legendre
c            polynomial P-sub-M using a cubically convergent
c            refinement of Newton's method.  Compute the
c            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
c            that Newton's method can very easily diverge; only a
c            very good initial guess can guarantee convergence.
c            The initial guess used here has never led to divergence
c            even for M up to 1000.

c   ACCURACY:  relative error no better than TOL or computer
c              precision (machine epsilon), whichever is larger

c   INTERNAL VARIABLES:

c    ITER      : number of Newton Method iterations
c    MAXIT     : maximum allowed iterations of Newton Method
c    PM2,PM1,P : 3 successive Legendre polynomials
c    PPR       : derivative of Legendre polynomial
c    P2PRI     : 2nd derivative of Legendre polynomial
c    TOL       : convergence criterion for Legendre poly root iteration
c    X,XI      : successive iterates in cubically-convergent version
c                of Newtons Method (seeking roots of Legendre poly.)

c   Called by- SETDIS, SURFAC
c   Calls- D1MACH, ERRMSG
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..

      INTEGER   M
c     ..
c     .. Array Arguments ..

      REAL      GMU( M ), GWT( M )
c     ..
c     .. Local Scalars ..

      INTEGER   ITER, K, LIM, MAXIT, NN, NP1
      REAL      CONA, PI, T
      DOUBLE PRECISION EN, NNP1, ONE, P, P2PRI, PM1, PM2, PPR, PROD,
     &                 TMP, TOL, TWO, X, XI
c     ..
c     .. External Functions ..

      DOUBLE PRECISION D1MACH
      EXTERNAL  D1MACH
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, ASIN, COS, FLOAT, MOD, TAN
c     ..
      SAVE      PI, TOL

      DATA      PI / 0.0 / , MAXIT / 1000 / , ONE / 1.D0 / ,
     &          TWO / 2.D0 /


      IF( PI.EQ.0.0 ) THEN

         PI   = 2.*ASIN( 1.0 )
         TOL  = 10.*D1MACH( 4 )

      END IF


      IF( M.LT.1 ) CALL ERRMSG( 'QGAUSN--Bad value of M',.True.)

      IF( M.EQ.1 ) THEN

         GMU( 1 ) = 0.5
         GWT( 1 ) = 1.0
         RETURN

      END IF

      EN   = M
      NP1  = M + 1
      NNP1 = M*NP1
      CONA = FLOAT( M - 1 ) / ( 8*M**3 )

      LIM  = M / 2

      DO 30 K = 1, LIM
c                                        ** Initial guess for k-th root
c                                           of Legendre polynomial, from
c                                           Davis/Rabinowitz (2.7.3.3a)
         T  = ( 4*K - 1 )*PI / ( 4*M + 2 )
         X  = COS( T + CONA / TAN( T ) )
         ITER = 0
c                                        ** Upward recurrence for
c                                           Legendre polynomials
   10    CONTINUE
         ITER   = ITER + 1
         PM2    = ONE
         PM1    = X

         DO 20 NN = 2, M
            P    = ( ( 2*NN - 1 )*X*PM1 - ( NN - 1 )*PM2 ) / NN
            PM2  = PM1
            PM1  = P
   20    CONTINUE
c                                              ** Newton Method
         TMP    = ONE / ( ONE - X**2 )
         PPR    = EN*( PM2 - X*P )*TMP
         P2PRI  = ( TWO*X*PPR - NNP1*P )*TMP
         XI     = X - ( P / PPR )*( ONE +
     &            ( P / PPR )*P2PRI / ( TWO*PPR ) )

c                                              ** Check for convergence
         IF( ABS( XI - X ).GT.TOL ) THEN

            IF( ITER.GT.MAXIT )
     &          CALL ERRMSG( 'QGAUSN--max iteration count',.True.)

            X  = XI
            GO TO  10

         END IF
c                             ** Iteration finished--calculate weights,
c                                abscissae for (-1,1)
         GMU( K ) = -X
         GWT( K ) = TWO / ( TMP*( EN*PM2 )**2 )
         GMU( NP1 - K ) = -GMU( K )
         GWT( NP1 - K ) = GWT( K )
   30 CONTINUE
c                                    ** Set middle abscissa and weight
c                                       for rules of odd order
      IF( MOD( M,2 ).NE.0 ) THEN

         GMU( LIM + 1 ) = 0.0
         PROD   = ONE

         DO 40 K = 3, M, 2
            PROD   = PROD * K / ( K - 1 )
   40    CONTINUE

         GWT( LIM + 1 ) = TWO / PROD**2
      END IF

c                                        ** Convert from (-1,1) to (0,1)
      DO 50 K = 1, M
         GMU( K ) = 0.5*GMU( K ) + 0.5
         GWT( K ) = 0.5*GWT( K )
   50 CONTINUE


      END

      SUBROUTINE SETDIS( dsdh, nid, tausla, tauslau, mu2,
     &                   CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM,
     &                   FLYR, GL, HL, HLPR, IBCND, LAMBER, LAYRU,
     &                   LYRCUT, MAXUMU, MAXCMU, MXCMU, NCUT, NLYR,
     &                   NTAU, NN, NSTR, PLANK, NUMU, ONLYFL, OPRIM,
     &                   PMOM, SSALB, TAUC, TAUCPR, UTAU, UTAUPR, UMU,
     &                   UMU0, USRTAU, USRANG )

c          Perform miscellaneous setting-up operations

c       INPUT :  all are DISORT input variables (see DOC file)

c       OUTPUT:  NTAU,UTAU   if USRTAU = FALSE
c                NUMU,UMU    if USRANG = FALSE
c                CMU,CWT     computational polar angles and
c                               corresponding quadrature weights
c                EXPBEA      transmission of direct beam
c                FLYR        truncated fraction in delta-M method
c                GL          phase function Legendre coefficients multi-
c                              plied by (2L+1) and single-scatter albedo
c                HLPR        Legendre moments of surface bidirectional
c                              reflectivity, times 2K+1
c                LAYRU       Computational layer in which UTAU falls
c                LYRCUT      flag as to whether radiation will be zeroed
c                              below layer NCUT
c                NCUT        computational layer where absorption
c                              optical depth first exceeds  ABSCUT
c                NN          NSTR / 2
c                OPRIM       delta-M-scaled single-scatter albedo
c                TAUCPR      delta-M-scaled optical depth
c                UTAUPR      delta-M-scaled version of  UTAU

c   Called by- DISORT
c   Calls- QGAUSN, ERRMSG
c ----------------------------------------------------------------------

      INCLUDE 'params'

c     .. Scalar Arguments ..

      LOGICAL   DELTAM, LAMBER, LYRCUT, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND, MAXCMU, MAXUMU, MXCMU, NCUT, NLYR, NN, NSTR,
     &          NTAU, NUMU
      REAL      FBEAM, UMU0

c geometry
      REAL dsdh(0:kz,kz)
      INTEGER nid(0:kz)
      REAL tausla(0:kz), tauslau(0:kz), mu2(0:kz)
      REAL sum, sumu
c     ..
c     .. Array Arguments ..

      INTEGER   LAYRU( * )
      REAL      CMU( MXCMU ), CWT( MXCMU ), DTAUC( * ), DTAUCP( * ),
     &          EXPBEA( 0:* ), FLYR( * ), GL( 0:MXCMU, * ),
     &          HL( 0:MAXCMU ), HLPR( 0:MXCMU ), OPRIM( * ),
     &          PMOM( 0:MAXCMU, * ), SSALB( * ), TAUC( 0:* ),
     &          TAUCPR( 0:* ), UMU( MAXUMU ), UTAU( * ), UTAUPR( * )
c     ..
c     .. Local Scalars ..

      INTEGER   IQ, IU, K, LC, LU
      REAL      ABSCUT, ABSTAU, F

      REAL      R1MACH
      EXTERNAL  R1MACH
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, QGAUSN
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, EXP
c     ..
      DATA      ABSCUT / 10. /


      IF( .NOT.USRTAU ) THEN
c                              ** Set output levels at computational
c                                 layer boundaries
         NTAU  = NLYR + 1

         DO 10 LC = 0, NTAU - 1
            UTAU( LC + 1 ) = TAUC( LC )
   10    CONTINUE

      END IF
c                        ** Apply delta-M scaling and move description
c                           of computational layers to local variables
      EXPBEA( 0 ) = 1.0
      TAUCPR( 0 ) = 0.0
      ABSTAU      = 0.0
      do i = 0, kz
       tausla( i ) = 0.0
       tauslau( i ) = 0.0
       mu2(i) = 1./largest
      end do

      DO 40 LC = 1, NLYR

         PMOM( 0, LC ) = 1.0

         IF( ABSTAU.LT.ABSCUT ) NCUT  = LC

         ABSTAU = ABSTAU + ( 1.- SSALB( LC ) )*DTAUC( LC )

         IF( .NOT.DELTAM ) THEN

            OPRIM( LC )  = SSALB( LC )
            DTAUCP( LC ) = DTAUC( LC )
            TAUCPR( LC ) = TAUC( LC )

            DO 20 K = 0, NSTR - 1
               GL( K, LC ) = ( 2*K + 1 )*OPRIM( LC )*PMOM( K, LC )
   20       CONTINUE

            F  = 0.0


         ELSE
c                                    ** Do delta-M transformation

            F  = PMOM( NSTR, LC )
            OPRIM(LC) = SSALB(LC) * ( 1.- F ) / ( 1.- F * SSALB(LC) )
            DTAUCP( LC ) = ( 1.- F*SSALB( LC ) )*DTAUC( LC )
            TAUCPR( LC ) = TAUCPR( LC-1 ) + DTAUCP( LC )

            DO 30 K = 0, NSTR - 1
               GL( K, LC ) = ( 2*K + 1 ) * OPRIM( LC ) *
     &                       ( PMOM( K,LC ) - F ) / ( 1.- F )
   30       CONTINUE

         END IF

         FLYR( LC )   = F
         EXPBEA( LC ) = 0.0

   40 CONTINUE
c 
* calculate slant optical depth
*              
         IF(umu0 .LT. 0.0) THEN
           IF(nid(0) .LT. 0) THEN
             tausla(0) = largest
             tauslau(0) = largest
           ELSE
             sum = 0.0
             sumu = 0.0
             DO lc = 1, nid(0)
               sum = sum + 2.*dtaucp(lc)*dsdh(0,lc)
               sumu = sumu + 2.*dtauc(lc)*dsdh(0,lc)
             END DO
             tausla(0) = sum 
             tauslau(0) = sumu 
           END IF
         END IF

         expbea( 0 ) = EXP( -tausla( 0 ) )

*
         DO 41, lc = 1, nlyr
          IF(nid(lc) .LT. 0) THEN
            tausla(lc) = largest
            tauslau(lc) = largest
          ELSE
            sum = 0.0
            sumu = 0.0
            DO lu = 1, MIN(nid(lc),lc)
               sum = sum + dtaucp(lu)*dsdh(lc,lu)
               sumu = sumu + dtauc(lu)*dsdh(lc,lu)
            ENDDO
            DO lu = MIN(nid(lc),lc)+1,nid(lc)
               sum = sum + 2.*dtaucp(lu)*dsdh(lc,lu)
               sumu = sumu + 2.*dtauc(lu)*dsdh(lc,lu)
            ENDDO
            tausla(lc) = sum 
            tauslau(lc) = sumu 
            IF(tausla(lc) .EQ. tausla(lc-1)) THEN
              mu2(lc) = largest
            ELSE
              mu2(lc) = (taucpr(lc)-taucpr(lc-1))
     $                         /(tausla(lc)-tausla(lc-1))
              mu2(lc) = SIGN( AMAX1(ABS(mu2(lc)),1./largest),
     $                     mu2(lc) )
            END IF
          END IF
          expbea(lc) = EXP( -tausla( lc ) )
 41      CONTINUE

c                      ** If no thermal emission, cut off medium below
c                         absorption optical depth = ABSCUT ( note that
c                         delta-M transformation leaves absorption
c                         optical depth invariant ).  Not worth the
c                         trouble for one-layer problems, though.
      LYRCUT = .FALSE.

      IF( ABSTAU.GE.ABSCUT .AND. .NOT.PLANK .AND. IBCND.NE.1 .AND.
     &    NLYR.GT.1 ) LYRCUT = .TRUE.

      IF( .NOT.LYRCUT ) NCUT   = NLYR

c                             ** Set arrays defining location of user
c                             ** output levels within delta-M-scaled
c                             ** computational mesh
      DO 70 LU = 1, NTAU

         DO 50 LC = 1, NLYR

            IF( UTAU( LU ).GE.TAUC( LC - 1 ) .AND.
     &          UTAU( LU ).LE.TAUC( LC ) ) GO TO  60

   50    CONTINUE
         LC   = NLYR

   60    CONTINUE
         UTAUPR( LU ) = UTAU( LU )
         IF( DELTAM ) UTAUPR( LU ) = TAUCPR( LC - 1 ) +
     &                               ( 1.- SSALB( LC )*FLYR( LC ) )*
     &                               ( UTAU( LU ) - TAUC( LC-1 ) )
         LAYRU( LU ) = LC

   70 CONTINUE
c                      ** Calculate computational polar angle cosines
c                         and associated quadrature weights for Gaussian
c                         quadrature on the interval (0,1) (upward)
      NN   = NSTR / 2

      CALL QGAUSN( NN, CMU, CWT )
c                                  ** Downward (neg) angles and weights
      DO 80 IQ = 1, NN
         CMU( IQ + NN ) = - CMU( IQ )
         CWT( IQ + NN ) = CWT( IQ )
   80 CONTINUE


c     IF( FBEAM.GT.0.0 ) THEN
c                               ** Compare beam angle to comput. angles
         DO 90 IQ = 1, NN

C                      ** Dither mu2 if it is close to one of the 
C                         quadrature angles.

         DO  lc = 1, nlyr
          IF (  ABS(mu2(lc)) .lt. 1.E5 ) THEN
            IF( ABS( 1. - ABS(mu2(lc))/CMU( IQ ) ) .LT. 0.05 ) 
     &           mu2(lc) = mu2(lc)*0.999
          ENDIF
         END DO

   90    CONTINUE

c     END IF

      IF( .NOT.USRANG .OR. ( ONLYFL .AND. MAXUMU.GE.NSTR ) ) THEN

c                                   ** Set output polar angles to
c                                      computational polar angles
         NUMU   = NSTR

         DO 100 IU = 1, NN
            UMU( IU ) = - CMU( NN + 1 - IU )
  100    CONTINUE

         DO 110 IU = NN + 1, NSTR
            UMU( IU ) = CMU( IU - NN )
  110    CONTINUE

      END IF


      IF( USRANG .AND. IBCND.EQ.1 ) THEN

c                               ** Shift positive user angle cosines to
c                                  upper locations and put negatives
c                                  in lower locations
         DO 120 IU = 1, NUMU
            UMU( IU + NUMU ) = UMU( IU )
  120    CONTINUE

         DO 130 IU = 1, NUMU
            UMU( IU ) = -UMU( 2*NUMU + 1 - IU )
  130    CONTINUE

         NUMU   = 2*NUMU

      END IF


      IF( .NOT.LYRCUT .AND. .NOT.LAMBER ) THEN

         DO 140 K = 0, NSTR
            HLPR( K ) = ( 2*K + 1 )*HL( K )
  140    CONTINUE

      END IF


      END

      SUBROUTINE SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,
     &                   LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT,
     &                   NNLYRI, NN, NSTR, TAUCPR, WK )

c        Calculate coefficient matrix for the set of equations
c        obtained from the boundary conditions and the continuity-
c        of-intensity-at-layer-interface equations;  store in the
c        special banded-matrix format required by LINPACK routines

c     I N P U T      V A R I A B L E S:

c       BDR      :  Surface bidirectional reflectivity
c       CMU      :  Abscissae for Gauss quadrature over angle cosine
c       CWT      :  Weights for Gauss quadrature over angle cosine
c       DELM0    :  Kronecker delta, delta-sub-m0
c       GC       :  Eigenvectors at polar quadrature angles, SC(1)
c       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
c       LYRCUT   :  Logical flag for truncation of comput. layer
c       NN       :  Number of streams in a hemisphere (NSTR/2)
c       NCUT     :  Total number of computational layers considered
c       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
c       (remainder are DISORT input variables)

c   O U T P U T     V A R I A B L E S:

c       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
c                      scaled by Eq. SC(12); in banded form required
c                      by LINPACK solution routines
c       NCOL     :  Counts of columns in CBAND

c   I N T E R N A L    V A R I A B L E S:

c       IROW     :  Points to row in CBAND
c       JCOL     :  Points to position in layer block
c       LDA      :  Row dimension of CBAND
c       NCD      :  Number of diagonals below or above main diagonal
c       NSHIFT   :  For positioning number of rows in band storage
c       WK       :  Temporary storage for EXP evaluations

c   Called by- DISORT, ALBTRN
c   Calls- ZEROIT
c +--------------------------------------------------------------------+


c     .. Scalar Arguments ..

      LOGICAL   LAMBER, LYRCUT
      INTEGER   MI, MI9M2, MXCMU, NCOL, NCUT, NN, NNLYRI, NSTR
      REAL      DELM0
c     ..
c     .. Array Arguments ..

      REAL      BDR( MI, 0:MI ), CBAND( MI9M2, NNLYRI ), CMU( MXCMU ),
     &          CWT( MXCMU ), DTAUCP( * ), GC( MXCMU, MXCMU, * ),
     &          KK( MXCMU, * ), TAUCPR( 0:* ), WK( MXCMU )
c     ..
c     .. Local Scalars ..

      INTEGER   IQ, IROW, JCOL, JQ, K, LC, LDA, NCD, NNCOL, NSHIFT
      REAL      EXPA, SUM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ZEROIT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC EXP
c     ..


      CALL ZEROIT( CBAND, MI9M2*NNLYRI )

      NCD    = 3*NN - 1
      LDA    = 3*NCD + 1
      NSHIFT = LDA - 2*NSTR + 1
      NCOL   = 0
c                         ** Use continuity conditions of Eq. STWJ(17)
c                            to form coefficient matrix in STWJ(20);
c                            employ scaling transformation STWJ(22)
      DO 60 LC = 1, NCUT

         DO 10 IQ = 1, NN
            WK( IQ ) = EXP( KK( IQ,LC )*DTAUCP( LC ) )
   10    CONTINUE

         JCOL  = 0

         DO 30 IQ = 1, NN

            NCOL  = NCOL + 1
            IROW  = NSHIFT - JCOL

            DO 20 JQ = 1, NSTR
               CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )
               CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )*WK( IQ )
               IROW  = IROW + 1
   20       CONTINUE

            JCOL  = JCOL + 1

   30    CONTINUE


         DO 50 IQ = NN + 1, NSTR

            NCOL  = NCOL + 1
            IROW  = NSHIFT - JCOL

            DO 40 JQ = 1, NSTR
               CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )*
     &                                          WK( NSTR + 1 - IQ )
               CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )
               IROW  = IROW + 1
   40       CONTINUE

            JCOL  = JCOL + 1

   50    CONTINUE

   60 CONTINUE
c                  ** Use top boundary condition of STWJ(20a) for
c                     first layer

      JCOL  = 0

      DO 80 IQ = 1, NN

         EXPA  = EXP( KK( IQ,1 )*TAUCPR( 1 ) )
         IROW  = NSHIFT - JCOL + NN

         DO 70 JQ = NN, 1, -1
            CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )*EXPA
            IROW  = IROW + 1
   70    CONTINUE

         JCOL  = JCOL + 1

   80 CONTINUE


      DO 100 IQ = NN + 1, NSTR

         IROW  = NSHIFT - JCOL + NN

         DO 90 JQ = NN, 1, -1
            CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )
            IROW  = IROW + 1
   90    CONTINUE

         JCOL  = JCOL + 1

  100 CONTINUE
c                           ** Use bottom boundary condition of
c                              STWJ(20c) for last layer

      NNCOL = NCOL - NSTR
      JCOL  = 0

      DO 130 IQ = 1, NN

         NNCOL  = NNCOL + 1
         IROW   = NSHIFT - JCOL + NSTR

         DO 120 JQ = NN + 1, NSTR

            IF( LYRCUT .OR. ( LAMBER .AND. DELM0.EQ.0 ) ) THEN

c                          ** No azimuthal-dependent intensity if Lam-
c                             bert surface; no intensity component if
c                             truncated bottom layer

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )

            ELSE

               SUM  = 0.0

               DO 110 K = 1, NN
                  SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )*
     &                     GC( NN + 1 - K, IQ, NCUT )
  110          CONTINUE

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT ) -
     &                                ( 1.+ DELM0 )*SUM
            END IF

            IROW  = IROW + 1

  120    CONTINUE

         JCOL  = JCOL + 1

  130 CONTINUE


      DO 160 IQ = NN + 1, NSTR

         NNCOL  = NNCOL + 1
         IROW   = NSHIFT - JCOL + NSTR
         EXPA   = WK( NSTR + 1 - IQ )

         DO 150 JQ = NN + 1, NSTR

            IF( LYRCUT .OR. ( LAMBER .AND. DELM0.EQ.0 ) ) THEN

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )*EXPA

            ELSE

               SUM  = 0.0

               DO 140 K = 1, NN
                  SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )*
     &                         GC( NN + 1 - K, IQ, NCUT )
  140          CONTINUE

               CBAND( IROW, NNCOL ) = ( GC( JQ,IQ,NCUT ) -
     &                                ( 1.+ DELM0 )*SUM )*EXPA
            END IF

            IROW  = IROW + 1

  150    CONTINUE

         JCOL  = JCOL + 1

  160 CONTINUE

      END


      SUBROUTINE SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MI, MAZIM,
     &                   MXCMU, NN, NSTR, YLMC, CC, EVECC, EVAL, KK, GC,
     &                   AAD, EVECCD, EVALD, WKD )

c         Solves eigenvalue/vector problem necessary to construct
c         homogeneous part of discrete ordinate solution; STWJ(8b)
c         ** NOTE ** Eigenvalue problem is degenerate when single
c                    scattering albedo = 1;  present way of doing it
c                    seems numerically more stable than alternative
c                    methods that we tried

c   I N P U T     V A R I A B L E S:

c       GL     :  Delta-M scaled Legendre coefficients of phase function
c                    (including factors 2l+1 and single-scatter albedo)
c       CMU    :  Computational polar angle cosines
c       CWT    :  Weights for quadrature over polar angle cosine
c       MAZIM  :  Order of azimuthal component
c       NN     :  Half the total number of streams
c       YLMC   :  Normalized associated Legendre polynomial
c                    at the quadrature angles CMU
c       (remainder are DISORT input variables)

c   O U T P U T    V A R I A B L E S:

c       CC     :  C-sub-ij in Eq. SS(5); needed in SS(15&18)
c       EVAL   :  NN eigenvalues of Eq. SS(12) on return from ASYMTX
c                    but then square roots taken
c       EVECC  :  NN eigenvectors  (G+) - (G-)  on return
c                    from ASYMTX ( column j corresponds to EVAL(j) )
c                    but then  (G+) + (G-)  is calculated from SS(10),
c                    G+  and  G-  are separated, and  G+  is stacked on
c                    top of  G-  to form NSTR eigenvectors of SS(7)
c       GC     :  Permanent storage for all NSTR eigenvectors, but
c                    in an order corresponding to KK
c       KK     :  Permanent storage for all NSTR eigenvalues of SS(7),
c                    but re-ordered with negative values first ( square
c                    roots of EVAL taken and negatives added )

c   I N T E R N A L   V A R I A B L E S:

c       AMB,APB :  Matrices (alpha-beta), (alpha+beta) in reduced
c                    eigenvalue problem
c       ARRAY   :  Complete coefficient matrix of reduced eigenvalue
c                    problem: (alfa+beta)*(alfa-beta)
c       GPPLGM  :  (G+) + (G-) (cf. Eqs. SS(10-11))
c       GPMIGM  :  (G+) - (G-) (cf. Eqs. SS(10-11))
c       WKD     :  Scratch array required by ASYMTX

c   Called by- DISORT, ALBTRN
c   Calls- ASYMTX, ERRMSG
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      INTEGER   MAZIM, MI, MXCMU, NN, NSTR
c     ..
c     .. Array Arguments ..

      REAL      AMB( MI, MI ), APB( MI, MI ), ARRAY( MI, * ),
     &          CC( MXCMU, MXCMU ), CMU( MXCMU ), CWT( MXCMU ),
     &          EVAL( MI ), EVECC( MXCMU, MXCMU ), GC( MXCMU, MXCMU ),
     &          GL( 0:MXCMU ), KK( MXCMU ), YLMC( 0:MXCMU, MXCMU )
      DOUBLE PRECISION AAD( MI, MI ), EVALD( MI ), EVECCD( MI, MI ),
     &                 WKD( MXCMU )
c     ..
c     .. Local Scalars ..

      INTEGER   IER, IQ, JQ, KQ, L
      REAL      ALPHA, BETA, GPMIGM, GPPLGM, SUM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ASYMTX, ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, SQRT
c     ..

c                             ** Calculate quantities in Eqs. SS(5-6)
      DO 40 IQ = 1, NN

         DO 20 JQ = 1, NSTR

            SUM  = 0.0
            DO 10 L = MAZIM, NSTR - 1
               SUM  = SUM + GL( L )*YLMC( L, IQ )*YLMC( L, JQ )
   10       CONTINUE

            CC( IQ, JQ ) = 0.5*SUM*CWT( JQ )

   20    CONTINUE

         DO 30 JQ = 1, NN
c                             ** Fill remainder of array using symmetry
c                                relations  C(-mui,muj) = C(mui,-muj)
c                                and        C(-mui,-muj) = C(mui,muj)

            CC( IQ + NN, JQ ) = CC( IQ, JQ + NN )
            CC( IQ + NN, JQ + NN ) = CC( IQ, JQ )

c                                       ** Get factors of coeff. matrix
c                                          of reduced eigenvalue problem

            ALPHA  = CC( IQ, JQ ) / CMU( IQ )
            BETA   = CC( IQ, JQ + NN ) / CMU( IQ )
            AMB( IQ, JQ ) = ALPHA - BETA
            APB( IQ, JQ ) = ALPHA + BETA

   30    CONTINUE

         AMB( IQ, IQ ) = AMB( IQ, IQ ) - 1.0 / CMU( IQ )
         APB( IQ, IQ ) = APB( IQ, IQ ) - 1.0 / CMU( IQ )

   40 CONTINUE
c                      ** Finish calculation of coefficient matrix of
c                         reduced eigenvalue problem:  get matrix
c                         product (alfa+beta)*(alfa-beta); SS(12)
      DO 70 IQ = 1, NN

         DO 60 JQ = 1, NN

            SUM  = 0.
            DO 50 KQ = 1, NN
               SUM  = SUM + APB( IQ, KQ )*AMB( KQ, JQ )
   50       CONTINUE

            ARRAY( IQ, JQ ) = SUM

   60    CONTINUE

   70 CONTINUE
c                      ** Find (real) eigenvalues and eigenvectors

      CALL ASYMTX( ARRAY, EVECC, EVAL, NN, MI, MXCMU, IER, WKD, AAD,
     &             EVECCD, EVALD )

      IF( IER.GT.0 ) THEN

         WRITE( *, FMT = '(//,A,I4,A)' ) ' ASYMTX--eigenvalue no. ',
     &      IER, '  didnt converge.  Lower-numbered eigenvalues wrong.'

         CALL ERRMSG( 'ASYMTX--convergence problems',.True.)

      END IF

CDIR$ IVDEP
      DO 80 IQ = 1, NN
         EVAL( IQ )    = SQRT( ABS( EVAL( IQ ) ) )
         KK( IQ + NN ) = EVAL( IQ )
c                                      ** Add negative eigenvalue
         KK( NN + 1 - IQ ) = -EVAL( IQ )
   80 CONTINUE

c                          ** Find eigenvectors (G+) + (G-) from SS(10)
c                             and store temporarily in APB array
      DO 110 JQ = 1, NN

         DO 100 IQ = 1, NN

            SUM  = 0.
            DO 90 KQ = 1, NN
               SUM  = SUM + AMB( IQ, KQ )*EVECC( KQ, JQ )
   90       CONTINUE

            APB( IQ, JQ ) = SUM / EVAL( JQ )

  100    CONTINUE

  110 CONTINUE


      DO 130 JQ = 1, NN
CDIR$ IVDEP
         DO 120 IQ = 1, NN

            GPPLGM = APB( IQ, JQ )
            GPMIGM = EVECC( IQ, JQ )
c                                ** Recover eigenvectors G+,G- from
c                                   their sum and difference; stack them
c                                   to get eigenvectors of full system
c                                   SS(7) (JQ = eigenvector number)

            EVECC( IQ,      JQ ) = 0.5*( GPPLGM + GPMIGM )
            EVECC( IQ + NN, JQ ) = 0.5*( GPPLGM - GPMIGM )

c                                ** Eigenvectors corresponding to
c                                   negative eigenvalues (corresp. to
c                                   reversing sign of 'k' in SS(10) )
            GPPLGM = - GPPLGM
            EVECC(IQ,   JQ+NN) = 0.5 * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,JQ+NN) = 0.5 * ( GPPLGM - GPMIGM )
            GC( IQ+NN,   JQ+NN )   = EVECC( IQ,    JQ )
            GC( NN+1-IQ, JQ+NN )   = EVECC( IQ+NN, JQ )
            GC( IQ+NN,   NN+1-JQ ) = EVECC( IQ,    JQ+NN )
            GC( NN+1-IQ, NN+1-JQ ) = EVECC( IQ+NN, JQ+NN )

  120    CONTINUE

  130 CONTINUE


      END

      SUBROUTINE SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,
     &                   FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM,
     &                   MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, NNLYRI,
     &                   PI, TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )

c        Construct right-hand side vector B for general boundary
c        conditions STWJ(17) and solve system of equations obtained
c        from the boundary conditions and the continuity-of-
c        intensity-at-layer-interface equations.
c        Thermal emission contributes only in azimuthal independence.

c     I N P U T      V A R I A B L E S:

c       BDR      :  Surface bidirectional reflectivity
c       BEM      :  Surface bidirectional emissivity
c       BPLANK   :  Bottom boundary thermal emission
c       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
c                   scaled by Eq. SC(12); in banded form required
c                   by LINPACK solution routines
c       CMU      :  Abscissae for Gauss quadrature over angle cosine
c       CWT      :  Weights for Gauss quadrature over angle cosine
c       EXPBEA   :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
c       LYRCUT   :  Logical flag for truncation of comput. layer
c       MAZIM    :  Order of azimuthal component
c       ncol     :  Counts of columns in CBAND
c       NN       :  Order of double-Gauss quadrature (NSTR/2)
c       NCUT     :  Total number of computational layers considered
c       TPLANK   :  Top boundary thermal emission
c       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
c       ZZ       :  Beam source vectors in Eq. SS(19)
c       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
c       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
c       (remainder are DISORT input variables)

c   O U T P U T     V A R I A B L E S:

c       B        :  Right-hand side vector of Eq. SC(5) going into
c                   SGBSL; returns as solution vector of Eq. SC(12),
c                   constants of integration without exponential term
c
c      LL        :  Permanent storage for B, but re-ordered

c   I N T E R N A L    V A R I A B L E S:

c       IPVT     :  Integer vector of pivot indices
c       IT       :  Pointer for position in  B
c       NCD      :  Number of diagonals below or above main diagonal
c       RCOND    :  Indicator of singularity for CBAND
c       Z        :  Scratch array required by SGBCO

c   Called by- DISORT
c   Calls- ZEROIT, SGBCO, ERRMSG, SGBSL
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      LOGICAL   LAMBER, LYRCUT
      INTEGER   MAZIM, MI, MI9M2, MXCMU, NCOL, NCUT, NN, NNLYRI, NSTR
      REAL      BPLANK, FBEAM, FISOT, PI, TPLANK, UMU0
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      B( NNLYRI ), BDR( MI, 0:MI ), BEM( MI ),
     &          CBAND( MI9M2, NNLYRI ), CMU( MXCMU ), CWT( MXCMU ),
     &          EXPBEA( 0:* ), LL( MXCMU, * ), TAUCPR( 0:* ),
     &          Z( NNLYRI ), ZPLK0( MXCMU, * ), ZPLK1( MXCMU, * ),
     &          ZZ( MXCMU, * )
c     ..
c     .. Local Scalars ..

      INTEGER   IPNT, IQ, IT, JQ, LC, NCD
      REAL      RCOND, SUM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, SGBCO, SGBSL, ZEROIT
c     ..


      CALL ZEROIT( B, NNLYRI )
c                              ** Construct B,  STWJ(20a,c) for
c                                 parallel beam + bottom reflection +
c                                 thermal emission at top and/or bottom

      IF( MAZIM.GT.0 .AND. FBEAM.GT.0.0 ) THEN

c                                         ** Azimuth-dependent case
c                                            (never called if FBEAM = 0)
         IF( LYRCUT .OR. LAMBER ) THEN

c               ** No azimuthal-dependent intensity for Lambert surface;
c                  no intensity component for truncated bottom layer

            DO 10 IQ = 1, NN
c                                                  ** Top boundary
               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 )
c                                                  ** Bottom boundary

               B( NCOL - NN + IQ ) = -ZZ( IQ + NN, NCUT )*EXPBEA( NCUT )

   10       CONTINUE


         ELSE

            DO 30 IQ = 1, NN

               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 )

               SUM  = 0.
               DO 20 JQ = 1, NN
                  SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ )*
     &                         ZZ( NN + 1 - JQ, NCUT )*EXPBEA( NCUT )
   20          CONTINUE

               B( NCOL - NN + IQ ) = SUM
               IF( FBEAM.GT.0.0 ) B( NCOL - NN + IQ ) = SUM +
     &             ( BDR( IQ,0 )*UMU0*FBEAM / PI - ZZ( IQ + NN,NCUT ) )*
     &             EXPBEA( NCUT )

   30       CONTINUE

         END IF
c                             ** Continuity condition for layer
c                                interfaces of Eq. STWJ(20b)
         IT   = NN

         DO 50 LC = 1, NCUT - 1

            DO 40 IQ = 1, NSTR
               IT   = IT + 1
               B( IT ) = ( ZZ( IQ, LC+1 ) - ZZ( IQ, LC ) )*EXPBEA( LC )
   40       CONTINUE

   50    CONTINUE


      ELSE
c                                   ** Azimuth-independent case

         IF( FBEAM.EQ.0.0 ) THEN

            DO 60 IQ = 1, NN
c                                      ** Top boundary

               B( IQ ) = -ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK

   60       CONTINUE


            IF( LYRCUT ) THEN
c                               ** No intensity component for truncated
c                                  bottom layer
               DO 70 IQ = 1, NN
c                                      ** Bottom boundary

                  B( NCOL - NN + IQ ) = - ZPLK0( IQ + NN, NCUT ) -
     &                                    ZPLK1( IQ + NN, NCUT )*
     &                                    TAUCPR( NCUT )
   70          CONTINUE


            ELSE

               DO 90 IQ = 1, NN

                  SUM  = 0.
                  DO 80 JQ = 1, NN
                     SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ )*
     &                            ( ZPLK0( NN + 1 - JQ,NCUT ) +
     &                        ZPLK1( NN + 1 - JQ,NCUT )*TAUCPR( NCUT ) )
   80             CONTINUE

                  B( NCOL - NN + IQ ) = 2.*SUM + BEM( IQ )*BPLANK -
     &                                  ZPLK0( IQ + NN, NCUT ) -
     &                                  ZPLK1( IQ + NN, NCUT )*
     &                                  TAUCPR( NCUT )
   90          CONTINUE

            END IF
c                             ** Continuity condition for layer
c                                interfaces, STWJ(20b)
            IT   = NN
            DO 110 LC = 1, NCUT - 1

               DO 100 IQ = 1, NSTR
                  IT   = IT + 1
                  B( IT ) =   ZPLK0( IQ, LC + 1 ) - ZPLK0( IQ, LC ) +
     &                      ( ZPLK1( IQ, LC + 1 ) - ZPLK1( IQ, LC ) )*
     &                      TAUCPR( LC )
  100          CONTINUE

  110       CONTINUE


         ELSE

            DO 120 IQ = 1, NN
               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 ) -
     &                   ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK
  120       CONTINUE

            IF( LYRCUT ) THEN

               DO 130 IQ = 1, NN
                  B(NCOL-NN+IQ) = - ZZ(IQ+NN, NCUT) * EXPBEA(NCUT)
     &                            - ZPLK0(IQ+NN, NCUT)
     &                            - ZPLK1(IQ+NN, NCUT) * TAUCPR(NCUT)
  130          CONTINUE


            ELSE

               DO 150 IQ = 1, NN

                  SUM  = 0.
                  DO 140 JQ = 1, NN
                     SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)
     &                          * ( ZZ(NN+1-JQ, NCUT) * EXPBEA(NCUT)
     &                            + ZPLK0(NN+1-JQ, NCUT)
     &                            + ZPLK1(NN+1-JQ, NCUT) * TAUCPR(NCUT))
  140             CONTINUE

                  B(NCOL-NN+IQ) = 2.*SUM + ( BDR(IQ,0) * UMU0*FBEAM/PI
     &                                - ZZ(IQ+NN, NCUT) ) * EXPBEA(NCUT)
     &                            + BEM(IQ) * BPLANK
     &                            - ZPLK0(IQ+NN, NCUT)
     &                            - ZPLK1(IQ+NN, NCUT) * TAUCPR(NCUT)
  150          CONTINUE

            END IF


            IT   = NN

            DO 170 LC = 1, NCUT - 1

               DO 160 IQ = 1, NSTR

                  IT   = IT + 1
                  B(IT) = ( ZZ(IQ,LC+1) - ZZ(IQ,LC) ) * EXPBEA(LC)
     &                    + ZPLK0(IQ,LC+1) - ZPLK0(IQ,LC) +
     &                    ( ZPLK1(IQ,LC+1) - ZPLK1(IQ,LC) ) * TAUCPR(LC)
  160          CONTINUE

  170       CONTINUE

         END IF

      END IF
c                     ** Find L-U (lower/upper triangular) decomposition
c                        of band matrix CBAND and test if it is nearly
c                        singular (note: CBAND is destroyed)
c                        (CBAND is in LINPACK packed format)
      RCOND  = 0.0
      NCD    = 3*NN - 1

      CALL SGBCO( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z )

      IF( 1.0 + RCOND.EQ.1.0 )
     &    CALL ERRMSG('SOLVE0--SGBCO says matrix near singular',.FALSE.)

c                   ** Solve linear system with coeff matrix CBAND
c                      and R.H. side(s) B after CBAND has been L-U
c                      decomposed.  Solution is returned in B.

      CALL SGBSL( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, B, 0 )

c                   ** Zero CBAND (it may contain 'foreign'
c                      elements upon returning from LINPACK);
c                      necessary to prevent errors

      CALL ZEROIT( CBAND, MI9M2*NNLYRI )

      DO 190 LC = 1, NCUT

         IPNT  = LC*NSTR - NN

         DO 180 IQ = 1, NN
            LL( NN + 1 - IQ, LC ) = B( IPNT + 1 - IQ )
            LL( IQ + NN,     LC ) = B( IQ + IPNT )
  180    CONTINUE

  190 CONTINUE

      RETURN
      END

      SUBROUTINE SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER, MI, MAZIM,
     &                   MXCMU, MXUMU, NN, NUMU, NSTR, ONLYFL, UMU,
     &                   USRANG, YLM0, YLMC, YLMU, BDR, EMU, BEM, RMU )

c       Specifies user's surface bidirectional properties, STWJ(21)

c   I N P U T     V A R I A B L E S:

c       DELM0  :  Kronecker delta, delta-sub-m0
c       HLPR   :  Legendre moments of surface bidirectional reflectivity
c                    (with 2K+1 factor included)
c       MAZIM  :  Order of azimuthal component
c       NN     :  Order of double-Gauss quadrature (NSTR/2)
c       YLM0   :  Normalized associated Legendre polynomial
c                 at the beam angle
c       YLMC   :  Normalized associated Legendre polynomials
c                 at the quadrature angles
c       YLMU   :  Normalized associated Legendre polynomials
c                 at the user angles
c       (remainder are DISORT input variables)

c    O U T P U T     V A R I A B L E S:

c       BDR :  Surface bidirectional reflectivity (computational angles)
c       RMU :  Surface bidirectional reflectivity (user angles)
c       BEM :  Surface directional emissivity (computational angles)
c       EMU :  Surface directional emissivity (user angles)

c    I N T E R N A L     V A R I A B L E S:

c       DREF      Directional reflectivity
c       NMUG   :  Number of angle cosine quadrature points on (0,1) for
c                   integrating bidirectional reflectivity to get
c                   directional emissivity (it is necessary to use a
c                   quadrature set distinct from the computational
c                   angles, because the computational angles may not be
c                   dense enough--NSTR may be too small--to give an
c                   accurate approximation for the integration).
c       GMU    :  The NMUG angle cosine quadrature points on (0,1)
c       GWT    :  The NMUG angle cosine quadrature weights on (0,1)
c       YLMG   :  Normalized associated Legendre polynomials
c                   at the NMUG quadrature angles

c   Called by- DISORT
c   Calls- QGAUSN, LEPOLY, ZEROIT, ERRMSG
c +-------------------------------------------------------------------+

c     .. Parameters ..

      INTEGER   NMUG, MAXSTR
      PARAMETER ( NMUG = 10, MAXSTR = 100 )
c     ..
c     .. Scalar Arguments ..

      LOGICAL   LAMBER, ONLYFL, USRANG
      INTEGER   MAZIM, MI, MXCMU, MXUMU, NN, NSTR, NUMU
      REAL      ALBEDO, DELM0, FBEAM
c     ..
c     .. Array Arguments ..

      REAL      BDR( MI, 0:MI ), BEM( MI ), EMU( MXUMU ),
     &          HLPR( 0:MXCMU ), RMU( MXUMU, 0:MI ), UMU( * ),
     &          YLM0( 0:MXCMU ), YLMC( 0:MXCMU, MXCMU ),
     &          YLMU( 0:MXCMU, MXUMU )
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      INTEGER   IQ, IU, JG, JQ, K
      REAL      DREF, SGN, SUM
c     ..
c     .. Local Arrays ..

      REAL      GMU( NMUG ), GWT( NMUG ), YLMG( 0:MAXSTR, NMUG )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, LEPOLY, QGAUSN, ZEROIT
c     ..
      SAVE      PASS1, GMU, GWT, YLMG
      DATA      PASS1 / .TRUE. /


      IF( PASS1 ) THEN

         PASS1  = .FALSE.

         CALL QGAUSN( NMUG, GMU, GWT )

         CALL LEPOLY( NMUG, 0, MAXSTR, MAXSTR, GMU, YLMG )

c                       ** Convert Legendre polys. to negative GMU
         SGN  = - 1.0

         DO 20 K = 0, MAXSTR

            SGN  = - SGN

            DO 10 JG = 1, NMUG
               YLMG( K, JG ) = SGN*YLMG( K, JG )
   10       CONTINUE

   20    CONTINUE

      END IF


      CALL ZEROIT( BDR, MI*( MI + 1 ) )
      CALL ZEROIT( BEM, MI )

      IF( LAMBER .AND. MAZIM.EQ.0 ) THEN

         DO 40 IQ = 1, NN

            BEM( IQ ) = 1.- ALBEDO

            DO 30 JQ = 0, NN
               BDR( IQ, JQ ) = ALBEDO
   30       CONTINUE

   40    CONTINUE


      ELSE IF( .NOT.LAMBER ) THEN
c                                  ** Compute surface bidirectional
c                                     properties at computational angles
         DO 80 IQ = 1, NN

            DO 60 JQ = 1, NN

               SUM  = 0.0
               DO 50 K = MAZIM, NSTR - 1
                  SUM  = SUM + HLPR( K )*YLMC( K, IQ )*
     &                         YLMC( K, JQ + NN )
   50          CONTINUE

               BDR( IQ, JQ ) = ( 2.- DELM0 )*SUM

   60       CONTINUE


            IF( FBEAM.GT.0.0 ) THEN

               SUM  = 0.0
               DO 70 K = MAZIM, NSTR - 1
                  SUM  = SUM + HLPR( K )*YLMC( K, IQ )*YLM0( K )
   70          CONTINUE

               BDR( IQ, 0 ) = ( 2.- DELM0 )*SUM

            END IF

   80    CONTINUE


         IF( MAZIM.EQ.0 ) THEN

            IF( NSTR.GT.MAXSTR )
     &          CALL ERRMSG('SURFAC--parameter MAXSTR too small',.True.)

c                              ** Integrate bidirectional reflectivity
c                                 at reflection polar angles CMU and
c                                 incident angles GMU to get
c                                 directional emissivity at
c                                 computational angles CMU.
            DO 110 IQ = 1, NN

               DREF  = 0.0

               DO 100 JG = 1, NMUG

                  SUM  = 0.0
                  DO 90 K = 0, NSTR - 1
                     SUM  = SUM + HLPR( K )*YLMC( K, IQ )*
     &                            YLMG( K, JG )
   90             CONTINUE

                  DREF  = DREF + 2.*GWT( JG )*GMU( JG )*SUM

  100          CONTINUE

               BEM( IQ ) = 1.- DREF

  110       CONTINUE

         END IF

      END IF
c                                       ** Compute surface bidirectional
c                                          properties at user angles

      IF( .NOT.ONLYFL .AND. USRANG ) THEN

         CALL ZEROIT( EMU, MXUMU )
         CALL ZEROIT( RMU, MXUMU*( MI + 1 ) )

         DO 180 IU = 1, NUMU

            IF( UMU( IU ).GT.0.0 ) THEN

               IF( LAMBER .AND. MAZIM.EQ.0 ) THEN

                  DO 120 IQ = 0, NN
                     RMU( IU, IQ ) = ALBEDO
  120             CONTINUE

                  EMU( IU ) = 1.- ALBEDO


               ELSE IF( .NOT.LAMBER ) THEN

                  DO 140 IQ = 1, NN

                     SUM  = 0.0
                     DO 130 K = MAZIM, NSTR - 1
                        SUM  = SUM + HLPR( K )*YLMU( K, IU )*
     &                               YLMC( K, IQ + NN )
  130                CONTINUE

                     RMU( IU, IQ ) = ( 2.- DELM0 )*SUM

  140             CONTINUE


                  IF( FBEAM.GT.0.0 ) THEN

                     SUM  = 0.0
                     DO 150 K = MAZIM, NSTR - 1
                        SUM  = SUM + HLPR( K )*YLMU( K, IU )*YLM0( K )
  150                CONTINUE

                     RMU( IU, 0 ) = ( 2.- DELM0 )*SUM

                  END IF


                  IF( MAZIM.EQ.0 ) THEN

c                               ** Integrate bidirectional reflectivity
c                                  at reflection angles UMU and
c                                  incident angles GMU to get
c                                  directional emissivity at
c                                  user angles UMU.
                     DREF  = 0.0

                     DO 170 JG = 1, NMUG

                        SUM  = 0.0
                        DO 160 K = 0, NSTR - 1
                           SUM  = SUM + HLPR( K )*YLMU( K, IU )*
     &                                  YLMG( K, JG )
  160                   CONTINUE

                        DREF  = DREF + 2.*GWT( JG )*GMU( JG )*SUM

  170                CONTINUE

                     EMU( IU ) = 1.- DREF

                  END IF

               END IF

            END IF

  180    CONTINUE

      END IF

      END


*bm  SOLVEC calls SOLEIG and UPBEAM; if UPBEAM reports a potenially 
*bm  unstable solution, the calculation is repeated with a slightly 
*bm  changed single scattering albedo; this process is iterates 
*bm  until a stable solution is found; as stable solutions may be 
*bm  reached either by increasing or by decreasing the single 
*bm  scattering albedo, both directions are explored ('upward' and
*bm  'downward' iteration); the solution which required the smaller 
*bm  change in the single scattering albedo is finally returned 
*bm  by SOLVEC.

      SUBROUTINE SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL, MI,
     &     MAZIM, MXCMU, NN, NSTR, YLM0, YLMC, CC, 
     &     EVECC, EVAL, KK, GC, AAD, EVECCD, EVALD,
     &     WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0, ZJ, ZZ,
     &     OPRIM, LC, DITHER, mu2, glsave, dgl)

cgy added glsave and dgl to call to allow adjustable dimensioning


c     .. Scalar Arguments ..

      INTEGER   MAZIM, MI, MXCMU, NN, NSTR, LC
      REAL      DELM0, FBEAM, PI, UMU0, OPRIM, DITHER
      REAL      mu2

c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      
      REAL      AMB( MI, MI ), APB( MI, MI ), ARRAY( MI, * ),
     &     CC( MXCMU, MXCMU ), CMU( MXCMU ), CWT( MXCMU ),
     &     EVAL( MI ), EVECC( MXCMU, MXCMU ), GC( MXCMU, MXCMU ),
     &     GL( 0:MXCMU ), KK( MXCMU ), 
     &     YLM0( 0:MXCMU ), YLMC( 0:MXCMU, MXCMU ),
     &     WK( MXCMU ), ZJ( MXCMU ), ZZ( MXCMU )

      DOUBLE PRECISION AAD( MI, MI ), EVALD( MI ), EVECCD( MI, MI ),
     &                 WKD( MXCMU )

*bm   Variables for instability fix
      
      INTEGER UAGAIN, DAGAIN
      REAL MINRCOND, ADD, UADD, DADD, SSA, DSSA, FACTOR
      REAL GLSAVE( 0:MXCMU ), DGL( 0:MXCMU )
      
      LOGICAL  DONE, NOUP, NODN, DEBUG, INSTAB
      
*bm   reset parameters

      DONE = .FALSE.
      NOUP = .FALSE.
      NODN = .FALSE.


*bm   flag for printing debugging output      
*      DEBUG  = .TRUE.
      DEBUG  = .FALSE.

*bm   instability parameter; the solution is considered 
*bm   unstable, if the RCOND reported by SGECO is smaller 
*bm   than MINRCOND
      MINRCOND = 5000. * R1MACH(4)

*bm   if an instability is detected, the single scattering albedo
*bm   is iterated downwards in steps of DADD and upwards in steps 
*bm   of UADD; in practice, MINRCOND and -MINRCOND should 
*bm   be reasonable choices for these parameters
      DADD    = -MINRCOND
      UADD    = MINRCOND

      UAGAIN = 0
      DAGAIN = 0
      ADD   = DADD
      

*bm   save array GL( ) because it will be 
*bm   changed if an iteration should be neccessary
      DO K = MAZIM, NSTR - 1
         GLSAVE( K ) =  GL( K )
      ENDDO
      
      SSA = OPRIM


*bm   in case of an instability reported by UPBEAM (INSTAB)
*bm   the single scattering albedo will be changed by a small 
*bm   amount (ADD); this is indicated by DAGAIN or UAGAIN 
*bm   being larger than 0; a change in the single scattering 
*bm   albedo is equivalent to scaling the array GL( )

 666  IF ( DAGAIN .GT. 0 .OR. UAGAIN .GT. 0)  THEN
         FACTOR = (SSA + ADD) / SSA
         DO K = MAZIM, NSTR - 1
            GL( K ) =  GL( K ) * FACTOR
         ENDDO

         SSA = SSA + ADD
         
*bm   if the single scattering albedo is now smaller than 0
*bm   the downward iteration is stopped and upward iteration 
*bm   is forced instead

         IF( SSA .LT. DITHER) THEN
            NODN = .TRUE.
            DAGAIN = -1
            goto 778
         ENDIF

*bm   if the single scattering albedo is now larger than its maximum 
*bm   allowed value (1.0 - DITHER), the upward iteration is 
*bm   stopped and downward iteration is forced instead

         IF( SSA .GT. 1.0 - DITHER) THEN
            NOUP = .TRUE.
            UAGAIN = -1
            goto 888
         ENDIF
      ENDIF


c     ** Solve eigenfunction problem in Eq. STWJ(8B);
c        return eigenvalues and eigenvectors

 777     CALL SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MI,
     &     MAZIM, MXCMU, NN, NSTR, YLMC, CC, EVECC, EVAL,
     &     KK, GC, AAD, EVECCD, EVALD,
     &     WKD )

c     ** Calculate particular solutions of
c        q.SS(18) for incident beam source

      IF ( FBEAM.GT.0.0 ) THEN
         CALL  UPBEAM( mu2,
     $        ARRAY, CC, CMU, DELM0, FBEAM, GL,
     $        IPVT, MAZIM, MXCMU, NN, NSTR, PI, UMU0, WK,
     $        YLM0, YLMC, ZJ, ZZ, MINRCOND, INSTAB)
      ENDIF
      
c     ** Calculate particular solutions of
c        Eq. SS(15) for thermal emission source
c        (not available in psndo.f)
      
*bm   finished if the result is stable on the first try
      IF ( (.NOT. INSTAB) .AND. 
     $     (UAGAIN .EQ. 0) .AND. (DAGAIN .EQ. 0)) THEN
         goto 999
      ENDIF

*bm   downward iteration
      IF( INSTAB .AND. UAGAIN .EQ. 0 )  THEN
         DAGAIN = DAGAIN + 1
         GOTO 666
      ENDIF
      
*bm   upward iteration
      IF( INSTAB .AND. UAGAIN .GT. 0 )  THEN
         UAGAIN = UAGAIN + 1
         GOTO 666
      ENDIF


*bm   ( DAGAIN .NE. 0 ) at this place means that the downward
*bm   iteration is finished 

 778  IF (DAGAIN .NE. 0 .AND. UAGAIN .EQ. 0) THEN
         
*bm   save downward iteration data for later use and 
*bm   restore original input data
         DO K = MAZIM, NSTR - 1
            DGL( K ) =  GL( K )
            GL( K ) =  GLSAVE( K )
         ENDDO

         DSSA = SSA
         SSA = OPRIM

*bm   start upward iteration
         ADD = UADD
         UAGAIN = UAGAIN + 1
         GOTO 666
      ENDIF

*bm   both iterations finished
 888  IF (DONE) THEN
         goto 998
      ENDIF


*bm  if neither upward nor downward iteration converged, the 
*bm  original conditions are restored and SOLEIG/UPBEAM 
*bm  is called for the last time 
         
      IF (NOUP .AND. NODN) THEN
         
         DO K = MAZIM, NSTR - 1
            GL( K ) =  GLSAVE( K )
         ENDDO
         
         SSA = OPRIM
         
         IF (DEBUG) THEN
            write (*,*) '! *** Neither upward nor downward iteration'
            write (*,*) '! *** converged; using original result.'
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ENDIF

*bm  if upward iteration did not converge, the stable downward conditions
*bm  are restored and SOLEIG/UPBEAM is called for the last time
      IF (NOUP) THEN
         DO K = MAZIM, NSTR - 1
            GL( K ) =  DGL( K )
         ENDDO
         
         SSA = DSSA
         
         IF (DEBUG) THEN
            write (*,*) '! *** The upward iteration did not converge.'
            write (*,*) '! *** Had to iterate ', DAGAIN,
     $           ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ENDIF

*bm  if downward iteration did not converge, we are done 
*bm  (the result of the upward iteration will be used)
      IF (NODN) THEN
         IF (DEBUG) THEN
            write (*,*) '! *** The downward iteration did not converge.'
            write (*,*) '! *** Had to iterate ', UAGAIN,
     $           ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF
         
         DONE = .TRUE.
         GOTO 998
      ENDIF

      
*bm   if both iterations converged, and if the upward iteration 
*bm   required more steps than the downward iteration, the stable 
*bm   downward conditions are restored and SOLEIG/UPBEAM is 
*bm   called for the last time 
         
      IF (UAGAIN .GT. DAGAIN) THEN
         DO K = MAZIM, NSTR - 1
            GL( K ) =  DGL( K )
         ENDDO
         
         SSA = DSSA
         
         IF (DEBUG) THEN
            write (*,*) '! *** Both iterations converged;',
     $           ' using downward.'
            write (*,*) '! *** Had to iterate ', DAGAIN,
     $        ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ELSE
         
         IF (DEBUG) THEN
            write (*,*) '! *** Both iterations converged;',
     $           ' using upward.'
            write (*,*) '! *** Had to iterate ', UAGAIN,
     $        ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         goto 998
      ENDIF
      
*bm   finally restore original input data
 998  DO K = MAZIM, NSTR - 1
         GL( K ) =  GLSAVE( K )
      ENDDO
      
 999  CONTINUE
      END



      SUBROUTINE UPBEAM( mu2,
     &                   ARRAY, CC, CMU, DELM0, FBEAM, GL, IPVT, MAZIM,
     &                   MXCMU, NN, NSTR, PI, UMU0, WK, YLM0, YLMC, ZJ,
     &                   ZZ, MINRCOND, INSTAB )

c         Finds the incident-beam particular solution of SS(18)

c   I N P U T    V A R I A B L E S:

c       CC     :  C-sub-ij in Eq. SS(5)
c       CMU    :  Abscissae for Gauss quadrature over angle cosine
c       DELM0  :  Kronecker delta, delta-sub-m0
c       GL     :  Delta-M scaled Legendre coefficients of phase function
c                    (including factors 2L+1 and single-scatter albedo)
c       MAZIM  :  Order of azimuthal component
c       YLM0   :  Normalized associated Legendre polynomial
c                    at the beam angle
c       YLMC   :  Normalized associated Legendre polynomial
c                    at the quadrature angles
c       (remainder are DISORT input variables)

c   O U T P U T    V A R I A B L E S:

c       ZJ     :  Right-hand side vector X-sub-zero in SS(19); also the
c                 solution vector Z-sub-zero after solving that system

c       ZZ     :  Permanent storage for ZJ, but re-ordered

c   I N T E R N A L    V A R I A B L E S:

c       ARRAY  :  Coefficient matrix in left-hand side of Eq. SS(19)
c       IPVT   :  Integer vector of pivot indices required by LINPACK
c       WK     :  Scratch array required by LINPACK

c   Called by- DISORT
c   Calls- SGECO, ERRMSG, SGESL
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      INTEGER   MAZIM, MXCMU, NN, NSTR
      LOGICAL   INSTAB
      REAL      MINRCOND
      REAL      DELM0, FBEAM, PI, UMU0
      REAL mu2
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ARRAY( MXCMU, MXCMU ), CC( MXCMU, MXCMU ), CMU( MXCMU ),
     &          GL( 0:MXCMU ), WK( MXCMU ), YLM0( 0:MXCMU ),
     &          YLMC( 0:MXCMU, * ), ZJ( MXCMU ), ZZ( MXCMU )
c     ..
c     .. Local Scalars ..

      INTEGER   IQ, JOB, JQ, K
      REAL      RCOND, SUM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, SGECO, SGESL
c     ..


      DO 30 IQ = 1, NSTR

         DO 10 JQ = 1, NSTR
            ARRAY( IQ, JQ ) = -CC( IQ, JQ )
   10    CONTINUE

         ARRAY( IQ, IQ ) = 1.+ CMU( IQ ) / mu2 + ARRAY( IQ, IQ )

         SUM  = 0.
         DO 20 K = MAZIM, NSTR - 1
            SUM  = SUM + GL( K )*YLMC( K, IQ )*YLM0( K )
   20    CONTINUE

         ZJ( IQ ) = ( 2.- DELM0 )*FBEAM*SUM / ( 4.*PI )
   30 CONTINUE

c                  ** Find L-U (lower/upper triangular) decomposition
c                     of ARRAY and see if it is nearly singular
c                     (NOTE:  ARRAY is destroyed)
      RCOND  = 0.0

      CALL SGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )

*bm      IF( 1.0 + RCOND.EQ.1.0 )
*bm     &    CALL ERRMSG('UPBEAM--SGECO says matrix near singular',.FALSE.)
*bm
*bm   replaced original check of RCOND by the following:

      INSTAB = .FALSE.
      IF( ABS(RCOND) .LT. MINRCOND )  THEN
         INSTAB = .TRUE.
         RETURN
      ENDIF

c                ** Solve linear system with coeff matrix ARRAY
c                   (assumed already L-U decomposed) and R.H. side(s)
c                   ZJ;  return solution(s) in ZJ
      JOB  = 0

      CALL SGESL( ARRAY, MXCMU, NSTR, IPVT, ZJ, JOB )

CDIR$ IVDEP
      DO 40 IQ = 1, NN
         ZZ( IQ + NN )     = ZJ( IQ )
         ZZ( NN + 1 - IQ ) = ZJ( IQ + NN )
   40 CONTINUE

      END


      SUBROUTINE ZEROAL( ND1, EXPBEA, FLYR, OPRIM, TAUCPR, XR0, XR1,
     &                    ND2, CMU, CWT, PSI, WK, Z0, Z1, ZJ,
     &                    ND3, HLPR, YLM0,
     &                    ND4, ARRAY, CC, EVECC,
     &                    ND5, GL,
     &                    ND6, YLMC,
     &                    ND7, YLMU,
     &                    ND8, KK, LL, ZZ, ZPLK0, ZPLK1,
     &                    ND9, GC,
     &                    ND10, LAYRU, UTAUPR,
     &                    ND11, GU,
     &                    ND12, Z0U, Z1U, ZBEAM,
     &                    ND13, EVAL,
     &                    ND14, AMB, APB,
     &                    ND15, IPVT, Z,
     &                    ND16, RFLDIR, RFLDN, FLUP, UAVG, DFDT,
     &                    ND17, ALBMED, TRNMED,
     &                    ND18, U0U,
     &                    ND19, UU )

c         ZERO ARRAYS; NDn is dimension of all arrays following
c         it in the argument list

c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   ND1, ND10, ND11, ND12, ND13, ND14, ND15, ND16, ND17,
     &          ND18, ND19, ND2, ND3, ND4, ND5, ND6, ND7, ND8, ND9
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * ), LAYRU( * )
      REAL      ALBMED( * ), AMB( * ), APB( * ), ARRAY( * ), CC( * ),
     &          CMU( * ), CWT( * ), DFDT( * ), EVAL( * ), EVECC( * ),
     &          EXPBEA( * ), FLUP( * ), FLYR( * ), GC( * ), GL( * ),
     &          GU( * ), HLPR( * ), KK( * ), LL( * ), OPRIM( * ),
     &          PSI( * ), RFLDIR( * ), RFLDN( * ), TAUCPR( * ),
     &          TRNMED( * ), U0U( * ), UAVG( * ), UTAUPR( * ), UU( * ),
     &          WK( * ), XR0( * ), XR1( * ), YLM0( * ), YLMC( * ),
     &          YLMU( * ), Z( * ), Z0( * ), Z0U( * ), Z1( * ), Z1U( * ),
     &          ZBEAM( * ), ZJ( * ), ZPLK0( * ), ZPLK1( * ), ZZ( * )
c     ..
c     .. Local Scalars ..

      INTEGER   N
c     ..


      DO 10 N = 1, ND1
         EXPBEA( N ) = 0.0
         FLYR( N )   = 0.0
         OPRIM( N )  = 0.0
         TAUCPR( N ) = 0.0
         XR0( N )    = 0.0
         XR1( N )    = 0.0
   10 CONTINUE

      DO 20 N = 1, ND2
         CMU( N ) = 0.0
         CWT( N ) = 0.0
         PSI( N ) = 0.0
         WK( N )  = 0.0
         Z0( N )  = 0.0
         Z1( N )  = 0.0
         ZJ( N )  = 0.0
   20 CONTINUE

      DO 30 N = 1, ND3
         HLPR( N ) = 0.0
         YLM0( N ) = 0.0
   30 CONTINUE

      DO 40 N = 1, ND4
         ARRAY( N ) = 0.0
         CC( N )    = 0.0
         EVECC( N ) = 0.0
   40 CONTINUE

      DO 50 N = 1, ND5
         GL( N ) = 0.0
   50 CONTINUE

      DO 60 N = 1, ND6
         YLMC( N ) = 0.0
   60 CONTINUE

      DO 70 N = 1, ND7
         YLMU( N ) = 0.0
   70 CONTINUE

      DO 80 N = 1, ND8
         KK( N )    = 0.0
         LL( N )    = 0.0
         ZZ( N )    = 0.0
         ZPLK0( N ) = 0.0
         ZPLK1( N ) = 0.0
   80 CONTINUE

      DO 90 N = 1, ND9
         GC( N ) = 0.0
   90 CONTINUE

      DO 100 N = 1, ND10
         LAYRU( N )  = 0
         UTAUPR( N ) = 0.0
  100 CONTINUE

      DO 110 N = 1, ND11
         GU( N ) = 0.0
  110 CONTINUE

      DO 120 N = 1, ND12
         Z0U( N )   = 0.0
         Z1U( N )   = 0.0
         ZBEAM( N ) = 0.0
  120 CONTINUE

      DO 130 N = 1, ND13
         EVAL( N ) = 0.0
  130 CONTINUE

      DO 140 N = 1, ND14
         AMB( N ) = 0.0
         APB( N ) = 0.0
  140 CONTINUE

      DO 150 N = 1, ND15
         IPVT( N ) = 0
         Z( N )    = 0.0
  150 CONTINUE

      DO 160 N = 1, ND16
         RFLDIR( N ) = 0.
         RFLDN( N )  = 0.
         FLUP( N )   = 0.
         UAVG( N )   = 0.
         DFDT( N )   = 0.
  160 CONTINUE

      DO 170 N = 1, ND17
         ALBMED( N ) = 0.
         TRNMED( N ) = 0.
  170 CONTINUE

      DO 180 N = 1, ND18
         U0U( N ) = 0.
  180 CONTINUE

      DO 190 N = 1, ND19
         UU( N ) = 0.
  190 CONTINUE


      END

      SUBROUTINE ZEROIT( A, LENGTH )

c         Zeros a real array A having LENGTH elements
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   LENGTH
c     ..
c     .. Array Arguments ..

      REAL      A( LENGTH )
c     ..
c     .. Local Scalars ..

      INTEGER   L
c     ..

      DO 10 L = 1, LENGTH
         A( L ) = 0.0
   10 CONTINUE

      END

      REAL FUNCTION DREF( MU, HL, NSTR )

c        Exact flux albedo for given angle of incidence, given
c        a bidirectional reflectivity characterized by its
c        Legendre coefficients ( NOTE** these will only agree
c        with bottom-boundary albedos calculated by DISORT in
c        the limit as number of streams go to infinity, because
c        DISORT evaluates the integral 'CL' only approximately,
c        by quadrature, while this routine calculates it exactly.)

c  INPUT :   MU     Cosine of incidence angle
c            HL     Legendre coefficients of bidirectional reflectivity
c          NSTR     Number of elements of HL to consider

c  INTERNAL VARIABLES (P-sub-L is the L-th Legendre polynomial) :

c       CL      Integral from 0 to 1 of  MU * P-sub-L(MU)
c                   (vanishes for  L = 3, 5, 7, ... )
c       PL      P-sub-L
c       PLM1    P-sub-(L-1)
c       PLM2    P-sub-(L-2)

c   Called by- CHEKIN
c   Calls- ERRMSG
c +-------------------------------------------------------------------+

c     .. Parameters ..

      INTEGER   MAXTRM
      PARAMETER ( MAXTRM = 100 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   NSTR
      REAL      MU
c     ..
c     .. Array Arguments ..

      REAL      HL( 0:NSTR )
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      INTEGER   L
      REAL      CL, PL, PLM1, PLM2
c     ..
c     .. Local Arrays ..

      REAL      C( MAXTRM )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..
      SAVE      PASS1, C
      DATA      PASS1 / .TRUE. /
c     ..


      IF( PASS1 ) THEN

         PASS1  = .FALSE.
         CL     = 0.125
         C( 2 ) = 10.*CL

         DO 10 L = 4, MAXTRM, 2
            CL     = - CL*( L - 3 ) / ( L + 2 )
            C( L ) = 2.*( 2*L + 1 )*CL
   10    CONTINUE

      END IF


      IF( NSTR.LT.2 .OR. ABS(MU).GT.1.0 )
     &    CALL ERRMSG( 'DREF--input argument error(s)',.True. )

      IF( NSTR.GT.MAXTRM )
     &    CALL ERRMSG( 'DREF--parameter MAXTRM too small',.True. )


      DREF  = HL( 0 ) - 2.*HL( 1 )*MU
      PLM2  = 1.0
      PLM1  = - MU

      DO 20 L = 2, NSTR - 1
c                                ** Legendre polynomial recurrence

         PL = ( ( 2*L - 1 )*( -MU )*PLM1 - ( L-1 )*PLM2 ) / L

         IF( MOD( L,2 ).EQ.0 ) DREF   = DREF + C( L )*HL( L )*PL

         PLM2  = PLM1
         PLM1  = PL

   20 CONTINUE

      IF( DREF.LT.0.0 .OR. DREF.GT.1.0 )
     &    CALL ERRMSG( 'DREF--albedo value not in (0,1)',.False. )

      END

      REAL FUNCTION RATIO( A, B )

c        Calculate ratio  A/B  with over- and under-flow protection
c        (thanks to Prof. Jeff Dozier for some suggestions here).
c        Since this routine takes two logs, it is no speed demon,
c        but it is invaluable for comparing results from two runs
c        of a program under development.

c        NOTE:  In Fortran90, built-in functions TINY and HUGE
c               can replace the R1MACH calls.
c ---------------------------------------------------------------

c     .. Scalar Arguments ..

      REAL      A, B
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      REAL      ABSA, ABSB, HUGE, POWA, POWB, POWMAX, POWMIN, TINY
c     ..
c     .. External Functions ..

      REAL      R1MACH
      EXTERNAL  R1MACH
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, LOG10, SIGN
c     ..
      SAVE      PASS1, TINY, HUGE, POWMAX, POWMIN
      DATA      PASS1 / .TRUE. /
c     ..


      IF( PASS1 ) THEN

         TINY   = R1MACH( 1 )
         HUGE   = R1MACH( 2 )
         POWMAX = LOG10( HUGE )
         POWMIN = LOG10( TINY )
         PASS1  = .FALSE.

      END IF


      IF( A.EQ.0.0 ) THEN

         IF( B.EQ.0.0 ) THEN

            RATIO  = 1.0

         ELSE

            RATIO  = 0.0

         END IF


      ELSE IF( B.EQ.0.0 ) THEN

         RATIO  = SIGN( HUGE, A )

      ELSE

         ABSA   = ABS( A )
         ABSB   = ABS( B )
         POWA   = LOG10( ABSA )
         POWB   = LOG10( ABSB )

         IF( ABSA.LT.TINY .AND. ABSB.LT.TINY ) THEN

            RATIO  = 1.0

         ELSE IF( POWA - POWB.GE.POWMAX ) THEN

            RATIO  = HUGE

         ELSE IF( POWA - POWB.LE.POWMIN ) THEN

            RATIO  = TINY

         ELSE

            RATIO  = ABSA / ABSB

         END IF
c                      ** DONT use old trick of determining sign
c                      ** from A*B because A*B may (over/under)flow

         IF( ( A.GT.0.0 .AND. B.LT.0.0 ) .OR.
     &       ( A.LT.0.0 .AND. B.GT.0.0 ) ) RATIO = -RATIO

      END IF

      END
      SUBROUTINE  ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error
c        after making symbolic dump (machine-specific)

      LOGICAL       FATAL, MsgLim, Cray
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /


      IF ( FATAL )  THEN
         WRITE ( *, '(//,2A,//)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
      ELSE
         WRITE ( *,99 )
         MsgLim = .True.
      ENDIF

      RETURN

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',
     $   'They will no longer be printed  <<<<<<<', // )
      END

      LOGICAL FUNCTION  WrtBad ( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

      CHARACTER*(*)  VarNam
      INTEGER        MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 50 /


      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( *, '(3A)' )  ' ****  Input variable  ', VarNam,
     $                     '  in error  ****'
      IF ( NumMsg.EQ.MaxMsg )
     $   CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )
      RETURN
      END

      LOGICAL FUNCTION  WrtDim ( DimNam, MinVal )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

      CHARACTER*(*)  DimNam
      INTEGER        MinVal


      WRITE ( *, '(3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,
     $                     '  should be increased to at least ', MinVal
      WrtDim = .TRUE.
      RETURN
      END

      LOGICAL FUNCTION  TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

      CHARACTER*(*)  VarNam
      REAL           RelErr


      TstBad = .FALSE.
      WRITE( *, '(/,3A,1P,E11.2,A)' )
     $       ' Output variable ', VarNam,' differed by ', 100.*RelErr,
     $       ' per cent from correct value.  Self-test failed.'
      RETURN
      END
	SUBROUTINE  SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )

C         FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION 
C         AND ESTIMATES THE CONDITION OF THE MATRIX.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C     IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW SBGCO BY SGBSL.

C     INPUT:

C        ABD     REAL(LDA, N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.

C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. 2*ML + MU + 1 .

C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.

C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0 .LE. ML .LT. N .

C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. MU .LT. N .
C                MORE EFFICIENT IF  ML .LE. MU .

C     ON RETURN

C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.

C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.

C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

C     BAND STORAGE

C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.

C                   ML = (BAND WIDTH BELOW THE DIAGONAL)
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE

C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.

C     EXAMPLE:  IF THE ORIGINAL MATRIX IS

C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66

C      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN

C            *  *  *  +  +  +  , * = NOT USED
C            *  * 13 24 35 46  , + = USED FOR PIVOTING
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C           21 32 43 54 65  *


C     ROUTINES CALLED:  FROM LINPACK: SGBFA
C                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
C                       FROM FORTRAN: ABS, AMAX1, MAX0, MIN0, SIGN

	INTEGER  LDA, N, ML, MU, IPVT(*)
	REAL     ABD(LDA,*), Z(*)
	REAL     RCOND

	REAL     SDOT, EK, T, WK, WKM
	REAL     ANORM, S, SASUM, SM, YNORM
	INTEGER  IS, INFO, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM


C                       ** COMPUTE 1-NORM OF A
	ANORM = 0.0E0
	L = ML + 1
	IS = L + MU
	DO 10 J = 1, N
	   ANORM = AMAX1(ANORM, SASUM(L,ABD(IS,J), 1))
	   IF (IS .GT. ML + 1) IS = IS - 1
	   IF (J .LE. MU) L = L + 1
	   IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
C                                               ** FACTOR
	CALL SGBFA(ABD, LDA, N, ML, MU, IPVT, INFO)

C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.

C                     ** SOLVE TRANS(U)*W = E
	EK = 1.0E0
	DO 20 J = 1, N
	   Z(J) = 0.0E0
   20 CONTINUE

	M = ML + MU + 1
	JU = 0
	DO 100 K = 1, N
	   IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
	   IF (ABS(EK-Z(K)) .GT. ABS(ABD(M,K))) THEN
	      S = ABS(ABD(M,K))/ABS(EK-Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      EK = S*EK
	   ENDIF
	   WK = EK - Z(K)
	   WKM = -EK - Z(K)
	   S = ABS(WK)
	   SM = ABS(WKM)
	   IF (ABD(M,K) .NE. 0.0E0) THEN
	      WK  = WK /ABD(M,K)
	      WKM = WKM/ABD(M,K)
	   ELSE
	      WK  = 1.0E0
	      WKM = 1.0E0
	   ENDIF
	   KP1 = K + 1
	   JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
	   MM = M
	   IF (KP1 .LE. JU) THEN
	      DO 60 J = KP1, JU
	         MM = MM - 1
	         SM = SM + ABS(Z(J)+WKM*ABD(MM,J))
	         Z(J) = Z(J) + WK*ABD(MM,J)
	         S = S + ABS(Z(J))
   60       CONTINUE
	      IF (S .LT. SM) THEN
	         T = WKM - WK
	         WK = WKM
	         MM = M
	         DO 70 J = KP1, JU
	            MM = MM - 1
	            Z(J) = Z(J) + T*ABD(MM,J)
   70          CONTINUE
	      ENDIF
	   ENDIF
	   Z(K) = WK
  100 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)

C                         ** SOLVE TRANS(L)*Y = W
	DO 120 KB = 1, N
	   K = N + 1 - KB
	   LM = MIN0(ML, N-K)
	   IF (K .LT. N) Z(K) = Z(K) + SDOT(LM, ABD(M+1,K), 1, Z(K+1), 1)
	   IF (ABS(Z(K)) .GT. 1.0E0) THEN
	      S = 1.0E0 / ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	   ENDIF
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
  120 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)

	YNORM = 1.0E0
C                         ** SOLVE L*V = Y
	DO 140 K = 1, N
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
	   LM = MIN0(ML, N-K)
	   IF (K .LT. N) CALL SAXPY(LM, T, ABD(M+1,K), 1, Z(K+1), 1)
	   IF (ABS(Z(K)) .GT. 1.0E0) THEN
	      S = 1.0E0 / ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
  140 CONTINUE

	S = 1.0E0/SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
	YNORM = S*YNORM
C                           ** SOLVE  U*Z = W
	DO 160 KB = 1, N
	   K = N + 1 - KB
	   IF (ABS(Z(K)) .GT. ABS(ABD(M,K))) THEN
	      S = ABS(ABD(M,K)) / ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
	   IF (ABD(M,K) .NE. 0.0E0) Z(K) = Z(K)/ABD(M,K)
	   IF (ABD(M,K) .EQ. 0.0E0) Z(K) = 1.0E0
	   LM = MIN0(K, M) - 1
	   LA = M - LM
	   LZ = K - LM
	   T = -Z(K)
	   CALL SAXPY(LM, T, ABD(LA,K), 1, Z(LZ), 1)
  160 CONTINUE
C                              ** MAKE ZNORM = 1.0
	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
	YNORM = S*YNORM

	IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
	IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
	RETURN
	END
	SUBROUTINE  SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

C         FACTORS A REAL BAND MATRIX BY ELIMINATION.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.

C     INPUT:  SAME AS 'SGBCO'

C     ON RETURN:

C        ABD,IPVT    SAME AS 'SGBCO'

C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.

C     (SEE 'SGBCO' FOR DESCRIPTION OF BAND STORAGE MODE)

C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
C                       FROM FORTRAN: MAX0, MIN0

	INTEGER  LDA, N, ML, MU, IPVT(*), INFO
	REAL     ABD(LDA,*)

	REAL     T
	INTEGER  I,ISAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1


	M = ML + MU + 1
	INFO = 0
C                        ** ZERO INITIAL FILL-IN COLUMNS
	J0 = MU + 2
	J1 = MIN0(N, M) - 1
	DO 20 JZ = J0, J1
	   I0 = M + 1 - JZ
	   DO 10 I = I0, ML
	      ABD(I,JZ) = 0.0E0
   10    CONTINUE
   20 CONTINUE
	JZ = J1
	JU = 0

C                       ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
	NM1 = N - 1
	DO 120 K = 1, NM1
	   KP1 = K + 1
C                                  ** ZERO NEXT FILL-IN COLUMN
	   JZ = JZ + 1
	   IF (JZ .LE. N) THEN
	      DO 40 I = 1, ML
	         ABD(I,JZ) = 0.0E0
   40       CONTINUE
	   ENDIF
C                                  ** FIND L = PIVOT INDEX
	   LM = MIN0(ML, N-K)
	   L = ISAMAX(LM+1, ABD(M,K), 1) + M - 1
	   IPVT(K) = L + K - M

	   IF (ABD(L,K) .EQ. 0.0E0) THEN
C                                      ** ZERO PIVOT IMPLIES THIS COLUMN 
C                                      ** ALREADY TRIANGULARIZED
	      INFO = K
	   ELSE
C                                ** INTERCHANGE IF NECESSARY
	      IF (L .NE. M) THEN
	         T = ABD(L,K)
	         ABD(L,K) = ABD(M,K)
	         ABD(M,K) = T
	      ENDIF
C                                   ** COMPUTE MULTIPLIERS
	      T = -1.0E0 / ABD(M,K)
	      CALL SSCAL(LM, T, ABD(M+1,K), 1)

C                               ** ROW ELIMINATION WITH COLUMN INDEXING

	      JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
	      MM = M
	      DO 80 J = KP1, JU
	         L = L - 1
	         MM = MM - 1
	         T = ABD(L,J)
	         IF (L .NE. MM) THEN
	            ABD(L,J) = ABD(MM,J)
	            ABD(MM,J) = T
	         ENDIF
	         CALL SAXPY(LM, T, ABD(M+1,K), 1, ABD(MM+1,J), 1)
   80       CONTINUE

	   ENDIF

  120 CONTINUE

	IPVT(N) = N
	IF (ABD(M,N) .EQ. 0.0E0) INFO = N
	RETURN
	END
	SUBROUTINE  SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )

C         SOLVES THE REAL BAND SYSTEM
C            A * X = B  OR  TRANSPOSE(A) * X = B
C         USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C     INPUT:

C        ABD     REAL(LDA, N)
C                THE OUTPUT FROM SBGCO OR SGBFA.

C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .

C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.

C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.

C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.

C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SBGCO OR SGBFA.

C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.

C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
C                            TRANS(A)  IS THE TRANSPOSE.

C     ON RETURN

C        B       THE SOLUTION VECTOR  X .

C     ERROR CONDITION

C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0
C        OR SGBFA HAS SET INFO .EQ. 0 .

C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE

C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
C                       FROM FORTRAN: MIN0

	INTEGER  LDA, N, ML, MU, IPVT(*), JOB
	REAL     ABD(LDA,*), B(*)

	REAL     SDOT,T
	INTEGER  K,KB,L,LA,LB,LM,M,NM1


	M = MU + ML + 1
	NM1 = N - 1
	IF (JOB .EQ. 0) THEN
C                               ** JOB = 0 , SOLVE  A * X = B
C                               ** FIRST SOLVE L*Y = B
	   IF (ML .NE. 0) THEN
	      DO 20 K = 1, NM1
	         LM = MIN0(ML, N-K)
	         L = IPVT(K)
	         T = B(L)
	         IF (L .NE. K) THEN
	            B(L) = B(K)
	            B(K) = T
	         ENDIF
	         CALL SAXPY( LM, T, ABD(M+1,K), 1, B(K+1), 1 )
   20       CONTINUE
	   ENDIF
C                           ** NOW SOLVE  U*X = Y
	   DO 40 KB = 1, N
	      K = N + 1 - KB
	      B(K) = B(K) / ABD(M,K)
	      LM = MIN0(K, M) - 1
	      LA = M - LM
	      LB = K - LM
	      T = -B(K)
	      CALL SAXPY(LM, T, ABD(LA,K), 1, B(LB), 1)
   40    CONTINUE

	ELSE
C                          ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
C                                  ** FIRST SOLVE  TRANS(U)*Y = B
	   DO 60 K = 1, N
	      LM = MIN0(K, M) - 1
	      LA = M - LM
	      LB = K - LM
	      T = SDOT(LM, ABD(LA,K), 1, B(LB), 1)
	      B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C                                  ** NOW SOLVE TRANS(L)*X = Y
	   IF (ML .NE. 0) THEN
	      DO 80 KB = 1, NM1
	         K = N - KB
	         LM = MIN0(ML, N-K)
	         B(K) = B(K) + SDOT(LM, ABD(M+1,K), 1, B(K+1), 1)
	         L = IPVT(K)
	         IF (L .NE. K) THEN
	            T = B(L)
	            B(L) = B(K)
	            B(K) = T
	         ENDIF
   80       CONTINUE
	   ENDIF

	ENDIF

	RETURN
	END
	SUBROUTINE  SGECO( A, LDA, N,IPVT, RCOND, Z )

C         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
C         AND ESTIMATES THE CONDITION OF THE MATRIX.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C         IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
C         TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.

C     ON ENTRY

C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.

C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .

C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .

C     ON RETURN

C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.

C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.

C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

C     ROUTINES CALLED:  FROM LINPACK: SGEFA
C                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
C                       FROM FORTRAN: ABS, AMAX1, SIGN

	INTEGER  LDA, N, IPVT(*)
	REAL     A(LDA,*), Z(*)
	REAL     RCOND

	REAL     SDOT,EK,T,WK,WKM
	REAL     ANORM,S,SASUM,SM,YNORM
	INTEGER  INFO,J,K,KB,KP1,L


C                        ** COMPUTE 1-NORM OF A
	ANORM = 0.0E0
	DO 10 J = 1, N
	   ANORM = AMAX1( ANORM, SASUM(N,A(1,J),1) )
   10 CONTINUE
C                                      ** FACTOR
	CALL SGEFA(A,LDA,N,IPVT,INFO)

C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.

C                        ** SOLVE TRANS(U)*W = E
	EK = 1.0E0
	DO 20 J = 1, N
	   Z(J) = 0.0E0
   20 CONTINUE

	DO 100 K = 1, N
	   IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
	   IF (ABS(EK-Z(K)) .GT. ABS(A(K,K))) THEN
	      S = ABS(A(K,K)) / ABS(EK-Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      EK = S*EK
	   ENDIF
	   WK = EK - Z(K)
	   WKM = -EK - Z(K)
	   S = ABS(WK)
	   SM = ABS(WKM)
	   IF (A(K,K) .NE. 0.0E0) THEN
	      WK  = WK  / A(K,K)
	      WKM = WKM / A(K,K)
	   ELSE
	      WK  = 1.0E0
	      WKM = 1.0E0
	   ENDIF
	   KP1 = K + 1
	   IF (KP1 .LE. N) THEN
	      DO 60 J = KP1, N
	         SM = SM + ABS(Z(J)+WKM*A(K,J))
	         Z(J) = Z(J) + WK*A(K,J)
	         S = S + ABS(Z(J))
   60       CONTINUE
	      IF (S .LT. SM) THEN
	         T = WKM - WK
	         WK = WKM
	         DO 70 J = KP1, N
	            Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
	      ENDIF
	   ENDIF
	   Z(K) = WK
  100 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
C                                ** SOLVE TRANS(L)*Y = W
	DO 120 KB = 1, N
	   K = N + 1 - KB
	   IF (K .LT. N) Z(K) = Z(K) + SDOT(N-K, A(K+1,K), 1, Z(K+1), 1)
	   IF (ABS(Z(K)) .GT. 1.0E0) THEN
	      S = 1.0E0/ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	   ENDIF
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
  120 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
C                                 ** SOLVE L*V = Y
	YNORM = 1.0E0
	DO 140 K = 1, N
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
	   IF (K .LT. N) CALL SAXPY(N-K, T, A(K+1,K), 1, Z(K+1), 1)
	   IF (ABS(Z(K)) .GT. 1.0E0) THEN
	      S = 1.0E0/ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
  140 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
C                                  ** SOLVE  U*Z = V
	YNORM = S*YNORM
	DO 160 KB = 1, N
	   K = N + 1 - KB
	   IF (ABS(Z(K)) .GT. ABS(A(K,K))) THEN
	      S = ABS(A(K,K))/ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
	   IF (A(K,K) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
	   IF (A(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
	   T = -Z(K)
	   CALL SAXPY(K-1, T, A(1,K), 1, Z(1), 1)
  160 CONTINUE
C                                   ** MAKE ZNORM = 1.0
	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
	YNORM = S*YNORM

	IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
	IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
	RETURN
	END
	SUBROUTINE  SGEFA( A, LDA, N, IPVT, INFO )

C         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .

C     INPUT:  SAME AS 'SGECO'

C     ON RETURN:

C        A,IPVT  SAME AS 'SGECO'

C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.

C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX

	INTEGER  LDA, N, IPVT(*), INFO
	REAL     A(LDA,*)

	REAL     T
	INTEGER  ISAMAX,J,K,KP1,L,NM1


C                      ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
	INFO = 0
	NM1 = N - 1
	DO 60 K = 1, NM1
	   KP1 = K + 1
C                                            ** FIND L = PIVOT INDEX
	   L = ISAMAX( N-K+1, A(K,K), 1) + K-1
	   IPVT(K) = L

	   IF (A(L,K) .EQ. 0.0E0) THEN
C                                     ** ZERO PIVOT IMPLIES THIS COLUMN 
C                                     ** ALREADY TRIANGULARIZED
	      INFO = K
	   ELSE
C                                     ** INTERCHANGE IF NECESSARY
	      IF (L .NE. K) THEN
	         T = A(L,K)
	         A(L,K) = A(K,K)
	         A(K,K) = T
	      ENDIF
C                                     ** COMPUTE MULTIPLIERS
	      T = -1.0E0 / A(K,K)
	      CALL SSCAL( N-K, T, A(K+1,K), 1 )

C                              ** ROW ELIMINATION WITH COLUMN INDEXING
	      DO 30 J = KP1, N
	         T = A(L,J)
	         IF (L .NE. K) THEN
	            A(L,J) = A(K,J)
	            A(K,J) = T
	         ENDIF
	         CALL SAXPY( N-K, T, A(K+1,K), 1, A(K+1,J), 1 )
   30       CONTINUE

	   ENDIF

   60 CONTINUE

	IPVT(N) = N
	IF (A(N,N) .EQ. 0.0E0) INFO = N
	RETURN
	END
	SUBROUTINE  SGESL( A, LDA, N,IPVT, B, JOB )

C         SOLVES THE REAL SYSTEM
C            A * X = B  OR  TRANS(A) * X = B
C         USING THE FACTORS COMPUTED BY SGECO OR SGEFA.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C     ON ENTRY

C        A       REAL(LDA, N)
C                THE OUTPUT FROM SGECO OR SGEFA.

C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .

C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .

C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA.

C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.

C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.

C     ON RETURN

C        B       THE SOLUTION VECTOR  X .

C     ERROR CONDITION

C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
C        OR SGEFA HAS SET INFO .EQ. 0 .

C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE


C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT

	INTEGER  LDA, N, IPVT(*), JOB
	REAL     A(LDA,*), B(*)

	REAL     SDOT,T
	INTEGER  K,KB,L,NM1


	NM1 = N - 1
	IF (JOB .EQ. 0) THEN
C                                 ** JOB = 0 , SOLVE  A * X = B
C                                     ** FIRST SOLVE  L*Y = B
	   DO 20 K = 1, NM1
	      L = IPVT(K)
	      T = B(L)
	      IF (L .NE. K) THEN
	         B(L) = B(K)
	         B(K) = T
	      ENDIF
	      CALL SAXPY( N-K, T, A(K+1,K), 1, B(K+1), 1 )
   20    CONTINUE
C                                    ** NOW SOLVE  U*X = Y
	   DO 40 KB = 1, N
	      K = N + 1 - KB
	      B(K) = B(K) / A(K,K)
	      T = -B(K)
	      CALL SAXPY( K-1, T, A(1,K), 1, B(1), 1 )
   40    CONTINUE

	ELSE
C                         ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
C                                    ** FIRST SOLVE  TRANS(U)*Y = B
	   DO 60 K = 1, N
	      T = SDOT( K-1, A(1,K), 1, B(1), 1 )
	      B(K) = (B(K) - T) / A(K,K)
   60    CONTINUE
C                                    ** NOW SOLVE  TRANS(L)*X = Y
	   DO 80 KB = 1, NM1
	      K = N - KB
	      B(K) = B(K) + SDOT( N-K, A(K+1,K), 1, B(K+1), 1 )
	      L = IPVT(K)
	      IF (L .NE. K) THEN
	         T = B(L)
	         B(L) = B(K)
	         B(K) = T
	      ENDIF
   80    CONTINUE

	ENDIF

	RETURN
	END
	REAL FUNCTION  SASUM( N, SX, INCX )

C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR TO BE SUMMED
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

C --OUTPUT-- SASUM   SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))

	REAL SX(*)


	SASUM = 0.0
	IF( N.LE.0 )  RETURN
	IF( INCX.NE.1 ) THEN
C                                          ** NON-UNIT INCREMENTS
	    DO 10 I = 1, 1+(N-1)*INCX, INCX
	       SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
	ELSE
C                                          ** UNIT INCREMENTS
	   M = MOD(N,6)
	   IF( M.NE.0 ) THEN
C                             ** CLEAN-UP LOOP SO REMAINING VECTOR 
C                             ** LENGTH IS A MULTIPLE OF 6.
	      DO 30  I = 1, M
	        SASUM = SASUM + ABS(SX(I))
   30       CONTINUE
	   ENDIF
C                              ** UNROLL LOOP FOR SPEED
	   DO 50  I = M+1, N, 6
	     SASUM = SASUM + ABS(SX(I))   + ABS(SX(I+1)) + ABS(SX(I+2))
     $                   + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
   50    CONTINUE
	ENDIF

	RETURN
	END
	SUBROUTINE     SAXPY( N, SA, SX, INCX, SY, INCY )

C          Y = A*X + Y  (X, Y = VECTORS, A = SCALAR)

C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SA  SINGLE PRECISION SCALAR MULTIPLIER 'A'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

C --OUTPUT--
C       SY   FOR I = 0 TO N-1, OVERWRITE  SY(LY+I*INCY) WITH 
C                 SA*SX(LX+I*INCX) + SY(LY+I*INCY), 
C            WHERE LX = 1          IF INCX .GE. 0,
C                     = (-INCX)*N  IF INCX .LT. 0
C            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.

	REAL SX(*), SY(*), SA


	IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN

	IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       SY(I) = SY(I) + SA * SX(I)
   10     CONTINUE

	ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN

C                                        ** EQUAL, UNIT INCREMENTS
	   M = MOD(N,4)
	   IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 4.
	      DO 20  I = 1, M
	        SY(I) = SY(I) + SA * SX(I)
   20       CONTINUE
	   ENDIF
C                              ** UNROLL LOOP FOR SPEED
	   DO 30  I = M+1, N, 4
	      SY(I)   = SY(I)   + SA * SX(I)
	      SY(I+1) = SY(I+1) + SA * SX(I+1)
	      SY(I+2) = SY(I+2) + SA * SX(I+2)
	      SY(I+3) = SY(I+3) + SA * SX(I+3)
   30    CONTINUE

	ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
	   IX = 1
	   IY = 1
	   IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
	   IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
	   DO 40  I = 1, N
	      SY(IY) = SY(IY) + SA*SX(IX)
	      IX = IX + INCX
	      IY = IY + INCY
   40    CONTINUE

	ENDIF

	RETURN
	END
	REAL FUNCTION  SDOT( N, SX, INCX, SY, INCY )

C          S.P. DOT PRODUCT OF VECTORS  'X'  AND  'Y'

C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

C --OUTPUT--
C     SDOT   SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
C            WHERE  LX = 1          IF INCX .GE. 0, 
C                      = (-INCX)*N  IF INCX .LT. 0,
C            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.

	REAL SX(*), SY(*)


	SDOT = 0.0
	IF( N.LE.0 )  RETURN

	IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       SDOT = SDOT + SX(I) * SY(I)
   10     CONTINUE

	ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN

C                                        ** EQUAL, UNIT INCREMENTS
	   M = MOD(N,5)
	   IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 4.
	      DO 20  I = 1, M
	         SDOT = SDOT + SX(I) * SY(I)
   20       CONTINUE
	   ENDIF
C                              ** UNROLL LOOP FOR SPEED
	   DO 30  I = M+1, N, 5
	      SDOT = SDOT + SX(I)*SY(I)     + SX(I+1)*SY(I+1)
     $                  + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3)
     $                  + SX(I+4)*SY(I+4)
   30    CONTINUE

	ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
	   IX = 1
	   IY = 1
	   IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
	   IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
	   DO 40  I = 1, N
	      SDOT = SDOT + SX(IX) * SY(IY)
	      IX = IX + INCX
	      IY = IY + INCY
   40    CONTINUE

	ENDIF

	RETURN
	END
	SUBROUTINE     SSCAL( N, SA, SX, INCX )

C         CALCULATE  X = A*X  (X = VECTOR, A = SCALAR)

C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR
C            SA  SINGLE PRECISION SCALE FACTOR
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

C --OUTPUT-- SX  REPLACE  SX(1+I*INCX)  WITH  SA * SX(1+I*INCX) 
C                FOR I = 0 TO N-1

	REAL SA, SX(*)


	IF( N.LE.0 ) RETURN

	IF( INCX.NE.1 ) THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       SX(I) = SA * SX(I)
   10     CONTINUE

	ELSE

	   M = MOD(N,5)
	   IF( M.NE.0 ) THEN
C                           ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                           ** IS A MULTIPLE OF 5.
	      DO 30  I = 1, M
	         SX(I) = SA * SX(I)
   30       CONTINUE
	   ENDIF
C                             ** UNROLL LOOP FOR SPEED
	   DO 50  I = M+1, N, 5
	      SX(I)   = SA * SX(I)
	      SX(I+1) = SA * SX(I+1)
	      SX(I+2) = SA * SX(I+2)
	      SX(I+3) = SA * SX(I+3)
	      SX(I+4) = SA * SX(I+4)
   50    CONTINUE

	ENDIF

	RETURN
	END
	SUBROUTINE     SSWAP( N, SX, INCX, SY, INCY )

C          INTERCHANGE S.P VECTORS  X  AND  Y

C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

C --OUTPUT--
C       SX  INPUT VECTOR SY (UNCHANGED IF N .LE. 0)
C       SY  INPUT VECTOR SX (UNCHANGED IF N .LE. 0)

C     FOR I = 0 TO N-1, INTERCHANGE  SX(LX+I*INCX) AND SY(LY+I*INCY),
C     WHERE LX = 1          IF INCX .GE. 0, 
C              = (-INCX)*N  IF INCX .LT. 0
C     AND LY IS DEFINED IN A SIMILAR WAY USING INCY.

	REAL SX(*), SY(*), STEMP1, STEMP2, STEMP3


	IF( N.LE.0 ) RETURN

	IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       STEMP1 = SX(I)
	       SX(I) = SY(I)
	       SY(I) = STEMP1
   10     CONTINUE

	ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN

C                                        ** EQUAL, UNIT INCREMENTS
	   M = MOD(N,3)
	   IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 3.
	      DO 20  I = 1, M
	         STEMP1 = SX(I)
	         SX(I) = SY(I)
	         SY(I) = STEMP1
   20       CONTINUE
	   ENDIF
C                              ** UNROLL LOOP FOR SPEED
	   DO 30  I = M+1, N, 3
	      STEMP1  = SX(I)
	      STEMP2  = SX(I+1)
	      STEMP3  = SX(I+2)
	      SX(I)   = SY(I)
	      SX(I+1) = SY(I+1)
	      SX(I+2) = SY(I+2)
	      SY(I)   = STEMP1
	      SY(I+1) = STEMP2
	      SY(I+2) = STEMP3
   30    CONTINUE

	ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
	   IX = 1
	   IY = 1
	   IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
	   IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
	   DO 40  I = 1, N
	      STEMP1 = SX(IX)
	      SX(IX) = SY(IY)
	      SY(IY) = STEMP1
	      IX = IX + INCX
	      IY = IY + INCY
   40    CONTINUE

	ENDIF

	RETURN
	END
	INTEGER FUNCTION  ISAMAX( N, SX, INCX )

C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR OF INTEREST
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

C --OUTPUT-- ISAMAX   FIRST I, I = 1 TO N, TO MAXIMIZE
C                         ABS(SX(1+(I-1)*INCX))

	REAL SX(*), SMAX, XMAG


	IF( N.LE.0 ) THEN
	   ISAMAX = 0
	ELSE IF( N.EQ.1 ) THEN
	   ISAMAX = 1
	ELSE
	   SMAX = 0.0
	   II = 1
	   DO 20  I = 1, 1+(N-1)*INCX, INCX
	      XMAG = ABS(SX(I))
	      IF( SMAX.LT.XMAG ) THEN
	         SMAX = XMAG
	         ISAMAX = II
	      ENDIF
	      II = II + 1
   20    CONTINUE
	ENDIF

	RETURN
	END
      FUNCTION D1MACH(i)

*-----------------------------------------------------------------------------*
*= PURPOSE:                                                                  =*
*= D1MACH calculates various machine constants in single precision.          =*
*-----------------------------------------------------------------------------*
*= PARAMETERS:                                                               =*
*=   I       -  INTEGER, identifies the machine constant (0<I<5)         (I) =*
*=   D1MACH  -  REAL, machine constant in single precision               (O) =*
*=      I=1     - the smallest non-vanishing normalized floating-point       =*
*=                power of the radix, i.e., D1MACH=FLOAT(IBETA)**MINEXP      =*
*=      I=2     - the largest finite floating-point number.  In              =*
*=                particular D1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*
*=                Note - on some machines D1MACH will be only the            =*
*=                second, or perhaps third, largest number, being            =*
*=                too small by 1 or 2 units in the last digit of             =*
*=                the significand.                                           =*
*=      I=3     - A small positive floating-point number such that           =*
*=                1.0-D1MACH .NE. 1.0. In particular, if IBETA = 2           =*
*=                or  IRND = 0, D1MACH = FLOAT(IBETA)**NEGEPS.               =*
*=                Otherwise,  D1MACH = (IBETA**NEGEPS)/2.  Because           =*
*=                NEGEPS is bounded below by -(IT+3), D1MACH may not         =*
*=                be the smallest number that can alter 1.0 by               =*
*=                subtraction.                                               =*
*=      I=4     - the smallest positive floating-point number such           =*
*=                that  1.0+D1MACH .NE. 1.0. In particular, if either        =*
*=                IBETA = 2  or  IRND = 0, D1MACH=FLOAT(IBETA)**MACHEP.      =*
*=                Otherwise, D1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*
*=  (see routine T665D for more information on different constants)          =*
*-----------------------------------------------------------------------------*

      DOUBLE PRECISION d1mach
      INTEGER i
   
      LOGICAL doinit
      DATA doinit/.TRUE./
      SAVE doinit

      DOUBLE PRECISION dmach(4) 
      SAVE dmach

      IF (( i .GE. 1 ) .AND. ( i .LE. 4 )) THEN
* compute constants at first call only
        IF (doinit) THEN
           CALL t665d(dmach)
           doinit = .FALSE.
        ENDIF
        d1mach = dmach(i)
      ELSE
        WRITE(0,*) '>>> ERROR (D1MACH) <<<  invalid argument'
        STOP
      ENDIF

*!csm
*!!! over-ride by sm on 5/26/03.  For some compilers than don't allow
* calculation of d1mach(4).  Use value found on ACD server.

c      if( i .eq. 4) d1mach = 2.22e-15

      END


C      ALGORITHM 665, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 14, NO. 4, PP. 303-311.
      SUBROUTINE T665D(DMACH)
C-----------------------------------------------------------------------
C This subroutine is a double precision version of subroutine T665R.
C See code of T665R for detailed comments and explanation
C-----------------------------------------------------------------------
      DOUBLE PRECISION DMACH(4)
      INTEGER I,IBETA,IEXP,IRND,IT,ITEMP,IZ,J,K,MACHEP,MAXEXP,
     1        MINEXP,MX,NEGEP,NGRD,NXRES
CS    REAL A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,T,TEMP,TEMPA,
CS   1     TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
      DOUBLE PRECISION A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,
     1                 T,TEMP,TEMPA,TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
C-----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      ONE = CONV(1)
      TWO = ONE + ONE
      ZERO = ONE - ONE
C-----------------------------------------------------------------------
C  Determine IBETA, BETA ala Malcolm.
C-----------------------------------------------------------------------
      A = ONE
   10 A = A + A
         TEMP = A+ONE
         TEMP1 = TEMP-A
         IF (TEMP1-ONE .EQ. ZERO) GO TO 10
      B = ONE
   20 B = B + B
         TEMP = A+B
         ITEMP = INT(TEMP-A)
         IF (ITEMP .EQ. 0) GO TO 20
      IBETA = ITEMP
      BETA = CONV(IBETA)
C-----------------------------------------------------------------------
C  Determine IT, IRND.
C-----------------------------------------------------------------------
      IT = 0
      B = ONE
  100 IT = IT + 1
         B = B * BETA
         TEMP = B+ONE
         TEMP1 = TEMP-B
         IF (TEMP1-ONE .EQ. ZERO) GO TO 100
      IRND = 0
      BETAH = BETA / TWO
      TEMP = A+BETAH
      IF (TEMP-A .NE. ZERO) IRND = 1
      TEMPA = A + BETA
      TEMP = TEMPA+BETAH
      IF ((IRND .EQ. 0) .AND. (TEMP-TEMPA .NE. ZERO)) IRND = 2
C-----------------------------------------------------------------------
C  Determine NEGEP, EPSNEG.
C-----------------------------------------------------------------------
      NEGEP = IT + 3
      BETAIN = ONE / BETA
      A = ONE
      DO 200 I = 1, NEGEP
         A = A * BETAIN
  200 CONTINUE
      B = A
  210 TEMP = ONE-A
         IF (TEMP-ONE .NE. ZERO) GO TO 220
         A = A * BETA
         NEGEP = NEGEP - 1
      GO TO 210
  220 NEGEP = -NEGEP
      EPSNEG = A
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 300
      A = (A*(ONE+A)) / TWO
      TEMP = ONE-A
      IF (TEMP-ONE .NE. ZERO) EPSNEG = A
C-----------------------------------------------------------------------
C  Determine MACHEP, EPS.
C-----------------------------------------------------------------------
  300 MACHEP = -IT - 3
      A = B
  310 TEMP = ONE+A
         IF (TEMP-ONE .NE. ZERO) GO TO 320
         A = A * BETA
         MACHEP = MACHEP + 1
      GO TO 310
  320 EPS = A
      TEMP = TEMPA+BETA*(ONE+EPS)
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 350
      A = (A*(ONE+A)) / TWO
      TEMP = ONE+A
      IF (TEMP-ONE .NE. ZERO) EPS = A
C-----------------------------------------------------------------------
C  Determine NGRD.
C-----------------------------------------------------------------------
  350 NGRD = 0
      TEMP = ONE+EPS
      IF ((IRND .EQ. 0) .AND. (TEMP*ONE-ONE .NE. ZERO)) NGRD = 1
C-----------------------------------------------------------------------
C  Determine IEXP, MINEXP, XMIN.
C
C  Loop to determine largest I and K = 2**I such that
C         (1/BETA) ** (2**(I))
C  does not underflow.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
      I = 0
      K = 1
      Z = BETAIN
      T = ONE + EPS
      NXRES = 0
  400 Y = Z
         Z = Y * Y
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Z * ONE
         TEMP = Z * T
         IF ((A+A .EQ. ZERO) .OR. (ABS(Z) .GE. Y)) GO TO 410
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .EQ. Z) GO TO 410
         I = I + 1
         K = K + K
      GO TO 400
  410 IF (IBETA .EQ. 10) GO TO 420
      IEXP = I + 1
      MX = K + K
      GO TO 450
C-----------------------------------------------------------------------
C  This segment is for decimal machines only.
C-----------------------------------------------------------------------
  420 IEXP = 2
      IZ = IBETA
  430 IF (K .LT. IZ) GO TO 440
         IZ = IZ * IBETA
         IEXP = IEXP + 1
      GO TO 430
  440 MX = IZ + IZ - 1
C-----------------------------------------------------------------------
C  Loop to determine MINEXP, XMIN.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
  450 XMIN = Y
         Y = Y * BETAIN
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Y * ONE
         TEMP = Y * T
         IF (((A+A) .EQ. ZERO) .OR. (ABS(Y) .GE. XMIN)) GO TO 460
         K = K + 1
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .NE. Y) GO TO 450
      NXRES = 3
      XMIN = Y
  460 MINEXP = -K
C-----------------------------------------------------------------------
C  Determine MAXEXP, XMAX.
C-----------------------------------------------------------------------
      IF ((MX .GT. K+K-3) .OR. (IBETA .EQ. 10)) GO TO 500
      MX = MX + MX
      IEXP = IEXP + 1
  500 MAXEXP = MX + MINEXP
C-----------------------------------------------------------------
C  Adjust IRND to reflect partial underflow.
C-----------------------------------------------------------------
      IRND = IRND + NXRES
C-----------------------------------------------------------------
C  Adjust for IEEE-style machines.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 2) .OR. (IRND .EQ. 5)) MAXEXP = MAXEXP - 2
C-----------------------------------------------------------------
C  Adjust for non-IEEE machines with partial underflow.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 3) .OR. (IRND .EQ. 4)) MAXEXP = MAXEXP - IT
C-----------------------------------------------------------------
C  Adjust for machines with implicit leading bit in binary
C  significand, and machines with radix point at extreme
C  right of significand.
C-----------------------------------------------------------------
      I = MAXEXP + MINEXP
      IF ((IBETA .EQ. 2) .AND. (I .EQ. 0)) MAXEXP = MAXEXP - 1
      IF (I .GT. 20) MAXEXP = MAXEXP - 1
      IF (A .NE. Y) MAXEXP = MAXEXP - 2
      XMAX = ONE - EPSNEG
      IF (XMAX*ONE .NE. XMAX) XMAX = ONE - BETA * EPSNEG
      XMAX = XMAX / (BETA * BETA * BETA * XMIN)
      I = MAXEXP + MINEXP + 3
      IF (I .LE. 0) GO TO 520
      DO 510 J = 1, I
          IF (IBETA .EQ. 2) XMAX = XMAX + XMAX
          IF (IBETA .NE. 2) XMAX = XMAX * BETA
  510 CONTINUE
      DMACH(1) = XMIN
      DMACH(2) = XMAX
      DMACH(3) = EPSNEG
      DMACH(4) = EPS
  520 RETURN
C---------- LAST CARD OF T665D ----------
      END


      FUNCTION R1MACH(i)

*-----------------------------------------------------------------------------*
*= PURPOSE:                                                                  =*
*= R1MACH calculates various machine constants in single precision.          =*
*-----------------------------------------------------------------------------*
*= PARAMETERS:                                                               =*
*=   I       -  INTEGER, identifies the machine constant (0<I<5)         (I) =*
*=   R1MACH  -  REAL, machine constant in single precision               (O) =*
*=      I=1     - the smallest non-vanishing normalized floating-point       =*
*=                power of the radix, i.e., R1MACH=FLOAT(IBETA)**MINEXP      =*
*=      I=2     - the largest finite floating-point number.  In              =*
*=                particular R1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*
*=                Note - on some machines R1MACH will be only the            =*
*=                second, or perhaps third, largest number, being            =*
*=                too small by 1 or 2 units in the last digit of             =*
*=                the significand.                                           =*
*=      I=3     - A small positive floating-point number such that           =*
*=                1.0-R1MACH .NE. 1.0. In particular, if IBETA = 2           =*
*=                or  IRND = 0, R1MACH = FLOAT(IBETA)**NEGEPS.               =*
*=                Otherwise,  R1MACH = (IBETA**NEGEPS)/2.  Because           =*
*=                NEGEPS is bounded below by -(IT+3), R1MACH may not         =*
*=                be the smallest number that can alter 1.0 by               =*
*=                subtraction.                                               =*
*=      I=4     - the smallest positive floating-point number such           =*
*=                that  1.0+R1MACH .NE. 1.0. In particular, if either        =*
*=                IBETA = 2  or  IRND = 0, R1MACH=FLOAT(IBETA)**MACHEP.      =*
*=                Otherwise, R1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*
*=  (see routine T665R for more information on different constants)          =*
*-----------------------------------------------------------------------------*

      REAL r1mach
      INTEGER i
   
      LOGICAL doinit
      DATA doinit/.TRUE./
      SAVE doinit

      REAL rmach(4) 
      SAVE rmach

      IF (( i .GE. 1 ) .AND. ( i .LE. 4 )) THEN
* compute constants at first call only
        IF (doinit) THEN
           CALL t665r(rmach)
           doinit = .FALSE.
        ENDIF
        r1mach = rmach(i)
      ELSE
        WRITE(0,*) '>>> ERROR (R1MACH) <<<  invalid argument'
        STOP
      ENDIF

      END


C      ALGORITHM 665, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 14, NO. 4, PP. 303-311.
      SUBROUTINE T665R(RMACH)
C-----------------------------------------------------------------------
C  This Fortran 77 subroutine is intended to determine the parameters
C   of the floating-point arithmetic system specified below.  The
C   determination of the first three uses an extension of an algorithm
C   due to M. Malcolm, CACM 15 (1972), pp. 949-951, incorporating some,
C   but not all, of the improvements suggested by M. Gentleman and S.
C   Marovich, CACM 17 (1974), pp. 276-277.  An earlier version of this
C   program was published in the book Software Manual for the
C   Elementary Functions by W. J. Cody and W. Waite, Prentice-Hall,
C   Englewood Cliffs, NJ, 1980.
C
C  The program as given here must be modified before compiling.  If
C   a single (double) precision version is desired, change all
C   occurrences of CS (CD) in columns 1 and 2 to blanks.
C
C  Parameter values reported are as follows:
C
C       IBETA   - the radix for the floating-point representation
C       IT      - the number of base IBETA digits in the floating-point
C                 significand
C       IRND    - 0 if floating-point addition chops
C                 1 if floating-point addition rounds, but not in the
C                   IEEE style
C                 2 if floating-point addition rounds in the IEEE style
C                 3 if floating-point addition chops, and there is
C                   partial underflow
C                 4 if floating-point addition rounds, but not in the
C                   IEEE style, and there is partial underflow
C                 5 if floating-point addition rounds in the IEEE style,
C                   and there is partial underflow
C       NGRD    - the number of guard digits for multiplication with
C                 truncating arithmetic.  It is
C                 0 if floating-point arithmetic rounds, or if it
C                   truncates and only  IT  base  IBETA digits
C                   participate in the post-normalization shift of the
C                   floating-point significand in multiplication;
C                 1 if floating-point arithmetic truncates and more
C                   than  IT  base  IBETA  digits participate in the
C                   post-normalization shift of the floating-point
C                   significand in multiplication.
C       MACHEP  - the largest negative integer such that
C                 1.0+FLOAT(IBETA)**MACHEP .NE. 1.0, except that
C                 MACHEP is bounded below by  -(IT+3)
C       NEGEPS  - the largest negative integer such that
C                 1.0-FLOAT(IBETA)**NEGEPS .NE. 1.0, except that
C                 NEGEPS is bounded below by  -(IT+3)
C       IEXP    - the number of bits (decimal places if IBETA = 10)
C                 reserved for the representation of the exponent
C                 (including the bias or sign) of a floating-point
C                 number
C       MINEXP  - the largest in magnitude negative integer such that
C                 FLOAT(IBETA)**MINEXP is positive and normalized
C       MAXEXP  - the smallest positive power of  BETA  that overflows
C       EPS     - the smallest positive floating-point number such
C                 that  1.0+EPS .NE. 1.0. In particular, if either
C                 IBETA = 2  or  IRND = 0, EPS = FLOAT(IBETA)**MACHEP.
C                 Otherwise,  EPS = (FLOAT(IBETA)**MACHEP)/2
C       EPSNEG  - A small positive floating-point number such that
C                 1.0-EPSNEG .NE. 1.0. In particular, if IBETA = 2
C                 or  IRND = 0, EPSNEG = FLOAT(IBETA)**NEGEPS.
C                 Otherwise,  EPSNEG = (IBETA**NEGEPS)/2.  Because
C                 NEGEPS is bounded below by -(IT+3), EPSNEG may not
C                 be the smallest number that can alter 1.0 by
C                 subtraction.
C       XMIN    - the smallest non-vanishing normalized floating-point
C                 power of the radix, i.e.,  XMIN = FLOAT(IBETA)**MINEXP
C       XMAX    - the largest finite floating-point number.  In
C                 particular  XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP
C                 Note - on some machines  XMAX  will be only the
C                 second, or perhaps third, largest number, being
C                 too small by 1 or 2 units in the last digit of
C                 the significand.
C
C     Latest revision - April 20, 1987
C
C     Author - W. J. Cody
C              Argonne National Laboratory
C
C-----------------------------------------------------------------------
      REAL rmach(4)
      INTEGER I,IBETA,IEXP,IRND,IT,ITEMP,IZ,J,K,MACHEP,MAXEXP,
     1        MINEXP,MX,NEGEP,NGRD,NXRES
      REAL A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,T,TEMP,TEMPA,
     1     TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
CD    DOUBLE PRECISION A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,
CD   1                 T,TEMP,TEMPA,TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
C-----------------------------------------------------------------------
      CONV(I) = REAL(I)
CD    CONV(I) = DBLE(I)
      ONE = CONV(1)
      TWO = ONE + ONE
      ZERO = ONE - ONE
C-----------------------------------------------------------------------
C  Determine IBETA, BETA ala Malcolm.
C-----------------------------------------------------------------------
      A = ONE
   10 A = A + A
         TEMP = A+ONE
         TEMP1 = TEMP-A
         IF (TEMP1-ONE .EQ. ZERO) GO TO 10
      B = ONE
   20 B = B + B
         TEMP = A+B
         ITEMP = INT(TEMP-A)
         IF (ITEMP .EQ. 0) GO TO 20
      IBETA = ITEMP
      BETA = CONV(IBETA)
C-----------------------------------------------------------------------
C  Determine IT, IRND.
C-----------------------------------------------------------------------
      IT = 0
      B = ONE
  100 IT = IT + 1
         B = B * BETA
         TEMP = B+ONE
         TEMP1 = TEMP-B
         IF (TEMP1-ONE .EQ. ZERO) GO TO 100
      IRND = 0
      BETAH = BETA / TWO
      TEMP = A+BETAH
      IF (TEMP-A .NE. ZERO) IRND = 1
      TEMPA = A + BETA
      TEMP = TEMPA+BETAH
      IF ((IRND .EQ. 0) .AND. (TEMP-TEMPA .NE. ZERO)) IRND = 2
C-----------------------------------------------------------------------
C  Determine NEGEP, EPSNEG.
C-----------------------------------------------------------------------
      NEGEP = IT + 3
      BETAIN = ONE / BETA
      A = ONE
      DO 200 I = 1, NEGEP
         A = A * BETAIN
  200 CONTINUE
      B = A
  210 TEMP = ONE-A
         IF (TEMP-ONE .NE. ZERO) GO TO 220
         A = A * BETA
         NEGEP = NEGEP - 1
      GO TO 210
  220 NEGEP = -NEGEP
      EPSNEG = A
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 300
      A = (A*(ONE+A)) / TWO
      TEMP = ONE-A
      IF (TEMP-ONE .NE. ZERO) EPSNEG = A
C-----------------------------------------------------------------------
C  Determine MACHEP, EPS.
C-----------------------------------------------------------------------
  300 MACHEP = -IT - 3
      A = B
  310 TEMP = ONE+A
         IF (TEMP-ONE .NE. ZERO) GO TO 320
         A = A * BETA
         MACHEP = MACHEP + 1
      GO TO 310
  320 EPS = A
      TEMP = TEMPA+BETA*(ONE+EPS)
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 350
      A = (A*(ONE+A)) / TWO
      TEMP = ONE+A
      IF (TEMP-ONE .NE. ZERO) EPS = A
C-----------------------------------------------------------------------
C  Determine NGRD.
C-----------------------------------------------------------------------
  350 NGRD = 0
      TEMP = ONE+EPS
      IF ((IRND .EQ. 0) .AND. (TEMP*ONE-ONE .NE. ZERO)) NGRD = 1
C-----------------------------------------------------------------------
C  Determine IEXP, MINEXP, XMIN.
C
C  Loop to determine largest I and K = 2**I such that
C         (1/BETA) ** (2**(I))
C  does not underflow.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
      I = 0
      K = 1
      Z = BETAIN
      T = ONE + EPS
      NXRES = 0
  400 Y = Z
         Z = Y * Y
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Z * ONE
         TEMP = Z * T
         IF ((A+A .EQ. ZERO) .OR. (ABS(Z) .GE. Y)) GO TO 410
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .EQ. Z) GO TO 410
         I = I + 1
         K = K + K
      GO TO 400
  410 IF (IBETA .EQ. 10) GO TO 420
      IEXP = I + 1
      MX = K + K
      GO TO 450
C-----------------------------------------------------------------------
C  This segment is for decimal machines only.
C-----------------------------------------------------------------------
  420 IEXP = 2
      IZ = IBETA
  430 IF (K .LT. IZ) GO TO 440
         IZ = IZ * IBETA
         IEXP = IEXP + 1
      GO TO 430
  440 MX = IZ + IZ - 1
C-----------------------------------------------------------------------
C  Loop to determine MINEXP, XMIN.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
  450 XMIN = Y
         Y = Y * BETAIN
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Y * ONE
         TEMP = Y * T
         IF (((A+A) .EQ. ZERO) .OR. (ABS(Y) .GE. XMIN)) GO TO 460
         K = K + 1
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .NE. Y) GO TO 450
      NXRES = 3
      XMIN = Y
  460 MINEXP = -K
C-----------------------------------------------------------------------
C  Determine MAXEXP, XMAX.
C-----------------------------------------------------------------------
      IF ((MX .GT. K+K-3) .OR. (IBETA .EQ. 10)) GO TO 500
      MX = MX + MX
      IEXP = IEXP + 1
  500 MAXEXP = MX + MINEXP
C-----------------------------------------------------------------
C  Adjust IRND to reflect partial underflow.
C-----------------------------------------------------------------
      IRND = IRND + NXRES
C-----------------------------------------------------------------
C  Adjust for IEEE-style machines.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 2) .OR. (IRND .EQ. 5)) MAXEXP = MAXEXP - 2
C-----------------------------------------------------------------
C  Adjust for non-IEEE machines with partial underflow.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 3) .OR. (IRND .EQ. 4)) MAXEXP = MAXEXP - IT
C-----------------------------------------------------------------
C  Adjust for machines with implicit leading bit in binary
C  significand, and machines with radix point at extreme
C  right of significand.
C-----------------------------------------------------------------
      I = MAXEXP + MINEXP
      IF ((IBETA .EQ. 2) .AND. (I .EQ. 0)) MAXEXP = MAXEXP - 1
      IF (I .GT. 20) MAXEXP = MAXEXP - 1
      IF (A .NE. Y) MAXEXP = MAXEXP - 2
      XMAX = ONE - EPSNEG
      IF (XMAX*ONE .NE. XMAX) XMAX = ONE - BETA * EPSNEG
      XMAX = XMAX / (BETA * BETA * BETA * XMIN)
      I = MAXEXP + MINEXP + 3
      IF (I .LE. 0) GO TO 520
      DO 510 J = 1, I
          IF (IBETA .EQ. 2) XMAX = XMAX + XMAX
          IF (IBETA .NE. 2) XMAX = XMAX * BETA
  510 CONTINUE
      RMACH(1) = XMIN
      RMACH(2) = XMAX
      RMACH(3) = EPSNEG
      RMACH(4) = EPS
  520 RETURN
C---------- LAST CARD OF T665R ----------
      END

