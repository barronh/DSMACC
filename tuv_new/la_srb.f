* This file contains the following subroutines, related to the calculation
* of radiation at Lyman-alpha and Schumann-Runge wavelengths:
*     la_srb
*     lymana
*     schum
*     effxs
*     calc_params
*     init_xs
*     sjo2   
* and the following functions
*     chebev
*=============================================================================*

      SUBROUTINE la_srb(nz,z,tlev,nw,wl,vcol,scol,o2xs1,dto2,o2xs)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Compute equivalent optical depths for O2 absorption, and O2 effective    =*
*=  absorption cross sections, parameterized in the Lyman-alpha and SR bands =*
*-----------------------------------------------------------------------------* 
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL      - REAL, vector of lxower limits of wavelength intervals in    (I)=*
*=            working wavelength grid                                        =*
*=  CZ      - REAL, number of air molecules per cm^2 at each specified    (I)=*
*=            altitude layer                                                 =*
*=  ZEN     - REAL, solar zenith angle                                    (I)=*
*=                                                                           =*
*=  O2XS1   - REAL, O2 cross section from rdo2xs                          (I)=*
*=                                                                           =*
*=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer at each specified wavelength                    =*
*=  O2XS    - REAL, molecular absorption cross section in SR bands at     (O)=*
*=            each specified altitude and wavelength.  Includes Herzberg     =*
*=            continuum.                                                     =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      REAL wl(kw)
      REAL z(kz)
      INTEGER nz, nw, iz, iw

      REAL vcol(kz), scol(kz)
      REAL o2col(kz)
      REAL o2xs1(kw)
      REAL dto2(kz,kw), o2xs(kz,kw)
      REAL secchi(kz)
      REAL tlev(kz)

* Lyman-alpha variables
* O2 optical depth and equivalent cross section in the Lyman-alpha region

      INTEGER ila, nla, kla
      PARAMETER (kla = 2)
      REAL wlla(kla)
      REAL dto2la(kz, kla-1), o2xsla(kz, kla-1)
      SAVE ila

* grid on which Koppers' parameterization is defined
* O2 optical depth and equivalent cross section on Koppers' grid

      INTEGER isrb, nsrb, ksrb
      PARAMETER(ksrb = 18)
      REAL wlsrb(ksrb)
      REAL dto2k(kz, ksrb-1), o2xsk(kz, ksrb-1)
      SAVE isrb

      INTEGER i

      LOGICAL call1
      DATA call1/.TRUE./
      SAVE call1

* Wavelengths for Lyman alpha and SRB parameterizations:

      DATA nla /1/
      DATA wlla/ 121.4, 121.9/

      DATA nsrb /17/
      DATA wlsrb/174.4, 177.0, 178.6, 180.2, 181.8, 183.5, 185.2, 186.9,
     $     188.7, 190.5, 192.3, 194.2, 196.1, 198.0, 200.0, 202.0, 
     $     204.1, 205.8/

*----------------------------------------------------------------------
* initalize O2 cross sections 
*----------------------------------------------------------------------

      DO iz = 1, nz
         DO iw =1, nw - 1   
            o2xs(iz,iw) = o2xs1(iw)
         ENDDO  
      ENDDO

      IF(wl(1) .GT. wlsrb(nsrb)) RETURN


*----------------------------------------------------------------------
* Slant O2 column and x-sections.
*----------------------------------------------------------------------

      DO iz = 1, nz
         o2col(iz) = 0.2095 * scol(iz)
      ENDDO

*----------------------------------------------------------------------
* On first call, check that the user wavelength grid, WL(IW), is compatible 
* with the wavelengths for the parameterizations of the Lyman-alpha and SRB.
* Also compute and save corresponding grid indices (ILA, ISRB)
*----------------------------------------------------------------------

      IF (call1) THEN

** locate Lyman-alpha wavelengths on grid

         ila = 0
         DO iw = 1, nw
            IF(ABS(wl(iw) - wlla(1)) .LT. 10.*precis) THEN
               ila = iw
               GO TO 5
            ENDIF
         ENDDO
 5       CONTINUE

* check 

         IF(ila .EQ. 0) THEN
            WRITE(*,*) 'For wavelengths below 205.8 nm, only the'
            WRITE(*,*) 'pre-specified wavelength grid is permitted'
            WRITE(*,*) 'Use nwint=-156, or edit subroutine gridw.f'
            STOP ' Lyman alpha grid mis-match - 1'
         ENDIF
         DO i = 2, nla + 1
            IF(ABS(wl(ila + i - 1) - wlla(i)) .GT. 10.*precis) THEN
               WRITE(*,*) 'Lyman alpha grid mis-match - 2'
               STOP
            ENDIF
         ENDDO

** locate Schumann-Runge wavelengths on grid

         isrb = 0
         DO iw = 1, nw
            IF(ABS(wl(iw) - wlsrb(1)) .LT. 10.*precis) THEN
               isrb = iw
               GO TO 6
            ENDIF
         ENDDO
 6       CONTINUE

* check

         IF(isrb .EQ. 0) THEN
            WRITE(*,*) 'For wavelengths below 205.8 nm, only the'
            WRITE(*,*) 'pre-specified wavelength grid is permitted'
            WRITE(*,*) 'Use nwint=-156, or edit subroutine gridw.f'
            STOP ' SRB grid mis-match - 1'
         ENDIF
         DO i = 2, nsrb + 1
            IF(ABS(wl(isrb + i - 1) - wlsrb(i)) .GT. 10.* precis) THEN
               WRITE(*,*) ' SRB grid mismatch - w'
               STOP
            ENDIF
         ENDDO

         IF (call1) call1 = .FALSE.
      ENDIF

*----------------------------------------------------------------------
* Effective secant of solar zenith angle.  
* Use 2.0 if no direct sun (value for isotropic radiation)
* For nz, use value at nz-1
*----------------------------------------------------------------------

      DO i = 1, nz - 1
         secchi(i) = scol(i)/vcol(i)
         IF(scol(i) .GT. largest/10.) secchi(i) = 2.
      ENDDO
      secchi(nz) = secchi(nz-1)

*---------------------------------------------------------------------
* Lyman-Alpha parameterization, output values of O2 optical depth
* and O2 effective (equivalent) cross section
*----------------------------------------------------------------------

      CALL lymana(nz,o2col,secchi,dto2la,o2xsla)
      DO iw = ila, ila + nla - 1
         DO iz = 1, nz
            dto2(iz,iw) = dto2la(iz, iw - ila + 1)
            o2xs(iz,iw) = o2xsla(iz, iw - ila + 1)
         ENDDO
      ENDDO

*------------------------------------------------------------------------------
* Koppers' parameterization of the SR bands, output values of O2
* optical depth and O2 equivalent cross section 
*------------------------------------------------------------------------------

      CALL schum(nz,o2col,tlev,secchi,dto2k,o2xsk)
      DO iw = isrb, isrb + nsrb - 1
         DO iz = 1, nz
            dto2(iz,iw) = dto2k(iz, iw - isrb + 1)
            o2xs(iz,iw) = o2xsk(iz, iw - isrb + 1)
         ENDDO
      ENDDO

      RETURN
      END

*=============================================================================*

      SUBROUTINE lymana(nz,o2col,secchi,dto2la,o2xsla)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the effective absorption cross section of O2 in the Lyman-Alpha=*
*=  bands and an effective O2 optical depth at all altitudes.  Parameterized =*
*=  after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the  =*
*=  absorption of the solar Lyman-Alpha line, Geophysical Research Letters,  =*
*=  Vol.24, No.21, pp 2659-2662, 1997.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
*=            altitude                                                       =*
*=  DTO2LA  - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer                                                 =*
*=  O2XSLA  - REAL, molecular absorption cross section in LA bands        (O)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:

      INCLUDE 'params'
      INTEGER nz
      REAL o2col(kz)
      REAL secchi(kz)

* output

      REAL dto2la(kz,*), o2xsla(kz,*)

* local variables

      DOUBLE PRECISION rm(kz), ro2(kz)
      DOUBLE PRECISION b(3), c(3), d(3), e(3)
      DATA b/ 6.8431D-01, 2.29841D-01,  8.65412D-02/,
     >     c/8.22114D-21, 1.77556D-20,  8.22112D-21/,
     >     d/ 6.0073D-21, 4.28569D-21,  1.28059D-20/,
     >     e/8.21666D-21, 1.63296D-20,  4.85121D-17/

      INTEGER iz, i
      REAL xsmin

*------------------------------------------------------------------------------*
*sm:  set minimum cross section

      xsmin = 1.e-20

* calculate reduction factors at every altitude

      DO iz = 1, nz
        rm(iz) = 0.D+00
        ro2(iz) = 0.D+00
        DO i = 1, 3
          rm(iz) = rm(iz) + b(i) * DEXP(-c(i) * DBLE(o2col(iz)))
          ro2(iz) = ro2(iz) + d(i) * DEXP(-e(i) * DBLE(o2col(iz)))
        ENDDO
      ENDDO

* calculate effective O2 optical depths and effective O2 cross sections

      DO iz = 1, nz-1

         IF (rm(iz) .GT. 1.0D-100) THEN
            IF (ro2(iz) .GT. 1.D-100) THEN
               o2xsla(iz,1) = ro2(iz)/rm(iz)
            ELSE
               o2xsla(iz,1) = xsmin
            ENDIF

            IF (rm(iz+1) .GT. 0.) THEN

               dto2la(iz,1) = LOG(rm(iz+1)) / secchi(iz+1) 
     $                      - LOG(rm(iz))   / secchi(iz)

            ELSE
               dto2la(iz,1) = 1000.
            ENDIF
         ELSE
            dto2la(iz,1) = 1000.
            o2xsla(iz,1) = xsmin
         ENDIF

      ENDDO

* do top layer separately

      dto2la(nz,1) = 0.
      IF(rm(nz) .GT. 1.D-100) THEN
         o2xsla(nz,1) = ro2(nz)/rm(nz)
      ELSE
         o2xsla(nz,1) = xsmin
      ENDIF

*----------------------------------------------------------------------------*

      END

*=============================================================================*

      SUBROUTINE schum(nz, o2col, tlev, secchi, dto2, o2xsk)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the equivalent absorption cross section of O2 in the SR bands. =*
*=  The algorithm is based on parameterization of G.A. Koppers, and          =*
*=  D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]                         =*
*=  Final values do include effects from the Herzberg continuum.             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
*=            altitude                                                       =*
*=  TLEV    - tmeperature at each level                                   (I)=*
*=  SECCHI  - ratio of slant to vertical o2 columns                       (I)=*
*=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer at each specified wavelength                    =*
*=  O2XSK  - REAL, molecular absorption cross section in SR bands at     (O)=*
*=            each specified wavelength.  Includes Herzberg continuum        =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER nz
      REAL o2col(kz), o2col1(kz)
      REAL tlev(kz), secchi(kz)

      REAL dto2(kz,17), o2xsk(kz,17)

      INTEGER i, k, ktop, ktop1, kbot

      REAL XS(17), X
      REAL xslod(17)
      LOGICAL firstcall
      SAVE firstcall
      DATA firstcall /.TRUE./

      DATA xslod  /6.2180730E-21, 5.8473627E-22, 5.6996334E-22,
     $             4.5627094E-22, 1.7668250E-22, 1.1178808E-22,
     $             1.2040544E-22, 4.0994668E-23, 1.8450616E-23,
     $             1.5639540E-23, 8.7961075E-24, 7.6475608E-24,
     $             7.6260556E-24, 7.5565696E-24, 7.6334338E-24,
     $             7.4371992E-24, 7.3642966E-24/

c------------------------------------------
*sm	 Initialize cross sections to values
*sm	 at large optical depth
c------------------------------------------

      DO k = 1, nz
         DO i = 1, 17
            o2xsk(k,i) = xslod(i)
         ENDDO	
      ENDDO

c------------------------------------------
c      Loads Chebyshev polynomial Coeff.
c------------------------------------------

      IF (firstcall) THEN 
        CALL INIT_XS
	firstcall = .FALSE.
      ENDIF

c------------------------------------------
c     Calculate cross sections
*sm:  Set smallest O2col = exp(38.) molec cm-2
*sm     to stay in range of parameterization
*sm     given by Koppers et al. at top of atm.
c------------------------------------------

      ktop = 121
      kbot = 0

      DO k=1,nz    !! loop for alt
         o2col1(k) = MAX(o2col(k),EXP(38.))

         x  = ALOG(o2col1(k))
         
         IF (x .LT. 38.0) THEN
            ktop1 = k-1
            ktop  = MIN(ktop1,ktop)
         ELSE IF (x .GT. 56.0) THEN
            kbot = k
         ELSE
            CALL effxs( x, tlev(k), xs )
            DO i=1,17
               o2xsk(k,i) = xs(i)
            END DO
         ENDIF

      END DO                    !! finish loop for alt

c------------------------------------------
c  fill in cross section where X is out of range 
c  by repeating edge table values
c------------------------------------------

*sm do not allow kbot = nz to avoid division by zero in
*   no light case.
       
      IF(kbot .EQ. nz) kbot = nz - 1

      DO k=1,kbot
         DO i=1,17
            o2xsk(k,i) = o2xsk(kbot+1,i)
         END DO
      END DO
      
      DO k=ktop+1,nz
         DO i=1,17
            o2xsk(k,i) = o2xsk(ktop,i)
         END DO
      END DO

c------------------------------------------
c  Calculate incremental optical depths 
c------------------------------------------

      DO i=1,17                   ! loop over wavelength

         DO k=1,nz-1            ! loop for alt

c... calculate an optical depth weighted by density
*sm:  put in mean value estimate, if in shade

            IF (ABS(1. - o2col1(k+1)/o2col1(k)) .LE. 2.*precis) THEN

               dto2(k,i) = o2xsk(k+1,i)*o2col1(k+1)/(nz-1)

            ELSE

            dto2(k,i) = ABS(
     $           ( o2xsk(k+1,i)*o2col1(k+1) - o2xsk(k,i)*o2col1(k) )
     $           / ( 1. + ALOG(o2xsk(k+1,i)/o2xsk(k,i)) 
     $           / ALOG(o2col1(k+1)/o2col1(k)) ) )

c... change to vertical optical depth

            dto2(k,i) = 2. * dto2(k,i) / (secchi(k)+secchi(k+1))

            ENDIF

         END DO
         dto2(nz,i) = 0.0       ! set optical depth to zero at top


      END DO 

      RETURN
      END

*=============================================================================*

      SUBROUTINE EFFXS( X, T, XS )


C     Subroutine for evaluating the effective cross section
C     of O2 in the Schumann-Runge bands using parameterization
C     of G.A. Koppers, and D.P. Murtagh [ref. Ann.Geophys., 14
C     68-79, 1996]
C      
C     method:
C     ln(xs) = A(X)[T-220]+B(X)
C     X = log of slant column of O2
C     A,B calculated from Chebyshev polynomial coeffs
C     AC and BC using NR routine chebev.  Assume interval
C     is 38<ln(NO2)<56.
C
C     Revision History:
C
C     drm 2/97  initial coding
C
C-------------------------------------------------------------

	IMPLICIT NONE

	REAL*4 NO2, T, X
	REAL*4 XS(17)
	REAL*4 A(17), B(17) 
	INTEGER I

	CALL CALC_PARAMS( X, A, B )

	DO I = 1,17
	  XS(I) = EXP( A(I)*( T - 220.) + B(I) )
	ENDDO

        RETURN

	END

*=============================================================================*

	SUBROUTINE CALC_PARAMS( X, A, B )

C-------------------------------------------------------------
C
C       calculates coefficients (A,B), used in calculating the
C	effective cross section, for 17 wavelength intervals
C       as a function of log O2 column density (X)
C       Wavelength intervals are defined in WMO1985
C
C-------------------------------------------------------------

	IMPLICIT NONE

	REAL*4 X
	REAL*4 A(17), B(17)

	REAL*4   CHEBEV

	REAL*8 AC(20,17)
        REAL*8 BC(20,17) ! Chebyshev polynomial coeffs
	REAL*4 WAVE_NUM(17)
	COMMON /XS_COEFFS/ AC, BC, WAVE_NUM

	INTEGER I

C       call Chebyshev Evaluation routine to calc A and B from
C	set of 20 coeficients for each wavelength

	DO I=1,17
	  A(I) = CHEBEV(38.0 , 56.0, AC(1,I), 20, X)
	  B(I) = CHEBEV(38.0 , 56.0, BC(1,I), 20, X)
	ENDDO

	RETURN

	END

*=============================================================================*

	SUBROUTINE INIT_XS

C-------------------------------------------------------------
C       loads COMMON block XS_COEFFS containing the Chebyshev
C	polynomial coeffs necessary to calculate O2 effective
C       cross-sections
C
C-------------------------------------------------------------
	REAL*8 AC(20,17)
	REAL*8 BC(20,17) ! Chebyshev polynomial coeffs
	REAL*4 WAVE_NUM(17)
	COMMON /XS_COEFFS/ AC, BC, WAVE_NUM
	

C       locals
	INTEGER*4 IN_LUN	! file unit number
	INTEGER*4 IOST		! i/o status
	INTEGER*4 I, J

        IN_LUN = 11

	OPEN (UNIT=IN_LUN, FILE=
     $       'DATAE1/O2/effxstex.txt',FORM='FORMATTED')

	READ( IN_LUN, 901 )
	DO I = 1,20
	  READ( IN_LUN, 903 ) ( AC(I,J), J=1,17 )
	ENDDO
	READ( IN_LUN, 901 )
	DO I = 1,20
	  READ( IN_LUN, 903 ) ( BC(I,J), J=1,17 )
	ENDDO

 901    FORMAT( / )
 903    FORMAT( 17(E23.14,1x))

 998	CLOSE (IN_LUN)
	
	DO I=1,17
	  WAVE_NUM(18-I) = 48250. + (500.*I)
	ENDDO

        END

*=============================================================================*

	FUNCTION chebev(a,b,c,m,x)

C-------------------------------------------------------------
C
C     Chebyshev evaluation algorithm
C     See Numerical recipes p193
C
C-------------------------------------------------------------
      
	INTEGER M
        REAL*4 CHEBEV,A,B,X
	REAL*8 C(M)
        INTEGER J
        REAL D,DD,SV,Y,Y2

        IF ((X-A)*(X-B).GT.0.) THEN
	  WRITE(6,*) 'X NOT IN RANGE IN CHEBEV', X
	  CHEBEV = 0.0
	  RETURN
        ENDIF

	D=0.
        DD=0.
        Y=(2.*X-A-B)/(B-A)
        Y2=2.*Y
        DO 11 J=M,2,-1
          SV=D
          D=Y2*D-DD+C(J)
          DD=SV
 11     CONTINUE
        CHEBEV=Y*D-DD+0.5*C(1)
      
	RETURN
        END

*=============================================================================*

       SUBROUTINE sjo2(nz,nw,xso2,nj,sq)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Update the weighting function (cross section x quantum yield) for O2     =*
*=  photolysis.  The strong spectral variations in the O2 cross sections are =*
*=  parameterized into a few bands for Lyman-alpha (121.4-121.9 nm, one band)=*
*=  and Schumann-Runge (174.4-205.8, 17 bands) regions. The parameterizations=*
*=  depend on the overhead O2 column, and therefore on altitude and solar    =*
*=  zenith angle, so they need to be updated at each time/zenith step.       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  XSO2   - REAL, molecular absorption cross section in SR bands at      (I)=*
*=           each specified altitude and wavelength.  Includes Herzberg      =*
*=            continuum.                                                     =*
*=  NJ     - INTEGER, index of O2 photolysis in array SQ                  (I)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction, at each wavelength and each altitude level =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* calling parameters

      INTEGER nz, nw, nj
      REAL xso2(kz,kw)
      REAL sq(kj,kz,kw)

* local

      INTEGER iw, iz
*______________________________________________________________________________

* O2 + hv -> O + O
* quantum yield assumed to be unity
* assign cross section values at all wavelengths and at all altitudes
*      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(nj,iz,iw) = xso2(iz,iw)
        ENDDO
      ENDDO
*______________________________________________________________________________


      RETURN
      END
