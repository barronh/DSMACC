         SUBROUTINE r01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product of (cross section) x (quantum yield) for the two     =*
*=  O3 photolysis reactions:                                                 =*
*=             (a) O3 + hv -> O2 + O(1D)                                     =*
*=             (b) O3 + hv -> O2 + O(3P)                                     =*
*=  Cross section:  Combined data from WMO 85 Ozone Assessment (use 273K     =*
*=                  value from 175.439-847.5 nm) and data from Molina and    =*
*=                  Molina (use in Hartley and Huggins bans (240.5-350 nm)   =*
*=  Quantum yield:  Choice between                                           =*
*=                   (1) data from Michelsen et al, 1994                     =*
*=                   (2) JPL 87 recommendation                               =*
*=                   (3) JPL 90/92 recommendation (no "tail")                =*
*=                   (4) data from Shetter et al., 1996                      =*
*=                   (5) JPL 97 recommendation                               =*
*=                   (6) JPL 00 recommendation                               =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER n1, n2, n3, n4, n5
      INTEGER kdata
      PARAMETER (kdata = 500)
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata)

* local

      INTEGER mabs
      REAL xs(kz,kw)

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy1d, qy3p
      REAL tau, tau2, tau3
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL xl, xl0
      REAL so3
      REAL dum
      INTEGER myld
      INTEGER kmich, kjpl87, kjpl92, kshet, kjpl97, kjpl00, kmats
      INTEGER i, iw, n, idum
      INTEGER ierr

      REAL fo3qy, fo3qy2
      EXTERNAL fo3qy, fo3qy2

****************************************************************

*************       jlabel(j) = 'O3 -> O2 + O(1D)'
*************       jlabel(j) = 'O3 -> O2 + O(3P)'

      j = j + 1
      jlabel(j) = 'O3 -> O2 + O(1D)'
      
      j = j + 1
      jlabel(j) = 'O3 -> O2 + O(3P)'

* call cross section read/interpolate routine
* cross sections from WMO 1985 Ozone Assessment
* from 175.439 to 847.500 nm. Using value at 273 K.
* Values are over-written in Hartly and Huggins bands, using different
* options depending on value of mopt:
*     mabs = 1 = Molina and Molina
*     mabs = 2 = Malicet et al.
*     mabs = 3 = Bass et al.

      mabs = 1
      IF(mabs. EQ. 1) CALL o3xs_mm(nw,wl,nz,tlev, xs)
      IF(mabs. EQ. 2) CALL o3xs_mal(nw,wl,nz,tlev, xs)
      IF(mabs. EQ. 3) CALL o3xs_bass(nw,wl,nz,tlev, xs)

******* quantum yield:

      kmich = 1
      kjpl87 = 2
      kjpl92 = 3
      kshet = 4
      kjpl97 = 5
      kjpl00 = 6
      kmats = 7

* choose quantum yield recommendation:
*    kjpl87:  JPL recommendation 1987                - JPL 87, 90, 92 do not "tail"
*    kjpl92:  JPL recommendations 1990/92 (identical) - still with no "tail"
*    kjpl97:  JPL recommendation 1997, includes tail, similar to Shetter et al.
*    kmich :  Michelsen et al., 1994
*    kshet :  Shetter et al., 1996
*    kjpl00:  JPL 2000
*    kmats:  Matsumi et al., 2002

c      myld = kjpl87
c      myld = kjpl92
c      myld = kshet
c      myld = kmich
c      myld = kjpl97
c      myld = kjpl00

      myld = kmats

* read parameters from JPL'97

      IF (myld .EQ. kjpl97) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/YLD/O3.param_jpl97.yld',STATUS='old')
        READ(kin,*)
        READ(kin,*)
        READ(kin,*)
        n1 = 21
        n2 = n1
        DO i = 1, n1
           READ(kin,*) x1(i), y1(i), y2(i)
           x2(i) = x1(i)
        ENDDO
        CLOSE(kin)

        CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
        CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
        CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),y1(n1))
        CALL addpnt(x1,y1,kdata,n1,            1.e+38,y1(n1))
        CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF

        CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
        CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
        CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
        CALL addpnt(x2,y2,kdata,n2,            1.e+38,y2(n2))
        CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF
      ENDIF

* read parameters from Michelsen, H. A., R.J. Salawitch, P. O. Wennber, 
* and J. G. Anderson, Geophys. Res. Lett., 21, 2227-2230, 1994.

      IF (myld .EQ. kmich) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/YLD/O3.param.yld',STATUS='old')
        READ(kin,*)
        READ(kin,*)
        READ(kin,*)
        n1 = 21
        n2 = n1
        DO i = 1, n1
           READ(kin,*) x1(i), y1(i), y2(i)
           x2(i) = x1(i)
        ENDDO
        CLOSE(kin)

        CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
        CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
        CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),y1(n1))
        CALL addpnt(x1,y1,kdata,n1,            1.e+38,y1(n1))
        CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF

        CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
        CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
        CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
        CALL addpnt(x2,y2,kdata,n2,            1.e+38,y2(n2))
        CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF
      ENDIF

* quantum yield data from 
* Shetter et al, J.Geophys.Res., v 101 (D9), pg. 14,631-14,641, June 20, 1996

      IF (myld .EQ. kshet) THEN
        OPEN(UNIT=kin,FILE='DATAJ1/YLD/O3_shetter.yld',STATUS='OLD')
        READ(kin,*) idum, n
        DO i = 1, idum-2
          READ(kin,*)
        ENDDO
        n = n-2
        DO i = 1, n
          READ(kin,*) x1(i),y3(i),y4(i),y1(i),y2(i)
          x2(i) = x1(i)
          x3(i) = x1(i)
          x4(i) = x1(i)
        ENDDO
        DO i = n+1, n+2
           READ(kin,*) x3(i),y3(i),y4(i)
           x4(i) = x3(i)
        ENDDO
        CLOSE(kin)

        n1 = n
        n2 = n
        n3 = n+2
        n4 = n+2

* coefficients for exponential fit:

        CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax), y1(1))
        CALL addpnt(x1,y1,kdata,n1,                0., y1(1))
        CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
        CALL addpnt(x1,y1,kdata,n1,              1E38,0.)

        CALL inter2(nw,wl,yg1, n1,x1,y1, ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF

        CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
        CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
        CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
        CALL addpnt(x2,y2,kdata,n2,              1E38,0.)

        CALL inter2(nw,wl,yg2, n2,x2,y2, ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF

* phi data at 298 and 230 K

        CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),y3(1))
        CALL addpnt(x3,y3,kdata,n3,               0.,y3(1))
        CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
        CALL addpnt(x3,y3,kdata,n3,              1E38,0.)

        CALL inter2(nw,wl,yg3, n3,x3,y3, ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr,jlabel(j)
           STOP
        ENDIF

        CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),y4(1))
        CALL addpnt(x4,y4,kdata,n4,               0.,y4(1))
        CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
        CALL addpnt(x4,y4,kdata,n4,              1E38,0.)

        CALL inter2(nw,wl,yg4, n4,x4,y4, ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr,jlabel(j)
           STOP
        ENDIF
      ENDIF

* compute cross sections and yields at different wavelengths, altitudes:

      DO 10 iw = 1, nw-1

         DO 20 i = 1, nz

* quantum yields
* coefficients from jpl 87:

             IF (myld .EQ. kjpl87) THEN
               tau = tlev(i) - 230.
               tau2 = tau*tau
               tau3 = tau2*tau
               xl = wc(iw)
               xl0 = 308.2 + 4.4871e-2*tau + 6.938e-5*tau2 -
     >               2.5452e-6*tau3
               a = 0.9*(0.369 + 2.85e-4*tau + 1.28e-5*tau2 + 
     >                  2.57e-8*tau3)
               b     = -0.575 + 5.59e-3*tau - 1.439e-5*tau2 - 
     >                  3.27e-8*tau3
               c = 0.9*(0.518 + 9.87e-4*tau - 3.94e-5*tau2 + 
     >                  3.91e-7*tau3)
               qy1d = a*atan(b*(xl-xl0)) + c
               qy1d = amax1(0.,qy1d)
               qy1d = amin1(0.9,qy1d)
             ENDIF

* from jpl90, jpl92:
* (caution: error in JPL92 for first term of a3)

             IF (myld .EQ. kjpl92) THEN
               tau = 298. - tlev(i)
               tau2 = tau*tau
               xl0 = wc(iw) - 305.
               a0 =   .94932   - 1.7039e-4*tau + 1.4072E-6*tau2
               a1 = -2.4052e-2 + 1.0479e-3*tau - 1.0655e-5*tau2
               a2 =  1.8771e-2 - 3.6401e-4*tau - 1.8587e-5*tau2
               a3 = -1.4540e-2 - 4.7787e-5*tau + 8.1277e-6*tau2
               a4 =  2.3287e-3 + 1.9891e-5*tau - 1.1801e-6*tau2
               a5 = -1.4471e-4 - 1.7188e-6*tau + 7.2661e-8*tau2
               a6 =  3.1830e-6 + 4.6209e-8*tau - 1.6266e-9*tau2
               qy1d = a0 + a1*xl0 + a2*(xl0)**2 + a3*(xl0)**3 +
     >                a4*(xl0)**4 + a5*(xl0)**5 + a6*(xl0)**6
               IF (wc(iw) .LT. 305.) qy1d = 0.95
               IF (wc(iw) .GT. 320.) qy1d = 0.
               IF (qy1d .LT. 0.02) qy1d = 0.
             ENDIF

* from JPL'97

           IF (myld .EQ. kjpl97) THEN
             IF (wc(iw) .LT. 271.) THEN
                qy1d = 0.87
             ELSE IF (wc(iw) .GE. 271. .AND. wc(iw) .LT. 290.) THEN
                qy1d = 0.87 + (wc(iw)-271.)*(.95-.87)/(290.-271.)
             ELSE IF (wc(iw) .GE. 290. .AND. wc(iw) .LT. 305.) THEN
                qy1d = 0.95
             ELSE IF (wc(iw) .GE. 305. .AND. wc(iw) .LE. 325.) THEN
                qy1d = yg1(iw) * EXP ( -yg2(iw) /tlev(i) )
             ELSE
                qy1d = 0.
             ENDIF
           ENDIF
 
* from Michelsen, H. A., R.J. Salawitch, P. O. Wennber, and J. G. Anderson
* Geophys. Res. Lett., 21, 2227-2230, 1994.

           IF (myld .EQ. kmich) THEN
             IF (wc(iw) .LT. 271.) THEN
                qy1d = 0.87
             ELSE IF (wc(iw) .GE. 271. .AND. wc(iw) .LT. 305.) THEN
                qy1d = 1.98 - 301./wc(iw)
             ELSE IF (wc(iw) .GE. 305. .AND. wc(iw) .LE. 325.) THEN
                qy1d = yg1(iw) * EXP (-yg2(iw) /(0.6951*tlev(i)))
             ELSE
                qy1d = 0.
             ENDIF
           ENDIF
 
* Shetter et al.:
* phi = A * exp(-B/T), A and B are based on meas. at 298 and 230 K
* do linear interpolation between phi(298) and phi(230) for wavelengths > 321
* as phi(230)=0. for those wavelengths, so there are no A and B factors

           IF (myld .EQ. kshet) THEN
             IF (wl(iw+1) .LE. 321.) THEN
               qy1d = yg1(iw) * EXP(-1. * yg2(iw)/tlev(i))
             ELSE
               qy1d = (yg3(iw) - yg4(iw))/(298.-230.) * (tlev(i)-230.) +
     >                 yg4(iw)
             ENDIF
           ENDIF

* JPL 2000:

           IF (myld .EQ. kjpl00) THEN
              qy1d = fo3qy(wc(iw),tlev(i))
           ENDIF

* Matsumi et al.

           IF (myld .EQ. kmats) THEN
              qy1d = fo3qy2(wc(iw),tlev(i))
           ENDIF

* compute product

           sq(j-1,i,iw) = qy1d*xs(i,iw)
           qy3p = 1.0 - qy1d
           sq(j,i,iw) = qy3p*xs(i,iw)

 20     CONTINUE
 10   CONTINUE

      RETURN
      END

* This file contains the following subroutines, related to reading/loading
* the product (cross section) x (quantum yield) for photo-reactions:
*     r01 through r47
*     r101 through r110
*=============================================================================*

*=============================================================================*

      SUBROUTINE r02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for NO2            =*
*=  photolysis:                                                              =*
*=         NO2 + hv -> NO + O(3P)                                            =*
*=  Cross section from JPL94 (can also have Davidson et al.)                 =*
*=  IUPAC 2003 (data of Merienne, Coquart and Jenouvrier)                    =*
*=    added by Judit Zador												   =*
*=  Quantum yield from Gardiner, Sperry, and Calvert                         =*
*=  IUPAC 2003, evaluation of Troe, added by Judit Zador					   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*


      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER n1
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL no2xs(kz,kw)
      REAL dum
      INTEGER i, iw, n, idum, ierr
      integer mabs


**************** NO2 photodissociation

      j = j + 1
      jlabel(j) = 'NO2 -> NO + O(3P)'

* cross section

*------------NEED TO CHANGE kdata = 1000 FOR DAVIDSON ET AL. DATA---------
* measurements by:
* Davidson, J. A., C. A. Cantrell, A. H. McDaniel, R. E. Shetter,
* S. Madronich, and J. G. Calvert, Visible-ultraviolet absorption
* cross sections for NO2 as a function of temperature, J. Geophys.
* Res., 93, 7105-7112, 1988.
*     from 263.8 to 648.8 nm in approximately 0.5 nm intervals
C     OPEN(UNIT=kin,FILE='DATAE1/NO2/NO2_ncar_00.abs',STATUS='old')
C     n = 750
C     DO i = 1, n
C        READ(kin,*) x1(i), y1(i), dum, dum, idum
C     ENDDO
C     CLOSE(kin)

C     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,               0.,0.)
C     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
C     CALL inter2(nw,wl,yg,n,x1,y1,ierr)
C     IF (ierr .NE. 0) THEN
C        WRITE(*,*) ierr, jlabel(j)
C        STOP
C     ENDIF

* options:
*     mabs = 1    NO2_jpl94.abs  (same as JPL97)
*     mabs = 2    HarNO2CS.rxt  from Harder et al.
*     mabs = 3    IUPAC 2003 November (Zador)


      mabs = 3

* cross section data from JPL 94 recommendation
* JPL 97 recommendation is identical

      if (mabs .eq. 1) then

         OPEN(UNIT=kin,FILE='DATAE1/NO2/NO2_jpl94.abs',STATUS='old')
         READ(kin,*) idum, n
         DO i = 1, idum-2
            READ(kin,*)
         ENDDO 

* read in wavelength bins, cross section at T0 and temperature correction
* coefficient a;  see input file for details.
* data need to be scaled to total area per bin so that they can be used with
* inter3

         DO i = 1, n
            READ(kin,*) x1(i), x3(i), y1(i), dum, y2(i)
            y1(i) = (x3(i)-x1(i)) * y1(i)*1.E-20
            y2(i) = (x3(i)-x1(i)) * y2(i)*1.E-22
            x2(i) = x1(i) 
         ENDDO
         CLOSE(kin)

         x1(n+1) = x3(n)
         x2(n+1) = x3(n)
         n = n+1
         n1 = n

         CALL inter3(nw,wl,yg1,n,x1,y1,0)
         CALL inter3(nw,wl,yg2,n1,x2,y2,0)

* yg1, yg2 are per nm, so rescale by bin widths

         DO iw = 1, nw-1
            yg1(iw) = yg1(iw)/(wl(iw+1)-wl(iw))
            yg2(iw) = yg2(iw)/(wl(iw+1)-wl(iw))
         ENDDO

         DO iw = 1, nw-1
            DO i = 1, nz
               no2xs(i,iw) = yg1(iw) + yg2(iw)*(tlev(i)-273.15)
            ENDDO
         ENDDO 

      elseif (mabs .eq. 2) then

         OPEN(UNIT=kin,FILE='DATAE1/NO2/NO2_Har.abs',status='old')
         DO i = 1, 9
            READ(kin,*)
         ENDDO
         n = 135
         DO i = 1, n
            READ(kin,*) idum, y1(i)
            x1(i) = FLOAT(idum)
         enddo

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,   0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)

         DO iw = 1, nw-1
            DO i = 1, nz
               no2xs(i,iw) = yg1(iw)
            ENDDO
         ENDDO 

      elseif (mabs .eq. 3) then

         OPEN(UNIT=kin,FILE='DATAE1/NO2/NO2_IUPAC2003.abs',status='old')
         DO i = 1, 5
            READ(kin,*)
         ENDDO
         n = 59
         DO i = 1, n
            READ(kin,*) idum, y1(i)
            x1(i) = FLOAT(idum)
	      y1(i) = y1(i)/1e20
         enddo

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,   0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)

         DO iw = 1, nw-1
            DO i = 1, nz
               no2xs(i,iw) = yg1(iw)
            ENDDO
         ENDDO 

      ENDIF



* quantum yield
* from Gardiner, Sperry, and Calvert

      OPEN(UNIT=kin,FILE='DATAJ1/YLD/NO2_calvert.yld',STATUS='old')
      DO i = 1, 8
         READ(kin,*) 
      ENDDO
      n = 66
      DO i = 1, n
         READ(kin,*) x1(i),y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,   0.)
      CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* combine

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = no2xs(i,iw)*yg1(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (absorptioon cross section) x (quantum yield) for    =*
*=  both channels of NO3 photolysis:                                         =*
*=          (a) NO3 + hv -> NO2 + O(3P)                                      =*
*=          (b) NO3 + hv -> NO + O2                                          =*
*=  Cross section combined from Graham and Johnston (<600 nm) and JPL 94     =*
*=  Quantum yield from Madronich (1988)                                      =*
*=             or from Johnston (1996) added by Judit Zador                  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=350)

      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw)
      REAL qy
      INTEGER irow, icol
      INTEGER i, iw, n, idum
      INTEGER ierr
	INTEGER myld, n1, n2

****************      jlabel(j) = 'NO3 -> NO2 + O(3P)'
****************      jlabel(j) = 'NO3 -> NO + O2'

	myld = 2

* cross section
*     measurements of Graham and Johnston 1978

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/NO3_gj78.abs',STATUS='old')
      DO i = 1, 9
         READ(kin,*)
      ENDDO
      n = 305
      DO irow = 1, 30
         READ(kin,*) ( y1(10*(irow-1) + icol), icol =  1, 10 )
      ENDDO
      READ(kin,*) ( y1(300 + icol), icol = 1, 5 )
      CLOSE (kin)
      DO i = 1, n
         y1(i) =  y1(i) * 1.E-19
         x1(i) = 400. + 1.*FLOAT(i-1)
      ENDDO

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

*     cross section from JPL94:

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/NO3_jpl94.abs',STATUS='old')
      READ(kin,*) idum, n
      DO i = 1, idum-2
         READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i)*1E-20
      ENDDO 
      CLOSE (kin)
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* use JPL94 for wavelengths longer than 600 nm

      DO iw = 1, nw-1
         IF(wl(iw) .GT. 600.) yg(iw) = yg1(iw)
      ENDDO

* quantum yield:
* from Madronich (1988) see CEC NO3 book.

	IF (myld.EQ.1) THEN

* for   NO3 ->NO+O2

		j = j + 1
		jlabel(j) = 'NO3 -> NO + O2'
		DO iw = 1, nw - 1
		   IF (wc(iw).LT.584.) THEN 
			  qy = 0.
		   ELSEIF (wc(iw).GE.640.) THEN
			  qy = 0.
		   ELSEIF (wc(iw).GE.595.) THEN 
			  qy = 0.35*(1.-(wc(iw)-595.)/45.)
		   ELSE
			  qy = 0.35*(wc(iw)-584.)/11.
		   ENDIF
		   DO i = 1, nz
			  sq(j,i,iw) = yg(iw)*qy
		   ENDDO
		ENDDO

* for  NO3 ->NO2+O
		j = j + 1
		jlabel(j) = 'NO3 -> NO2 + O(3P)'
		DO iw = 1, nw - 1
		   IF (wc(iw).LT.584.) THEN
			  qy = 1.
		   ELSEIF (wc(iw).GT.640.) THEN
			  qy = 0.
		   ELSEIF (wc(iw).GT.595.) THEN
			  qy = 0.65*(1-(wc(iw)-595.)/45.)
		   ELSE
			  qy = 1.-0.35*(wc(iw)-584.)/11.
		   ENDIF
		   DO i = 1, nz
			  sq(j,i,iw) = yg(iw)*qy
		   ENDDO
		ENDDO
	
	ELSEIF (myld.EQ.2) THEN

	   OPEN(UNIT=kin,FILE='DATAJ1/YLD/NO3_Johnston96.yld',
     &     STATUS='old') 
         READ(kin,*) 
         READ(kin,*) 
         READ(kin,*) 
         n1 = 56 
         n2 = n1 
         DO i = 1, n1 
            READ(kin,*) x1(i), y1(i), y2(i) 
            x2(i) = x1(i) 
         ENDDO 
         CLOSE(kin) 

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1)) 
         CALL addpnt(x1,y1,kdata,n1,               0.,y1(1)) 
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),y1(n1))
         CALL addpnt(x1,y1,kdata,n1,            1.e+38,y1(n1))
         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr) 
         IF (ierr .NE. 0) THEN 
            WRITE(*,*) ierr, jlabel(j) 
            STOP 
         ENDIF 

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1)) 
         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1)) 
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
         CALL addpnt(x2,y2,kdata,n2,            1.e+38,y2(n2))
         CALL inter2(nw,wl,yg2,n2,x2,y2,ierr) 
         IF (ierr .NE. 0) THEN 
            WRITE(*,*) ierr, jlabel(j) 
            STOP 
         ENDIF 

* combine
	   j = j + 1 
	   jlabel(j) = 'NO3 -> NO + O2' 
         DO iw = 1, nw - 1 
            DO i = 1, nz 
	         IF (wc(iw).LT.584.) THEN 
	            qy = 0. 
               ELSEIF (wc(iw).GT.640.) THEN 
                  qy = 0. 
               ELSE 
                  qy = yg1(iw)/1000
               ENDIF
                  sq(j,i,iw) = yg(iw)*qy
            ENDDO
         ENDDO 

	   j = j + 1 
	   jlabel(j) = 'NO3 -> NO2 + O(3P)' 
         DO iw = 1, nw - 1 
            DO i = 1, nz 
	         IF (wc(iw).LT.584.) THEN 
	            qy = 1. 
               ELSEIF (wc(iw).GT.640.) THEN 
                  qy = 0. 
               ELSE 
                  qy = yg2(iw)/1000
               ENDIF
                  sq(j,i,iw) = yg(iw)*qy
            ENDDO
         ENDDO 

	ENDIF
	
      END

*=============================================================================*

      SUBROUTINE r04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yiels) for N2O5 photolysis =*
*=  reactions:                                                               =*
*=       (a) N2O5 + hv -> NO3 + NO + O(3P)                                   =*
*=       (b) N2O5 + hv -> NO3 + NO2                                          =*
*=  Cross section from JPL97: use tabulated values up to 280 nm, use expon.  =*
*=                            expression for >285nm, linearly interpolate    =*
*=                            between s(280) and s(285,T) in between         =*
*=  Quantum yield: Analysis of data in JPL94 (->DATAJ1/YLD/N2O5.qy)          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      REAL xs, xs270, xs280, xst290, xst300
      REAL dum1, dum2
      REAL t
      INTEGER i, iw, n, idum
      INTEGER ierr

**************** N2O5 photodissociation

      j = j + 1
      jlabel(j) = 'N2O5 -> NO3 + NO + O(3P)'

      j = j + 1
      jlabel(j) = 'N2O5 -> NO3 + NO2'

* cross section from jpl97, table up to 280 nm

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/N2O5_jpl97.abs',STATUS='old')
      READ(kin,*) idum, n
      DO i = 1, idum-2
         READ(kin,*)
      ENDDO
      DO i = 1,  n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      xs270 = y1(n-2)
      xs280 = y1(n)

      CLOSE(kin)

      CALL addpnt(x1,y1,kdata, n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata, n,               0.,0.)
      CALL addpnt(x1,y1,kdata, n,x1(n)*(1.+deltax),y1(n))
      CALL addpnt(x1,y1,kdata, n,            1.E36,y1(n))

      CALL inter2(nw,wl,yg, n,x1,y1, ierr)
      IF (ierr .NE. 0) THEN
         WRITE(0,*) ierr,jlabel(j)
         STOP
      ENDIF 

* quantum yield : see DATAJ1/YLD/N2O5.qy for explanation
* correct for T-dependence of cross section

      DO iw = 1, nw - 1

         qy = MIN( 1., 3.832441 - 0.012809638 * wc(iw) )
         qy = MAX( 0., qy )

         DO i = 1, nz

* temperature dependence only valid for 225 - 300 K.

            t = MAX(225.,MIN(tlev(i),300.))

* evaluation of exponential

            IF (wc(iw) .GE. 285. .AND. wc(iw) .LE. 380.) THEN

               sq(j-1,i,iw) = qy *
     $            1.E-20*EXP( 2.735 + (4728.5-17.127*wc(iw)) / t )
               sq(j,i,iw) = (1.-qy) * 
     $            1.E-20*EXP( 2.735 + (4728.5-17.127*wc(iw)) / t )

* between 280 and 285 nm:  Extrapolate from both sides, then average.

            ELSEIF (wc(iw) .GE. 280. .AND. wc(iw) .LT. 285.) THEN

               xst290 = 1.E-20*
     >              EXP( 2.735 + (4728.5-17.127*290.) / t )
               xst300 = 1.E-20*
     >              EXP( 2.735 + (4728.5-17.127*300.) / t )
               
               dum1 = xs270 + (wc(iw) - 270.)*(xs280 - xs270)/10.
               dum2 = xst290 + (wc(iw) - 290.)*(xst300 - xst290)/10.
               xs = 0.5*(dum1 + dum2)

               sq(j-1,i,iw) = qy * xs
               sq(j,i,iw) = (1.-qy) * xs

* for less than 280 nm, use tabulated values

            ELSEIF (wc(iw) .LT. 280.) THEN

               sq(j-1,i,iw) = qy * yg(iw)
               sq(j,i,iw) = (1.-qy) * yg(iw)

* beyond 380 nm, set to zero

            ELSE
               sq(j-1,i,iw) = 0.
               sq(j,i,iw) = 0.
            ENDIF
         ENDDO
      ENDDO
      END

*=============================================================================*

      SUBROUTINE r05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for HNO2 photolysis=*
*=     HNO2 + hv -> NO + OH                                                  =*
*=  Cross section:  (1) from JPL97 (Bongartz 1991)                           =*
*=                  (2) from IUPAC 2003 (Bongartz 1994) added by Judit Zador =*
*=  Quantum yield:  assumed to be unity                                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=1010)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n
      INTEGER ierr
	INTEGER mabs

**************** HNO2 photodissociation
* cross section from JPL92
* mabs = 1 from Bongartz et al., identical to JPL94, JPL97 recommendation
* mabs = 2 from the original paper, recommendation of IUPAC 2003

	mabs = 2

	IF (mabs.eq.1) THEN

		j = j + 1
		jlabel(j) = 'HNO2 -> OH + NO'
		OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO2_jpl92.abs',STATUS='old')
		DO i = 1, 13
		   READ(kin,*)
		ENDDO
		n = 91
		DO i = 1, n
		   READ(kin,*) x1(i), y1(i)
		   y1(i) = y1(i) * 1.E-20
		ENDDO
		CLOSE (kin)

		CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
		CALL addpnt(x1,y1,kdata,n,               0.,0.)
		CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
		CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
		CALL inter2(nw,wl,yg,n,x1,y1,ierr)
		IF (ierr .NE. 0) THEN
		   WRITE(*,*) ierr, jlabel(j)
		   STOP
		ENDIF

	ELSE

		j = j + 1
		jlabel(j) = 'HNO2 -> OH + NO'
		OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO2_bongartz94.abs',
     >		STATUS='old')
		DO i = 1, 7
		   READ(kin,*)
		ENDDO
		n = 998
		DO i = 1, n
		   READ(kin,*) x1(i), y1(i)
		ENDDO
		CLOSE (kin)

		CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
		CALL addpnt(x1,y1,kdata,n,               0.,0.)
		CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
		CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
		CALL inter2(nw,wl,yg,n,x1,y1,ierr)
		IF (ierr .NE. 0) THEN
		   WRITE(*,*) ierr, jlabel(j)
		   STOP
		ENDIF


	ENDIF

* quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HNO3 photolysis =*
*=        HNO3 + hv -> OH + NO2                                              =*
*=  Cross section: Burkholder et al., 1993                                   =*
*=  Quantum yield: Assumed to be unity                                       =*
*=              or combined qy from IUPAC 2003 added by Judit Zador          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      INTEGER i, iw
      INTEGER ierr

**************** HNO3 photodissociation

       j = j + 1
       jlabel(j) = 'HNO3 -> OH + NO2'

C* cross section from JPL85
C
C      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO3.abs',STATUS='old')
C      DO i = 1, 9
C         READ(kin,*)
C      ENDDO
C      n = 29
C      DO i = 1, n
C         READ(kin,*) x1(i), y1(i)
C         y1(i) = y1(i) * 1.E-20
C      ENDDO
C      CLOSE (kin)
C
C      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C      CALL addpnt(x1,y1,kdata,n,               0.,0.)
C      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
C      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
C      IF (ierr .NE. 0) THEN
C         WRITE(*,*) ierr, jlabel(j)
C         STOP
C      ENDIF
C
C* quantum yield = 1
C
C      qy = 1.
C      DO iw = 1, nw - 1
C         DO i = 1, nz
C            sq(j,i,iw) = yg(iw)*qy
C         ENDDO
C      ENDDO


* HNO3 cross section parameters from Burkholder et al. 1993

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO3_burk.abs',STATUS='old')
      DO i = 1, 6
         READ(kin,*)
      END DO
      n1 =  83
      n2 = n1
      DO i = 1, n1
         READ(kin,*) y1(i), y2(i)
         x1(i) = 184. + i*2.
         x2(i) = x1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,            1.e+38,0.)
      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF


      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
      CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,y2(n2))
      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield = 1
* correct for temperature dependence

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg1(iw) * 1.E-20
     $           * exp( yg2(iw)/1.e3*(tlev(i)-298.) )
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HNO4 photolysis =*
*=       HNO4 + hv -> HO2 + NO2                                              =*
*=  Cross section:  from JPL97                                               =*
*=  Quantum yield:  Assumed to be unity                                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)
 
C* local
 
      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n
      INTEGER ierr

**************** HNO4 photodissociation

* cross section from JPL85 (identical to JPL92 and JPL94 and JPL97)

      j = j + 1
      jlabel(j) = 'HNO4 -> HO2 + NO2'
C      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO4.abs',STATUS='old')
      OPEN(UNIT=kin,FILE='DATAJ1/ABS/HNO4_jpl92.abs',STATUS='old')
      DO i = 1, 4
         READ(kin,*)
      ENDDO
      n = 31
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for H2O2 photolysis =*
*=         H2O2 + hv -> 2 OH                                                 =*
*=  Cross section:  From JPL97, tabulated values @ 298K for <260nm, T-depend.=*
*=                  parameterization for 260-350nm                           =*
*=  Quantum yield:  Assumed to be unity                                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

C     INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL xs
      REAL t
      INTEGER i, iw, n, idum
      INTEGER ierr
      REAL lambda
      REAL sumA, sumB, chi

**************** H2O2 photodissociation

* cross section from Lin et al. 1978

      j = j + 1
      jlabel(j) = 'H2O2 -> 2 OH'
C     OPEN(UNIT=kin,FILE='DATAJ1/ABS/H2O2_lin.abs',STATUS='old')
C     DO i = 1, 7
C        READ(kin,*)
C     ENDDO
C     n = 32
C     DO i = 1, n
C        READ(kin,*) x1(i), y1(i)
C        y1(i) = y1(i) * 1.E-20
C     ENDDO
C     CLOSE (kin)
C
C      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C      CALL addpnt(x1,y1,kdata,n,               0.,0.)
C      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
C      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
C      IF (ierr .NE. 0) THEN
C         WRITE(*,*) ierr, jlabel(j)
C         STOP
C      ENDIF

* cross section from JPL94 (identical to JPL97)
* tabulated data up to 260 nm

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/H2O2_jpl94.abs',STATUS='old')
      READ(kin,*) idum,n
      DO i = 1, idum-2
         READ(kin,*)
      ENDDO
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE (kin)


      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      A0 = 6.4761E+04            
      A1 = -9.2170972E+02        
      A2 = 4.535649              
      A3 = -4.4589016E-03        
      A4 = -4.035101E-05         
      A5 = 1.6878206E-07
      A6 = -2.652014E-10
      A7 = 1.5534675E-13

      B0 = 6.8123E+03
      B1 = -5.1351E+01
      B2 = 1.1522E-01
      B3 = -3.0493E-05
      B4 = -1.0924E-07

* quantum yield = 1

      qy = 1.

      DO iw = 1, nw - 1

* Parameterization (JPL94)
* Range 260-350 nm; 200-400 K

         IF ((wl(iw) .GE. 260.) .AND. (wl(iw) .LT. 350.)) THEN

           lambda = wc(iw)
           sumA = ((((((A7*lambda + A6)*lambda + A5)*lambda + 
     >                  A4)*lambda +A3)*lambda + A2)*lambda + 
     >                  A1)*lambda + A0
           sumB = (((B4*lambda + B3)*lambda + B2)*lambda + 
     >               B1)*lambda + B0

           DO i = 1, nz
              t = MIN(MAX(tlev(i),200.),400.)            
              chi = 1./(1.+EXP(-1265./t))
              xs = (chi * sumA + (1.-chi)*sumB)*1E-21
              sq(j,i,iw) = xs*qy
           ENDDO
         ELSE
           DO i = 1, nz
              sq(j,i,iw) = yg(iw)*qy
           ENDDO
         ENDIF

      ENDDO

      END

*=============================================================================*

      SUBROUTINE r09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CHBr3 photolysis=*
*=          CHBr3 + hv -> Products                                           =*
*=  Cross section: Choice of data from Atlas (?Talukdar???) or JPL97         =*
*=  Quantum yield: Assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)

      real t
      real qy

      INTEGER i, iw, n
      INTEGER ierr
      INTEGER iz

      integer kopt


*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE



**************** CHBr3 photodissociation

      j = j + 1
      jlabel(j) = 'CHBr3 -> Products'

* option:

* kopt = 1:  cross section from Elliot Atlas, 1997
* kopt = 2:  cross section from JPL 1997

      kopt = 2
      if (kopt .eq. 1) then

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/CHBr3.abs',STATUS='old')
      DO i = 1, 5
         READ(kin,*)
      ENDDO

      n5 = 25
      n4 = 27
      n3 = 29
      n2 = 31
      n1 = 39
      DO i = 1, n5
         READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i), y5(i)
      ENDDO
      do i = n5 + 1, n4
         READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i)
      enddo
      do i = n4 + 1, n3
         READ(kin,*) x1(i), y1(i), y2(i), y3(i)
      enddo
      do i = n3 + 1, n2
         READ(kin,*) x1(i), y1(i), y2(i)
      enddo
      do i = n2 + 1, n1
         READ(kin,*) x1(i), y1(i)
      enddo
      CLOSE (kin)

      do i = 1, n1
         y1(i) = y1(i) * 1.e-23
      enddo
      do i = 1, n2
         x2(i) = x1(i)
         y2(i) = y2(i) * 1.e-23
      enddo
      do i = 1, n3
         x3(i) = x1(i)
         y3(i) = y3(i) * 1.e-23
      enddo
      do i = 1, n4
         x4(i) = x1(i)
         y4(i) = y4(i) * 1.e-23
      enddo
      do i = 1, n5
         x5(i) = x1(i)
         y5(i) = y5(i) * 1.e-23
      enddo

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
      CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),y3(1))
      CALL addpnt(x3,y3,kdata,n3,               0.,y3(1))
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,           1.e+38,0.)
      CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)

      ENDIF

      CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),y4(1))
      CALL addpnt(x4,y4,kdata,n4,               0.,y4(1))
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,           1.e+38,0.)
      CALL inter2(nw,wl,yg4,n4,x4,y4,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      CALL addpnt(x5,y5,kdata,n5,x5(1)*(1.-deltax),y5(1))
      CALL addpnt(x5,y5,kdata,n5,               0.,y5(1))
      CALL addpnt(x5,y5,kdata,n5,x5(n5)*(1.+deltax),0.)
      CALL addpnt(x5,y5,kdata,n5,           1.e+38,0.)
      CALL inter2(nw,wl,yg5,n5,x5,y5,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF


* quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO iz = 1, nz

            t = tlev(iz)

            if (t .ge. 296.) then
               yg(iw) = yg1(iw)

            else if(t .ge. 286.) then
               yg(iw) = yg1(iw) + (t-286.)*(yg2(iw)-yg1(iw))/10.

            else if(t .ge. 276.) then
               yg(iw) = yg2(iw) + (t-276.)*(yg3(iw)-yg2(iw))/10.

            else if(t .ge. 266.) then
               yg(iw) = yg3(iw) + (t-266.)*(yg4(iw)-yg3(iw))/10.

            else if(t .ge. 256.) then
               yg(iw) = yg4(iw) + (t-256.)*(yg5(iw)-yg4(iw))/10.

            else if(t .lt. 256.) then
               yg(iw) = yg5(iw)

            endif

            sq(j,iz,iw) = yg(iw)*qy

         ENDDO
      ENDDO

* jpl97, with temperature dependence formula,
*w = 290 nm to 340 nm, 
*T = 210K to 300 K
*sigma, cm2 = exp((0.06183-0.000241*w)*(273.-T)-(2.376+0.14757*w))

      ELSEIF (kopt .EQ. 2) THEN

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/CHBr3.jpl97',STATUS='old')
      DO i = 1, 6
         READ(kin,*)
      ENDDO
      n1 = 87
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO iz = 1, nz

            t = tlev(iz)
            yg(iw) = yg1(iw)

            IF (wc(iw) .GT. 290. .AND. wc(iw) .LT. 340. 
     $           .AND. t .GT. 210 .AND. t .LT. 300) THEN
               yg(iw) = EXP((0.06183-0.000241*wc(iw))*(273.-T)-
     $              (2.376+0.14757*wc(iw)))
            ENDIF

            sq(j,iz,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      ENDIF

      END

*=============================================================================*

      SUBROUTINE r10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH2O photolysis =*
*=        (a) CH2O + hv -> H + HCO                                           =*
*=        (b) CH2O + hv -> H2 + CO                                           =*
*=  Cross section: Choice between                                            =*
*=                 1) Bass et al., 1980 (resolution: 0.025 nm)               =*
*=                 2) Moortgat and Schneider (resolution: 1 nm)              =*
*=                 3) Cantrell et al. (orig res.) for > 301 nm,              =*
*=                    IUPAC 92, 97 elsewhere                                 =*
*=                 4) Cantrell et al. (2.5 nm res.) for > 301 nm,            =*
*=                    IUPAC 92, 97 elsewhere                                 =*
*=                 5) Rogers et al., 1990                                    =*
*=                 6) new NCAR recommendation, based on averages of          =*
*=                    Cantrell et al., Moortgat and Schneider, and Rogers    =*
*=                    et al.                                                 =*
*=                 7) Meller and Moortgat (2000) added by Judit Zador        =*
*=  Quantum yield: Choice between                                            =*
*=                 1) Evaluation by Madronich 1991 (unpublished)             =*
*=                 2) IUPAC 89, 92, 97                                       =*
*=                 3) Madronich, based on 1), updated 1998.                  =*
*=                 4) IUPAC 2003, added by Judit Zador					   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=16000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mopt1, mopt2

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** CH2O photodissociatation

      j = j+1
      jlabel(j) = 'CH2O -> H + HCO' 

      j = j+1
      jlabel(j) = 'CH2O -> H2 + CO'

* working grid arrays:
*     yg1 = cross section at a specific temperature
*     yg2, yg3 = cross sections at different temp or slope, for calculating
*                temperature depedence
*     yg4 = quantum yield data for radical channel
*     yg5 = quantum yield data for molecular channel

* Input data options:
* mopt1 for absorption:
* 1:  DATAJ1/CH2O/CH2O_nbs.abs'
*     from Bass et al., Planet. Space. Sci. 28, 675, 1980.
*     over 258.750-359.525 in 0.025 nm steps
* 2:  DATAJ1/CH2O_iupac1.abs 
*     Moortgat and Schneider, personal communication as reported in IUPAC 89, 92, 97
*     at 285K.  Over 240-360 nm in 1 nm bins (note that IUPAC 89,92,97 incorectly 
*     claims 0.5 nm intervals in footnote)
* 3:  DATAJ1/CH2O/ch2o_can_hr.abs for wc > 301 nm, temperature dependent
*     DATAJ1/CH2O/ch2o_iupac1.abs elsewhere
*     from Cantrell et al. 1990 for wc > 301 nm.  Original data from Cantrell,
*     at high resolution
* 4:  DATAJ1/CH2O/CH2O_can_lr.abs for wc > 301 nm, temperature dependent
*     DATAJ1/CH2O/CH2O_iupac1.abs elsewhere
*     from Cantrell et al. 1990 for wc > 301 nm.  Data from Cantrell et al., as
*     reported by IUPAC'92,'97.  On 2.5 nm intervals.
* 5:  DATAJ1/CH2O/CH2O_rog.abs'
*     from Rogers et al., J. Phys. Chem. 94, 4011, 1990.
* 6:  DATAJ2/CH2O_ncar.abs
*     new NCAR recommendation, based on averages of Moortgat and Schneider, Cantrell et al.,
*     and Rogers.
* 7:  DATAJ2/CH2O/CH2O_meller_moortgat.abs
*     new IUPAC 2002 recommendation, 1 nm resolution, wavelenght is calibrated for air
*     data of Meller and Moortgat 2000, added by Judit Zador
* mopt2 for quantum yields:
* 1:  DATAJ1/CH2O/CH2O_i_mad.yld and 
*     DATAJ1/CH2O/CH2O_ii_mad.yld
*     evaluated by Madronich, 1991, unpublished
* 2:  DATAJ1/CH2O/CH2O_iupac.yld
*     from IUPAC'89, '92, '97
* 3:  DATAJ1/CH2O/CH2O_jpl97.dat'
*     based on Madronich 1991 unpublished evaluation, updated Jan 1998.
* 4:  DATAJ1/CH2O/CH2O_iupac2003.yld'
*     based on Smith et al., 2002, added by Judit Zador
      mopt1 = 7
      mopt2 = 4

      IF (mopt1 .EQ. 1) THEN

* read NBS/Bass data

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_nbs.abs'
     $        ,STATUS='old')
         n = 4032
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)

         CALL inter2(nw,wl,yg1,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j-1)
            STOP
         ENDIF

      ELSEIF (mopt1 .EQ. 2 .OR. mopt1 .EQ. 3 .OR. mopt1 .EQ. 4) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O_iupac1.abs',STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 121
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)
         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j-1)
            STOP
         ENDIF

         IF(mopt1 .EQ. 3) THEN

* data are on wavenumber grid (cm-1), so convert to wavelength in nm:
* grid was on increasing wavenumbers, so need to reverse to get increasing
* wavelengths
* cross section assumed to be zero for wavelengths longer than 360 nm
* if y1 < 0, then make = 0 (some negative cross sections, actually 273 K intercepts
* are in the original data,  Here, make equal to zero)

         OPEN(kin,FILE='DATAJ1/CH2O/CH2O_can_hr.abs',STATUS='old')
         READ(kin,*) idum, n
         DO i = 1, idum-2
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x1(i) = 1./x1(i) * 1E7
            IF (x1(i) .GT. 360.) THEN
               y1(i) = 0.
               y2(i) = 0.
            ENDIF
         ENDDO
         CLOSE(kin)

         DO i = 1, n/2
            irev = n+1-i
            dum = x1(i)
            x1(i) = x1(irev)
            x1(irev) = dum
            dum = y1(i)
            y1(i) = y1(irev)
            y1(irev) = dum
            dum = y2(i)
            y2(i) = y2(irev)
            y2(irev) = dum
         ENDDO
         DO i = 1, n
            x2(i) = x1(i)
            y1(i) = max(y1(i),0.)
         ENDDO
         n1 = n
         n2 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,             1E38,0.)
         CALL inter2(nw,wl,yg2,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,               0.,0.)
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,              1E38,0.)
         CALL inter2(nw,wl,yg3,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mopt1 .eq. 4) THEN

            OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_can_lr.abs',
     $        STATUS='old')
            DO i = 1, 4
               READ(kin,*)
            ENDDO
            n = 23
            DO i = 1, n
               READ(kin,*) x2(i), y2(i), y3(i), dum, dum
               x3(i) = x2(i)
            ENDDO
            CLOSE(kin)
            n2 = n
            n3 = n

            CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
            CALL addpnt(x2,y2,kdata,n2,               0.,0.)
            CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
            CALL addpnt(x2,y2,kdata,n2,             1E38,0.)
            CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF

            CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
            CALL addpnt(x3,y3,kdata,n3,               0.,0.)
            CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
            CALL addpnt(x3,y3,kdata,n3,              1E38,0.)
            CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF

         ENDIF

      ELSEIF (mopt1 .EQ. 5) THEN

* read Rodgers data

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_rog.abs'
     $        ,STATUS='old')
         DO i = 1, 10
            READ(kin,*)
         ENDDO
         n = 261
         DO i = 1, n
            READ(kin,*) x(i), y(i), dum
            y(i) = y(i) * 1.e-20
         ENDDO
         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j-1)
            STOP
         ENDIF

      ELSEIF(mopt1 .EQ. 6) THEN

            OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_ncar.abs',STATUS='old')
            DO i = 1, 3
               READ(kin,*)
            ENDDO
            n = 126
            DO i = 1, n
               READ(kin,*) x2(i), y2(i), y3(i)
               x3(i) = x2(i)
            ENDDO
            CLOSE(kin)
            n2 = n
            n3 = n

            CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
            CALL addpnt(x2,y2,kdata,n2,               0.,0.)
            CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
            CALL addpnt(x2,y2,kdata,n2,             1E38,0.)
            CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF

            CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
            CALL addpnt(x3,y3,kdata,n3,               0.,0.)
            CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
            CALL addpnt(x3,y3,kdata,n3,              1E38,0.)
            CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF

* Meller and Moortgat (2000), with temperature dependence
	ELSEIF(mopt1 .EQ. 7) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/CH2O/CH2O_meller_moortgat.abs',STATUS='old')
            DO i = 1, 8
               READ(kin,*)
            ENDDO
            n = 150
            DO i = 1, n
               READ(kin,*) x2(i), y2(i), y3(i)
			 y2(i) = y2(i)*1.e-21
	         y3(i) = y3(i)*1.e-24
	         x3(i) = x2(i)
            ENDDO
            CLOSE(kin)
            n2 = n		
	      n3 = n

            CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
            CALL addpnt(x2,y2,kdata,n2,               0.,0.)
            CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
            CALL addpnt(x2,y2,kdata,n2,             1E38,0.)
            CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF

            CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
            CALL addpnt(x3,y3,kdata,n3,               0.,0.)
            CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
            CALL addpnt(x3,y3,kdata,n3,              1E38,0.)
            CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF

      ENDIF
      
* quantum yield

      IF (mopt2 .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_i_mad.yld',STATUS='old')
         DO i = 1, 11
            READ(kin,*)
         ENDDO
         n = 20
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)
         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),y(1))
         CALL addpnt(x,y,kdata,n,               0.,y(1))
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg4,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j-1)
            STOP
         ENDIF

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_ii_mad.yld',STATUS='old')
         DO i = 1, 9
            READ(kin,*)
         ENDDO
         n = 33
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)
         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),y(1))
         CALL addpnt(x,y,kdata,n,               0.,y(1))
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg5,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mopt2 .EQ. 2) then

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_iupac.yld',STATUS='old')
         DO i = 1, 7
            READ(kin,*) 
         ENDDO
         n = 13
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)
         n1 = n
         n2 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg4,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
         CALL inter2(nw,wl,yg5,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

* box-filling interpolation.  
c         DO i = 1, n
c            READ(kin,*) x1(i), y1(i), y2(i)
c            x1(i) = x1(i) - 5.0
c            x2(i) = x1(i)
c         ENDDO
c         n = n + 1
c         x1(n) = x1(n-1) + 5.0
c         x2(n) = x1(n)
c         CLOSE(kin)
c         DO i = 1, n-1
c            y1(i) = y1(i) * (x1(i+1)-x1(i))
c         ENDDO
c         CALL inter3(nw,wl,yg4,n,x1,y1,0)
c         DO iw = 1, nw-1
c            yg4(iw) = yg4(iw)/(wl(iw+1)-wl(iw))
c         ENDDO
c         DO i = 1, n-1
c            y2(i) = y2(i) * (x2(i+1)-x2(i))
c         ENDDO
c         CALL inter3(nw,wl,yg5,n,x2,y2,0)
c         DO iw = 1, nw-1
c            yg5(iw) = yg5(iw)/(wl(iw+1)-wl(iw))
c         ENDDO

      ELSE IF(mopt2 .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_jpl97.dat',STATUS='old')
         DO i = 1, 4
            READ(kin,*) 
         ENDDO
         n = 23
         DO i = 1, n
            READ(kin,*) x1(i), dum, dum, dum, dum, y1(i), y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)
         n1 = n
         n2 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg4,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
         CALL inter2(nw,wl,yg5,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mopt2 .EQ. 4) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH2O/CH2O_iupac2003.yld',
     >	   STATUS='old')
         DO i = 1, 4
            READ(kin,*) 
         ENDDO
         n = 26
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)
         n1 = n
         n2 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg4,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
         CALL inter2(nw,wl,yg5,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* combine
* y1 = xsect
* y2 = xsect(223), Cantrell et al.
* y3 = xsect(293), Cantrell et al.
* y4 = qy for radical channel
* y5 = qy for molecular channel
* pressure and temperature dependent for w > 330.

      DO iw = 1, nw - 1

         IF (mopt1 .eq. 6) THEN
            sig = yg2(iw)
         ELSE
            sig = yg1(iw)
         ENDIF

         DO i = 1, nz

* correct cross section for temperature dependence for > 301. nm
         
            IF (mopt1 .LE. 6) THEN !everything except for Meller & Moortgat
			  IF (wl(iw) .GE. 301.) THEN 
				 t = MAX(223.15, MIN(tlev(i), 293.15))
				 IF (mopt1 .EQ. 3 .OR. mopt1 .EQ. 6) THEN
					sig = yg2(iw) + yg3(iw) * (t - 273.15)

				 ELSEIF (mopt1 .EQ. 4) THEN
					slope = (yg3(iw) - yg2(iw)) / (293. - 223.)
					sig = yg2(iw) + slope * (t - 223.)

				 ENDIF

			  ENDIF
			  sig = MAX(sig, 0.)

            ELSE 
	         sig = yg2(iw) + yg3(iw) * (tlev(i) - 298.)
            ENDIF   

* quantum yields:
* temperature and pressure dependence beyond 330 nm

            qy1 = yg4(iw)
            IF ( (wc(iw) .GE. 330.) .AND. (yg5(iw) .GT. 0.) ) THEN
               phi1 = yg4(iw)
               phi2 = yg5(iw)
               phi20 = 1. - phi1
               ak300=((1./phi2)-(1./phi20))/2.54E+19
               akt=ak300*(1.+61.69*(1.-tlev(i)/300.)*(wc(iw)/329.-1.))
               qy2 = 1. / ( (1./phi20) + airden(i)*akt)

            ELSE
               qy2 = yg5(iw)
            ENDIF
            qy2 = MAX(0.,qy2)
            qy2 = MIN(1.,qy2)
            
            sq(j-1,i,iw) = sig * qy1
            sq(j  ,i,iw) = sig * qy2

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3CHO photolysis: =*
*=      (a)  CH3CHO + hv -> CH3 + HCO                                        =*
*=      (b)  CH3CHO + hv -> CH4 + CO                                         =*
*=      (c)  CH3CHO + hv -> CH3CO + H                                        =*
*=  Cross section:  Choice between                                           =*
*=                   (1) IUPAC 97 data, from Martinez et al.                 =*
*=                   (2) Calvert and Pitts                                   =*
*=                   (3) Martinez et al., Table 1 scanned from paper         =*
*=                   (4) KFA tabulations                                     =*
*=  Quantum yields: Choice between                                           =*
*=                   (1) IUPAC 97, pressure correction using Horowith and    =*
*=                                 Calvert, 1982                             =*
*=                   (2) NCAR data file, from Moortgat, 1986                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=150)

      INTEGER i, n
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy1, qy2, qy3
      REAL sig
      REAL dum
      INTEGER ierr
      INTEGER  iz, iw

      INTEGER mabs, myld

****************************************************************
************************* CH3CHO photolysis
* 1:  CH3 + HCO
* 2:  CH4 + CO
* 3:  CH3CO + H


      j = j+1
      jlabel(j) = 'CH3CHO -> CH3 + HCO'
      j = j+1
      jlabel(j) = 'CH3CHO -> CH4 + CO'
      j = j+1
      jlabel(j) = 'CH3CHO -> CH3CO + H'

* options
* mabs for cross sections
* myld for quantum yields

* Absorption:
* 1:  IUPAC-97 data, from Martinez et al.
* 2:  Calvert and Pitts
* 3:  Martinez et al., Table 1 scanned from paper
* 4:  KFA tabulations, 6 choices, see file OPEN statements

* Quantum yield
* 1:  DATAJ1/CH3CHO/CH3CHO_iup.yld
* pressure correction using Horowitz and Calvert 1982, based on slope/intercept
* of Stern-Volmer plots

* 2:  ncar data file, from Moortgat 1986.
*     DATAJ1/CH3CHO/d021_i.yld
*     DATAJ1/CH3CHO/d021_i.yld
*     DATAJ1/CH3CHO/d021_i.yld

      mabs = 3
      myld = 1

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/CH3CHO_iup.abs',STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 106
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 2) THEN

* cross section from Calvert and  Pitts
         
         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/d021_cp.abs',STATUS='old')
         DO i = 1, 14
            READ(kin,*)
         ENDDO
         n = 54
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            x1(i) = x1(i)/10.
            y1(i) = y1(i) * 3.82E-21
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/CH3CHO_mar.abs',STATUS='old')
         DO i = 1, 3
            READ(kin,*)
         ENDDO
         n = 106
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 4) THEN

* cross section from KFA tables
* ch3cho.001 - Calvert and Pitts 1966
* ch3cho.002 - Meyrahn thesis 1984
* ch3cho.003 - Schneider and Moortgat, priv comm. MPI Mainz 1989, 0.012 nm resol.
* ch3cho.004 - Schneider and Moortgat, priv comm. MPI Mainz 1989, 0.08  nm resol.
* ch3cho.005 - IUPAC'92
* ch3cho.006 - Libuda, thesis Wuppertal 1992
         
c         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.001',STATUS='old')
C         n = 217
c         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.002',STATUS='old')
c         n = 63
c         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.003',STATUS='old')
c         n = 13738
c         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.004',STATUS='old')
c         n = 2053
         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.005',STATUS='old')
         n = 18
c         OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cho.006',STATUS='old')
c         n = 1705

         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yields

      IF (myld .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/CH3CHO_iup.yld',STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 12
         DO i = 1, n
            READ(kin,*) x1(i), y2(i), y1(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)
         n1 = n
         n2 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,               0.,0.)
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
         CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         DO iw = 1, nw-1
            yg3(iw) = 0.
         ENDDO

      ELSEIF (myld .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/d021_i.yld',STATUS='old')
         DO i = 1, 18
            READ(kin,*)
         ENDDO
         n = 10
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF
      
         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/d021_ii.yld',STATUS='old')
         DO i = 1, 10
            READ(kin,*)
         ENDDO
         n = 9
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/d021_iii.yld',STATUS='old')
         DO i = 1, 10
            READ(kin,*)
         ENDDO
         n = 9
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg3,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* pressure-dependence parameters
      
         OPEN(UNIT=kin,FILE='DATAJ1/CH3CHO/CH3CHO_press.yld',
     $     STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 5
         DO i = 1, n
            READ(kin,*) x1(i), dum, dum, y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg4,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

* combine:

      DO iw = 1, nw - 1
         DO i = 1, nz

            sig = yg(iw)

* quantum yields:

            qy1 = yg1(iw)
            qy2 = yg2(iw)
            qy3 = yg3(iw)

* pressure correction for channel 1, CH3 + CHO
* based on Horowitz and Calvert 1982.

            qy1 = qy1 * (1. + yg4(iw))/(1. + yg4(iw)*airden(i)/2.465E19)
            qy1 = MIN(1., qy1)
            qy1 = MAX(0., qy1)

            sq(j-2,i,iw) = sig * qy1
            sq(j-1,i,iw) = sig * qy2
            sq(j  ,i,iw) = sig * qy3

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for C2H5CHO        =*
*=  photolysis:                                                              =*
*=         C2H5CHO + hv -> C2H5 + HCO                                        =*
*=                                                                           =*
*=  Cross section:  Choice between                                           =*
*=                   (1) IUPAC 97 data, from Martinez et al.                 =*
*=                   (2) Calvert and Pitts, as tabulated by KFA              =*
*=  Quantum yield:  Choice between 	         	                     =*
*=		     (1) IUPAC 97 recommendation                             =*
*=		     (2) IUPAC 05 recommendation (extended out to 340 nm)    =*
*=		     (3) Chen and Zhu (2001) pressure dependent              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=150)

      INTEGER i, n
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      REAL qy1
      REAL sig
      INTEGER ierr
      INTEGER iw

      INTEGER mabs, myld

************************* C2H5CHO photolysis
* 1:  C2H5 + HCO

      j = j+1
      jlabel(j) = 'C2H5CHO -> C2H5 + HCO'

* options
* mabs for cross sections
* myld for quantum yields

* Absorption:
* 1:  IUPAC-97 data, from Martinez et al.
* 2:  Calvert and Pitts, as tabulated by KFA.

* Quantum yield
* 1:  IUPAC-97 data
* 2:  IUPAC-05 data (extended out to 340 nm)
* 3:  Chen and Zhu, J. Phys. Chem. A., 105, 9689 (2001)
*     Pressure dependent

      mabs = 1
      myld = 2

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/C2H5CHO/C2H5CHO_iup.abs',
     $        STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 106
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 2) THEN

* cross section from KFA tables
* c2h5cho.001 - Calvert and Pitts 1966
         
         OPEN(UNIT=kin,FILE='DATAJ2/KFA/c2h5cho.001',STATUS='old')
         n = 83

         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yields

      IF (myld .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/C2H5CHO/C2H5CHO_iup97.yld',
     $        STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 5
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)
         n1 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,340.,0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

       ELSEIF (myld .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/C2H5CHO/C2H5CHO_iup05.yld',
     $        STATUS='old')
         do i = 1, 6
            read(kin,*)
         enddo
         n = 7
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)
         n1 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,340.,0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF (myld .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/C2H5CHO/C2H5CHO_Chen01.yld',
     $        STATUS='old')
         DO i = 1, 11
            read(kin,*)
         ENDDO

         n = 11

         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)
	   n1 = n
	   n2 = n
         
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,               0.,0.)
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
         CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF   
   
      ENDIF

* combine:

      DO iw = 1, nw - 1
         DO i = 1, nz

            sig = yg(iw)

* quantum yields:
* use Stern-Volmer pressure dependence:

         IF (myld .EQ. 1) THEN

            IF (yg1(iw) .LT. pzero) THEN
               qy1 = 0.
            ELSE
               qy1 = 1./(1. + (1./yg1(iw) - 1.)*airden(i)/2.45e19)
            ENDIF
            
         ELSEIF (myld .EQ. 2) THEN 
         
            qy1 = yg1(iw)

         ELSEIF (myld .EQ. 3) THEN

               qy1 = 1/((1/yg1(iw))+(yg2(iw)*airden(i)))
         
         ENDIF

            qy1 = MIN(qy1,1.)
            sq(j,i,iw) = sig * qy1
         
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CHOCHO         =*
*=  photolysis:                                                              =*
*=              CHOCHO + hv -> Products                                      =*
*=                                                                           =*
*=  Cross section: Choice between                                            =*
*=                  (1) Plum et al., as tabulated by IUPAC 97                =*
*=                  (2) Plum et al., as tabulated by KFA.                    =*
*=                  (3) Orlando et al.                                       =*
*=                  (4) Horowitz et al., 2001                                =*
*=                  (5) Volkamer et al., 2005                                =*
*=  Quantum yield: (1) IUPAC 97 recommendation (Calvert, 2000)               =*
*=                 (2) RADICAL REPORT vI (Moortgat et al., 2002)             =*
*=                 (3) RADICAL REPORT vII (Moortgat et al., 2002)            =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=500)

      INTEGER i, n, n2, n3
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      REAL qyI, qyII, qyIII
      REAL sig
      INTEGER ierr
      INTEGER iw

      INTEGER mabs, myld

************************* CHOCHO photolysis
* see review by Madronich, Chapter VII in "The Mechansims of 
*  Atmospheric Oxidation of the Alkanes, Calvert et al, Oxford U.
*  Press, 2000.
* Four possible channels:
*     I     H2 + 2 CO
*     II    2 HCO
*     III   HCHO + CO
*     IV    HCO + H + CO
*
*  Based on that review, the following quantum yield assignments are made:
*
*     qy_I = 0
*     qy_II = 0.63 for radiation between 280 and 380 nm
*     qy_III = 0.2  for radiation between 280 and 380 nm
*     qy_IV = 0
* The yields for channels II and III were determined by Bauerle et al. (personal
* communication from G. Moortgat, still unpublished as of Dec 2000).
* Bauerle et al. used broad-band irradiation 280-380 nm.
* According to Zhu et al., the energetic threshold (for II) is 417 nm.  Therefore,
* here the quantum yields were set to zero for wc > 417.  Furthermore, the
* qys of Bauerle et al. were reduced to give the same J values when using full solar
* spectrum.  The reduction factor was calculated by comparing the J-values (for 
* high sun) using the 380 and 417 cut offs.  The reduction factor is 7.1

      j = j + 1
      jlabel(j) = 'CHOCHO -> 2CO + H2'
      
      j = j + 1
      jlabel(j) = 'CHOCHO -> HCO + HCO'

      j = j + 1
      jlabel(j) = 'CHOCHO -> CH2O + CO'

* options
* mabs for cross sections
* myld for quantum yields

* Absorption:
* 1:  Plum et al., as tabulated by IUPAC-97
* 2:  Plum et al., as tabulated by KFA.
* 3:  Orlando, J. J.; G. S. Tyndall, 2001:  The atmospheric chemistry of the
*        HC(O)CO radical. Int. J. Chem. Kinet., 33, 149-156.
* 4:  Horowitz, A., R. Meller, and G. K. Moortgat, 
*       The UV-VIS absorption cross sectiono of the a-dicarbonyl compounds:
*       pyruvic acid, biacetyl, and glyoxal.
*       J. Photochem. Photobiol. A:Chemistry, v.146, pp.19-27, 2001.
* 5:  Volkamer, R., Spietz, P., Burrows, J. and Platt, U: 
*       High-resolution absorption cross-section of glyoxal 
*       in the uv-vis and ir spectral ranges.
*       J. Photochem. Photobiol. A:Chemistry, v.172, pp.35-46, 2005.
*
* Quantum yield
* 1:  IUPAC-97 data (Calvert (2000))
* 2:  RADICAL REPORT vI:  Moortgat, G.K. (Ed.), RADICAL:  Evaluation of radical sources 
*        in atmospheric chemistry through chamber and laboratory studies. 
*        e-mail moo@mpch-mainz.mpg.de for copies (EUPHORE Studies)
* 3:  RADICAL REPORT vII (model yields) 

      mabs = 5
      myld = 3

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CHOCHO/CHOCHO_iup.abs',
     $        STATUS='old')
         DO i = 1, 4
            read(kin,*)
         ENDDO
         n = 110
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE(kin)


         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF


      ELSEIF(mabs .EQ. 2) THEN

* cross section from KFA tables
* chocho.001 - Plum et al. 1983
         
         OPEN(UNIT=kin,FILE='DATAJ2/KFA/chocho.001',STATUS='old')
         n = 219

         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 3) THEN

* cross section from Orlando et la.
* Orlando, J. J.; G. S. Tyndall, 2001:  The atmospheric chemistry of the
* HC(O)CO radical. Int. J. Chem. Kinet., 33, 149-156.

         OPEN(UNIT=kin,
     $        FILE='DATAJ1/CHOCHO/glyoxal_orl.abs',STATUS='old')

         do i = 1, 6
            read(kin,*)
         enddo
         n = 481
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 4) THEN

         OPEN(UNIT=kin,
     $        FILE='DATAJ1/CHOCHO/glyoxal_horowitz.abs',STATUS='old')

         DO i = 1, 8
            read(kin,*)
         ENDDO
         n = 270
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 5) THEN

         OPEN(UNIT=kin,
     $        FILE='DATAJ1/CHOCHO/glyoxal_volkamer.abs',STATUS='old')

         DO i = 1, 15
            read(kin,*)
         ENDDO
         n = 276
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yields

* combine:

      DO iw = 1, nw - 1

              sig = yg(iw)

* quantum yields:
* Use values from Bauerle, but corrected to cutoff at 417 rather than 380.
* this correction is a reduction by 7.1.
* so that qyII = 0.63/7.1  and qyIII = 0.2/7.1

              IF (myld .EQ. 1) THEN

                if(wc(iw) .lt. 417. ) then
                   qyII = 0.089
                   qyIII = 0.028
                else
                   qyII = 0.
                   qyIII = 0.
                endif
               
                DO i = 1, nz
                   yg1(iw) = qyII
                   yg2(iw) = qyIII
                   
                   sq(j-1,i,iw) = sig * qyII
                   sq(j,i, iw) = sig * qyIII
                ENDDO

              ELSEIF(myld .EQ. 2) THEN

* 2:  RADICAL REPORT vI:  Moortgat, G.K. (Ed.), RADICAL:  Evaluation of radical sources 
*     in atmospheric chemistry through chamber and laboratory studies.
*     Using the relative qy as measured @EUPHORE multiplied by the effective qy of 
*     0.0035 (Volkamer et al. 2005) 

                OPEN(UNIT=kin, FILE
     $               ='DATAJ1/CHOCHO/glyoxal_RADICALa.yld',STATUS='old')

                DO i = 1, 14
                   read(kin,*)
                ENDDO
                n = 48
                DO i = 1, n
                   READ(kin,*) x1(i), y1(i), y2(i), y3(i)
                   x2(i) = x1(i)
                   x3(i) = x1(i)
                ENDDO
                n2 = n
                n3 = n
                CLOSE (kin)

                 CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
                 CALL addpnt(x1,y1,kdata,n,               0.,0.)
                 CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
                 CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
                 CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
                 IF (ierr .NE. 0) THEN
                    WRITE(*,*) ierr, jlabel(j)
                    STOP
                 ENDIF   

                 CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
                 CALL addpnt(x2,y2,kdata,n2,               0.,0.)
                 CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
                 CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
                 CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
                 IF (ierr .NE. 0) THEN
                    WRITE(*,*) ierr, jlabel(j)
                    STOP
                 ENDIF   

                 CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
                 CALL addpnt(x3,y3,kdata,n3,               0.,0.)
                 CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
                 CALL addpnt(x3,y3,kdata,n3,           1.e+38,0.)
                 CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)
                 IF (ierr .NE. 0) THEN
                    WRITE(*,*) ierr, jlabel(j)
                    STOP
                 ENDIF   

                qyI = yg1(iw)
		        qyII = yg2(iw)
		        qyIII = yg3(iw)

                DO i = 1, nz
                   sq(j-2,i,iw) = sig * qyI
                   sq(j-1,i,iw) = sig * qyII
                   sq(j,i, iw) = sig * qyIII
                ENDDO

             ELSEIF (myld .EQ. 3) THEN
             
* 3:  RADICAL REPORT vII (model yields) 
*     Yeilds multiplied by effective qy of 0.0035 (Volkamer et al. (2005))
*     Cutt off at 418 nm (see RADICAL report)
              
               IF(wc(iw) .LT. 418. ) THEN
        
                   qyI = 0.0044
                   qyII = 0.026
                   qyIII = 0.0044
                ELSE
                   qyI = 0.
                   qyII = 0.
                   qyIII = 0.
                   
               ENDIF
               
                DO i = 1, nz
                
                   yg1(iw) = qyI
                   yg2(iw) = qyII
                   yg3(iw) = qyIII 
                   
                   sq(j-2,i,iw) = sig * qyI
                   sq(j-1,i,iw) = sig * qyII
                   sq(j,i, iw) = sig * qyIII
                ENDDO

             ENDIF

      ENDDO

      END

*=============================================================================*

      SUBROUTINE r14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3COCHO       =*
*=  photolysis:                                                              =*
*=           CH3COCHO + hv -> CH3CO + HCO                                    =*
*=                                                                           =*
*=  Cross section: Choice between                                            =*
*=                  (1) from Meller et al., 1991, as tabulated by IUPAC 97   =*
*=                         5 nm resolution (table 1) for < 402 nm            =*
*=                         2 nm resolution (table 2) for > 402 nm            =*
*=                  (2) average at 1 nm of Staffelbach et al., 1995, and     =*
*=                      Meller et al., 1991                                  =*
*=                  (3) Plum et al., 1983, as tabulated by KFA	             =*
*=                  (4) Meller et al., 1991 (0.033 nm res.), as tab. by KFA  =*
*=                  (5) Meller et al., 1991 (1.0 nm res.), as tab. by KFA    =*
*=                  (6) Staffelbach et al., 1995, as tabulated by KFA        =*
*=                  (7) synthetic spectrum (average)                         =*
*=                  (8) JPL 2005, average Staffelbach et al., 1995, and      =*
*=                      Meller et al., 1991                                  =*
*=  Quantum yield: Choice between                                            =*
*=                  (1) Plum et al., fixed at 0.107                          =*
*=                  (2) Plum et al., divided by 2, fixed at 0.0535           =*
*=                  (3) Staffelbach et al., 0.45 for < 300 nm, 0 for > 430 nm=*
*=                      linear interp. in between                            =*
*=                  (4) Koch and Moortgat, prv. comm., 1997                  =*
*=                  (5) Chen et al., J. Phys. Chem. A., 104, 11126, 2000     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=500)

      INTEGER i, n
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)
      real x(kdata), y(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw)
      REAL qy
      REAL sig
      INTEGER ierr
      INTEGER iw

      INTEGER mabs, myld

      REAL phi0, kq, pa, pt


************************* CH3COCHO photolysis
* 1:  CH3COCHO

      j = j+1
      jlabel(j) = 'CH3COCHO -> CH3CO + HCO'

* options
* mabs for cross sections
* myld for quantum yields

* Absorption:
* 1:  from Meller et al. (1991), as tabulated by IUPAC-97
*     for wc < 402, use coarse data (5 nm, table 1)
*     for wc > 402, use finer data (2 nm, table 2)
* 2: average at 1nm of  Staffelbach et al. 1995 and Meller et al. 1991
*     Cross section from KFA tables:
* 3: ch3cocho.001 - Plum et al. 1983
* 4: ch3cocho.002 - Meller et al. 1991, 0.033 nm resolution
* 5: ch3cocho.003 - Meller et al. 1991, 1.0   nm resolution
* 6: ch3cocho.004 - Staffelbach et al. 1995
* 7: use synthetic spectrum, average of CHOCHO and CH3COCOCH3:
* 8: JPL 2005, average Staffelbach et al., 1995, and Meller et al. 1991



* Quantum yield
* 1:  Plum et al., 0.107
* 2:  Plum et al., divided by two = 0.0535
* 3:  Staffelbach et al., 0.45 at wc .le. 300, 0 for wc .gt. 430, linear 
*     interpl in between
* 4:  Koch and Moortgat, prv. comm. 1997. - pressure-dependent
* 5:  Chen, Y., W. Wang, and L. Zhu, Wavelength-dependent photolysis of methylglyoxal
*      in the 290-440 nm region, J Phys Chem A, 104, 11126-11131, 2000.
* 6:  Chen et al. (2000), new calculation

      mabs = 8
      myld = 4

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCHO/CH3COCHO_iup1.abs',
     $        STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 38
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCHO/CH3COCHO_iup2.abs',
     $        STATUS='old')
         do i = 1, 4
            read(kin,*)
         enddo
         n = 75
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         DO iw = 1, nw-1 
            IF(wc(iw) .LT. 402.) THEN
               yg(iw) = yg1(iw)
            ELSE
               yg(iw) = yg2(iw)
            ENDIF               
         ENDDO

      ELSEIF(mabs .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCHO/CH3COCHO_ncar.abs',
     $        STATUS='old')
         n = 271
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .GT. 2 .and. mabs .lt. 7) THEN

* cross section from KFA tables
* ch3cocho.001 - Plum et al. 1983
* ch3cocho.002 - Meller et al. 1991, 0.033 nm resolution
* ch3cocho.003 - Meller et al. 1991, 1.0   nm resolution
* ch3cocho.004 - Staffelbach et al. 1995
         
         IF(mabs .EQ. 3) THEN
            OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cocho.001',STATUS='old')
            n = 136
         ELSEIF(mabs .EQ. 4) THEN
            OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cocho.002',STATUS='old')
            n = 8251
         ELSEIF(mabs .EQ. 5) THEN
            OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cocho.003',STATUS='old')
            n = 275
         ELSEIF(mabs .EQ. 6) THEN
            OPEN(UNIT=kin,FILE='DATAJ2/KFA/ch3cocho.004',STATUS='old')
            n = 162
         ENDIF
         
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         ELSEIF(mabs .EQ. 7) THEN

      OPEN(UNIT=kin,FILE='DATAJ1/CH3COCOCH3/biacetyl_plum.abs',
     $     STATUS='old')
      DO i = 1, 7
         READ(kin,*)
      ENDDO
      n = 55
      DO i = 1, n
         READ(kin,*) x(i), y(i)
         y(i) = y(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,               0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,x,y,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF


         OPEN(UNIT=kin,
     $        FILE='DATAJ1/CHOCHO/glyoxal_orl.abs',STATUS='old')
         do i = 1, 6
            read(kin,*)
         enddo
         n = 481
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         do iw = 1, nw-1
            yg(iw) = 0.5*(yg1(iw) + yg2(iw))
         enddo

      ELSEIF(mabs .EQ. 8) THEN

        OPEN(UNIT=kin,FILE='DATAJ1/CH3COCHO/ch3cocho_jpl05.abs',
     $     STATUS='old')
        DO i = 1, 16
           READ(kin,*)
        ENDDO
        n = 294
        DO i = 1, n
           READ(kin,*) x(i), y(i)
           y(i) = y(i) * 1.e-20
        ENDDO
        CLOSE(kin)

        CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
        CALL addpnt(x,y,kdata,n,              0.,0.)
        CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
        CALL addpnt(x,y,kdata,n,          1.e+38,0.)
        CALL inter2(nw,wl,yg,n,x,y,ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF

      ENDIF

* quantum yields

         IF(myld .EQ. 4) THEN
            OPEN(UNIT=kin,FILE='DATAJ1/CH3COCHO/CH3COCHO_km.yld',
     $           STATUS='old')
            DO i = 1, 5
               READ(kin,*)
            ENDDO
            n = 5
            DO i = 1, n
               READ(kin,*) x1(i), y1(i), y2(i)
               x2(i) = x1(i)
            ENDDO
            CLOSE (kin)
            n1 = n
            n2 = n

            CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),1.)
            CALL addpnt(x1,y1,kdata,n1,               0.,1.)
            CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
            CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF

            CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),1.)
            CALL addpnt(x2,y2,kdata,n2,               0.,1.)
            CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
            CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
            CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF


         ELSEIF(myld .EQ. 6) THEN

            OPEN(UNIT=kin,FILE='DATAJ1/CH3COCHO/CH3COCHO_chen00.yld',
     $           STATUS='old')

            DO i = 1, 10
               READ(kin,*)
            ENDDO

            n = 16

            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
               x2(i) = x1(i)
            ENDDO

            CLOSE (kin)

            n1 = n
           
            CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),1.)
            CALL addpnt(x1,y1,kdata,n1,               0.,1.)
            CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
            CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF

         ENDIF


* combine:

      DO iw = 1, nw - 1

         sig = yg(iw)

         DO i = 1, nz

* quantum yields:

            IF (myld .EQ. 1) THEN
               qy = 0.107

            ELSEIF(myld .EQ. 2) THEN
               qy = 0.107/2.

            ELSEIF(myld .EQ. 3) THEN
               IF(wc(iw) .LE. 300.) THEN
                  qy = 0.45
               ELSE IF (wc(iw) .GE. 430.) THEN 
                  qy = 0.
               ELSE
                  qy = 0.45 + (0-0.45)*(wc(iw)-300.)/(430.-300.)
               ENDIF

            ELSEIF(myld .EQ. 4) THEN

               IF (yg1(iw) .GT. 0.) THEN

                  qy = yg2(iw)/( 1. + (airden(i)/2.465E19) 
     $                 * ( (yg2(iw)/yg1(iw)) - 1.))

               ELSE
                  qy = 0.
               ENDIF
                
            ELSEIF(myld .EQ. 5) THEN
               
* zero pressure yield:
* 1.0 for wc < 380 nm
* 0.0 for wc > 440 nm
* linear in between:

               phi0 = 1. - (wc(iw) - 380.)/60.
               phi0 = MIN(phi0,1.)
               phi0 = MAX(phi0,0.)

* Pressure correction: quenching coefficient, torr-1
* in air, Koch and Moortgat:

               kq = 1.36e8 * EXP(-8793/wc(iw))

* in N2, Chen et al:

c               kq = 1.93e4 * EXP(-5639/wc(iw))

               IF(phi0 .GT. 0.) THEN
                  IF (wc(iw) .GE. 380. .AND. wc(iw) .LE. 440.) THEN
                     qy = phi0 / (phi0 + kq * airden(i) * 760./2.456E19)
                  ELSE
                     qy = phi0
                  ENDIF
               ELSE
                  qy = 0.
               ENDIF


            ELSEIF(myld .EQ. 6) THEN

* Pressure correction: quenching coefficient, torr-1
* in air, Koch and Moortgat:

               kq = 1.36e8 * EXP(-8793/wc(iw))
               
* Conversionn of air density (molec cm-3) air density (Torr):
* pa is pressure in Pa (Pascal)
* pt is pressure in Torr 

              pa = (airden(i)/7.243E16)*tlev(i)
              
              pt = pa/133.322

               IF (wc(iw) .GE. 290. .AND. wc(iw) .LE. 440.) THEN
                  
                  qy = 1/((1/yg1(iw)) + (kq * pt))

               ELSEIF (wc(iw) .LT. 290. ) THEN 

                  qy = 1.

               ELSE
                  qy = 0.
               ENDIF
                  yg2(iw) = qy
            ENDIF

            sq(j,i,iw) = sig * qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r15(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3COCH3 photolysis=*
*=          CH3COCH3 + hv -> Products                                        =*
*=                                                                           =*
*=  Cross section:  Choice between                                           =*
*=                   (1) Calvert and Pitts                                   =*
*=                   (2) Martinez et al., 1991, also in IUPAC 97             =*
*=                   (3) NOAA, 1998, unpublished as of 01/98		         =*
*=		             (4) Gierczak et al., 1998, also in IUPAC 05             =*
*=  Quantum yield:  Choice between                                           =*
*=                   (1) Gardiner et al, 1984                                =*
*=                   (2) IUPAC 97                                            =*
*=                   (3) McKeen et al., 1997 				                 =*
*=                   (4) Warneck (2001)                     		         =*
*=                   (5) Blitz et al. (2004)				                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=150)

      INTEGER i, n
      INTEGER n1, n2, n3
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy
      REAL sig
      REAL qs, qco, qch3co
      INTEGER ierr
      INTEGER iw

      REAL a, b, t
      REAL a0, a1, a2, a3, a4, a5
      REAL b0, b1, b2, b3, b4, b5
      REAL ca, cb, cc, cd, ce, cf
      REAL cg, ch, ci, cj, ck, cl
      REAL c1, c2, c3
      REAL AA0, AA1, AA2, AA3, AA4
      
      INTEGER mabs, myld
      
**************** CH3COCH3 photodissociation

      j = j + 1
      jlabel(j) = 'CH3COCH3 -> CH3CO + CH3'

* options
* mabs for cross sections
* myld for quantum yields

* Absorption:
* 1:  cross section from Calvert and  Pitts
* 2:  Martinez et al. 1991, also in IUPAC'97
* 3:  NOAA 1998, unpublished as of Jan 98.
* 4:  Gierczak et al. 1998, also in IUPAC'05

* Quantum yield
* 1:  Gardiner et al. 1984
* 2:  IUPAC 97
* 3:  McKeen, S. A., T. Gierczak, J. B. Burkholder, P. O. Wennberg, T. F. Hanisco,
*       E. R. Keim, R.-S. Gao, S. C. Liu, A. R. Ravishankara, and D. W. Fahey, 
*       The photochemistry of acetone in the upper troposphere:  a source of 
*       odd-hydrogen radicals, Geophys. Res. Lett., 24, 3177-3180, 1997.
* 4:  Warneck 2001, also in IUPAC'05
* 5:  Blitz et al. 2004

      mabs = 4
      myld = 5

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCH3/CH3COCH3_cp.abs',
     $        STATUS='old')
         DO i = 1, 6
            READ(kin,*)
         ENDDO
         n = 35
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 3.82E-21
         ENDDO
         CLOSE (kin)
         
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCH3/CH3COCH3_iup.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 96
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE (kin)
         
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCH3/CH3COCH3_noaa.abs',
     $        STATUS='old')
         DO i = 1, 12
            READ(kin,*)
         ENDDO
         n = 135
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i), y3(i)
            x2(i) = x1(i)
            x3(i) = x1(i)
         ENDDO
         CLOSE (kin)
         n1 = n
         n2 = n
         n3 = n
         
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,               0.,0.)
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
         CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF


         CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
         CALL addpnt(x3,y3,kdata,n3,               0.,0.)
         CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
         CALL addpnt(x3,y3,kdata,n3,           1.e+38,0.)
         CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 4) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCH3/CH3COCH3_gierczak98.abs',
     $        STATUS='old')
         DO i = 1, 12
            READ(kin,*)
         ENDDO
         
         n = 135

         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i)*1.e-20
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF
      
      ENDIF

      IF (myld .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCH3/CH3COCH3_iup.yld',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 9
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

      DO iw = 1, nw - 1

         DO i = 1, nz

            sig = yg(iw)

            IF(mabs .EQ. 3) THEN
               t = 298. - tlev(i)
               t = MIN(t, 298.-235.)
               t = MAX(t, 0.)

               sig = yg(iw)*(1. + yg2(iw)*t + yg3(iw)*t*t)

	    ELSEIF (mabs .EQ. 4) THEN
	       
	       a0 = 13.1852
	       a1 = -0.252309
	       a2 = 0.00192398
	       a3 = -7.30847e-6
	       a4 = 1.3832e-8
	       a5 = -1.04375e-11
               
             b0 = -0.0445238
             b1 = 0.000851385
             b2 = -6.48766e-6
             b3 = 2.46273e-8
             b4 = -4.65788e-11
	       b5 = 3.51253e-14
          
             ca = a0
	       cb = a1*wc(iw)
             cc = a2*(wc(iw)**2)
             cd = a3*(wc(iw)**3)
             ce = a4*(wc(iw)**4) 
             cf = a5*(wc(iw)**5)

	       cg = b0
             ch = b1*wc(iw)
             ci = b2*(wc(iw)**2)
             cj = b3*(wc(iw)**3)
             ck = b4*(wc(iw)**4)
             cl = b5*(wc(iw)**5)

             c1 = ca+cb+cc+cd+ce+cf
             c2 = cg+ch+ci+cj+ck+cl

	       yg2(iw) = yg1(iw)*(1+(c1*tlev(i))+(c2*(tlev(i)**2)))
	       
	       sig = yg1(iw)*(1+(c1*tlev(i))+(c2*(tlev(i)**2)))
	     
            ENDIF

            IF (myld .EQ. 1) THEN
               qy = 0.0766 + 0.09415*EXP(-airden(i)/3.222e18)

            ELSEIF (myld .EQ. 2) THEN
               qy = yg1(iw)

            ELSEIF (myld .EQ. 3) THEN
               IF (wc(iw) .LE. 292.) THEN
                  qy = 1.
               ELSEIF (wc(iw) .GE. 292.  .AND. wc(iw) .LT. 308. ) THEN
                  a = -15.696 + 0.05707*wc(iw)
                  b = EXP(-88.81+0.15161*wc(iw))
                  qy = 1./(a + b*airden(i))
               ELSEIF (wc(iw) .GE. 308.  .AND. wc(iw) .LT. 337. ) THEN
                  a = -130.2 + 0.42884*wc(iw)
                  b = EXP(-55.947+0.044913*wc(iw))
                  qy = 1./(a + b*airden(i))
               ELSEIF (wc(iw) .GE. 337.) THEN
                  qy = 0.
               ENDIF

               qy = max(0., qy)
               qy = min(1., qy)

            ELSEIF (myld .EQ. 4) THEN
               
	       yg3(iw) = (0.887)*(1/(1 + EXP((wc(iw)-307.5)/3))) + 0.113

	       yg4(iw) = 1/((1/yg3(iw))+(7.145e-7*airden(i))*
     $                 EXP(-(8780.6/wc(iw)))) 
     
	       qs = yg3(iw)
	       qy = yg4(iw)

            ELSEIF (myld .EQ. 5) THEN

	       a0 = 0.35*((tlev(i)/295)**-1.28)
	       b0 = 0.068*((tlev(i)/295)**-2.65)
             a1 = 1.6e-19*((tlev(i)/295)**-2.38)   
	       b1 = 0.55e-3*((tlev(i)/295)**-3.19)
	       a2 = 1.62e-17*((tlev(i)/295)**-10.03)
	       b2 = 1.79e-3*((tlev(i)/295)**-1.364)
	       a3 = 26.29*((tlev(i)/295)**-6.59)
	       b3 = 5.72e-7*((tlev(i)/295)**-2.93)
      	     c3 = 30006*((tlev(i)/295)**-0.064)
	       a4 = 1.67e-15*((tlev(i)/295)**-7.25)
	       b4 = 2.08e-3*((tlev(i)/295)**-1.16)

	       AA0 = (a0/(1 - a0))*EXP(b0*(wc(iw) - 248))
	       AA1 = a1*EXP((-b1)*((1.0e+7/wc(iw)) - 33113))
	       AA2 = a2*EXP((-b2)*((1.0e+7/wc(iw)) - 30488))
	       AA3 = a3*EXP((-b3)*(((1.0e+7/wc(iw)) - c3)**2))
	       AA4 = a4*EXP((-b4)*((1.0e+7/wc(iw)) - 30488))

	         IF (wc(iw) .GE. 279. .AND. wc(iw) .LT. 302.) THEN	       

	            yg3(iw) = 1/(1+AA0)
	            yg4(iw) = (1 - yg3(iw))/(1 + (AA1*airden(i)))

	         ELSEIF (wc(iw) .GE. 302. .AND. wc(iw) .LT. 327.5) THEN
	       
	            yg3(iw) = 1/(1+AA0)
       	  yg4(iw) = ((1 + (AA4*airden(i)) + AA3)/((1 + (AA2*airden(i))
     $                      + AA3)*(1 + (AA4*airden(i)))))*(1 - yg3(iw))

                  qco = yg3(iw)
                  qch3co = yg4(iw)
               
                  qy = qco + qch3co
               
	         ENDIF	         
            
               qy = max(0., qy)
               qy = min(1., qy)

          
            ENDIF

            sq(j,i,iw) = sig*qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r16(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3OOH photolysis: =*
*=         CH3OOH + hv -> CH3O + OH                                          =*
*=                                                                           =*
*=  Cross section: Choice between                                            =*
*=                  (1) JPL 97 recommendation (based on Vaghjiana and        =*
*=                      Ravishankara, 1989), 10 nm resolution                =*
*=                  (2) IUPAC 97 (from Vaghjiana and Ravishankara, 1989),    =*
*=                      5 nm resolution                                      =*
*=                  (3) Cox and Tyndall, 1978; only for wavelengths < 280 nm =*
*=                  (4) Molina and Arguello, 1979;  might be 40% too high    =*
*=  Quantum yield: Assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER ierr
      INTEGER idum
      INTEGER iw

      INTEGER mabs


**************** CH3OOH photodissociation
         j = j + 1
         jlabel(j) = 'CH3OOH -> CH3O + OH'

* mabs: Absorption cross section options:
* 1:  JPL data base (1985,92,94,97). 1997 is from  Vaghjiani and Ravishankara (1989), at
*     10 nm resolution
* 2:  IUPAC97 (from  Vaghjiani and Ravishankara (1989) at 5 nm resolution).
* 3:  Cox and Tyndall (1978), only for wavelengths < 280 nm
* 4:  Molina and Arguello (1979).  According to Vaghjiani and Ravishankara (1989), 
*     Molina and Arguello had a problem measuring CH3OOH, cross sections 40% too high.

      mabs = 2

      IF (mabs .EQ. 1) THEN

c         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_jpl85.abs',
c     $        STATUS='old')
c         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_jpl92.abs',
c     $        STATUS='old')
c         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_jpl94.abs',
c     $        STATUS='old')
         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_jpl94.abs',
     $        STATUS='old')
         READ(kin,*) idum, n
         DO i = 1, idum-2
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-20
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF (mabs .EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_iup.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 32
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-20
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF (mabs .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_ct.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 12
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF (mabs .EQ. 4) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3OOH/CH3OOH_ma.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 15
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1

         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r17(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3ONO2            =*
*=  photolysis:                                                              =*
*=          CH3ONO2 + hv -> CH3O + NO2                                       =*
*=                                                                           =*
*=  Cross section: Choice between                                            =*
*=                  (1) Calvert and Pitts, 1966                              =*
*=                  (2) Talukdar, Burkholder, Hunter, Gilles, Roberts,       =*
*=                      Ravishankara, 1997                                   =*
*=                  (3) IUPAC 97, table of values for 198K                   =*
*=                  (4) IUPAC 97, temperature-dependent equation             =*
*=                  (5) Taylor et al, 1980                                   =*
*=                  (6) fit from Roberts and Fajer, 1989                     =*
*=                  (7) Rattigan et al., 1992                                =*
*=                  (8) Libuda and Zabel, 1995                               =*
*=                  (9) IUPAC 02  temperature-dependent equation,            =* 
*=                      added by Judit Zador                                 =*
*=  Quantum yield: Assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER (kdata = 2000)

      INTEGER i, n
      INTEGER iw
      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg1(kw)
      REAL qy
      REAL sig
      INTEGER ierr

      INTEGER mabs, myld

**************** CH3ONO2 photodissociation

      j = j + 1
      jlabel(j) = 'CH3ONO2 -> CH3O + NO2'

* mabs: absorption cross section options:
* 1:  Calvert and  Pitts 1966
* 2:  Talukdar, Burkholder, Hunter, Gilles, Roberts, Ravishankara, 1997.
* 3:  IUPAC-97, table of values for 298K.
* 4:  IUPAC-97, temperature-dependent equation
* 5:  Taylor et al. 1980
* 6:  fit from Roberts and Fajer, 1989
* 7:  Rattigan et al. 1992
* 8:  Libuda and Zabel 1995
* 9:  IUPAC-05, temperature-dependent equation added by Judit Zador

      mabs = 9 ! original: 2

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_cp.abs',STATUS='old')
         DO i = 1, 3
            READ(kin,*)
         ENDDO
         n = 15
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 2) THEN

*        sigma(T,lambda) = sigma(298,lambda) * exp(B * (T-298))

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_tal.abs',STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 55
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
            y1(i) = y1(i) * 1.e-20
         ENDDO
         CLOSE (kin)

         n1 = n
         n2 = n
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
         CALL addpnt(x2,y2,kdata,n2,            1.e+38,y2(n2))
         CALL inter2(nw,wl,yg1,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF (mabs .EQ. 3) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_iup1.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 13
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i)*1e-20
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 4) THEN

*        sigma(T,lambda) = sigma(298,lambda) * 10**(B * T)

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_iup2.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 7
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
            y1(i) = y1(i) * 1.e-21
            y2(i) = y2(i) * 1.e-3
         ENDDO
         CLOSE (kin)

         n1 = n
         n2 = n
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),-36.)
         CALL addpnt(x1,y1,kdata,n1,               0.,-36.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),-36.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,-36.)
         CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
         CALL addpnt(x2,y2,kdata,n2,            1.e+38,y2(n2))
         CALL inter2(nw,wl,yg1,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF (mabs .EQ. 5) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_tay.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 13
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF (mabs .EQ. 6) THEN

         DO iw = 1, nw-1
            IF(wc(iw) .GT. 284.) THEN
               yg(iw) = EXP(-1.044e-3*wc(iw)*wc(iw) + 
     $              0.5309*wc(iw) - 112.4)
            ELSE
               yg(iw) = 0.
            ENDIF
         ENDDO

      ELSEIF (mabs .EQ. 7) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_rat.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 24
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF (mabs .EQ. 8) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_lib.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 1638
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)

         n1 = n
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF



      ELSEIF(mabs .EQ. 9) THEN

*        sigma(T,lambda) = sigma(298,lambda) * EXP(B * (T - 298))

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3ONO2_iup_zj.abs',
     $        STATUS='old')
         DO i = 1, 4
            READ(kin,*)
         ENDDO
         n = 19
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
            y1(i) = y1(i) * 1.e-20
            y2(i) = y2(i) * 1.e-3
         ENDDO
         CLOSE (kin)

         n1 = n
         n2 = n
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
         CALL addpnt(x2,y2,kdata,n2,            1.e+38,y2(n2))
         CALL inter2(nw,wl,yg1,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield = 1

      qy = 1.

      DO iw = 1, nw - 1
         sig = yg(iw)

         DO i = 1, nz
            
            IF(mabs .EQ. 2) THEN
               sig = yg(iw) * exp (yg1(iw) * (tlev(i)-298.))

            ELSEIF (mabs .EQ. 4) THEN
               sig = yg(iw)*10.**(yg1(iw)*tlev(i))

            ELSEIF (mabs .EQ. 9) THEN
            
               IF(wc(iw) .GE. 240. .AND. wc(iw) .LE. 340.)THEN ! arr

                   sig = yg(iw)*exp(yg1(iw)*(tlev(i)-298.))
                   
               ELSE !arr
               
                   sig = 0 !arr
               
               ENDIF    
            ENDIF

            sq(j,i,iw) = qy * sig

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r18(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for PAN photolysis:    =*
*=       PAN + hv -> Products                                                =*
*=                                                                           =*
*=  Cross section: from Talukdar et al., 1995                                =*
*=  Quantum yield: Assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      INTEGER iw
      INTEGER i, n
      INTEGER n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg(kw), yg2(kw)
      REAL qy
      REAL sig
      INTEGER ierr

**************** PAN photodissociation

      j = j+1
      jlabel(j) = 'CH3CO(OONO2) -> Products'

* cross section from Senum et al., 1984, J.Phys.Chem. 88/7, 1269-1270

C     OPEN(UNIT=kin,FILE='DATAJ1/RONO2/PAN_senum.abs',STATUS='OLD')
C     DO i = 1, 14
C        READ(kin,*)
C     ENDDO
C     n = 21
C     DO i = 1, n
C        READ(kin,*) x1(i), y1(i)
C        y1(i) = y1(i) * 1.E-20
C     ENDDO
C     CLOSE(kin)

C      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C      CALL addpnt(x1,y1,kdata,n,               0.,0.)
C      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
C      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
C      IF (ierr .NE. 0) THEN
C         WRITE(*,*) ierr, jlabel(j)
C         STOP
C      ENDIF

* cross section from 
*      Talukdar et al., 1995, J.Geophys.Res. 100/D7, 14163-14174
*      JPL-2005 recommended data 
*      Agree well with H.G. Libuda and F. Zabel, "UV absorption cross section 
*      of acetyl peroxynitrate and trifluoroacetyl peroxynitrate," Ber. Bunsenges. 
*      Phys. Chem. 99, 1205-1213 (1995) 
*      NB: IUPAC 2005 use an average of the two datasets

      OPEN(UNIT=kin,FILE='DATAJ1/RONO2/PAN_talukdar.abs',STATUS='OLD')
      DO i = 1, 14
         READ(kin,*)
      ENDDO
      n = 78
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), y2(i)
         y1(i) = y1(i) * 1.E-20
         y2(i) = y2(i) * 1E-3
         x2(i) = x1(i)
      ENDDO
      n2 = n
      CLOSE(kin)
 
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,          0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,      1.e+38,0.)
      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield
* yet unknown, but assumed to be 1.0 (Talukdar et al., 1995)

      qy = 1.0

      DO iw = 1, nw-1
        DO i = 1, nz

          sig = yg(iw) * EXP(yg2(iw)*(tlev(i)-298.))

          sq(j,i,iw) = qy * sig

        ENDDO
      ENDDO 

      END

*=============================================================================*

      SUBROUTINE r19(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CCl2O photolysis:  =*
*=        CCl2O + hv -> Products                                             =*
*=                                                                           =*
*=  Cross section: JPL 94 recommendation                                     =*
*=  Quantum yield: Unity (Calvert and Pitts)                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

************* CCl2O photodissociation

      j = j+1
      jlabel(j) = 'CCl2O -> Products'

*** cross sections from JPL94 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/CCl2O_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

*** quantum yield unity (Calvert and Pitts)
      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r20(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CCl4 photolysis:   =*
*=      CCl4 + hv -> Products                                                =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CCl4 photodissociation
      
      j = j+1
      jlabel(j) = 'CCl4 -> Products'

*** cross sections from JPL97 recommendation (identical to 94 data)

      OPEN(kin,FILE='DATAJ1/ABS/CCl4_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

*** quantum yield assumed to be unity
      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r21(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CClFO photolysis:  =*
*=         CClFO + hv -> Products                                            =*
*=  Cross section: from JPL 97                                               =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CClFO photodissociation

      j = j+1
      jlabel(j) = 'CClFO -> Products'

*** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CClFO_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

*** quantum yield unity
      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r22(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CF2O photolysis:   =*
*=        CF2O + hv -> Products                                              =*
*=  Cross section:  from JPL 97 recommendation                               =*
*=  Quantum yield:  unity (Nolle et al.)                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CF2O photodissociation

      j = j+1
      jlabel(j) = 'CF2O -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CF2O_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

*** quantum yield unity (Nolle et al.)
      qy = 1.
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r23(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CFC-113 photolysis:=*
*=          CF2ClCFCl2 + hv -> Products                                      =*
*=  Cross section:  from JPL 97 recommendation, linear interp. between       =*
*=                  values at 210 and 295K                                   =*
*=  Quantum yield:  assumed to be unity                                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER iz
      INTEGER ierr
      REAL slope

**************************************************************
************* CF2ClCFCl2 (CFC-113) photodissociation

      j = j+1
      jlabel(j) = 'CF2ClCFCl2 (CFC-113) -> Products'

*** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CFC-113_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n

** sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,         1E38,0.)

      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

* sigma @ 210 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,         1E38,0.)

      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
 
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

*** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = MAX(210.,MIN(tlev(iz),295.))
        slope = (t-210.)/(295.-210.)
        DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg2(iw) + slope*(yg1(iw)-yg2(iw)))
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r24(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CFC-144 photolysis:=*
*=              CF2ClCF2Cl + hv -> Products                                  =*
*=  Cross section: from JPL 97 recommendation, linear interp. between values =*
*=                 at 210 and 295K                                           =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz
      REAL slope

**************************************************************
************* CF2ClCF2Cl (CFC-114) photodissociation

      j = j+1
      jlabel(j) = 'CF2ClCF2Cl (CFC-114) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CFC-114_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        x2(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n

** sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,         1E38,0.)

      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

* sigma @ 210 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,         1E38,0.)

      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

*** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = MAX(210.,MIN(tlev(iz),295.))
        slope = (t-210.)/(295.-210.)
        DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg2(iw) + slope*(yg1(iw)-yg2(iw)))
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r25(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CFC-115 photolysis =*
*=             CF3CF2Cl + hv -> Products                                     =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CF3CF2Cl (CFC-115) photodissociation
      
      j = j+1
      jlabel(j) = 'CF3CF2Cl (CFC-115) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CFC-115_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
    
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r26(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CFC-111 photolysis =*
*=          CCl3F + hv -> Products                                           =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CCl3F (CFC-11) photodissociation
      
      j = j+1
      jlabel(j) = 'CCl3F (CFC-11) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CFC-11_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

** sigma @ 298 K

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF 

**** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = 1E-04 * (tlev(iz)-298.)
        DO iw = 1, nw-1
          sq(j,iz,iw) = qy * yg(iw) * EXP((wc(iw)-184.9) * t)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r27(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CFC-112 photolysis:=*
*=         CCl2F2 + hv -> Products                                           =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CCl2F2 (CFC-12) photodissociation
      
      j = j+1
      jlabel(j) = 'CCl2F2 (CFC-12) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CFC-12_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

** sigma @ 298 K

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = 1E-04 * (tlev(iz)-298.) 
        DO iw = 1, nw-1
          sq(j,iz,iw) = qy * yg(iw) * EXP((wc(iw)-184.9) * t)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r28(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3Br photolysis:  =*
*=         CH3Br + hv -> Products                                            =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CH3Br photodissociation

* data from JPL97 (identical to 94 recommendation)
      
      j = j+1
      jlabel(j) = 'CH3Br -> Products'
      OPEN(kin,FILE='DATAJ1/ABS/CH3Br_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r29(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3CCl3 photolysis =*
*=           CH3CCl3 + hv -> Products                                        =*
*=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
*=                 of data at 210, 250, and 295K                             =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      INTEGER n1, n2, n3
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz
      REAL slope

**************************************************************
************* CH3CCl3 photodissociation
      
      j = j+1
      jlabel(j) = 'CH3CCl3 -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CH3CCl3_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        y3(i) = y3(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n

** sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,         1E38,0.)

      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

** sigma @ 250 K
      
      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,         1E38,0.)

      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
      
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

** sigma @ 210 K

      CALL addpnt(x3,y3,kdata,n3, x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,           0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,         1E38,0.)

      CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = MIN(295.,MAX(tlev(iz),210.))
        IF (t .LE. 250.) THEN
          slope = (t-210.)/(250.-210.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg3(iw) + slope*(yg2(iw)-yg3(iw)))
          ENDDO
        ELSE
          slope = (t-250.)/(295.-250.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg2(iw) + slope*(yg1(iw)-yg2(iw)))
          ENDDO
        ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r30(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3Cl photolysis:  =*
*=            CH3Cl + hv -> Products                                         =*
*=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
*=                 from values at 255, 279, and 296K                         =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      INTEGER n1, n2, n3
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz
      REAL slope

**************************************************************
************* CH3Cl photodissociation

      j = j+1
      jlabel(j) = 'CH3Cl -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/CH3Cl_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        y3(i) = y3(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n

** sigma @ 296 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,         1E38,0.)

      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

** sigma @ 279 K
  
      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,         1E38,0.)

      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

** sigma @ 255 K

      CALL addpnt(x3,y3,kdata,n3, x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,           0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,         1E38,0.)

      CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)

      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
        t = MAX(255.,MIN(tlev(i),296.))
        IF (t .LE. 279.) THEN
          slope = (t-255.)/(279.-255.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg3(iw)+slope*(yg2(iw)-yg3(iw)))
          ENDDO
        ELSE
          slope = (t-279.)/(296.-279.)
          DO iw = 1, nw-1
            sq(j,iz,iw) = qy * (yg2(iw)+slope*(yg1(iw)-yg2(iw)))
          ENDDO
        ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r31(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for ClOO photolysis:   =*
*=          ClOO + hv -> Products                                            =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

C     INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* ClOO photodissociation

      j = j+1
      jlabel(j) = 'ClOO -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/ClOO_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
 
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r32(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-123 photolysis=*
*=       CF3CHCl2 + hv -> Products                                           =*
*=  Cross section: from Orlando et al., 1991                                 =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* local

      REAL qy
      REAL t
      INTEGER i, iw, idum
      INTEGER iz, k
      REAL lambda, sum
      CHARACTER*120 inline

      REAL coeff(4,3), TBar, LBar

**************************************************************
************* CF3CHCl2 (HCFC-123) photodissociation
      
      j = j+1
      jlabel(j) = 'CF3CHCl2 (HCFC-123) -> Products'

**** cross sections from JPL94 recommendation

C     OPEN(kin,FILE='DATAJ1/ABS/HCFC-123_jpl94.abs',STATUS='OLD')
C     READ(kin,*) idum, n
C     DO i = 1, idum-2
C       READ(kin,*)
C     ENDDO
C     DO i = 1, n
C       READ(kin,*) x1(i), y1(i)
C       y1(i) = y1(i) * 1E-20
C     ENDDO
C     CLOSE(kin)

C     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,          0.,0.)
C     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

C     CALL inter2(nw,wl,yg,n,x1,y1,ierr)

C     IF (ierr .NE. 0) THEN
C        WRITE(*,*) ierr, jlabel(j)
C        STOP
C     ENDIF

**** quantum yield assumed to be unity
C     qy = 1.

C     DO iw = 1, nw-1
C       DO iz = 1, nz
C         sq(j,iz,iw) = qy * yg(iw)
C       ENDDO
C     ENDDO


**** cross section from Orlando et al., 1991

      OPEN(kin,FILE='DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,100) inline
 100  FORMAT(A120)
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      LBar = 206.214

**** quantum yield assumed to be unity

      qy = 1. 

      DO iw = 1, nw-1

        lambda = wc(iw)

C use parameterization only up to 220 nm, as the error bars associated with
C the measurements beyond 220 nm are very large (Orlando, priv.comm.)

        IF (lambda .GE. 190. .AND. lambda .LE. 220.) THEN
          DO iz = 1, nz
             t = MIN(295.,MAX(tlev(i),203.))-TBar
             sum = 0.
             DO i = 1, 4
                sum = (coeff(i,1)+t*(coeff(i,2)+t*coeff(i,3))) *
     >                (lambda-LBar)**(i-1) + sum
             ENDDO 
             sq(j,iz,iw) = qy * EXP(sum)
          ENDDO
        ELSE
          DO iz = 1, nz
            sq(j,iz,iw) = 0.
          ENDDO
        ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r33(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-124 photolysis=*
*=        CF3CHFCl + hv -> Products                                          =*
*=  Cross section: from Orlando et al., 1991                                 =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays


* local

      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER iz, k
      REAL lambda, sum
      CHARACTER*120 inline

      REAL coeff(4,3), TBar, LBar

**************************************************************
************* CF3CHFCl (HCFC-124) photodissociation
      
      j = j+1
      jlabel(j) = 'CF3CHFCl (HCFC-124) -> Products'

**** cross sections from JPL94 recommendation

C     OPEN(kin,FILE='DATAJ1/ABS/HCFC-124_jpl94.abs',STATUS='OLD')
C     READ(kin,*) idum, n
C     DO i = 1, idum-2
C       READ(kin,*)
C     ENDDO
C     DO i = 1, n
C       READ(kin,*) x1(i), y1(i)
C       y1(i) = y1(i) * 1E-20
C     ENDDO
C     CLOSE(kin)

C     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,          0.,0.)
C     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

C     CALL inter2(nw,wl,yg,n,x1,y1,ierr)

C     IF (ierr .NE. 0) THEN
C       WRITE(*,*) ierr, jlabel(j)
C       STOP
C     ENDIF

**** quantum yield assumed to be unity
C     qy = 1.

C     DO iw = 1, nw-1
C       DO iz = 1, nz
C         sq(j,iz,iw) = qy * yg(iw)
C       ENDDO
C     ENDDO

**** cross section from Orlando et al., 1991

      OPEN(kin,FILE='DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      idum = idum+5
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,100) inline
 100  FORMAT(A120)
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      LBar = 206.214

**** quantum yield assumed to be unity

      qy = 1. 

      DO iw = 1, nw-1
        lambda = wc(iw)
        IF (lambda .GE. 190. .AND. lambda .LE. 230.) THEN
          DO iz = 1, nz
             t = MIN(295.,MAX(tlev(i),203.))-TBar
             sum = 0.
             DO i = 1, 4
                sum = (coeff(i,1)+t*(coeff(i,2)+t*coeff(i,3))) *
     >                (lambda-LBar)**(i-1) + sum
             ENDDO
             sq(j,iz,iw) = qy * EXP(sum)
          ENDDO
        ELSE
          DO iz = 1, nz
            sq(j,iz,iw) = 0.
          ENDDO
        ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r34(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-141b          =*
*=  photolysis:                                                              =*
*=         CH3CFCl2 + hv -> Products                                         =*
*=  Cross section: from JPL97 recommendation                                 =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CH3CFCl2 (HCFC-141b) photodissociation

      j = j+1
      jlabel(j) = 'CH3CFCl2 (HCFC-141b) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-141b_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r35(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-142b          =*
*=  photolysis:                                                              =*
*=          CH3CF2Cl + hv -> Products                                        =*
*=  Cross section: from Orlando et al., 1991                                 =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* local

      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz, k
      REAL lambda, sum
      CHARACTER*80 inline

      REAL coeff(4,3), TBar, LBar

**************************************************************
************* CH3CF2Cl (HCFC-142b) photodissociation

      j = j+1
      jlabel(j) = 'CH3CF2Cl (HCFC-142b) -> Products'

**** cross sections from JPL94 recommendation

C     OPEN(kin,FILE='DATAJ1/ABS/HCFC-142b_jpl94.abs',STATUS='OLD')
C     READ(kin,*) idum, n
C     DO i = 1, idum-2
C       READ(kin,*)
C     ENDDO
C     DO i = 1, n
C       READ(kin,*) x1(i), y1(i)
C       y1(i) = y1(i) * 1E-20
C     ENDDO
C     CLOSE(kin)

C     CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,          0.,0.)
C     CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
C     CALL addpnt(x1,y1,kdata,n,        1E38,0.)

C     CALL inter2(nw,wl,yg,n,x1,y1,ierr)

C     IF (ierr .NE. 0) THEN
C       WRITE(*,*) ierr, jlabel(j)
C       STOP
C     ENDIF

**** quantum yield assumed to be unity
C     qy = 1.

C     DO iw = 1, nw-1
C       DO iz = 1, nz
C         sq(j,iz,iw) = qy * yg(iw)
C       ENDDO
C     ENDDO

**** cross section from Orlando et al., 1991

      OPEN(kin,FILE='DATAJ1/ABS/HCFCs_orl.abs',STATUS='OLD')
      READ(kin,*) idum
      idum = idum+10
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      READ(kin,100) inline
 100  FORMAT(A80)
      READ(inline(6:),*) TBar,i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      READ(kin,*)           i,(coeff(i,k),k=1,3)
      CLOSE(kin)

      LBar = 206.214

**** quantum yield assumed to be unity

      qy = 1.

      DO iw = 1, nw-1

        lambda = wc(iw)
        IF (lambda .GE. 190. .AND. lambda .LE. 230.) THEN

          DO iz = 1, nz
             t = MIN(295.,MAX(tlev(i),203.))-TBar
             sum = 0.
             DO i = 1, 4
                sum = (coeff(i,1)+t*(coeff(i,2)+t*coeff(i,3))) *
     >                (lambda-LBar)**(i-1) + sum
             ENDDO

* offeset exponent by 40 (exp(-40.) = 4.248e-18) to prevent exp. underflow errors
* on some machines.

c             sq(j,iz,iw) = qy * EXP(sum)
             sq(j,iz,iw) = qy * 4.248e-18 * EXP(sum + 40.)

          ENDDO

        ELSE
          DO iz = 1, nz
            sq(j,iz,iw) = 0.
          ENDDO
        ENDIF

      ENDDO


      END

*=============================================================================*

      SUBROUTINE r36(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-225ca         =*
*=  photolysis:                                                              =*
*=           CF3CF2CHCl2 + hv -> Products                                    =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CF3CF2CHCl2 (HCFC-225ca) photodissociation
       
      j = j+1
      jlabel(j) = 'CF3CF2CHCl2 (HCFC-225ca) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-225ca_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r37(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-225cb         =*
*=  photolysis:                                                              =*
*=          CF2ClCF2CHFCl + hv -> Products                                   =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CF2ClCF2CHFCl (HCFC-225cb) photodissociation

      j = j+1
      jlabel(j) = 'CF2ClCF2CHFCl (HCFC-225cb) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-225cb_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r38(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HCFC-22 photolysis =*
*=          CHClF2 + hv -> Products                                          =*
*=  Cross section: from JPL 97 recommendation, piecewise linear interp.      =*
*=                 from values at 210, 230, 250, 279, and 295 K              =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL qy
      REAL t
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz
      REAL slope

**************************************************************
************* CHClF2 (HCFC-22) photodissociation
       
      j = j+1
      jlabel(j) = 'CHClF2 (HCFC-22) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/HCFC-22_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i), y4(i), y5(i)
        y1(i) = y1(i) * 1E-20
        y2(i) = y2(i) * 1E-20
        y3(i) = y3(i) * 1E-20
        y4(i) = y4(i) * 1E-20
        y5(i) = y5(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
        x4(i) = x1(i)
        x5(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      n2 = n
      n3 = n
      n4 = n
      n5 = n

** sigma @ 295 K

      CALL addpnt(x1,y1,kdata,n1, x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,         1E38,0.)

      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

** sigma @ 270 K

      CALL addpnt(x2,y2,kdata,n2, x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,         1E38,0.)

      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

** sigma @ 250 K

      CALL addpnt(x3,y3,kdata,n3, x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,           0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,         1E38,0.)

      CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

** sigma @ 230 K

      CALL addpnt(x4,y4,kdata,n4, x4(1)*(1.-deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,           0.,0.)
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,         1E38,0.)

      CALL inter2(nw,wl,yg4,n4,x4,y4,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

** sigma @ 210 K

      CALL addpnt(x5,y5,kdata,n5, x5(1)*(1.-deltax),0.)
      CALL addpnt(x5,y5,kdata,n5,           0.,0.)
      CALL addpnt(x5,y5,kdata,n5,x5(n5)*(1.+deltax),0.)
      CALL addpnt(x5,y5,kdata,n5,         1E38,0.)

      CALL inter2(nw,wl,yg5,n5,x5,y5,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.

      DO iz = 1, nz
         t = MIN(295.,MAX(tlev(iz),210.))
         IF (t .LE. 230.) THEN
            slope = (t-210.)/(230.-210.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = qy * (yg5(iw)+slope*(yg4(iw)-yg5(iw)))
            ENDDO
         ELSEIF (t .LE. 250.) THEN
            slope = (t-230.)/(250.-230.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = qy * (yg4(iw)+slope*(yg3(iw)-yg4(iw)))
            ENDDO
         ELSEIF (t .LE. 270.) THEN
            slope = (t-250.)/(270.-250.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = qy * (yg3(iw)+slope*(yg2(iw)-yg3(iw)))
            ENDDO
         ELSE
            slope = (t-270.)/(295.-270.)
            DO iw = 1, nw-1
              sq(j,iz,iw) = qy * (yg2(iw)+slope*(yg1(iw)-yg2(iw)))
            ENDDO
         ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r39(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for HO2 photolysis:    =*
*=          HO2 + hv -> OH + O                                               =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed shape based on work by Lee, 1982; normalized      =*
*=                 to unity at 248 nm                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* HO2 photodissociation

      j = j+1
      jlabel(j) = 'HO2 -> OH + O'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/HO2_jpl94.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield:  absolute quantum yield has not been reported yet, but
****                 Lee measured a quantum yield for O(1D) production at 248
****                 nm that was 15 time larger than at 193 nm
**** here:  a quantum yield of unity is assumed at 248 nm and beyond, for
****        shorter wavelengths a linear decrease with lambda is assumed

      DO iw = 1, nw-1
         IF (wc(iw) .GE. 248.) THEN
            qy = 1.
         ELSE
            qy = 1./15. + (wc(iw)-193.)*(14./15.)/(248.-193.)
            qy = MAX(qy,0.)
         ENDIF
         DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r40(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) Halon-1202 photolysis: =*
*=         CF2Br2 + hv -> Products                                           =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: unity (Molina and Molina)                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CF2Br2 (Halon-1202) photodissociation
      
      j = j+1
      jlabel(j) = 'CF2Br2 (Halon-1202) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/Halon-1202_jpl97.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield unity (Molina and Molina)
      qy = 1.
     
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r41(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for Halon-1211         =*
*=  photolysis:                                                              =*
*=           CF2ClBr + hv -> Products                                        =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CF2BrCl (Halon-1211) photodissociation

      j = j+1
      jlabel(j) = 'CF2BrCl (Halon-1211) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/Halon-1211_jpl97.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)
     
      CALL inter2(nw,wl,yg,n,x1,y1,ierr) 

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.
     
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r42(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for Halon-1301         =*
*=  photolysis:                                                              =*
*=         CF3Br + hv -> Products                                            =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CF3Br (Halon-1301) photodissociation

      j = j+1
      jlabel(j) = 'CF3Br (Halon-1301) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/Halon-1301_jpl97.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)
    
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.
     
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r43(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for Halon-2402         =*
*=  photolysis:                                                              =*
*=           CF2BrCF2Br + hv -> Products                                     =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

**************************************************************
************* CF2BrCF2Br (Halon-2402) photodissociation

      j = j+1
      jlabel(j) = 'CF2BrCF2Br (Halon-2402) -> Products'

**** cross sections from JPL97 recommendation (identical to 94 recommendation)

      OPEN(kin,FILE='DATAJ1/ABS/Halon-2402_jpl97.abs',STATUS='OLD')
      READ(kin,*) idum, n
      DO i = 1, idum-2
        READ(kin,*)
      ENDDO
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)+deltax,0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)
    
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF

**** quantum yield assumed to be unity
      qy = 1.
     
      DO iw = 1, nw-1
        DO iz = 1, nz
           sq(j,iz,iw) = qy * yg(iw)
        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r44(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for N2O photolysis:    =*
*=              N2O + hv -> N2 + O(1D)                                       =*
*=  Cross section: from JPL 97 recommendation                                =*
*=  Quantum yield: assumed to be unity, based on Greenblatt and Ravishankara =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* local

      REAL qy
      REAL a, b, c
      REAL a0, a1, a2, a3, a4
      REAL b0, b1, b2, b3
      REAL t
      INTEGER iw, iz
      REAL lambda

**************************************************************
************* N2O photodissociation

      j = j+1
      jlabel(j) = 'N2O -> N2 + O(1D)'

**** cross sections according to JPL97 recommendation (identical to 94 rec.)
**** see file DATAJ1/ABS/N2O_jpl94.abs for detail

      A0 = 68.21023                
      A1 = -4.071805               
      A2 = 4.301146E-02            
      A3 = -1.777846E-04           
      A4 = 2.520672E-07

      B0 = 123.4014
      B1 = -2.116255
      B2 = 1.111572E-02
      B3 = -1.881058E-05

**** quantum yield of N(4s) and NO(2Pi) is less than 1% (Greenblatt and
**** Ravishankara), so quantum yield of O(1D) is assumed to be unity
      qy = 1.

      DO iw = 1, nw-1
         lambda = wc(iw)   
         IF (lambda .GE. 173. .AND. lambda .LE. 240.) THEN
           DO iz = 1, nz
             t = MAX(194.,MIN(tlev(iz),320.))
             A = (((A4*lambda+A3)*lambda+A2)*lambda+A1)*lambda+A0
             B = (((B3*lambda+B2)*lambda+B1)*lambda+B0)
             B = (t-300.)*EXP(B)
             sq(j,iz,iw) = qy * EXP(A+B)
           ENDDO
         ELSE
           DO iz = 1, nz
             sq(j,iz,iw) = 0.
           ENDDO 
         ENDIF
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r45(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for ClONO2 photolysis: =*
*=        ClONO2 + hv -> Products                                            =*
*=                                                                           =*
*=  Cross section: JPL 97 recommendation                                     =*
*=  Quantum yield: JPL 97 recommendation                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=150)

      REAL x1(kdata),x2(kdata),x3(kdata)
      REAL y1(kdata),y2(kdata),y3(kdata)
      INTEGER n1, n2, n3

* local

      REAL yg1(kw), yg2(kw), yg3(kw)
      REAL qy1, qy2
      REAL xs 
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

************* ClONO2 photodissociation

      j = j+1
      jlabel(j) = 'ClONO2 -> Cl + NO3'

*** cross sections from JPL97 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/ClONO2_jpl97.abs',STATUS='OLD')
      n = 119
      DO i = 1, n
        READ(kin,*) x1(i), y1(i), y2(i), y3(i)
        y1(i) = y1(i) * 1E-20
        x2(i) = x1(i)
        x3(i) = x1(i)
      ENDDO
      CLOSE(kin)

      n1 = n
      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,          0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,        1E38,0.)
      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      n2 = n
      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,          0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,        1E38,0.)
      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      n3 = n
      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,          0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,        1E38,0.)
      CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw-1

*** quantum yields (from jpl97)

         IF( wc(iw) .LT. 308.) THEN
            qy1 = 0.6
         ELSEIF( (wc(iw) .GE. 308) .AND. (wc(iw) .LE. 364.) ) THEN
            qy1 = 7.143e-3 * wc(iw) - 1.6
         ELSEIF( wc(iw) .GT. 364. ) THEN
            qy1 = 1.0
         ENDIF
         qy2 = 1. - qy1
         
* compute T-dependent cross section

         DO iz = 1, nz
            xs = yg1(iw)*( 1. + 
     $           yg2(iw)*(tlev(iz)-296) + 
     $           yg3(iw)*(tlev(iz)-296)*(tlev(iz)-296))
            sq(j,iz,iw) = qy1 * xs
            sq(j+1,iz,iw) = qy2 * xs

C            if(iz .eq. 1) write(33,333) wc(iw), xs
C 333        format(0pf8.3,1pe11.4)

         ENDDO
      ENDDO

      j = j+1
      jlabel(j) = 'ClONO2 -> ClO + NO2'

      END

*=============================================================================*

      SUBROUTINE r46(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for BrONO2 photolysis: =*
*=        BrONO2 + hv -> Products                                            =*
*=                                                                           =*
*=  Cross section: JPL 03 recommendation                                     =*
*=  Quantum yield: JPL 03 recommendation                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=100)

      REAL x1(kdata)
      REAL y1(kdata)
      INTEGER n1, n2, n3

* local

      REAL yg1(kw)
      REAL qy1, qy2
      INTEGER i, iw, n, idum
      INTEGER ierr
      INTEGER iz

************* BrONO2 photodissociation

      j = j+1
      jlabel(j) = 'BrONO2 -> BrO + NO2'
      j = j+1
      jlabel(j) = 'BrONO2 -> Br + NO3'


*** cross sections from JPL03 recommendation

      OPEN(kin,FILE='DATAJ1/ABS/BrONO2_jpl03.abs',STATUS='OLD')
      DO i = 1, 13
         READ(kin,*)
      ENDDO
      n = 61
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      n1 = n
      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,          0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,        1E38,0.)
      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

*** quantum yields (from jpl97)

      qy1 = 0.71
      qy2 = 0.29
      DO iw = 1, nw-1
         DO iz = 1, nz
            sq(j-1,iz,iw) = qy1 * yg1(iw)
            sq(j,iz,iw) = qy2 * yg1(iw)
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r47(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for Cl2 photolysis:    =*
*=        Cl2 + hv -> 2 Cl                                                   =*
*=                                                                           =*
*=  Cross section: JPL 97 recommendation                                     =*
*=  Quantum yield: 1     (Calvert and Pitts, 1966)                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=150)

      INTEGER i, n
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER iz, iw
      INTEGER ierr


************* CL2 photodissociation

      j = j+1
      jlabel(j) = 'Cl2 -> Cl + Cl'

*** cross sections from JPL97 recommendation (as tab by Finlayson-Pitts
* and Pitts, 1999.

      OPEN(kin,FILE='DATAJ1/ABS/CL2_fpp.abs',STATUS='OLD')
      do i = 1, 5
         read(kin,*)
      enddo
      n = 22
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

*** quantum yield = 1 (Calvert and Pitts, 1966)

      qy = 1.
      DO iw = 1, nw-1
         DO iz = 1, nz
            sq(j,iz,iw) = qy * yg(iw)
         ENDDO
      ENDDO

      RETURN
      END

*=============================================================================*

      SUBROUTINE r101(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH2(OH)CHO     =*
*=  (glycolaldehye, hydroxy acetaldehyde) photolysis:                        =*
*=           CH2(OH)CHO + hv -> Products                                     =*
*=                                                                           =*
*=  Cross section from                                                       =*
*= (1) The Atmospheric Chemistry of Glycolaldehyde, C. Bacher, G. S. Tyndall =*
*= and J. J. Orlando, J. Atmos. Chem., 39 (2001) 171-189. (IUPAC 2005)       =*
*= (2) I. Magneron, A. Mellouki, G. Le Bras, G.K. Moortgat, A. Horowitz,     =*
*= and K. Wirtz, "Photolysis and OH-initiated oxidation of glycolaldehyde    =*
*= under atmospheric conditions," J. Phys. Chem. A 109, 4552-4561 (2005).    =*
*= (JPL 2005)                                                                =*                                                                      
*=  Quantum yield from                                                       =*
*=  (1) about 0.75 (>0.5) from Bacher et al. (2001) (IUPAC 05)               =*
*=  (2) 1.0 from Magneron et al. (2005) (hinted at in Bacher et al. (2001))  =*
*=                                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=300)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

************************* CH2(OH)CHO photolysis
* 1:  CH2(OH)CHO

      j = j+1
      jlabel(j) = 'CH2(OH)CHO -> Products'

      mabs=1
      myld=1
      
      IF(mabs. EQ. 1)THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH2OHCHO/glycolaldehyde.abs',
     $      STATUS='old')
         DO i = 1, 15
            READ(kin,*)
         ENDDO
         n = 131
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x,y,ierr)
          IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
          ENDIF

      ELSEIF(mabs. EQ. 2)THEN
 
         OPEN(UNIT=kin,FILE='DATAJ1/CH2OHCHO/glycald_jpl05.abs',
     $      STATUS='old')
         DO i = 1, 9
            READ(kin,*)
         ENDDO
         n = 60
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x,y,ierr)
          IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
          ENDIF
      ENDIF       

* combine:
*
* NB:  IUPAC 2005 recommend a qy=0.75 (+-0.25) from Bacher et al. 2001.
*      Magneron et al. 2005 speculate a qy of 0.1 for the CH3OH channel
*      and that OH is formed in the photolysis of glycolaldehyde.

      DO iw = 1, nw - 1
        DO i = 1, nz

           IF(myld. EQ. 1)THEN

              qy = 0.75

              sq(j,i,iw) = yg(iw) * qy
          
           ELSEIF(myld. EQ. 2)THEN

              qy = 1.0

              sq(j,i,iw) = yg(iw) * qy

           ENDIF

        ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r102(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3COCOCH3     =*
*=  (biacetyl) photolysis:                                                   =*
*=           CH3COCOCH3 + hv -> Products                                     =*
*=                                                                           =*
*=  Cross section from either                                                =*
*= 1.  Plum et al., Environ. Sci. Technol., Vol. 17, No. 8, 1983, p.480      =*
*= 2.  Horowitz et al., J. Photochem Photobio A, 146, 19-27, 2001.           =*
*=                                                                           =*
*=  Quantum yield =0.158                                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=300)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER ierr
      INTEGER iw
      INTEGER mabs

************************* CH3COCOCH3 photolysis
* 1:  CH3COCOCH3

* Cross section data bases:
* mabs = 1 Plum et al.
* mabs = 2 Horowitz et al.

      mabs = 2

      j = j+1
      jlabel(j) = 'CH3COCOCH3 -> Products'

      IF( mabs. EQ. 1) THEN
         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCOCH3/biacetyl_plum.abs',
     $        STATUS='old')
         DO i = 1, 7
            READ(kin,*)
         ENDDO
         n = 55
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs. EQ. 2) THEN
         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCOCH3/biacetyl_horowitz.abs',
     $        STATUS='old')
         DO i = 1, 8
            READ(kin,*)
         ENDDO
         n = 287
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF
            
      ENDIF

* quantum yield from Plum et al.

      qy = 0.158

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw) * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r103(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3COCHCH2     =*
*=  Methyl vinyl ketone photolysis:                                          =*
*=           CH3COCHCH2 + hv -> Products                                     =*
*=                                                                           =*
*=  Cross section from                                                       =*
*= (1) W. Schneider and G. K. Moorgat, priv. comm, MPI Mainz 1989 as reported=*
*= by Roeth, E.-P., R. Ruhnke, G. Moortgat, R. Meller, and W. Schneider,     =*
*= UV/VIS-Absorption Cross Sections and QUantum Yields for Use in            =*
*= Photochemistry and Atmospheric Modeling, Part 2: Organic Substances,      =*
*= Forschungszentrum Julich, Report Jul-3341, 1997.                          =*
*= (2) T. Gierczak, J.B. Burkholder, R.K. Talukdar, A. Mellouki, S.B. Barone,=* 
*= and A.R. Ravishankara: Atmospheric fate of methyl vinyl ketone and        =*
*= methacrolein, J. Photochem. Photobiol. A: Chem. 110, 1-10 (1997)          =*
*=                                                                           =*                                                                          
*=  Quantum yield from                                                       =*
*= (2) T. Gierczak, J.B. Burkholder, R.K. Talukdar, A. Mellouki, S.B. Barone,=* 
*= and A.R. Ravishankara: Atmospheric fate of methyl vinyl ketone and        =*
*= methacrolein, J. Photochem. Photobiol. A: Chem. 110, 1-10 (1997)          =*
*= values recommended by IUPAC 2005                                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=20000)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw), yg1(kw)
      REAL qy
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

************************* CH3COCHCH2 photolysis

      j = j+1

      jlabel(j) = 'CH3COCHCH2 -> Products'

* Cross section:
* mabs = 1 Schneider and Moorgat
* mabs = 2 Gierczak et al (1997)

      mabs = 2

      IF( mabs. EQ. 1) THEN
         OPEN(UNIT=kin,FILE='DATAJ1/ABS/methylvinylketone.abs',
     $        STATUS='old')
         DO i = 1, 9
            READ(kin,*)
         ENDDO
         n = 19682
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs. EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/mvk_gierczak97.abs',
     $        STATUS='old')
         DO i = 1, 15
            READ(kin,*)
         ENDDO
         n = 146
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF
      ENDIF

* quantum yield from
* Gierczak, T., J. B. Burkholder, R. K. Talukdar, A. Mellouki, S. B. Barone,
* and A. R. Ravishankara, Atmospheric fate of methyl vinyl ketone and methacrolein,
* J. Photochem. Photobiol A: Chemistry, 110 1-10, 1997.
* values as recommended by IUPAC 2005
*
* depends on pressure and wavelength, set upper limit to 1.0
*
* NB: RADICAL project (2000) reported an overall upper limit to the quantum yield of <0.004.
*
* Pinho et al. show that the change in quantum yield to the pressure dependent Gierczak et al. 1997
* measurements/fit as the important factor when optimising the model/chamber measurement fits   

      DO iw = 1, nw - 1
         DO i = 1, nz
            qy = exp(-0.055*(wc(iw)-308.)) / 
     $           (5.5 + 9.2e-19*airden(i))
            qy = min(qy, 1.)
            sq(j,i,iw) = yg(iw) * qy
         ENDDO
      ENDDO
      
      END

*=============================================================================*

      SUBROUTINE r104(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH2C(CH3)CHO   =*
*=  methacrolein        photolysis:                                          =*
*=           CH2C(CH3)CHO + hv -> Products                                   =*
*=                                                                           =*
*=  Cross section from                                                       =*
*= (1)  R. Meller, priv. comm, MPI Mainz 1990 as reported by                 =*
*= Roeth, E.-P., R. Ruhnke, G. Moortgat, R. Meller, and W. Schneider,        =*
*= UV/VIS-Absorption Cross Sections and QUantum Yields for Use in            =*
*= Photochemistry and Atmospheric Modeling, Part 2: Organic Substances,      =*
*= Forschungszentrum Julich, Report Jul-3341, 1997.                          =*
*= (2) T. Gierczak, J.B. Burkholder, R.K. Talukdar, A. Mellouki, S.B. Barone,=* 
*= and A.R. Ravishankara: Atmospheric fate of methyl vinyl ketone and        =*
*= methacrolein, J. Photochem. Photobiol. A: Chem. 110, 1-10 (1997)          =*
*= Values recommended by JPL/IUPAC 2005                                      =*
*=                                                                           =*
*=  Quantum Yield from                                                       =*
*= (1) T. Gierczak, J.B. Burkholder, R.K. Talukdar, A. Mellouki, S.B. Barone,=* 
*= and A.R. Ravishankara: Atmospheric fate of methyl vinyl ketone and        =*
*= methacrolein, J. Photochem. Photobiol. A: Chem. 110, 1-10 (1997)          =*
*= (2) Average from the effective qy from the EU RADICAL consortium (2000),  =*
*= the optimised qy's from Carter (2000) and Pinho et al. (2005).            =*
*= See Pinho et al., Atmos Env, 39(7), 1303-1322 (2005) for details          =*
*=                                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=20000)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw), yg1(kw)
      REAL qy
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

************************* CH2C(CH3)CHO photolysis

      j = j+1
      jlabel(j) = 'CH2C(CH3)CHO -> Products'

* Cross section:
* mabs = 1 Meller et al. (1990) and Raber and Moortgat (1996)
* mabs = 2 Gierczak et al (1997)

      mabs = 2
      myld = 2

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/methacrolein.abs',
     $        STATUS='old')
         DO i = 1, 10
            READ(kin,*)
         ENDDO
         n = 15213
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs. EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/ABS/macr_gierczak97.abs',
     $        STATUS='old')
         DO i = 1, 15
            READ(kin,*)
         ENDDO
         n = 146
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF
      ENDIF

* quantum yield from
*
* (1)  Gierczak, T., J. B. Burkholder, R. K. Talukdar, A. Mellouki, S. B. Barone,
*      and A. R. Ravishankara, Atmospheric fate of methyl vinyl ketone and methacrolein,
*      J. Photochem. Photobiol A: Chemistry, 110 1-10, 1997.
*      Upper limit, quantum yield < 0.01

      IF( mabs. EQ. 1) THEN

         qy = 0.01

* (2)  Average from:
*      RADICAL (2000): Effective quantum yield of < 0.004
*      Carter (2000):  Optimised SAPRC-99 quantum yield of 0.0041 using chamber data
*      Pinho et al. (2005):  Optimised MCMv3 quantum yield of 0.0036 using chamber data
* 
*      NB: The fit of the mcmv3 model to the SAPRC chamber data is very sensitive to 
*      changes in MACR qy in the above work of Pinho et al. 2005.
         
      ELSEIF ( myld. EQ. 2) THEN
      
         qy = 0.0039   
         
      ENDIF
      
* combine
      
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw) * qy
         ENDDO
      ENDDO     

      END

*=============================================================================*

      SUBROUTINE r105(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3COCO(OH)    =*
*=  pyruvic acid        photolysis:                                          =*
*=           CH3COCO(OH) + hv -> Products                                    =*
*=                                                                           =*
*=  Cross section from                                                       =*
*= (1)  Horowitz, A., R. Meller, and G. K. Moortgat, The UV-VIS absorption   =*
*=      cross section of the a-dicarbonyl compounds: pyruvic acid, biacetyl, =*
*=      and glyoxal. J. Photochem. Photobiol. A:Chemistry, v.146,            =*
*=      pp.19-27, 2001.                                                      =*
*=                                                                           =*
*= (2)  MPI-Mainz Spectral Atlas: JPL-2005 recommendation                    =*
*=                                                                           =*
*=  Quantum yield                                                            =*
*=  (1)  Assumed Unity                                                       =*
*=                                                                           =*
*=  (2)  Average quantum yield from the RADICAL project (Moortgat (2001)     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=20000)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

************************* CH3COCO(OH) photolysis

      j = j+1
      jlabel(j) = 'CH3COCO(OH) -> Products'
      
* Cross section:
* mabs = 1 Horowitz et al. (2001)
* mabs = 2 JPL 2005

      mabs = 2
      myld = 2

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCOOH/pyruvic_horowitz.abs',
     $        STATUS='old')
         DO i = 1, 8
            READ(kin,*)
         ENDDO
         n = 148
         DO i = 1, n
            READ(kin,*) x(i), y(i)
            y(i) = y(i) * 1.e-20
         ENDDO
         CLOSE(kin)

         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF
      
      ELSEIF (mabs. EQ. 2) THEN
      
         OPEN(UNIT=kin,FILE='DATAJ1/CH3COCOOH/pyruvic_JPL05.abs',
     $        STATUS='old')
         DO i = 1, 10
            READ(kin,*)
         ENDDO
         n = 140
         DO i = 1, n
            READ(kin,*) x(i), y(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
         CALL addpnt(x,y,kdata,n,               0.,0.)
         CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
         CALL addpnt(x,y,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x,y,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF
      ENDIF 
      
* quantum yield  = 1
      
      IF( myld. EQ. 1) THEN

        qy = 1.

        DO iw = 1, nw - 1
           DO i = 1, nz
              sq(j,i,iw) = yg(iw) * qy
           ENDDO
        ENDDO
* Average quantum yield taken from the RADICAL project (Moortgat (2001)        
        
      ELSEIF (myld. EQ. 2) THEN
      
      qy = 0.43

        DO iw = 1, nw - 1
           DO i = 1, nz
              sq(j,i,iw) = yg(iw) * qy
           ENDDO
        ENDDO
        
      ENDIF    

      END

*=============================================================================*

      SUBROUTINE r106(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3CH2ONO2     =*
*=  ethyl nitrate       photolysis:                                          =*
*=           CH3CH2ONO2 + hv -> CH3CH2O + NO2                                =*
*=                                                                           =*
*= (1)  Absorption cross sections of several organic from                    =*
*= Talukdar, R. K., J. B. Burkholder, M. Hunter, M. K. Gilles,               =*
*= J. M Roberts, and A. R. Ravishankara, Atmospheric fate of several         =*
*= alkyl nitrates, J. Chem. Soc., Faraday Trans., 93(16) 2797-2805, 1997.    =*
*=                                                                           =*
*= (2)  IUPAC 2005                                                           =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy, sig
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

************************* CH3CH2ONO2 photolysis

      j = j+1
      jlabel(j) = 'CH3CH2ONO2 -> CH3CH2O + NO2'

      mabs = 2

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/RONO2_talukdar.abs',
     $        STATUS='old')
         DO i = 1, 10
            READ(kin,*)
         ENDDO
         n1 = 0
         n2 = 0
         DO i = 1, 63
            READ(kin,*) x1(i), dum, dum, y1(i), y2(i), dum, dum
            if (y1(i) .gt. 0.) n1 = n1 + 1
            if (y2(i) .gt. 0.) n2 = n2 + 1
            x2(i) = x1(i)
            y1(i) = y1(i) * 1.e-20
            y2(i) = y2(i) * 1.e-3
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,y2(n2))
         CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF( mabs. EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/RONO2_iupac05.abs',
     $        STATUS='old')
         DO i = 1, 10
            READ(kin,*) 
         ENDDO
         n = 32
         DO i = 1, n
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)
         n1 = n
         n2 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg3,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

*         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
*         CALL addpnt(x2,y2,kdata,n2,               0.,0.)
*         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
*         CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
*         CALL inter2(nw,wl,yg4,n2,x2,y2,ierr)
         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,y2(n2))
         CALL inter2(nw,wl,yg4,n2,x2,y2,ierr)
         
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield  = 1

      qy = 1.

      DO iw = 1, nw - 1
         DO i = 1, nz

            IF( mabs. EQ. 1) THEN

               sig = yg1(iw)*exp(yg2(iw)*(tlev(i)-298.))
        
            ELSEIF( mabs. EQ. 2) THEN
            
               IF (wc(iw) .GE. 185. .AND. wc(iw) .LE. 235.) THEN

                  sig = yg3(iw)

               ELSEIF (wc(iw) .GT. 235. .AND. wc(iw) .LE. 340.) THEN

                  sig = yg3(iw)*exp(yg4(iw)*(tlev(i)-298.))

               ELSE
                  
                  sig = 0.
               
               ENDIF
               
            ENDIF

            sq(j,i,iw) = sig * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r107(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3CHONO2CH3   =*
*=  isopropyl nitrate   photolysis:                                          =*
*=           CH3CHONO2CH3 + hv -> CH3CHOCH3 + NO2                            =*
*=                                                                           =*
*= (1)  Absorption cross sections of several organic from                    =*
*= Talukdar, R. K., J. B. Burkholder, M. Hunter, M. K. Gilles,               =*
*= J. M Roberts, and A. R. Ravishankara, Atmospheric fate of several         =*
*= alkyl nitrates, J. Chem. Soc., Faraday Trans., 93(16) 2797-2805, 1997.    =*
*=                                                                           =*
*= (2)  Absorption cross sections from IUPAC 2005                            =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL qy, sig
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

************************* CH3CHONO2CH3 photolysis

      j = j+1
      jlabel(j) = 'CH3CHONO2CH3 -> CH3CHOCH3 + NO2'

      mabs = 2

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/RONO2_talukdar.abs',
     $        STATUS='old')
         DO i = 1, 10
            READ(kin,*)
         ENDDO
         n1 = 0
         n2 = 0
         DO i = 1, 63
            READ(kin,*) x1(i), dum, dum, dum, dum, y1(i), y2(i)
            if (y1(i) .gt. 0.) n1 = n1 + 1
            if (y2(i) .gt. 0.) n2 = n2 + 1
            x2(i) = x1(i)
            y1(i) = y1(i) * 1.e-20
            y2(i) = y2(i) * 1.e-3
         ENDDO
         CLOSE(kin)

*         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
*         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
*         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
*         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
*         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),y1(n1))
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,y1(n1))
         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,y2(n2))
         CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF( mabs. EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/iC3H7ONO2_iupac05.abs',
     $        STATUS='old')
         DO i = 1, 8
            READ(kin,*) 
         ENDDO
         n = 37 
         DO i = 1, n 
            READ(kin,*) x1(i), y1(i), y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)
         n1 = n
*         n2 = n

*         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
*         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
*         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
*         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
*         CALL inter2(nw,wl,yg3,n1,x1,y1,ierr)
         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),y1(n1))
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,y1(n1))
         CALL inter2(nw,wl,yg3,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/iC3H7ONO2_B_iupac05.abs',
     $        STATUS='old')
         DO i = 1, 9
            READ(kin,*) 
         ENDDO
         n = 21 
         DO i = 1, n 
            READ(kin,*) x2(i), y2(i)
         ENDDO
         CLOSE(kin)
         n2 = n

*         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
*         CALL addpnt(x2,y2,kdata,n2,               0.,0.)
*         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
*         CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
*         CALL inter2(nw,wl,yg4,n2,x2,y2,ierr)
*         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),y2(1))
         CALL addpnt(x2,y2,kdata,n2,               0.,y2(1))
*         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),y2(n2))
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,y2(n2))
         CALL inter2(nw,wl,yg4,n2,x2,y2,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield  = 1

      qy = 1.

      DO iw = 1, nw - 1
         DO i = 1, nz

            IF( mabs. EQ. 1) THEN

               IF(wc(iw) .GE. 236. .AND. wc(iw) .LE. 344.)THEN ! arr

               sig = yg1(iw)*exp(yg2(iw)*(tlev(i)-298.))
               
               ELSEIF(wc(iw) .GT. 344. .AND. wc(iw) .LE. 360.)THEN ! arr
               
               sig = yg1(iw) ! arr
               
               ELSE ! arr
               
               sig = 0. ! arr
               
               ENDIF ! arr

            ELSEIF( mabs. EQ. 2) THEN

               IF(wc(iw) .GE. 185. .AND. wc(iw) .LE. 240.)THEN

                  sig = yg3(iw)

               ELSEIF(wc(iw) .GT. 240. .AND. wc(iw) .LE. 340.)THEN

                  sig = yg3(iw)*exp(yg4(iw)*(tlev(i)-298.))

               ELSEIF(wc(iw) .GT. 340. .AND. wc(iw) .LE. 360.)THEN

                  sig = yg3(iw)

               ELSE
                  
                  sig = 0.

               ENDIF
            ENDIF
            
            yg5(iw) = sig ! arr
            sq(j,i,iw) = sig * qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r108(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=   nitroxy ethanol CH2(OH)CH2(ONO2) + hv -> CH2(OH)CH2(O.) + NO2           =*
*=                                                                           =*
*=  Cross section from Roberts, J. R. and R. W. Fajer, UV absorption cross   =*
*=    sections of organic nitrates of potential atmospheric importance and   =*
*=    estimation of atmospheric lifetimes, Env. Sci. Tech., 23, 945-951,     =*
*=    1989.
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* local

      REAL qy, sig
      INTEGER iw, i
      REAL a, b, c

************************* CH2(OH)CH2(ONO2) photolysis

      j = j+1
      jlabel(j) = 'CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2'


* coefficients from Roberts and Fajer 1989, over 270-306 nm

      a = -2.359E-3
      b = 1.2478
      c = -210.4

* quantum yield  = 1

      qy = 1.

      DO iw = 1, nw - 1
         IF (wc(iw) .GE. 270. .AND. wc(iw) .LE. 306.) THEN
            sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
         ELSE
            sig = 0.
         ENDIF
         DO i = 1, nz
            sq(j,i,iw) = sig * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r109(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=   nitroxy acetone CH3COCH2(ONO2) + hv -> CH3COCH2(O.) + NO2               =*
*=                                                                           =*
*=  (1) Cross section from Roberts, J. R. and R. W. Fajer, UV absorption     =*
*=      cross sections of organic nitrates of potential atmospheric          =*
*=      importance and importance and estimation of atmospheric lifetimes,   =*
*=      Env. Sci. Tech., 23, 945-951, 1989.                                  =*
*=                                                                           =*
*= (2)  I. Barnes, K.H. Becker, and T. Zhu, "Near UV absorption spectra      =*
*=      and photolysis products of difunctional organic nitrates: Possible   =*
*=      importance as NOx reservoirs," J. Atmos. Chem. 17, 353-373 (1993).   =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j
      
* data arrays      
      
      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL qy, sig
      INTEGER iw, ierr
      REAL a, b, c
      INTEGER mabs, myld

************************* CH3COCH2(ONO2) photolysis

      j = j+1
      jlabel(j) = 'CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2'

      mabs = 1
      
      IF( mabs. EQ. 1) THEN 

* coefficients from Roberts and Fajer 1989, over 284-335 nm

      a = -1.365E-3
      b = 0.7834
      c = -156.8

* quantum yield  = 1

      qy = 1.

         DO iw = 1, nw - 1
            IF (wc(iw) .GE. 284. .AND. wc(iw) .LE. 335.) THEN
               sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            ELSE
               sig = 0.
            ENDIF
            DO i = 1, nz
               sq(j,i,iw) = sig * qy
            ENDDO
         ENDDO
      
      ELSEIF( mabs. EQ. 2) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/NOA_barnes93.abs',
     $        STATUS='old')
         DO i = 1, 10
            READ(kin,*) 
         ENDDO
         n = 20 
         DO i = 1, n 
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)
         n1 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),y1(1))
         CALL addpnt(x1,y1,kdata,n1,               0.,y1(1))
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),y1(n1))
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,y1(n1))
         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)  
*         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
*         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
*         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
*         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
*         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF
      
* quantum yield  = 1

         qy = 1.

         DO iw = 1, nw - 1
            
            IF(wc(iw) .GE. 245. .AND. wc(iw) .LE. 340.)THEN
            sig = yg1(iw)
            ELSE
            sig = 0
            ENDIF
               
            DO i = 1, nz      
            
            sq(j,i,iw) = sig * qy

            ENDDO
         ENDDO

      ENDIF

      END

*=============================================================================*

      SUBROUTINE r110(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  t-butyl nitrate C(CH3)3(ONO2) + hv -> C(CH3)(O.) + NO2                   =*
*=                                                                           =*
*=  Cross section from Roberts, J. R. and R. W. Fajer, UV absorption cross   =*
*=    sections of organic nitrates of potential atmospheric importance and   =*
*=    estimation of atmospheric lifetimes, Env. Sci. Tech., 23, 945-951,     =*
*=    1989.
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* local

      REAL qy, sig
      INTEGER iw, i
      REAL a, b, c

************************* C(CH3)3(ONO2) photolysis

      j = j+1
      jlabel(j) = 'C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2'


* coefficients from Roberts and Fajer 1989, over 270-330 nm

      a = -0.993E-3
      b = 0.5307
      c = -115.5

* quantum yield  = 1

      qy = 1.

      DO iw = 1, nw - 1
         IF (wc(iw) .GE. 270. .AND. wc(iw) .LE. 330.) THEN
            sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
         ELSE
            sig = 0.
         ENDIF
         DO i = 1, nz
            sq(j,i,iw) = sig * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r111(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for ClOOCl         =*
*=  ClO dimer           photolysis:                                          =*
*=           ClOOCl + hv -> Cl + ClOO                                        =*
*=                                                                           =*
*=  Cross section from  JPL2002                                              =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=20000)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER ierr
      INTEGER iw

************************* ClOOCl photolysis
* from JPL-2002

      j = j+1
      jlabel(j) = 'ClOOCl -> Cl + ClOO'

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/CLOOCL_jpl02.abs',
     $     STATUS='old')
      DO i = 1, 25
         READ(kin,*)
      ENDDO
      n = 131
      DO i = 1, n
         READ(kin,*) x(i), y(i)
         y(i) = y(i) * 1.e-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,               0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x,y,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield  = 1

      qy = 1.

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw) * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r112(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for hydroxyacetone =*
*=  CH2(OH)COCH3        photolysis:                                          =*
*=           CH2(OH)COCH3  -> CH3CO + CH2OH
*=                         -> CH2(OH)CO + CH3                                =*
*=                                                                           =*
*=  Cross section from Orlando et al. (1999)                                 =*
*=                                                                           =*
*=  Quantum yield assumed 0.325 for each channel (J. Orlando, priv.comm.2003)=*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=20000)

      INTEGER i, n
      REAL x(kdata), y(kdata)

* local

      REAL yg(kw)
      REAL qy
      INTEGER ierr
      INTEGER iw

************************* CH2(OH)COCH3 photolysis
* from Orlando et al. 1999

      j = j+1
      jlabel(j) = 'CH2(OH)COCH3 -> CH3CO + CH2(OH)'
      j = j+1
      jlabel(j) = 'CH2(OH)COCH3 -> CH2(OH)CO + CH3'

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/Hydroxyacetone.abs',
     $     STATUS='old')
      DO i = 1, 8
         READ(kin,*)
      ENDDO
      n = 101
      DO i = 1, n
         READ(kin,*) x(i), y(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,               0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x,y,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* Total quantum yield  = 0.65, equal for each of the two channels

      qy = 0.325

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j-1,i,iw) = yg(iw) * qy
            sq(j,i,iw) = yg(iw) * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r113(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for HOBr           =*
*=  HOBr -> OH + Br                                                          =*
*=  Cross section from JPL 2003                                              =*
*=  Quantum yield assumed unity as in JPL2003                                =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, i

* data arrays

* local

      REAL qy, sig
      INTEGER iw

************************* HOBr photolysis
* from JPL2003

      j = j+1
      jlabel(j) = 'HOBr -> OH + Br'

      qy = 1.
      DO iw = 1, nw - 1
         sig = 24.77 * EXP( -109.80*(LOG(284.01/wc(iw)))**2 ) + 
     $         12.22 * exp(  -93.63*(LOG(350.57/wc(iw)))**2 ) + 
     $         2.283 * exp(- 242.40*(LOG(457.38/wc(iw)))**2 )
         sig = sig * 1.e-20
         IF(wc(iw) .LT. 250. .OR. wc(iw) .GT. 550.) sig = 0.

         DO i = 1, nz
            sq(j,i,iw) = sig * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r114(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for BrO            =*
*=  BrO -> Br + O                                                            =*
*=  Cross section from JPL 2003                                              =*
*=  Quantum yield assumed unity as in JPL2003                                =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* output weighting functions

      INTEGER j
      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* data arrays

      INTEGER n, i
      REAL x(20), y(20)

* local

      INTEGER iw
      REAL qy, yg(kw), dum

************************* HOBr photolysis
* from JPL2003

      j = j+1
      jlabel(j) = 'BrO -> Br + O'

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/BrO.jpl03',
     $     STATUS='old')
      DO i = 1, 14
         READ(kin,*)
      ENDDO
      n = 15
      DO i = 1, n
         READ(kin,*) x(i), dum, y(i)
         y(i) = y(i) * 1.e-20
      ENDDO
      n = n + 1
      x(n) = dum
      CLOSE(kin)

* use bin-to-bin interpolation

      CALL inter4(nw,wl,yg,n,x,y,1)

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw) * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r115(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for BrO            =*
*=  Br2 -> Br + Br                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* output weighting functions

      INTEGER j
      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* data arrays

      INTEGER kdata
      PARAMETER(kdata=50)
      INTEGER n, i
      REAL x(kdata), y(kdata)

* local

      INTEGER iw, ierr
      REAL qy, yg(kw)

************************* Br2 photolysis

      j = j + 1
      jlabel(j) = 'Br2 -> Br + Br'

* Absorption cross section from:
* Seery, D.J. and D. Britton, The continuous absorption spectra of chlorine, 
* bromine, bromine chloride, iodine chloride, and iodine bromide, J. Phys. 
* Chem. 68, p. 2263 (1964).

      OPEN(UNIT=kin,FILE='DATAJ1/ABS/Br2.abs',
     $     STATUS='old')

      DO i = 1, 6
         READ(kin,*) 
      ENDDO
      n = 29
      DO i = 1, n
         READ(kin,*) x(i),  y(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,               0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x,y,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw) * qy
         ENDDO
      ENDDO

      END


*=============================================================================*

      SUBROUTINE r116(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  ARR ADDITION							     =*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for nC3H7CHO       =*
*=  photolysis:                                                              =*
*=         nC3H7CHO + hv -> nC3H7 + HCO (Norrish I)                          =*
*=                       -> C2H4 + CH3CHO (Norrish II)                       =*
*=  Cross section:  IUPAC 05 recommendation (Martinez et al. (1992)          =*
*=  Quantum yield:  IUPAC 05 recommendation (Tadic et al. (2001)             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** nC3H7CHO photodissociation ********************

      j = j+1
      jlabel(j) = 'nC3H7CHO -> nC3H7 + HCO'

      j = j+1
      jlabel(j) = 'nC3H7CHO -> C2H4 + CH3CHO'
                         
* working grid arrays:
*     yg1 = cross section at a specific temperature (298K)
*     yg2 = quantum yield data for radical channel (NI)
*     yg3 = quantum yield data for molecular channel (NII)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  IUPAC-05 data, from Martinez et al. (1992)

* Quantum yield:
* 1:  IUPAC-05 data, from Tadic et al. (2001)

      mabs = 1
      myld = 1

	IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/C3H7CHO/C3H7CHO_Martinez.abs',STATUS='old')
            DO i = 1, 10
               READ(kin,*)
            ENDDO
            n = 105
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	       y1(i) = y1(i)*1.e-21
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
	ENDIF

* quantum yield

      DO iw = 1, nw - 1
      
        IF (myld .EQ. 1) THEN

             IF (wc(iw) .GE. 290. .AND. wc(iw) .LE. 380.) THEN
                yg2(iw) = 0.19
		yg3(iw) = 0.09
             ELSE
                yg2(iw) = 0.
		yg3(iw) = 0.
             
	     ENDIF
        ENDIF
      ENDDO
      
* combine
* yg1 = xsect
* yg2 = qy for radical channel
* yg3 = qy for molecular channel


      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            sig = yg1(iw)
            sig = MAX(sig, 0.)
           
* quantum yields:

            qy1 = yg2(iw)
            qy2 = yg3(iw)

	      qy1 = MAX(0., qy1)
            qy1 = MIN(1., qy1)
            
            qy2 = MAX(0., qy2)
            qy2 = MIN(1., qy2)

*            sq(j-1,iz,iw) = yg1(iw) * yg2(iw)
*            sq(j-1,iz,iw) = yg1(iw) * yg3(iw)
            
            sq(j-1,iz,iw) = sig * qy1
            sq(j  ,iz,iw) = sig * qy2
            
          ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r117(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for iC3H7CHO       =*
*=  photolysis:                                                              =*
*=         iC3H7CHO + hv -> iC3H7 + HCO                                      =*
*=                                                                           =*
*=  Cross section:  Choice between                                           =*
*=                 (1) MPI-Mainz Spectral Atlas, from Martinez et al. (1992) =*
*=                 (2) MPI-Mainz Spectral Atlas, Chen et al. (2002)          =*
*=  Quantum yield:  Chen et al. (2002)                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz

* data arrays

      INTEGER kdata
      PARAMETER(kdata=150)

      INTEGER i, n
      INTEGER n1
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw), yg1(kw)
      REAL qy1
      REAL sig
      INTEGER ierr
      INTEGER iw

      INTEGER mabs, myld

************************* iC3H7CHO photolysis
* 1:  iC3H7 + HCO

      j = j+1
      jlabel(j) = 'iC3H7CHO -> iC3H7 + HCO'

* options
* mabs for cross sections
* myld for quantum yields

* Absorption:
* 1:  Martinez et al. (1992): mabs = 1
* 2:  Chen et al. (2002): mabs = 2

* Quantum yield:
* 1:  Chen et al. (2002): myld = 1

      mabs = 1
      myld = 1

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE=
     >    'DATAJ1/C3H7CHO/i_C3H7CHO_Martinez.abs', STATUS='old')
         do i = 1, 10
            read(kin,*)
         enddo
         n = 106
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-21
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ELSEIF(mabs .EQ. 2) THEN

         OPEN(UNIT=kin,FILE=
     >    'DATAJ1/C3H7CHO/i_C3H7CHO_Chen.abs', STATUS='old')
      
         DO i = 1, 10
            read(kin,*)
         ENDDO
        
         n = 11

         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
	    y1(i) = y1(i) * 1e-21
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yields

      IF (myld .EQ. 1) THEN

         OPEN(UNIT=kin,FILE=
     >     'DATAJ1/C3H7CHO/i_C3H7CHO_Chen.yld', STATUS='old')

         DO i = 1, 8
            read(kin,*)
         ENDDO

         n = 11
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO

         CLOSE(kin)

         n1 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF
       ENDIF   

* combine:

      DO iw = 1, nw - 1
         DO iz = 1, nz

            sig = yg(iw)

* quantum yields:

            qy1 = yg1(iw)
            qy1 = MIN(qy1,1.)

            sq(j,iz,iw) = sig * qy1

         ENDDO
      ENDDO

      END


*=============================================================================*

      SUBROUTINE r118(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for nC4H9CHO       =*
*=  photolysis:                                                              =*
*=         nC4H9CHO + hv -> nC4H9 + HCO (Norrish I)                          =*
*=                       -> CH3CHO + CH2=CHCH3 (Norrish II)                  =*
*=  Cross section:  Choice between                                           =*
*=                 (1) MPI-Mainz Spectral Atlas, from Tadic et al. (2001)    =*
*=                 						             =*
*=  Quantum yield: (1) Tadic et al. (2001)                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** nC4H9CHO photodissociation ********************

      j = j+1
      jlabel(j) = 'nC4H9CHO -> nC4H9 + HCO'

      j = j+1
      jlabel(j) = 'nC4H9CHO -> C3H6 + CH3CHO'
                         
* working grid arrays:
*     yg1 = cross section at a specific temperature (298K)
*     yg2 = quantum yield data for radical channel (NI)
*     yg3 = quantum yield data for molecular channel (NII)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Tadic et al. (2001)

* Quantum yield:
* 1:  Tadic et al. (2001) 

      mabs = 1
      myld = 1

	IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/C4H9CHO/C4H9CHO_Tadic.abs',STATUS='old')
            DO i = 1, 10
               READ(kin,*)
            ENDDO
            n = 121
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	       y1(i) = y1(i)*1.e-21
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
	ENDIF

* quantum yield
*
*      DO iw = 1, nw - 1
*      
*        IF (myld .EQ. 1) THEN
*
*             IF (wc(iw) .GE. 230. .AND. wc(iw) .LE. 350.) THEN
*               yg2(iw) = 0.068
*		yg3(iw) = 0.272
*             ELSE
*                yg2(iw) = 0.
*		yg3(iw) = 0.
*             
*	     ENDIF
*        ENDIF
*      ENDDO
      
* combine
* yg1 = xsect
* yg2 = qy for radical channel
* yg3 = qy for molecular channel


      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            sig = yg1(iw)
            sig = MAX(sig, 0.)
           
* quantum yields:
* use Stern-Volmer pressure dependence:

        IF(myld .EQ. 1) THEN 
 
            qy1 = 0.13*(1./(2.44 + (1.931E-3)*airden(i)/2.45e19))
            qy2 = 0.52*(1./(2.44 + (1.931E-3)*airden(i)/2.45e19))

            qy1 = MIN(qy1,1.)
	    qy2 = MIN(qy2,1.)
        
        ENDIF
 
            sq(j-1,iz,iw) = sig * qy1
            sq(j  ,iz,iw) = sig * qy2
            
          ENDDO
      ENDDO

      END


*=============================================================================*

      SUBROUTINE r119(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  ARR ADDITION							     =*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for iC4H9CHO       =*
*=  photolysis:                                                              =*
*=         iC4H9CHO + hv -> iC4H9 + HCO (Norrish I)                          =*
*=                       -> C3H6 + CH3CHO (Norrish II)                       =*
*=  Cross section:  MPI-Mainz Spectral Atlas, from Zhu et al. (1999)  	     =*
*=  Quantum yield:  Zhu et al. (1999) - roughly scaled to QYeff of           =*
*=  Moortgat (2001)             		                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** iC4H9CHO photodissociation ********************

      j = j+1
      jlabel(j) = 'iC4H9CHO -> iC4H9 + HCO'

      j = j+1
      jlabel(j) = 'iC4H9CHO -> C3H6 + CH3CHO'
                         
* working grid arrays:
*     yg1 = cross section at a specific temperature (298K)
*     yg2 = quantum yield data for radical channel (NI)
*     yg3 = quantum yield data for molecular channel (NII)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Zhu et al. (1999)

* Quantum yield:
* 1:  Taken from Tadic et al. (1999)
* Estimated ch3cho qys roughly scaled to Moortgat et al (2001) effective 
* qy of 0.27 (+_0.01) 

      mabs = 1
      myld = 1

	IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/C4H9CHO/i_C4H9CHO_Zhu.abs',STATUS='old')
            DO i = 1, 10
               READ(kin,*)
            ENDDO
            n = 11
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	       y1(i) = y1(i)*1.e-21
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
	ENDIF

* quantum yield

      DO iw = 1, nw - 1
      
        IF (myld .EQ. 1) THEN

             IF (wc(iw) .GE. 280. .AND. wc(iw) .LE. 330.) THEN
                yg2(iw) = 0.06
		yg3(iw) = 0.16
             ELSE
                yg2(iw) = 0.
		yg3(iw) = 0.
             
	     ENDIF
        ENDIF
      ENDDO
      
* combine
* yg1 = xsect
* yg2 = qy for radical channel
* yg3 = qy for molecular channel


      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            sig = yg1(iw)
            sig = MAX(sig, 0.)
           
* quantum yields:

            qy1 = yg2(iw)
            qy2 = yg3(iw)

	      qy1 = MAX(0., qy1)
            qy1 = MIN(1., qy1)
            
            qy2 = MAX(0., qy2)
            qy2 = MIN(1., qy2)
            
            sq(j-1,iz,iw) = sig * qy1
            sq(j  ,iz,iw) = sig * qy2
            
          ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r120(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for tC4H9CHO       =*
*=  photolysis:                                                              =*
*=         tC4H9CHO + hv -> tC4H9 + HCO                                      =*
*=                                                                           =*
*=  Cross section:   MPI-Mainz Spectral Atlas, from Zhu et al. (1999)	     =*
*=                           						     =*
*=  Quantum yield:   Zhu et al. (1999)                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz

* data arrays

      INTEGER kdata
      PARAMETER(kdata=150)

      INTEGER i, n
      INTEGER n1
      REAL x1(kdata)
      REAL y1(kdata)

* local

      REAL yg(kw), yg1(kw)
      REAL qy1
      REAL sig
      INTEGER ierr
      INTEGER iw

      INTEGER mabs, myld

************************* tC4H9CHO photolysis
* 1:  tC4H9 + HCO

      j = j+1
      jlabel(j) = 'tC4H9CHO -> tC4H9 + HCO'

* options
* mabs for cross sections
* myld for quantum yields

* Absorption:
* 1:  Zhu et al. (1999): mabs = 1

* Quantum yield:
* 1:  Zhu et al. (1999): myld = 1

      mabs = 1
      myld = 1

      IF (mabs .EQ. 1) THEN

         OPEN(UNIT=kin,FILE=
     >    'DATAJ1/C4H9CHO/t_C4H9CHO_Zhu.abs', STATUS='old')
         do i = 1, 10
            read(kin,*)
         enddo
         n = 11
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.e-21
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yields

      IF (myld .EQ. 1) THEN

         OPEN(UNIT=kin,FILE=
     >     'DATAJ1/C4H9CHO/t_C4H9CHO_Zhu.yld', STATUS='old')

         DO i = 1, 10
            read(kin,*)
         ENDDO

         n = 11
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO

         CLOSE(kin)

         n1 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF
       ENDIF   

* combine:

      DO iw = 1, nw - 1
         DO iz = 1, nz

            sig = yg(iw)

* quantum yields:

            qy1 = yg1(iw)
            qy1 = MIN(qy1,1.)

            sq(j,iz,iw) = sig * qy1

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r121(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for nC4H9CHO       =*
*=  photolysis:                                                              =*
*=         nC5H11CHO + hv -> nC5H11 + HCO (Norrish I)                        =*
*=                        -> CH3CHO + CH2=CHCH2CH3 (Norrish II)              =*
*=  Cross section:  Choice between                                           =*
*=                 (1) MPI-Mainz Spectral Atlas, from Plagens et al. (1998)  =*
*=                 (2) MPI-Mainz Spectral Atlas, from Tang and Zhu (2004)    =*
*=									     =*
*=  Quantum yield: (1) Tadic et al. (2001)                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** nC5H11CHO photodissociation ********************

      j = j+1
      jlabel(j) = 'nC5H11CHO -> nC5H11 + HCO'

      j = j+1
      jlabel(j) = 'nC5H11CHO -> C4H8 + CH3CHO'
                         
* working grid arrays:
*     yg1 = cross section at a specific temperature (298K)
*     yg2 = quantum yield data for radical channel (NI)
*     yg3 = quantum yield data for molecular channel (NII)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Plagens et al. (1998)
* 2:  MPI-Mainz Spectral Atlas, from Tang and Zhu (2004)

* Quantum yield:
* 1:  Tadic et al. (2001) 

      mabs = 1
      myld = 1

	IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/C5H11CHO/C5H11CHO_Plagens.abs',STATUS='old')
            DO i = 1, 14
               READ(kin,*)
            ENDDO

            n = 117

            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	       y1(i) = y1(i)*1.e-21
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
	
        ELSEIF(mabs .EQ. 2) THEN

           OPEN(UNIT=kin,FILE=
     >      'DATAJ1/C5H11CHO/C5H11CHO_Tang.abs', STATUS='old')
      
           DO i = 1, 12
            read(kin,*)
           ENDDO
        
           n = 11

           DO i = 1, n
            READ(kin,*) x1(i), y1(i)
	    y1(i) = y1(i) * 1e-21
           ENDDO
           CLOSE (kin)

           CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
           CALL addpnt(x1,y1,kdata,n,               0.,0.)
           CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
           CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
           CALL inter2(nw,wl,yg,n,x1,y1,ierr)
           IF (ierr .NE. 0) THEN
              WRITE(*,*) ierr, jlabel(j)
              STOP
           ENDIF

      ENDIF


* quantum yield
    
* combine
* yg1 = xsect
* qy1 = qy for radical channel
* qy2 = qy for molecular channel


      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            sig = yg1(iw)
            sig = MAX(sig, 0.)
           
* quantum yields:
* use Stern-Volmer pressure dependence:

        IF(myld .EQ. 1) THEN 
 
            qy1 = 0.23*(1./(2.26 + (4.758E-4)*airden(i)/2.45e19))
            qy2 = 0.59*(1./(2.26 + (4.758E-4)*airden(i)/2.45e19))

            qy1 = MIN(qy1,1.)
	    qy2 = MIN(qy2,1.)
        
        ENDIF
 
            sq(j-1,iz,iw) = sig * qy1
            sq(j  ,iz,iw) = sig * qy2
            
          ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r122(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for nC4H9CHO       =*
*=  photolysis:                                                              =*
*=         nC6H13CHO + hv -> nC6H13 + HCO (Norrish I)                        =*
*=                        -> CH3CHO + CH2=CHCH2CH2CH3 (Norrish II)           =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     from Tang and Zhu et al. (2004)                       =*
*=                 						             =*
*=  Quantum yield: (1) Tadic et al. (2002)                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** nC6H13CHO photodissociation *******************

      j = j+1
      jlabel(j) = 'nC6H13CHO -> nC6H13 + HCO'

      j = j+1
      jlabel(j) = 'nC6H13CHO -> C5H10 + CH3CHO'
                         
* working grid arrays:
*     yg1 = cross section at a specific temperature (298K)
*     qy1 = quantum yield data for radical channel (NI)
*     qy2 = quantum yield data for molecular channel (NII)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Tang and Zhu (2004)

* Quantum yield:
* 1:  Tadic et al. (2002) 

      mabs = 1
      myld = 1

	IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/C6H13CHO/C6H13CHO_Tang.abs',STATUS='old')
            DO i = 1, 12
               READ(kin,*)
            ENDDO
            n = 11
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	       y1(i) = y1(i)*1.e-21
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
	ENDIF
      
* combine
* yg1 = xsect
* qy1 = qy for radical channel
* qy2 = qy for molecular channel


      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            sig = yg1(iw)
            sig = MAX(sig, 0.)
           
* quantum yields:
* use Stern-Volmer pressure dependence:

        IF(myld .EQ. 1) THEN 
 
            qy1 = 0.101*(1./(2.408 + (1.169E-3)*airden(i)/2.45e19))
            qy2 = 0.379*(1./(2.408 + (1.169E-3)*airden(i)/2.45e19))

            qy1 = MIN(qy1,1.)
	    qy2 = MIN(qy2,1.)
        
        ENDIF
 
            sq(j-1,iz,iw) = sig * qy1
            sq(j  ,iz,iw) = sig * qy2
            
          ENDDO
      ENDDO

      END


*=============================================================================*

      SUBROUTINE r123(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3COC2H5      =*
*=  photolysis:                                                              =*
*=         CH3COC2H5 + hv -> C2H5 + CH3CO (from MCMv3.1)                     =*
*=                        -> Products?                                       =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     from Martinez et al. 1992 (IUPAC 2005)                =*
*=                 (2) MPI-Mainz Spectral Atlas,                             =*
*=                     from YujingMellouki(2000)                             =*
*=                 						                                     =*
*=  Quantum yield: (1) Raber and Moortgat 1996 (IUPAC 2005)                  =*
*=                 (2) Pinho et al. 2005 and Carter (2000)		     =*
*=                 						                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** nC6H13CHO photodissociation *******************

      j = j+1
      jlabel(j) = 'CH3COC2H5 -> C2H5 + CH3CO'
                         
* working grid arrays:
*     yg = cross section at a specific temperature (298K)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Martinez et al. 1992
* 2:  MPI-Mainz Spectral Atlas, from YujingMellouki(2000)

* Quantum yield:
* 1:  IUPAC 2005, Raber and Moortget 1996
* 2.  Pinho et al., Atmos Env., 39(7), 1303 (2005):
*     From the evaluation and optimisation of the 
*     MCMv3 mechanism using SAPRC chamber data, and
*     Carter (2000): From the evaluation and optimisation 
*     of the SAPRC-99 mechanism using chamber data.
*     An average of the two optimised qy's is used

      mabs = 1
      myld = 2

	  IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/CH3COC2H5/CH3COC2H5_Martinez.abs',STATUS='old')
            DO i = 1, 11
               READ(kin,*)
            ENDDO
            n = 84
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	       y1(i) = y1(i)*1.e-20
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF

         ELSEIF(mabs .EQ. 2) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/CH3COC2H5/CH3COC2H5_Yujing.abs',STATUS='old')
            DO i = 1, 11
               READ(kin,*)
            ENDDO
            n = 104
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
	   ENDIF

*combine

      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            sig = yg(iw)
            sig = MAX(sig, 0.)
           
* quantum yields:
* 1:  IUPAC 2005, Raber and Moortgat 1996 (0.34)
* 2:  Average of Pinho et al. 2005 (0.17) and 
*     Carter 2000 (0.15)

        IF(myld .EQ. 1) THEN
         
           IF(wc(iw) .GE. 275. .AND. wc(iw) .LE. 380.) THEN
              
              qy1 = 0.34
           ELSE
              qy1 = 0.
           ENDIF   
        
        ELSEIF(myld .EQ. 2) THEN
              
              qy1 = 0.16
             
        ENDIF
        
            sq(j  ,iz,iw) = sig * qy1
            
          ENDDO
      ENDDO

      END


*=============================================================================*

      SUBROUTINE r124(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for pinonaldehyde  =*
*=  photolysis:                                                              =*
*=         pinonaldehyde + hv -> products                                    =*
*=                                                                           =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     from Hallquist et al., Environ. Sci. Technol. 31,     =*
*=                     3166-3172 (1997). (IUPAC 2005)                        =*
*=                 						                                     =*
*=  Quantum yield: (1) RADICAL project 2002                                  =*
*=                 (2) Jaoui and Kamens Atmos. Env., 37, 1835 (2003)         =*
*=                 						                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** pinonaldehyde photodissociation ***************

      j = j+1
      jlabel(j) = 'pinonaldehyde -> products'
                         
* working grid arrays:
*     yg = cross section at a specific temperature (298K)         

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Hallquist et al. (1997) 

* Quantum yield:
* 1.  RADICAL project 2002 (IUPAC 2005)
* 2.  Jaoui and Kamens Atmos. Env., 37, 1835 (2003) 

      mabs = 1
      myld = 1

	  IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/pinonaldehyde/pinal_hallquist.abs',STATUS='old')
            DO i = 1, 11
               READ(kin,*)
            ENDDO
            n = 14
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)  
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
        ENDIF
         
*combine

      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            sig = yg(iw)
            sig = MAX(sig, 0.)
           
* quantum yields:
* 1:  RADICAL project 2002 (IUPAC 2005)
* 2:  Jaoui and Kamens Atmos. Env., 37, 1835 (2003)

        IF(myld .EQ. 1) THEN
              
              qy1 = 0.14 
        
        ELSEIF(myld .EQ. 2) THEN
              
              qy1 = 0.4
             
        ENDIF
        
            sq(j  ,iz,iw) = sig * qy1
            
          ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r125(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for acrolein       =*
*=  photolysis:                                                              =*
*=         CH2CHCHO + hv -> CH2CHCO + H                                      =*
*=                       -> CH2CH + HCO                                      =*
*=                                                                           =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     I. Magneron, R. Thévenet, A. Mellouki, G. Le Bras,    =*
*=                     G.K. Moortgat, and K. Wirtz,                          =*
*=                     "A study of the photolysis and OH-initiated           =*
*=                     oxidation of acrolein and trans-crotonaldehyde,"      =*
*=                     J. Phys. Chem. A 106, 2526-2537 (2002) (JPL (2005)    =*
*=                                                                           =*
*=                 (2) RADICAL REPORT:  Moortgat, G.K. (Ed.), RADICAL:       =*
*=                     Evaluation of radical sources in atmospheric          =*
*=                     chemistry through chamber and laboratory studies.     =*
*=                     1999                                                  =*
*=                                                                           =*
*=  Quantum yield: (1) Effective QY: RADICAL project 2002                    =*
*=                                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** acrolein photodissociation ********************

      j = j+1
      jlabel(j) = 'CH2CHCHO -> CH2CHCO + H'

      j = j+1
      jlabel(j) = 'CH2CHCHO -> CH2CH + HCO'
                         
* working grid arrays:
*     yg = cross section at a specific temperature (298K)         

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Magneron et al. (2002)
*
* 2:  RADICAL REPORT:  Moortgat, G.K. (Ed.), RADICAL:  Evaluation 
*     of radical sources in atmospheric chemistry through chamber 
*     and laboratory studies. 2002

* Quantum yield:
* 1.  RADICAL project 2002

      mabs = 1
      myld = 1

	  IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/acrolein/magneron.abs',STATUS='old')
            DO i = 1, 13
               READ(kin,*)
            ENDDO
            n = 54
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)  
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
        ELSEIF(mabs .EQ. 2) THEN
        
            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/acrolein/radical.abs',STATUS='old')
            DO i = 1, 14
               READ(kin,*)
            ENDDO
            n = 35
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)  
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
        ENDIF
         
*combine

      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            sig = yg(iw)
            sig = MAX(sig, 0.)
           
* quantum yields:
* 1:  RADICAL project 2002 

        IF(myld .EQ. 1) THEN
              
              qy1 = 0.005/2
              qy2 = 0.005/2 
             
        ENDIF
        
            sq(j-1,iz,iw) = sig * qy1
            sq(j  ,iz,iw) = sig * qy2
            
          ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r126(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CH3CH2CH2ONO2  =*
*=  1-propyl nitrate photolysis:                                             =*
*=           nC3H7ONO2 + hv -> nC3H7O + NO2                                  =*
*=                                                                           =*
*= Absorption cross sections from IUPAC 2005                                 =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy, sig
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

*****************CH3CH2CH2ONO2 photolysis*******************

      j = j+1
      jlabel(j) = 'nC3H7ONO2 -> nC3H7O + NO2'

      mabs = 1

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/nC3H7ONO2_iupac05.abs',
     $        STATUS='old')
         DO i = 1, 8
            READ(kin,*) 
         ENDDO
         n = 32
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)
         n1 = n

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg3,n1,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield  = 1

      
      DO iw = 1, nw - 1
         DO i = 1, nz
         
         qy = 1.
         sig = yg3(iw)

         sq(j,i,iw) = sig * qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r127(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for nC4H9ONO2      =*
*=  1-butyl nitrate photolysis:                                              =*
*=           nC4H9ONO2 + hv -> nC4H9O + NO2                                  =*
*=                                                                           =*
*=  Absorption cross sections: average from the following datasets (at 298K):=*
*=  (from MPI-Mainz Spectral Atlas)                                          =*
*=                                                                           =*
*=  J.M. Roberts and R.W. Fajer, "UV absorption cross sections of organic    =*
*=  nitrates of potential atmospheric importance and estimation of           =*
*=  atmospheric lifetimes," Environ. Sci. Technol. 23, 945-951 (1989).       =*
*=                                                                           =*
*=  K.C. Clemitshaw, J. Williams, O.V. Rattigan, D.E. Shallcross, K.S. Law,  =*
*=  and R.A. Cox, "Gas-phase ultraviolet absorption cross-sections and       =*
*=  atmospheric lifetimes of several C2-C5 alkyl nitrates,"                  =*
*=  J. Photochem. Photobiol. A: Chem. 102, 117-126 (1996).                   =*
*=                                                                           =*
*=  M.P. Turberg, D.M. Giolando, C. Tilt, T. Soper, S. Mason, M. Davies,     =* 
*=  P. Klingensmith, and G.A. Takacs, "Atmospheric chemistry of alkyl        =* 
*=  nitrates,"  J. Photochem. Photobiol. A. Chem. 51, 281-292 (1990).        =* 
*=                                                                           =*
*=  L. Zhu and D. Kellis, "Temperature dependence of the UV absorption       =*
*=  cross sections and photodissociation products of C3 - C5 alkyl           =*
*=  nitrates," Chem. Phys. Lett. 278, 41-48 (1997)                           =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy, sig
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

*****************C4H9ONO2 photolysis*****************************

      j = j+1
      jlabel(j) = 'nC4H9ONO2 -> nC4H9O + NO2'

      mabs = 1

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/nC4H9ONO2_average.abs',
     $        STATUS='old')
         DO i = 1, 13
            READ(kin,*) 
         ENDDO
         n = 33
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield  = 1

      
      DO iw = 1, nw - 1
         DO i = 1, nz
         
         qy = 1.
         sig = yg(iw)

         sq(j,i,iw) = sig * qy

         ENDDO
      ENDDO

      END


*=============================================================================*

      SUBROUTINE r128(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for iC4H9ONO2      =*
*=  2-methyl-1-propyl nitrate (i-butyl nitrate) photolysis:                  =*
*=           iC4H9ONO2 + hv -> iC4H9O + NO2                                  =*
*=                                                                           =*
*=  Absorption cross sections: (from MPI-Mainz Spectral Atlas)               =*
*=                                                                           =*
*=  K.C. Clemitshaw, J. Williams, O.V. Rattigan, D.E. Shallcross, K.S. Law,  =*
*=  and R.A. Cox, "Gas-phase ultraviolet absorption cross-sections and       =*
*=  atmospheric lifetimes of several C2-C5 alkyl nitrates,"                  =*
*=  J. Photochem. Photobiol. A: Chem. 102, 117-126 (1996).                   =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy, sig
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

*****************iC4H9ONO2 photolysis***************************

      j = j+1
      jlabel(j) = 'iC4H9ONO2 -> iC4H9O + NO2'

      mabs = 1

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/iC4H9ONO2_clem.abs',
     $        STATUS='old')
         DO i = 1, 11
            READ(kin,*) 
         ENDDO
         n = 25
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield  = 1

      
      DO iw = 1, nw - 1
         DO i = 1, nz
         
         qy = 1.
         sig = yg(iw)

         sq(j,i,iw) = sig * qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r129(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for 2-C4H9ONO2     =*
*=  2-butyl nitrate photolysis:                                              =*
*=           2-C4H9ONO2 + hv -> 2-C4H9O + NO2                                =*
*=                                                                           =*
*=  Absorption cross sections: (from MPI-Mainz Spectral Atlas)               =*
*=                                                                           =*
*=  from IUPAC 2005.  Taken from:                                            =*
*=  J.M. Roberts and R.W. Fajer, "UV absorption cross sections of organic    =*
*=  nitrates of potential atmospheric importance and estimation of           =*
*=  atmospheric lifetimes," Environ. Sci. Technol. 23, 945-951 (1989).       =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy, sig
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

*****************2-C4H9ONO2 photolysis**************************

      j = j+1
      jlabel(j) = '2-C4H9ONO2 -> 2-C4H9O + NO2'

      mabs = 1

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/2_C4H9ONO2_iupac05.abs',
     $        STATUS='old')
         DO i = 1, 8
            READ(kin,*) 
         ENDDO
         n = 15
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield  = 1

      
      DO iw = 1, nw - 1
         DO i = 1, nz
         
         qy = 1.
         sig = yg1(iw)

         sq(j,i,iw) = sig * qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r130(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for nC5H11ONO2     =*
*=  1-pentyl nitrate photolysis:                                             =*
*=           nC5H11ONO2 + hv -> nC5H11O + NO2                                =*
*=                                                                           =*
*=  Absorption cross sections: average values from the following:            =*
*=  (from MPI-Mainz Spectral Atlas)                                          =*
*=                                                                           =*
*=  K.C. Clemitshaw, J. Williams, O.V. Rattigan, D.E. Shallcross, K.S. Law,  =*
*=  and R.A. Cox, "Gas-phase ultraviolet absorption cross-sections and       =*
*=  atmospheric lifetimes of several C2-C5 alkyl nitrates,"                  =*
*=  J. Photochem. Photobiol. A: Chem. 102, 117-126 (1996).                   =*
*=                                                                           =*
*=  L. Zhu and D. Kellis, "Temperature dependence of the UV absorption cross =* 
*=  sections and photodissociation products of C3 - C5 alkyl nitrates,"      =*
*=  Chem. Phys. Lett. 278, 41-48 (1997).                                     =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy, sig
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

*****************nC5H11ONO2 photolysis*************************

      j = j+1
      jlabel(j) = 'nC5H11ONO2 -> nC5H11O + NO2'

      mabs = 1

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/nC5H11ONO2_average.abs',
     $        STATUS='old')
         DO i = 1, 11
            READ(kin,*) 
         ENDDO
         n = 16
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield  = 1

      
      DO iw = 1, nw - 1
         DO i = 1, nz
         
         qy = 1.
         sig = yg1(iw)

         sq(j,i,iw) = sig * qy

         ENDDO
      ENDDO

      END


*=============================================================================*

      SUBROUTINE r131(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for 2-C5H11ONO2    =*
*=  2-pentyl nitrate photolysis:                                             =*
*=           2-C5H11ONO2 + hv -> 2-C5H11O + NO2                              =*
*=                                                                           =*
*=  Absorption cross sections:                                               =*
*=  (from MPI-Mainz Spectral Atlas)                                          =*
*=                                                                           =*
*=  J.M. Roberts and R.W. Fajer, "UV absorption cross sections of organic    =*
*=  nitrates of potential atmospheric importance and estimation of           =*
*=  atmospheric lifetimes," Environ. Sci. Technol. 23, 945-951 (1989).       =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy, sig
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

*****************2-C5H11ONO2 photolysis************************

      j = j+1
      jlabel(j) = '2-C5H11ONO2 -> 2-C5H11O + NO2'

      mabs = 1

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/2_C5H11ONO2_roberts.abs',
     $        STATUS='old')
         DO i = 1, 8
            READ(kin,*) 
         ENDDO
         n = 23
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield  = 1

      
      DO iw = 1, nw - 1
         DO i = 1, nz
         
         qy = 1.
         sig = yg1(iw)

         sq(j,i,iw) = sig * qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r132(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for 3-C5H11ONO2    =*
*=  3-pentyl nitrate photolysis:                                             =*
*=           3-C5H11ONO2 + hv -> 3-C5H11O + NO2                              =*
*=                                                                           =*
*=  Absorption cross sections:                                               =*
*=  (from MPI-Mainz Spectral Atlas)                                          =*
*=  									     =*
*=  Calculated using formula in:                                             =*
*=                                                                           =*
*=  J.M. Roberts and R.W. Fajer, "UV absorption cross sections of organic    =*
*=  nitrates of potential atmospheric importance and estimation of           =*
*=  atmospheric lifetimes," Environ. Sci. Technol. 23, 945-951 (1989).       =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy, sig
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld
      REAL a, b, c

*****************3-C5H11ONO2 photolysis************************

      j = j+1
      jlabel(j) = '3-C5H11ONO2 -> 3-C5H11O + NO2'


* coefficients from Roberts and Fajer 1989, over 270-335 nm

      a = -1.446E-3
      b = 0.7712
      c = -147.4

* quantum yield  = 1

      qy = 1.

      DO iw = 1, nw - 1
         IF (wc(iw) .GE. 270. .AND. wc(iw) .LE. 335.) THEN
            sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
         ELSE
            sig = 0.
         ENDIF

         DO i = 1, nz
            sq(j,i,iw) = sig * qy
         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r133(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for i-C5H11ONO2    =*
*=  3-methyl-1-butyl nitrate, isopentyl nitrate photolysis:                  =*
*=           i-C5H11ONO2 + hv -> i-C5H11O + NO2                              =*
*=                                                                           =*
*=  Absorption cross sections:                                               =*
*=  (from MPI-Mainz Spectral Atlas)                                          =*
*=                                                                           =*
*=  M.P. Turberg, D.M. Giolando, C. Tilt, T. Soper, S. Mason, M. Davies,     =*
*=  P. Klingensmith, and G.A. Takacs, "Atmospheric chemistry of alkyl        =*
*=  nitrates," J. Photochem. Photobiol. A. Chem. 51, 281-292 (1990)          =*                                                                      =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL dum
      REAL yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy, sig
      INTEGER ierr
      INTEGER iw
      INTEGER mabs, myld

*****************i-C5H11ONO2 photolysis************************

      j = j+1
      jlabel(j) = 'i-C5H11ONO2 -> i-C5H11O + NO2'

      mabs = 1

      IF( mabs. EQ. 1) THEN

         OPEN(UNIT=kin,FILE='DATAJ1/RONO2/i_C5H11ONO2_turberg.abs',
     $        STATUS='old')
         DO i = 1, 11
            READ(kin,*) 
         ENDDO
         n = 25
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,               0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      ENDIF

* quantum yield  = 1

      
      DO iw = 1, nw - 1
         DO i = 1, nz
         
         qy = 1.
         sig = yg1(iw)

         sq(j,i,iw) = sig * qy

         ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r134(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  1-nitroxy-2-butanone:                                                    =*
*=  CH3CH2C(O)CH2(ONO2) + hv -> CH3CH2C(O)CH2(O.) + NO2                      =*
*=                                                                           =*
*=  Cross section from:                                                      =*
*=                                                                           =*
*=  (1)  I. Barnes, K.H. Becker, and T. Zhu, "Near UV absorption spectra     =*
*=       and photolysis products of difunctional organic nitrates: Possible  =*
*=       importance as NOx reservoirs," J. Atmos. Chem. 17, 353-373 (1993).  =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j
      
* data arrays      
      
      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL qy, sig
      INTEGER iw, ierr
      REAL a, b, c
      INTEGER mabs, myld

************************* CH3CH2C(O)CH2(ONO2) photolysis

      j = j+1
      jlabel(j) = 'CH3CH2C(O)CH2(ONO2) -> CH3CH2C(O)CH2(O.) + NO2'

      mabs = 1
      
      IF( mabs. EQ. 1) THEN 

* coefficients from Barnes et al. 1993, over 270-340 nm

      a = -1.011E-3
      b = 0.5437
      c = -116.9

* quantum yield  = 1

      qy = 1.

         DO iw = 1, nw - 1
            IF (wc(iw) .GE. 270. .AND. wc(iw) .LE. 340.) THEN
               sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            ELSE
               sig = 0.
            ENDIF
            DO i = 1, nz
               sq(j,i,iw) = sig * qy
            ENDDO
         ENDDO

      ENDIF

      END

*=============================================================================*

      SUBROUTINE r135(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  3-nitroxy-2-butanone:                                                    =*
*=  CH3CH(ONO2)C(O)CH3 + hv -> CH3CH(O.)C(O)CH3 + NO2                        =*
*=                                                                           =*
*=  Cross section from:                                                      =*
*=                                                                           =*
*=  (1)  I. Barnes, K.H. Becker, and T. Zhu, "Near UV absorption spectra     =*
*=       and photolysis products of difunctional organic nitrates: Possible  =*
*=       importance as NOx reservoirs," J. Atmos. Chem. 17, 353-373 (1993).  =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j
      
* data arrays      
      
      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL qy, sig
      INTEGER iw, ierr
      REAL a, b, c
      INTEGER mabs, myld

************************* CH3CH2C(O)CH2(ONO2) photolysis

      j = j+1
      jlabel(j) = 'CH3CH(ONO2)C(O)CH3 -> CH3CH(O.)C(O)CH3 + NO2'

      mabs = 1
      
      IF( mabs. EQ. 1) THEN 

* coefficients from Barnes et al. 1993, over 270-340 nm

      a = -1.044E-3
      b = 0.5780
      c = -123.5

* quantum yield  = 1

      qy = 1.

         DO iw = 1, nw - 1
            IF (wc(iw) .GE. 270. .AND. wc(iw) .LE. 340.) THEN
               sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            ELSE
               sig = 0.
            ENDIF
            DO i = 1, nz
               sq(j,i,iw) = sig * qy
            ENDDO
         ENDDO

      ENDIF

      END

*=============================================================================*

      SUBROUTINE r136(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  1,2-propandiol dinitrate:                                                =*
*=  CH3CH(ONO2)CH2(ONO2) + hv -> CH3CH(O.)CH2(ONO2) + NO2                    =*
*=  CH3CH(ONO2)CH2(ONO2) + hv -> CH3CH(ONO2)CH2(O.) + NO2                    =*
*=                                                                           =*
*=  Cross section from:                                                      =*
*=                                                                           =*
*=  (1)  I. Barnes, K.H. Becker, and T. Zhu, "Near UV absorption spectra     =*
*=       and photolysis products of difunctional organic nitrates: Possible  =*
*=       importance as NOx reservoirs," J. Atmos. Chem. 17, 353-373 (1993).  =*
*=                                                                           =*
*=  Quantum yield assumed unity (0.5 each process)                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j
      
* data arrays      
      
      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL qy1, qy2, sig
      INTEGER iw, ierr
      REAL a, b, c
      INTEGER mabs, myld

****************CH3CH(ONO2)CH2(ONO2) photolysis ***************

      j = j+1
      jlabel(j) = 'CH3CH(ONO2)CH2(ONO2) -> CH3CH(O.)CH2(ONO2) + NO2'
      
      j = j+1
      jlabel(j) = 'CH3CH(ONO2)CH2(ONO2) -> CH3CH(ONO2)CH2(O.) + NO2'

      mabs = 1
      
      IF( mabs. EQ. 1) THEN 

* coefficients from Barnes et al. 1993, over 250-320 nm

      a = -5.99E-4
      b = 0.2915
      c = -79.24

* quantum yield  = 1

      qy1 = 0.5
      qy2 = 0.5
      
         DO iw = 1, nw - 1
            IF (wc(iw) .GE. 250. .AND. wc(iw) .LE. 320.) THEN
               sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            ELSE
               sig = 0.
            ENDIF
            DO i = 1, nz
               sq(j-1,i,iw) = sig * qy1
               sq(j  ,i,iw) = sig * qy2
            ENDDO
         ENDDO

      ENDIF

      END

*=============================================================================*

      SUBROUTINE r137(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  1,2-butandiol dinitrate:                                                 =*
*=  CH3CH2CH(ONO2)CH2(ONO2) + hv -> CH3CH2CH(O.)CH2(ONO2) + NO2              =*
*=  CH3CH2CH(ONO2)CH2(ONO2) + hv -> CH3CH2CH(ONO2)CH2(O.) + NO2              =*
*=                                                                           =*
*=  Cross section from:                                                      =*
*=                                                                           =*
*=  (1)  I. Barnes, K.H. Becker, and T. Zhu, "Near UV absorption spectra     =*
*=       and photolysis products of difunctional organic nitrates: Possible  =*
*=       importance as NOx reservoirs," J. Atmos. Chem. 17, 353-373 (1993).  =*
*=                                                                           =*
*=  Quantum yield assumed unity (0.5 each process)                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j
      
* data arrays      
      
      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL qy1, qy2, sig
      INTEGER iw, ierr
      REAL a, b, c
      INTEGER mabs, myld

****************CH3CH2CH(ONO2)CH2(ONO2) photolysis ******************

      j = j+1
      jlabel(j) = 'C2H5CH(ONO2)CH2(ONO2) -> C2H5CH(O.)CH2(ONO2) + NO2'
      j = j+1
      jlabel(j) = 'C2H5CH(ONO2)CH2(ONO2) -> C2H5CH(ONO2)CH2(O.) + NO2'
      
      mabs = 1
      
      IF( mabs. EQ. 1) THEN 

* coefficients from Barnes et al. 1993, over 250-320 nm
* NB: use parameterization as for 1,2-propandiol dinitrate

      a = -5.99E-4
      b = 0.2915
      c = -79.24

* quantum yield  = 1

      qy1 = 0.5
      qy2 = 0.5
      
         DO iw = 1, nw - 1
            IF (wc(iw) .GE. 250. .AND. wc(iw) .LE. 320.) THEN
               sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            ELSE
               sig = 0.
            ENDIF
            DO i = 1, nz
               sq(j-1,i,iw) = sig * qy1
               sq(j  ,i,iw) = sig * qy2
            ENDDO
         ENDDO

      ENDIF

      END

*=============================================================================*

      SUBROUTINE r138(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  2,3-butandiol dinitrate:                                                 =*
*=  CH3CH(ONO2)CH(ONO2)CH3 + hv -> CH3CH(O.)CH(ONO2)CH3 + NO2                =*
*=                                                                           =*
*=  Cross section from:                                                      =*
*=                                                                           =*
*=  (1)  I. Barnes, K.H. Becker, and T. Zhu, "Near UV absorption spectra     =*
*=       and photolysis products of difunctional organic nitrates: Possible  =*
*=       importance as NOx reservoirs," J. Atmos. Chem. 17, 353-373 (1993).  =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j
      
* data arrays      
      
      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL qy1, qy2, sig
      INTEGER iw, ierr
      REAL a, b, c
      INTEGER mabs, myld

*********************CH3CH(ONO2)CH(ONO2)CH3 photolysis****************

      j = j+1
      jlabel(j) = 'CH3CH(ONO2)CH(ONO2)CH3 -> CH3CH(O.)CH(ONO2)CH3 + NO2'
      
      mabs = 1
      
      IF( mabs. EQ. 1) THEN 

* coefficients from Barnes et al. 1993, over 250-320 nm
* NB: using parametrization as for 1,2-butandiol dinitrate

      a = -6.217E-4
      b = 0.3025
      c = -80.41

* quantum yield  = 1

      qy1 = 1.
      
         DO iw = 1, nw - 1
            IF (wc(iw) .GE. 250. .AND. wc(iw) .LE. 320.) THEN
               sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            ELSE
               sig = 0.
            ENDIF
            DO i = 1, nz
               sq(j  ,i,iw) = sig * qy1
            ENDDO
         ENDDO

      ENDIF

      END

*=============================================================================*

      SUBROUTINE r139(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  3,4-dinitroxy-1-butene:                                                  =*
*=  CH2=CHCH(ONO2)CH2(ONO2) + hv -> CH2=CHCH(O.)CH2(ONO2) + NO2              =*
*=  CH2=CHCH(ONO2)CH2(ONO2) + hv -> CH2=CHCH(ONO2)CH2(O.) + NO2              =*
*=                                                                           =*
*=  Cross section from:                                                      =*
*=                                                                           =*
*=  (1)  I. Barnes, K.H. Becker, and T. Zhu, "Near UV absorption spectra     =*
*=       and photolysis products of difunctional organic nitrates: Possible  =*
*=       importance as NOx reservoirs," J. Atmos. Chem. 17, 353-373 (1993).  =*
*=                                                                           =*
*=  Quantum yield assumed unity (0.5 for each process)                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j
      
* data arrays      
      
      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL qy1, qy2, sig
      INTEGER iw, ierr
      REAL a, b, c
      INTEGER mabs, myld

*********************CH2=CHCH(ONO2)CH2(ONO2) photolysis****************

      j = j+1
      jlabel(j) = 'CH2=CHCH(ONO2)CH2(ONO2) -> C2H3CH(O.)CH2(ONO2) + NO2'
      j = j+1
      jlabel(j) = 'CH2=CHCH(ONO2)CH2(ONO2) -> C2H3CH(ONO2)CH2(O.) + NO2'
      
      mabs = 1
      
      IF( mabs. EQ. 1) THEN 

* coefficients from Barnes et al. 1993, over 250-335 nm
* NB: using parametrization as for 2,3-butandiol dinitrate

      a = -5.740E-4
      b = 0.2771
      c = -77.47

* quantum yield  = 1

      qy1 = 0.5
      qy2 = 0.5
      
         DO iw = 1, nw - 1
            IF (wc(iw) .GE. 250. .AND. wc(iw) .LE. 335.) THEN
               sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            ELSE
               sig = 0.
            ENDIF
            DO i = 1, nz
               sq(j-1,i,iw) = sig * qy1
               sq(j  ,i,iw) = sig * qy2
            ENDDO
         ENDDO

      ENDIF

      END      
*=============================================================================*

      SUBROUTINE r140(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  1,4-dinitroxy-2-butene:                                                  =*
*=  CH2(ONO2)CH=CHCH2(ONO2) + hv -> CH2(ONO2)CH=CHCH2(O.) + NO2              =*
*=                                                                           =*
*=  Cross section from:                                                      =*
*=                                                                           =*
*=  (1)  I. Barnes, K.H. Becker, and T. Zhu, "Near UV absorption spectra     =*
*=       and photolysis products of difunctional organic nitrates: Possible  =*
*=       importance as NOx reservoirs," J. Atmos. Chem. 17, 353-373 (1993).  =*
*=                                                                           =*
*=  Quantum yield assumed unity                                              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j
      
* data arrays      
      
      INTEGER kdata
      PARAMETER(kdata=200)

      INTEGER i, n, n1, n2
      REAL x1(kdata), y1(kdata)
      REAL x2(kdata), y2(kdata)

* local

      REAL yg1(kw), yg2(kw)
      REAL qy1, qy2, sig
      INTEGER iw, ierr
      REAL a, b, c
      INTEGER mabs, myld

*********************CH2(ONO2)CH=CHCH2(ONO2) photolysis****************

      j = j+1
      jlabel(j) = 'CH2(ONO2)C2H2CH2(ONO2) -> CH2(ONO2)C2H2CH2(O.) + NO2'
      
      mabs = 1
      
      IF( mabs. EQ. 1) THEN 

* coefficients from Barnes et al. 1993, over 250-325 nm

      a = -5.432E-4
      b = 0.2631
      c = -75.92

* quantum yield  = 1

      qy1 = 1.
      
         DO iw = 1, nw - 1
            IF (wc(iw) .GE. 250. .AND. wc(iw) .LE. 325.) THEN
               sig = EXP(a*wc(iw)*wc(iw) + b*wc(iw) + c)
            ELSE
               sig = 0.
            ENDIF
            DO i = 1, nz
               sq(j  ,i,iw) = sig * qy1
            ENDDO
         ENDDO

      ENDIF

      END      


*=============================================================================*

      SUBROUTINE r141(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for crotonaldehyde =*
*=  photolysis:                                                              =*
*=         CH3CH=CHCHO + hv -> CH3CH=CHCO + H                                =*
*=                          -> CH3CH=CH + HCO                                =*
*=                          -> CH3CH=CH2 + CO                                =*
*=                                                                           =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     I. Magneron, R. Thévenet, A. Mellouki, G. Le Bras,    =*
*=                     G.K. Moortgat, and K. Wirtz,                          =*
*=                     "A study of the photolysis and OH-initiated           =*
*=                     oxidation of acrolein and trans-crotonaldehyde,"      =*
*=                     J. Phys. Chem. A 106, 2526-2537 (2002)                =*
*=                                                                           =*
*=  NB:  RADICAL project measurements of crotonaldehysde cross-sections      =*
*=       are very similar to those of Magneron 02                            =*
*=                                                                           =*
*=  Quantum yield: (1) Effective QY: RADICAL project 2002                    =*
*=                                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** crotonaldehyde photodissociation **************

      j = j+1
      jlabel(j) = 'CH3CH=CHCHO -> CH3CH=CHCO + H'

      j = j+1
      jlabel(j) = 'CH3CH=CHCHO -> CH3CH=CH + HCO'
      
      j = j+1
      jlabel(j) = 'CH3CH=CHCHO -> CH3CH=CH2 + CO'
      
* working grid arrays:
*     yg = cross section at a specific temperature (298K)         

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Magneron et al. (2002)

* Quantum yield:
* 1.  RADICAL project 2002

      mabs = 1
      myld = 1

	  IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/crotonaldehyde/magneron.abs',STATUS='old')
            DO i = 1, 13
               READ(kin,*)
            ENDDO
            n = 68
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)  
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
        ENDIF
         
*combine

      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            sig = yg(iw)
            sig = MAX(sig, 0.)
           
* quantum yields:
* 1:  RADICAL project 2002 

        IF(myld .EQ. 1) THEN
              
              qy1 = 0.03*0.333
              qy2 = 0.03*0.333
              qy3 = 0.03*0.333 
             
        ENDIF
        
            sq(j-2,iz,iw) = sig * qy1
            sq(j-1,iz,iw) = sig * qy2
            sq(j  ,iz,iw) = sig * qy3
            
          ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r142(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for 2-pentanone    =*
*=  photolysis:                                                              =*
*=         CH3COC3H7 + hv -> C3H7 + CH3CO (from MCMv3.1)                     =*
*=                        -> Products (NORRISH)?                             =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     from Martinez et al. 1992                             =*
*=                     NB: Horowitz (1999) cs are very similar to            =*
*=                     Martinez et al. (1992)                                =*
*=                 						                                     =*
*=  Quantum yield: (1) Raber and Moortgat 1996 (2-butanone)                  =*
*=                     c.f. Laval-Szopa thesis (2003)                        =*
*=                 (2) Pinho et al. 2005 and Carter (2000) (2-butanone)		 =*
*=                 						                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** 2-pentanone photodissociation *****************

      j = j+1
      jlabel(j) = 'CH3COC3H7 -> C3H7 + CH3CO'
                         
* working grid arrays:
*     yg = cross section at a specific temperature (298K)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Martinez et al. 1992

* Quantum yield:
* 1:  IUPAC 2005, Raber and Moortget 1996 (for 2-butanone)
* 2.  Pinho et al., Atmos Env., 39(7), 1303 (2005):
*     From the evaluation and optimisation of the 
*     MCMv3 mechanism using SAPRC chamber data, and
*     Carter (2000): From the evaluation and optimisation 
*     of the SAPRC-99 mechanism using chamber data.
*     An average of the two optimised qy's is used.  
*     Again this is originally for 2-butanone

      mabs = 1
      myld = 1

	  IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/CH3COC3H7/CH3COC3H7_Martinez.abs',STATUS='old')
            DO i = 1, 11
               READ(kin,*)
            ENDDO
            n = 84
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	         y1(i) = y1(i)*1.e-20
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
	   ENDIF

*combine

      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            IF(wc(iw) .GE. 218. .AND. wc(iw) .LE. 354.) THEN
          
              sig = yg(iw)
              sig = MAX(sig, 0.)
            
            ELSE 
              
              sig = 0.  
              
            ENDIF  
           
* quantum yields:
* 1:  IUPAC 2005, Raber and Moortgat 1996 (0.34)
* 2:  Average of Pinho et al. 2005 (0.17) and 
*     Carter 2000 (0.15)

        IF(myld .EQ. 1) THEN
              
              qy1 = 0.34   
        
        ELSEIF(myld .EQ. 2) THEN
              
              qy1 = 0.16      
             
        ENDIF
        
            sq(j  ,iz,iw) = sig * qy1
            
          ENDDO
      ENDDO

      END

*=============================================================================*

      SUBROUTINE r143(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for 3-pentanone    =*
*=  photolysis:                                                              =*
*=         C2H5COC2H5 + hv -> C2H5 + C2H5CO (from MCMv3.1)                   =*
*=                                                                           =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     from Martinez et al. 1992                             =*
*=                     NB: Horowitz (1999) cs are very similar to            =*
*=                     Martinez et al. (1992)                                =*
*=                 						                                     =*
*=  Quantum yield: (1) Raber and Moortgat 1996 (2-butanone)                  =*
*=                     c.f. Laval-Szopa thesis (2003)                        =*
*=                 (2) Pinho et al. 2005 and Carter (2000) (2-butanone)		 =*
*=                 						                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
**************** 3-pentanone photodissociation *****************

      j = j+1
      jlabel(j) = 'C2H5COC2H5 -> C2H5 + C2H5CO'
                         
* working grid arrays:
*     yg = cross section at a specific temperature (298K)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Martinez et al. 1992

* Quantum yield:
* 1:  IUPAC 2005, Raber and Moortget 1996 (for 2-butanone)
* 2.  Pinho et al., Atmos Env., 39(7), 1303 (2005):
*     From the evaluation and optimisation of the 
*     MCMv3 mechanism using SAPRC chamber data, and
*     Carter (2000): From the evaluation and optimisation 
*     of the SAPRC-99 mechanism using chamber data.
*     An average of the two optimised qy's is used.  
*     Again this is originally for 2-butanone

      mabs = 1
      myld = 1

	  IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	     'DATAJ1/C2H5COC2H5/C2H5COC2H5_Martinez.abs',STATUS='old')
            DO i = 1, 11
               READ(kin,*)
            ENDDO
            n = 81
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	         y1(i) = y1(i)*1.e-20
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
	   ENDIF

*combine

      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            IF(wc(iw) .GE. 210. .AND. wc(iw) .LE. 342.) THEN
          
              sig = yg(iw)
              sig = MAX(sig, 0.)
              
            ELSE
              sig = 0.
            ENDIF   
           
* quantum yields:
* 1:  IUPAC 2005, Raber and Moortgat 1996 (0.34)
* 2:  Average of Pinho et al. 2005 (0.17) and 
*     Carter 2000 (0.15)

        IF(myld .EQ. 1) THEN
         
           IF(wc(iw) .GE. 210. .AND. wc(iw) .LE. 342.) THEN
              
              qy1 = 0.34
           ELSE
              qy1 = 0.
           ENDIF   
        
        ELSEIF(myld .EQ. 2) THEN
        
           IF(wc(iw) .GE. 210. .AND. wc(iw) .LE. 342.) THEN
              
              qy1 = 0.16
           ELSE
              qy1 = 0.
           ENDIF      
             
        ENDIF
        
            sq(j  ,iz,iw) = sig * qy1
            
          ENDDO
      ENDDO
       
      END      

*=============================================================================*

      SUBROUTINE r144 (nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  2,4-dimethyl-3-pentanone photolysis:                                     =*
*=         ((CH3)2CH)2C(O) + hv -> (CH3)2C + (CH3)2CHCO                      =*
*=                                                                           =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     from Yujing and Mellouki (2000)                       =*
*=                 						                                     =*
*=  Quantum yield: (1) Raber and Moortgat 1996 (2-butanone)                  =*
*=                     c.f. Laval-Szopa thesis (2003)                        =*
*=                 (2) Pinho et al. 2005 and Carter (2000) (2-butanone)		 =*
*=                 						                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
********** 2,4-dimethyl-3-pentanone photodissociation **********

      j = j+1
      jlabel(j) = '((CH3)2CH)2C(O) -> (CH3)2C + (CH3)2CHCO'
                         
* working grid arrays:
*     yg = cross section at a specific temperature (298K)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Yujing and Mellouki (2000)

* Quantum yield:
* 1:  Raber and Moortgat 1996 (for 2-butanone)
* 2.  Pinho et al., Atmos Env., 39(7), 1303 (2005):
*     From the evaluation and optimisation of the 
*     MCMv3 mechanism using SAPRC chamber data, and
*     Carter (2000): From the evaluation and optimisation 
*     of the SAPRC-99 mechanism using chamber data.
*     An average of the two optimised qy's is used.  
*     Again this is originally for 2-butanone

      mabs = 1
      myld = 1

	  IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE=
     >	    'DATAJ1/((CH3)2CH)2CO/((CH3)2CH)2CO_Yujing.abs',STATUS='old')
            DO i = 1, 11
               READ(kin,*)
            ENDDO
            n = 111
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	         y1(i) = y1(i)*1.e-20
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
	   ENDIF

*combine

      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            IF(wc(iw) .GE. 240. .AND. wc(iw) .LE. 350.) THEN
          
              sig = yg(iw)
              sig = MAX(sig, 0.)
              
            ELSE
            
              sig = 0.
              
            ENDIF   
           
* quantum yields:
* 1:  IUPAC 2005, Raber and Moortgat 1996 (0.34)
* 2:  Average of Pinho et al. 2005 (0.17) and 
*     Carter 2000 (0.15)

            IF(myld .EQ. 1) THEN
              
              qy1 = 0.34   
        
            ELSEIF(myld .EQ. 2) THEN
              
              qy1 = 0.16
              
            ELSE
           
              qy1 = 0.
              
            ENDIF      
        
            sq(j  ,iz,iw) = sig * qy1
            
          ENDDO
      ENDDO
       
      END      

*=============================================================================*

      SUBROUTINE r145 (nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  4-methyl-2-pentanone photolysis:                                         =*
*=         CH3C(O)CH2CH(CH3)2 + hv -> CH3CO + (CH3)2CHCH2                    =*
*=                                 -> CH3C(O)CH3 + CH2=CHCH3                 =*
*=                                                                           =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     from Yujing and Mellouki (2000)                       =*
*=                 						                                     =*
*=  Quantum yield: (1) Raber and Moortgat 1996 (2-butanone)                  =*
*=                     c.f. Laval-Szopa thesis (2003)                        =*
*=                 (2) Pinho et al. 2005 and Carter (2000) (2-butanone)		 =*
*=                 						                                     =*
*=              Norrish I (30%) and II (70%) branching taken from            =*
*=              Laval-Szopa thesis (2003), based on results of Carter (2000) =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
************ 4-methyl-2-pentanone photodissociation ************

      j = j+1
      jlabel(j) = 'CH3C(O)CH2CH(CH3)2 -> CH3CO + (CH3)2CHCH2'
      
      j = j+1
      jlabel(j) = 'CH3C(O)CH2CH(CH3)2 -> CH3C(O)CH3 + CH2=CHCH3'
                         
* working grid arrays:
*     yg = cross section at a specific temperature (298K)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Yujing and Mellouki (2000)

* Quantum yield:
* 1:  Raber and Moortgat 1996 (for 2-butanone)
* 2.  Pinho et al., Atmos Env., 39(7), 1303 (2005):
*     From the evaluation and optimisation of the 
*     MCMv3 mechanism using SAPRC chamber data, and
*     Carter (2000): From the evaluation and optimisation 
*     of the SAPRC-99 mechanism using chamber data.
*     An average of the two optimised qy's is used.  
*     Again this is originally for 2-butanone
*
*     Norrish I (30%) and Norrish II (70%) branching taken from
*     Laval-Szopa thesis (2003), based on results of Carter (2000)
 
      mabs = 1
      myld = 1

	  IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE='DATAJ1/CH3C(O)CH2CH(CH3)2' //
     $'/CH3C(O)CH2CH(CH3)2_Yujing.abs',STATUS='old')
            DO i = 1, 11
               READ(kin,*)
            ENDDO
            n = 111
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	         y1(i) = y1(i)*1.e-20
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
	   ENDIF

*combine

      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            IF(wc(iw) .GE. 240. .AND. wc(iw) .LE. 350.) THEN
          
              sig = yg(iw)
              yg5 = sig
              sig = MAX(sig, 0.)
              
            ELSE
            
              sig = 0.
              
            ENDIF   
           
* quantum yields:
*     Norrish I (30%) and Norrish II (70%) branching taken from
*     Laval-Szopa thesis (2003), based on results of Carter (2000)
*
* 1:  IUPAC 2005, Raber and Moortgat 1996 (0.34)
* 2:  Average of Pinho et al. 2005 (0.17) and 
*     Carter 2000 (0.15)

            IF(myld .EQ. 1) THEN
              
              qy1 = 0.34*0.3
              qy2 = 0.34*0.7   
        
            ELSEIF(myld .EQ. 2) THEN
              
              qy1 = 0.16*0.3
              qy2 = 0.16*0.7
              
            ELSE
           
              qy1 = 0.
              qy2 = 0.
              
            ENDIF      
        
            sq(j-1,iz,iw) = sig * qy1
            sq(j  ,iz,iw) = sig * qy2
            
          ENDDO
      ENDDO
       
      END      

*=============================================================================*

      SUBROUTINE r146 (nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  5-methyl-2-hexanone photolysis:                                          =*
*=         CH3C(O)CH2CH2CH(CH3)2 + hv -> CH3CO + (CH3)2CHCH2CH2              =*
*=                                    -> CH3C(O)CH3 + CH2=CHCH(CH3)2         =*
*=                                                                           =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     from Yujing and Mellouki (2000)                       =*
*=                 						                                     =*
*=  Quantum yield: (1) Raber and Moortgat 1996 (2-butanone)                  =*
*=                     c.f. Laval-Szopa thesis (2003)                        =*
*=                 (2) Pinho et al. 2005 and Carter (2000) (2-butanone)		 =*
*=                 						                                     =*
*=              Norrish I (30%) and II (70%) branching taken from            =*
*=              Laval-Szopa thesis (2003), based on results of Carter (2000) =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
************ 5-methyl-2-hexanone photodissociation *************

      j = j+1
      jlabel(j) = 'CH3C(O)C3H5(CH3)2 -> CH3CO + (CH3)2C3H5'
      
      j = j+1
      jlabel(j) = 'CH3C(O)C3H5(CH3)2 -> CH3C(O)CH3 + CH2=CHCH(CH3)2'
                         
* working grid arrays:
*     yg = cross section at a specific temperature (298K)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from Yujing and Mellouki (2000)

* Quantum yield:
* 1:  Raber and Moortgat 1996 (for 2-butanone)
* 2.  Pinho et al., Atmos Env., 39(7), 1303 (2005):
*     From the evaluation and optimisation of the 
*     MCMv3 mechanism using SAPRC chamber data, and
*     Carter (2000): From the evaluation and optimisation 
*     of the SAPRC-99 mechanism using chamber data.
*     An average of the two optimised qy's is used.  
*     Again this is originally for 2-butanone
*
*     Norrish I (30%) and Norrish II (70%) branching taken from
*     Laval-Szopa thesis (2003), based on results of Carter (2000)
 
      mabs = 1
      myld = 1

	  IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE='DATAJ1/CH3C(O)C3H5(CH3)2' //
     $'/CH3C(O)C3H5(CH3)2_Yujing.abs',STATUS='old')
            DO i = 1, 11
               READ(kin,*)
            ENDDO
            n = 111
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	         y1(i) = y1(i)*1.e-20
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
	   ENDIF

*combine

      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            IF(wc(iw) .GE. 240. .AND. wc(iw) .LE. 350.) THEN
          
              sig = yg(iw)
              yg5 = sig
              sig = MAX(sig, 0.)
              
            ELSE
            
              sig = 0.
              
            ENDIF   
           
* quantum yields:
*     Norrish I (30%) and Norrish II (70%) branching taken from
*     Laval-Szopa thesis (2003), based on results of Carter (2000)
*
* 1:  IUPAC 2005, Raber and Moortgat 1996 (0.34)
* 2:  Average of Pinho et al. 2005 (0.17) and 
*     Carter 2000 (0.15)

            IF(myld .EQ. 1) THEN
              
              qy1 = 0.34*0.3
              qy2 = 0.34*0.7   
        
            ELSEIF(myld .EQ. 2) THEN
              
              qy1 = 0.16*0.3
              qy2 = 0.16*0.7
              
            ELSE
           
              qy1 = 0.
              qy2 = 0.
              
            ENDIF      
        
            sq(j-1,iz,iw) = sig * qy1
            sq(j  ,iz,iw) = sig * qy2
            
          ENDDO
      ENDDO
       
      END      


*=============================================================================*

      SUBROUTINE r147 (nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for                =*
*=  methyl peroxynitrate photolysis:                                         =*
*=         CH3OONO2 + hv -> CH3O2 + NO2                                      =*
*=                       -> CH3O + NO3                                       =*
*=                                                                           =*
*=  Cross section: (1) MPI-Mainz Spectral Atlas,                             =*
*=                     from IUPAC 2005:                                      =*
*=  200-280 nm:  I. Bridier, R. Lesclaux, and B. Veyret,                     =*
*=               Chem. Phys. Lett. 191, 259 (1992).                          =*
*=  > 290 nm:    Comparison with HO2NO2 spectrum                             =*
*=                                                                           =*                
*=  Quantum yield: (1) qy1 + qy2 = 1; each channel is assumed to be 50%      =*
*=                                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n
      real x(kdata), y(kdata)
      real xl(kdata), xc(kdata), xu(kdata)
      INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata), x5(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata), y5(kdata)

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mabs, myld

*_______________________________________________________________________

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

****************************************************************
************ methyl peroxynitrate photodissociation ************

      j = j+1
      jlabel(j) = 'CH3OONO2 -> CH3O2 + NO2'
      
      j = j+1
      jlabel(j) = 'CH3OONO2 -> CH3O + NO3'
      
* working grid arrays:
* yg = cross section at a specific temperature (298K)

* options
* mabs for cross sections
* myld for quantum yields

* Absorption Cross Section:
* 1:  MPI-Mainz Spectral Atlas, from IUPAC 2005

* Quantum yield:
* 1:  IUPAC 2005: 0.5 each channel
 
      mabs = 1
      myld = 1

	  IF(mabs .EQ. 1) THEN

            OPEN(UNIT=kin,FILE='DATAJ1/RONO2/CH3O2NO2_iup05.abs',
     $           STATUS='old')
            DO i = 1, 7
               READ(kin,*)
            ENDDO
            n = 26
            DO i = 1, n
               READ(kin,*) x1(i), y1(i)
	         y1(i) = y1(i)*1.e-20
            ENDDO
            CLOSE(kin)		
	    
            CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
            CALL addpnt(x1,y1,kdata,n,               0.,0.)
            CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
            CALL addpnt(x1,y1,kdata,n,            1E+38,0.)
            CALL inter2(nw,wl,yg,n,x1,y1,ierr)
            IF (ierr .NE. 0) THEN
               WRITE(*,*) ierr, jlabel(j)
               STOP
            ENDIF
            
	   ENDIF

*combine

      DO iw = 1, nw - 1
          DO iz = 1, nz
          
            IF(wc(iw) .GE. 200. .AND. wc(iw) .LE. 325.) THEN
          
              sig = yg(iw)
              sig = MAX(sig, 0.)
              
            ELSE
            
              sig = 0.
              
            ENDIF   
           
* quantum yields:
*     qy1 + qy2 = 1.0; 0.5 each channel

            IF(myld .EQ. 1) THEN
              
              qy1 = 0.5
              qy2 = 0.5   
              
            ELSE
           
              qy1 = 0.
              qy2 = 0.
              
            ENDIF      
        
            sq(j-1,iz,iw) = sig * qy1
            sq(j  ,iz,iw) = sig * qy2
            
          ENDDO
      ENDDO
       
      END      
