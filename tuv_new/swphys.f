* This file contains the following subroutines, related to specifying 
* physical spectral weighting functions:
*     swphys

*=============================================================================*

      SUBROUTINE swphys(nw,wl,wc,j,s,label)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create or read various spectral weighting functions, physically-based    =*
*=  e.g. UV-B, UV-A, visible ranges, instrument responses, etc.              =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of central wavelength of wavelength intervals    I)=*
*=           in working wavelength grid                                      =*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  S      - REAL, value of each defined weighting function at each       (O)=*
*=           defined wavelength                                              =*
*=  LABEL  - CHARACTER*40, string identifier for each weighting function  (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input:
      REAL wl(kw), wc(kw)
      INTEGER nw

* input/output:
      INTEGER j

* output: (weighting functions and labels)
      REAL s(ks,kw)
      CHARACTER*50 label(ks)

* internal:
      REAL x1(kdata)
      REAL y1(kdata)
      REAL yg(kw)

      INTEGER i, iw, n

      INTEGER ierr

      INTEGER idum
      REAL dum1, dum2
      REAL em, a, b, c
      REAL sum

*_______________________________________________________________________

      j = 0

********* UV-B (280-315 nm)
 
      j = j + 1
      label(j) = 'UV-B, 280-315 nm'
      DO iw = 1, nw-1
         IF (wc(iw) .GT. 280. .AND. wc(iw) .LT. 315.) THEN
            s(j,iw) = 1.
         ELSE
            s(j,iw) = 0.
         ENDIF
      ENDDO

********* UV-B* (280-320 nm)
 
      j = j + 1
      label(j) = 'UV-B*, 280-320 nm'
      DO iw = 1, nw-1
         IF (wc(iw) .GT. 280. .AND. wc(iw) .LT. 320.) THEN
            s(j,iw) = 1.
         ELSE
            s(j,iw) = 0.
         ENDIF
      ENDDO

********* UV-A (315-400 nm)
 
      j = j + 1
      label(j) = 'UV-A, 315-400 nm'
      DO iw = 1, nw-1
         IF (wc(iw) .GT. 315. .AND. wc(iw) .LT. 400.) THEN
            s(j,iw) = 1.
         ELSE
            s(j,iw) = 0.
         ENDIF
      ENDDO

********* visible+ (> 400 nm)
 
      j = j + 1
      label(j) = 'vis+, > 400 nm'
      DO iw = 1, nw-1
         IF (wc(iw) .GT. 400.) THEN
            s(j,iw) = 1.
         ELSE
            s(j,iw) = 0.
         ENDIF
      ENDDO

**********  Gaussian transmission functions

      j = j + 1
      label(j) = 'Gaussian, 305 nm, 10 nm FWHM'
      sum = 0.
      DO iw = 1, nw-1
         s(j,iw) = exp(- ( log(2.) * ((wc(iw)-305.)/(5.))**2) )
         sum = sum + s(j,iw)
      ENDDO
      DO iw = 1, nw-1
         s(j,iw) = s(j,iw)/sum
      ENDDO

      j = j + 1
      label(j) = 'Gaussian, 320 nm, 10 nm FWHM'
      sum = 0.
      DO iw = 1, nw-1
         s(j,iw) = exp(- ( log(2.) * ((wc(iw)-320.)/(5.))**2) )
         sum = sum + s(j,iw)
      ENDDO
      DO iw = 1, nw-1
         s(j,iw) = s(j,iw)/sum
      ENDDO

      j = j + 1
      label(j) = 'Gaussian, 340 nm, 10 nm FWHM'
      sum = 0.
      DO iw = 1, nw-1
         s(j,iw) = exp(- ( log(2.) * ((wc(iw)-340.)/(5.))**2) )
         sum = sum + s(j,iw)
      ENDDO
      DO iw = 1, nw-1
         s(j,iw) = s(j,iw)/sum
      ENDDO

      j = j + 1
      label(j) = 'Gaussian, 380 nm, 10 nm FWHM'
      sum = 0.
      DO iw = 1, nw-1
         s(j,iw) = exp(- ( log(2.) * ((wc(iw)-380.)/(5.))**2) )
         sum = sum + s(j,iw)
      ENDDO
      DO iw = 1, nw-1
         s(j,iw) = s(j,iw)/sum
      ENDDO

********** RB Meter, model 501
*  private communication, M. Morys (Solar Light Co.), 1994.
* From: morys@omni.voicenet.com (Marian Morys)
* Received: from acd.ucar.edu by sasha.acd.ucar.edu (AIX 3.2/UCB 5.64/4.03)
*          id AA17274; Wed, 21 Sep 1994 11:35:44 -0600

      j = j + 1
      label(j) = 'RB Meter, model 501'
      OPEN(UNIT=kin,FILE='DATAS1/rbm.501',STATUS='old')
      n = 57
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)
      
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
            
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

********** Eppley UV Photometer
* CAUTION:
*    Contrary to the manufacturer's claims, the Eppley Total UV photometer does 
* not measure integrated uv radiation in any wavelength band.  It actually 
* measures the *average* spectral irradiance over 295-385 nm, with the 
* averaging kernel given by the instrument's spectral response function 
* (which is the product of the photodiode response and the filter transmission 
* function). 
*    The calibration factor provided by Eppley Laboratories Inc. is only valid 
* when the photometer is exposed to a spectrum that has the same shape as their
* calibration lamp.  The calibration factor is not accurate for exposure to the 
* solar spectrum.  The values of Eppley total UV irradiance for outdoor 
* observations, although reported as W m-2, should be used with caution and 
* cannot be associated to any specific physical quantity, such as UVA, 
* UVA + UVB, the 295-385 nm band, or the 295-400 nm band.
*    The value calculated here is the average spectral irradiance using the 
* Eppley repsonse kernel, then multiplied by the nominal bandwith of 90 nm 
* (295-385 nm). This value is about 10% higher than the un-weighted integral 
* over 295-385 nm, and about 10% lower than the UVA (315-400 nm).
* - Sasha Madronich, March 2009

      j = j + 1
      label(j) = 'Eppley UV Photometer'
      OPEN(UNIT=kin,FILE='DATAS1/eppley_uv',STATUS='old')
      DO i = 1, 6
         READ(kin,*)
      ENDDO
      n = 19
      DO i = 1, n
         READ(kin,*) x1(i), dum1, dum2
         y1(i) = dum1*dum2
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF

* compute normalization

      sum = 0.
      DO iw = 1, nw-1
         sum = sum + yg(iw)*(wl(iw+1) - wl(iw))
      ENDDO

      DO iw = 1, nw-1
         s(j,iw) = 90.*yg(iw)/sum
      ENDDO

****************************************************************
****************************************************************

*_______________________________________________________________________

      IF (j .GT. ks) STOP '1001'
*_______________________________________________________________________

      RETURN
      END
