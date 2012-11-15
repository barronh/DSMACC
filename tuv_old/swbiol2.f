      SUBROUTINE swbiol2(nw,wl,wc,j,s,label)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create or read various weighting functions, e.g. biological action       =*
*=  spectra, instrument responses etc.                                       =*
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
*=  LABEL  - CHARACTER*50, string identifier for each weighting function  (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Change offset for grid-end interpolation to relative number       =*
*=         (x * (1 +- deltax))                                               =*
*=  05/96  Rename from LOADW2 to WSPEC2                                      =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input:
      REAL wl(kw), wc(kw)
      INTEGER nw

* input/output
      INTEGER j      

* output: (weighting functions)
      REAL s(ks,kw)
      CHARACTER*50 label(ks)

* internal:
      REAL x1(kdata)
      REAL y1(kdata)
      REAL yg(kw)

      REAL yy1, yy2
      REAL dum
      INTEGER i, n, iw
      INTEGER ierr
      REAL afit, bfit
      REAL chl, chl0, xk0, xkc, xkc1, xkc2, xkt, asp

      j = j + 1
      label(j) = 'Melanoma in fish'
      OPEN(UNIT=kin,FILE='DATAS2/cmm.fish',STATUS='old')
      READ(kin,*)
      READ(kin,*)
      n = 5
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) *  x1(i)/300.
      ENDDO

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
      CLOSE (kin)

********** RB Meter, old Temple model
* spectrum received from J.Frederick  1989.
* presumably per unit energy, not unit photon
* added  point at 280.0 nm by extrapolating from 300-290

      j = j + 1
      label(j) = 'RB Meter, model UT'
      OPEN(UNIT=kin,FILE='DATAS2/rbm.jef',STATUS='old')
      n = 96
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO

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
      CLOSE (kin)

***********

      j = j + 1
      label(j) = 'Banrud, binding'
      OPEN(UNIT=kin,FILE='DATAS2/banrud',STATUS='old')
      n = 4
      DO i = 1, n
         READ(kin,*) x1(i), y1(i), dum
         y1(i) = 10.**y1(i)
      ENDDO

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
      CLOSE (kin)

      j = j + 1
      label(j) = 'Banrud, survival'
      OPEN(UNIT=kin,FILE='DATAS2/banrud',STATUS='old')
      n = 4
      DO i = 1, n
         READ(kin,*) x1(i), dum, y1(i)
         y1(i) = 10.**y1(i)
      ENDDO

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
      CLOSE (kin)
      
**********

      j = j + 1
      label(j) = 'phot.syn. inh., Nodularia spumigena'
      OPEN(UNIT=kin,FILE='DATAS2/Nodul_spum',STATUS='old')
      n = 19
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO

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
      CLOSE (kin)
      
**********  from R. ZEPP, Sept. 94.  All are per quantum

      j = j + 1
      label(j) = 'CO.suwannee'
      OPEN(UNIT=kin,FILE='DATAS2/CO.suwannee',STATUS='old')
      n = 33
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i)  *  x1(i)/300.
      ENDDO

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
      CLOSE (kin)
      
      j = j + 1
      label(j) = 'COS.gmexico'
      OPEN(UNIT=kin,FILE='DATAS2/COS.gmexico',STATUS='old')
      n = 17
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i)  *  x1(i)/300.
      ENDDO

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
      CLOSE (kin)
      
      j = j + 1
      label(j) = 'COS.northsea'
      OPEN(UNIT=kin,FILE='DATAS2/COS.northsea',STATUS='old')
      n = 12
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i)  *  x1(i)/300.
      ENDDO

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
      CLOSE (kin)
      
      j = j + 1
      label(j) = 'NO3-'
      OPEN(UNIT=kin,FILE='DATAS2/NO3-',STATUS='old')
      n = 8
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i)  *  x1(i)/300.
      ENDDO

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
      CLOSE (kin)

      j = j + 1
      label(j) = 'CH2O.biscayne'
      OPEN(UNIT=kin,FILE='DATAS2/CH2O.biscayne',STATUS='old')
      n = 7
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i)  *  x1(i)/300.
      ENDDO

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
      CLOSE(kin)

      j = j + 1
      label(j) = 'H2O2.cooper'
      OPEN(UNIT=kin,FILE='DATAS2/H2O2.cooper',STATUS='old')
      n = 19
      DO i = 1, n
         READ(kin,*) x1(i), yy1, yy2
         y1(i) = yy1 * yy2  *  x1(i)/300.
      ENDDO


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
      CLOSE (kin)

      j = j + 1
      label(j) = 'alfalfa.sut'
      OPEN(UNIT=kin,FILE='DATAS2/alfalfa.sut',STATUS='old')
      n = 9
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) *  x1(i)/300.
         y1(i) = ALOG(y1(i))
      ENDDO

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
         IF(yg(iw) .EQ. 0.) THEN
            s(j,iw) = 0.
         ELSE
            s(j,iw) = EXP(yg(iw))
         ENDIF
      ENDDO
      CLOSE(kin)

**********
* Coral spectrum from Michael Lesser

      j = j + 1
      label(j) = 'coral.mpl'
      OPEN(UNIT=kin,FILE='DATAS2/coral.mpl',STATUS='old')
      DO i = 1, 9
         READ(kin,*)
      ENDDO
      n = 51
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO

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
      CLOSE (kin)
      
* compute Zepp and Andreae 1990 COS action spectrum
* exp fit by R. Najjar

      j = j + 1
      label(j) = 'cos.zepp.nj fit'
      afit = 14.521
      bfit = -0.0484
      DO iw = 1, nw - 1
         s(j,iw) = 10.**(afit + bfit*wc(iw))
      ENDDO

**** Zepp COS spectrum modified by Najjar et al. to account for
* transmission of water
** obviously does not make sense:
* what is meaning of data in cos.zepp.nj (old "moddat")?
* had mismatched parens in eq. for xkc.

      j = j + 1
      label(j) = 'cos.zepp.nj'
      OPEN(UNIT=kin,FILE='DATAS2/cos.zepp.nj',STATUS='old')
      n = 125
      chl = 1.5
      chl0 = 0.5
      DO i = 1, n
         READ(kin,*) xk0, xkc1, xkc2
         xkc = xkc1 * chl * EXP( - xkc2**2 * (ALOG10(chl)-
     $        ALOG10(chl0))**2 + 0.001*(chl)**2)
         xkt = xkc + xk0
      ENDDO

      DO iw = 1, nw-1
         asp = (3.3161E+14)*(10.**((-4.84E-2)*wc(iw)))
         s(j,iw) = asp / xkt
      ENDDO
      CLOSE(kin)


      j = j + 1
      label(j) = 'rocky.dna'
      OPEN(UNIT=kin,FILE='DATAS2/rocky.data',STATUS='old')
      read(kin,*)
      n = 55
      DO i = 1, n
         READ(kin,*) x1(i), dum, y1(i)
         y1(i) = y1(i)/ 0.0329959
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
      
**************** read dna damage action spectrum
* from: Setlow, R. B., The wavelengths in sunlight effective in 
*       producing skin cancer: a theoretical analysis, Proceedings 
*       of the National Academy of Science, 71, 3363 -3366, 1974.
* normalize to unity at 300 nm
* orig. data per quantum.  put on energy basis

      j = j + 1
      label(j) = 'Setlow dna'
      OPEN(UNIT=kin,FILE='DATAS2/dna.setlow',STATUS='old')
      n = 48
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) / 2.946E-02  *  x1(i)/300.
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

***** dna damage from john frederick

      j = j + 1
      label(j) = 'Setlow dna.jef'
      OPEN(UNIT=kin,FILE='DATAS2/dna.jef',STATUS='old')
      DO i = 1, 7
         READ(kin,*)
      ENDDO
      n = 96
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) / 5.850E-02  *  x1(i)/300.
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

********** cie-erythema
* spectrum received from J.Frederick  1989.
* presumably per unit energy, not unit photon
* added  point at 280.0 nm by extrapolating from 300-290

      j = j + 1
      label(j) = 'ery.jef'
      OPEN(UNIT=kin,FILE='DATAS2/ery.jef',STATUS='old')
      READ(kin,*)
      READ(kin,*)
      n = 105
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO

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
      CLOSE (kin)


      j = j + 1
      label(j) = 'cataract.rat'
      OPEN(UNIT=kin,FILE='DATAS2/cataract.rat',STATUS='old')
      READ(kin,*)
      READ(kin,*)
      READ(kin,*)
      n = 5
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
c         write(20,*) wc(iw), s(j,iw)
      ENDDO
      CLOSE (kin)


****************************************************************
****************************************************************

*_______________________________________________________________________

      IF (j .GT. ks) STOP '1001'
*_______________________________________________________________________

      RETURN
      END
