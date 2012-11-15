      SUBROUTINE swbiol3(nw,wl,wc,j,s,label)

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
      integer n0, n1, n2, n3, n4, n5, n6, n7

      REAL x0(kdata), x1(kdata), x2(kdata), x3(kdata), x4(kdata), 
     $     x5(kdata), x6(kdata), x7(kdata)

      REAL y0(kdata), y1(kdata), y2(kdata), y3(kdata), y4(kdata), 
     $     y5(kdata), y6(kdata), y7(kdata)

      REAL yg(kw)

      real dum

      INTEGER i, iw
      INTEGER ierr

*********YANKEE Filter functions
* first column is wavelength
* second column is unknown
* third column is diode response function
* fourth-tenth columns are filter transmission functions, normalized to unity
*   at peak response.
* I assume that the total response for each filter is product of 
*   (diode response * filter transmission)
* Until absolute calibration factors are known, the computed signal is only
* relative, not absolute W m-2 nm-1.
* However, even these relative signals should show the correct:
*     - zenith (time) dependence for each filter
*     - diffuse/direct ratios for each filter
*     - sensitivity (relative) to ozone, aerosols, clouds, etc.

* The data are read from files (one for YES231, another for YES232), then
* end-points are added at both extremes (low, high wavelegths) in order to
* span grid.  I assumed zero response beyond data range.  Then they are
* interpolated to our working grid, wg, and finally assigned to s(j,iw) which
* is the weighting function for the irradiance.

* instrument #YES231

      OPEN(UNIT=kin,FILE='DATAS3/YES231_fil.mod',STATUS='old')
      do i = 1, 3
         read(kin,*)
      enddo
      n0 = 600
      DO i = 1, n0
         READ(kin,*) x0(i), dum, Y0(i),
     $        y1(i),y2(i),y3(i),y4(i),y5(i),y6(i),y7(i)
         x1(i) = x0(i)
         x2(i) = x0(i)
         x3(i) = x0(i)
         x4(i) = x0(i)
         x5(i) = x0(i)
         x6(i) = x0(i)
         x7(i) = x0(i)
         y1(i) = y1(i) * y0(i)
         y2(i) = y2(i) * y0(i)
         y3(i) = y3(i) * y0(i)
         y4(i) = y4(i) * y0(i)
         y5(i) = y5(i) * y0(i)
         y6(i) = y6(i) * y0(i)
         y7(i) = y7(i) * y0(i)
      ENDDO
      CLOSE (kin)
      n1 = n0
      n2 = n0
      n3 = n0
      n4 = n0
      n5 = n0
      n6 = n0
      n7 = n0

      j = j + 1
      label(j) = 'YES231-300'
      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,          0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n1,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES231-305'
      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,          0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),   0.)
      CALL addpnt(x2,y2,kdata,n2,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES231-311'
      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,          0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),   0.)
      CALL addpnt(x3,y3,kdata,n3,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES231-317'
      CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,          0.,0.)
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),   0.)
      CALL addpnt(x4,y4,kdata,n4,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n4,x4,y4,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES231-325'
      CALL addpnt(x5,y5,kdata,n5,x5(1)*(1.-deltax),0.)
      CALL addpnt(x5,y5,kdata,n5,          0.,0.)
      CALL addpnt(x5,y5,kdata,n5,x5(n5)*(1.+deltax),   0.)
      CALL addpnt(x5,y5,kdata,n5,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n5,x5,y5,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES231-332'
      CALL addpnt(x6,y6,kdata,n6,x6(1)*(1.-deltax),0.)
      CALL addpnt(x6,y6,kdata,n6,          0.,0.)
      CALL addpnt(x6,y6,kdata,n6,x6(n6)*(1.+deltax),   0.)
      CALL addpnt(x6,y6,kdata,n6,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n6,x6,y6,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES231-368'
      CALL addpnt(x7,y7,kdata,n7,x7(1)*(1.-deltax),0.)
      CALL addpnt(x7,y7,kdata,n7,          0.,0.)
      CALL addpnt(x7,y7,kdata,n7,x7(n7)*(1.+deltax),   0.)
      CALL addpnt(x7,y7,kdata,n7,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n7,x7,y7,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

* YES232:

      OPEN(UNIT=kin,FILE='DATAS3/YES232_fil.mod',STATUS='old')
      do i = 1, 3
         read(kin,*)
      enddo
      n0 = 600
      DO i = 1, n0
         READ(kin,*) x0(i), dum, Y0(i),
     $        y1(i),y2(i),y3(i),y4(i),y5(i),y6(i),y7(i)
         x1(i) = x0(i)
         x2(i) = x0(i)
         x3(i) = x0(i)
         x4(i) = x0(i)
         x5(i) = x0(i)
         x6(i) = x0(i)
         x7(i) = x0(i)
         y1(i) = y1(i) * y0(i)
         y2(i) = y2(i) * y0(i)
         y3(i) = y3(i) * y0(i)
         y4(i) = y4(i) * y0(i)
         y5(i) = y5(i) * y0(i)
         y6(i) = y6(i) * y0(i)
         y7(i) = y7(i) * y0(i)
      ENDDO
      CLOSE (kin)
      n1 = n0
      n2 = n0
      n3 = n0
      n4 = n0
      n5 = n0
      n6 = n0
      n7 = n0

      j = j + 1
      label(j) = 'YES232-300'
      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,          0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n1,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES232-305'
      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,          0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),   0.)
      CALL addpnt(x2,y2,kdata,n2,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES232-311'
      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,          0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),   0.)
      CALL addpnt(x3,y3,kdata,n3,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES232-317'
      CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,          0.,0.)
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),   0.)
      CALL addpnt(x4,y4,kdata,n4,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n4,x4,y4,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES232-325'
      CALL addpnt(x5,y5,kdata,n5,x5(1)*(1.-deltax),0.)
      CALL addpnt(x5,y5,kdata,n5,          0.,0.)
      CALL addpnt(x5,y5,kdata,n5,x5(n5)*(1.+deltax),   0.)
      CALL addpnt(x5,y5,kdata,n5,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n5,x5,y5,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES232-332'
      CALL addpnt(x6,y6,kdata,n6,x6(1)*(1.-deltax),0.)
      CALL addpnt(x6,y6,kdata,n6,          0.,0.)
      CALL addpnt(x6,y6,kdata,n6,x6(n6)*(1.+deltax),   0.)
      CALL addpnt(x6,y6,kdata,n6,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n6,x6,y6,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

      j = j + 1
      label(j) = 'YES232-368'
      CALL addpnt(x7,y7,kdata,n7,x7(1)*(1.-deltax),0.)
      CALL addpnt(x7,y7,kdata,n7,          0.,0.)
      CALL addpnt(x7,y7,kdata,n7,x7(n7)*(1.+deltax),   0.)
      CALL addpnt(x7,y7,kdata,n7,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n7,x7,y7,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

****************************************************************
****************************************************************

*_______________________________________________________________________

      IF (j .GT. ks) STOP '1001'
*_______________________________________________________________________

      RETURN
      END
