* This file contains the following subroutine, related to foil trasmittance
* in a smog chamber:
*     foil

      SUBROUTINE foil(trans,nw,wl)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read the transmission of the smog chamber foil from a file, regrid it    =*
*=  according to the current wavelength grid. Wavelengths are in nm.         =*
*=  No gaps are allowed within the wavelength grid.                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  TRANS  - REAL, vector carrying the rescaled transmission data       (I/O)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'params'

* input:

      INTEGER nw
	INTEGER wl(nw)

* output:

      REAL trans(nw)

* data arrays:

      INTEGER kdata
      PARAMETER (kdata = 550)
	REAL x(kdata), y(kdata)

* local:
	
	INTEGER iw, i, n
	INTEGER ierr

*_______________________________________________________________________

**** read in data, measured by Bill Bloss

* foil_bloss.tra

      OPEN(UNIT=kin,FILE='DATAE1/foil_bloss.tra',
     $     STATUS='old')

      DO i = 1, 3
         READ(kin,*) 
      ENDDO
      n = 501
      DO i = 1, n
         READ(kin,*) x(i),  y(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),y(1))
      CALL addpnt(x,y,kdata,n,               0.,y(1))
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),y(n))
      CALL addpnt(x,y,kdata,n,            1.e+38,y(n))
      CALL inter2(nw,wl,trans,n,x,y,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr
         STOP
      ENDIF

	END
