* This file contains the following subroutines, related to the
* spherical geometry of the Earth's atmosphere
*     sphers
*     airmas
*=============================================================================*

      SUBROUTINE sphers(nz, z, zen, dsdh, nid)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate slant path over vertical depth ds/dh in spherical geometry.    =*
*=  Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model =*
*=  for computing the radiation field available for photolysis and heating   =*
*=  at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
*=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
*=            when travelling from the top of the atmosphere to layer i;     =*
*=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
*=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
*=            travelling from the top of the atmosphere to layer i;          =*
*=            NID(i), i = 0..NZ-1                                            =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  double precision fix for shallow layers - Julia Lee-Taylor Dec 2000      =*
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

* input
      INTEGER nz
      REAL zen, z(kz)

* output
      INTEGER nid(0:kz)
      REAL dsdh(0:kz,kz)

* more program constants
      REAL re, ze(kz)
      REAL  dr
      PARAMETER ( dr = pi/180.)

* local 

      DOUBLE PRECISION zenrad, rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm
      INTEGER i, j, k
      INTEGER id

      INTEGER nlayer
      REAL zd(0:kz-1)

*-----------------------------------------------------------------------------

      zenrad = zen*dr

* number of layers:
      nlayer = nz - 1

* include the elevation above sea level to the radius of the earth:
      re = radius + z(1)
* correspondingly z changed to the elevation above earth surface:
      DO k = 1, nz
         ze(k) = z(k) - z(1)
      END DO

* inverse coordinate of z
      zd(0) = ze(nz)
      DO k = 1, nlayer
        zd(k) = ze(nz - k)
      END DO

* initialize dsdh(i,j), nid(i)
      DO i = 0, kz
       nid(i) = 0
       DO j = 1, kz
        dsdh(i,j) = 0.
       END DO
      END DO

* calculate ds/dh of every layer
      DO 100 i = 0, nlayer

        rpsinz = (re + zd(i)) * SIN(zenrad)
 
        IF ( (zen .GT. 90.0) .AND. (rpsinz .LT. re) ) THEN
           nid(i) = -1
        ELSE

*
* Find index of layer in which the screening height lies
*
           id = i 
           IF( zen .GT. 90.0 ) THEN
              DO 10 j = 1, nlayer
                 IF( (rpsinz .LT. ( zd(j-1) + re ) ) .AND.
     $               (rpsinz .GE. ( zd(j) + re )) ) id = j
 10           CONTINUE
           END IF
 
           DO 20 j = 1, id

             sm = 1.0
             IF(j .EQ. id .AND. id .EQ. i .AND. zen .GT. 90.0)
     $          sm = -1.0
 
             rj = re + zd(j-1)
             rjp1 = re + zd(j)
 
             dhj = zd(j-1) - zd(j)
 
             ga = rj*rj - rpsinz*rpsinz
             gb = rjp1*rjp1 - rpsinz*rpsinz
             IF (ga .LT. 0.0) ga = 0.0
             IF (gb .LT. 0.0) gb = 0.0
 
             IF(id.GT.i .AND. j.EQ.id) THEN
                dsj = SQRT( ga )
             ELSE
                dsj = SQRT( ga ) - sm*SQRT( gb )
             END IF
             dsdh(i,j) = dsj / dhj
 20        CONTINUE
 
           nid(i) = id
 
        END IF

 100  CONTINUE

*-----------------------------------------------------------------------------

      RETURN
      END

*=============================================================================*

      SUBROUTINE airmas(nz, dsdh, nid, cz,
     $      vcol, scol)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate vertical and slant air columns, in spherical geometry, as a    =*
*=  function of altitude.                                                    =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
*=            when travelling from the top of the atmosphere to layer i;     =*
*=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
*=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
*=            travelling from the top of the atmosphere to layer i;          =*
*=            NID(i), i = 0..NZ-1                                            =*
*=  VCOL    - REAL, output, vertical air column, molec cm-2, above level iz  =*
*=  SCOL    - REAL, output, slant air column in direction of sun, above iz   =*
*=            also in molec cm-2                                             =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* Input:

      INTEGER nz
      INTEGER nid(0:kz)
      REAL dsdh(0:kz,kz)
      REAL cz(kz)

* output: 

      REAL vcol(kz), scol(kz)

* internal:

      INTEGER id, j
      REAL sum, vsum

* calculate vertical and slant column from each level:
* work downward

      vsum = 0.
      DO id = 0, nz - 1
         vsum = vsum + cz(nz-id)
         vcol(nz-id) = vsum
         sum = 0.
         IF(nid(id) .LT. 0) THEN
            sum = largest
         ELSE

* single pass layers:

            DO j = 1, MIN(nid(id), id)
               sum = sum + cz(nz-j)*dsdh(id,j)
            ENDDO

* double pass layers:

            DO j = MIN(nid(id),id)+1, nid(id)
               sum = sum + 2.*cz(nz-j)*dsdh(id,j)
            ENDDO

         ENDIF
         scol(nz - id) = sum

      ENDDO

      RETURN
      END




