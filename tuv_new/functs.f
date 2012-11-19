* This file contains the following user-defined fortran functions:
*     fery
*     fo3qy
*     fo3qy2
*     fsum
*     futr
*=============================================================================*

      FUNCTION fery(w)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the action spectrum value for erythema at a given wavelength   =*
*=  according to: McKinlay, A.F and B.L.Diffey, A reference action spectrum  =*
*=  for ultraviolet induced erythema in human skin, CIE Journal, vol 6,      =*
*=  pp 17-22, 1987.                                                          =*
*=  Value at 300 nm = 0.6486                                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  W - REAL, wavelength (nm)                                             (I)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      REAL w 

* function value:
      REAL fery

      IF (w .LT. 250.) THEN
          fery = 1.
C outside the ery spectrum range
      ELSEIF ((w .GE. 250.) .AND. (w .LT. 298)) THEN
          fery = 1.
      ELSEIF ((w .GE. 298.) .AND. (w .LT. 328.)) THEN
          fery = 10.**( 0.094*(298.-w) )
      ELSEIF ((w .GE. 328.) .AND. (w .LT. 400.)) THEN
          fery = 10.**( 0.015*(139.-w) )
      ELSE
         fery = 1.E-36
C outside the ery spectrum range
      ENDIF

      RETURN
      END

*=============================================================================*

      FUNCTION fo3qy(w,t)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
* function to calculate the quantum yield O3 + hv -> O(1D) + O2,             =*
* according to JPL 2000 recommendation:                                      =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      REAL w, t, kt, fo3qy
      REAL a(3), w0(3), nu(3), om(3)

      DATA a/ 0.887, 2.35, 57./
      DATA w0/ 302., 311.1, 313.9/
      DATA nu/ 0., 820., 1190./
      DATA om/ 7.9, 2.2, 7.4/

      fo3qy = 0.
      kt = 0.695 * t
      
      IF(w .LE. 300.) THEN
         fo3qy = 0.95
      ELSEIF(w .GT. 300. .AND. w .LE. 330.) THEN
         fo3qy = 0.06 + 
     $  a(1)                           *EXP(-((w-w0(1))/om(1))**4)+ 
     $  a(2)*(T/300.)**4*EXP(-nu(2)/kT)*EXP(-((w-w0(2))/om(2))**2)+
     $  a(3)            *EXP(-nu(3)/kT)*EXP(-((w-w0(3))/om(3))**2)
      ELSEIF(w .GT. 330. .AND. w .LE. 345.) THEN
         fo3qy = 0.06
      ELSEIF(w .GT. 345.) THEN
         fo3qy = 0.
      ENDIF

      END

      FUNCTION fo3qy2(w,t)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
* function to calculate the quantum yield O3 + hv -> O(1D) + O2,             =*
* according to:                                                             
* Matsumi, Y., F. J. Comes, G. Hancock, A. Hofzumanhays, A. J. Hynes,
* M. Kawasaki, and A. R. Ravishankara, QUantum yields for production of O(1D)
* in the ultraviolet photolysis of ozone:  Recommendation based on evaluation
* of laboratory data, J. Geophys. Res., 107, 10.1029/2001JD000510, 2002.
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      REAL w, t, kt, fo3qy2
      REAL A(3), X(3), om(3)
      REAL q1, q2 

      DATA A/ 0.8036, 8.9061, 0.1192/
      DATA X/ 304.225, 314.957, 310.737/
      DATA om/ 5.576, 6.601, 2.187/
      
      fo3qy2 = 0.
      kt = 0.695 * t
      q1 = 1.
      q2 = exp(-825.518/kt)
      
      IF(w .LE. 305.) THEN
         fo3qy2 = 0.90
      ELSEIF(w .GT. 305. .AND. w .LE. 328.) THEN

         fo3qy2 = 0.0765 + 
     $  a(1)*             (q1/(q1+q2))*EXP(-((x(1)-w)/om(1))**4)+ 
     $  a(2)*(T/300.)**2 *(q2/(q1+q2))*EXP(-((x(2)-w)/om(2))**2)+
     $  a(3)*(T/300.)**1.5            *EXP(-((x(3)-w)/om(3))**2)

      ELSEIF(w .GT. 328. .AND. w .LE. 340.) THEN
         fo3qy2 = 0.08
      ELSEIF(w .GT. 340.) THEN
         fo3qy2 = 0.
      ENDIF

      END

*=============================================================================*

      FUNCTION fsum(n,x)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Compute the sum of the first N elements of a floating point vector.      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  N  - INTEGER, number of elements to sum                               (I)=*
*=  X  - REAL, vector whose components are to be summed                   (I)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      INTEGER n
      REAL x(n)

* function value:
      REAL fsum

* local:
      INTEGER i

      fsum = 0.
      DO 10, i = 1, n
         fsum=fsum+x(i)
   10 CONTINUE

      RETURN
      END

*=============================================================================*

      FUNCTION futr(w)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the action spectrum value for skin cancer of albino hairless   =*
*=  mice at a given wavelength according to:  deGRuijl, F.R., H.J.C.M.Steren-=*
*=  borg, P.D.Forbes, R.E.Davies, C.Colse, G.Kelfkens, H.vanWeelden,         =*
*=  and J.C.van der Leun, Wavelength dependence of skin cancer induction by  =*
*=  ultraviolet irradiation of albino hairless mice, Cancer Research, vol 53,=*
*=  pp. 53-60, 1993                                                          =*
*=  (Action spectrum for carcinomas)                                         =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  W  - REAL, wavelength (nm)                                            (I)=*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      REAL w

* function value:
      REAL futr

* local:
      REAL a1, a2, a3, a4, a5,
     >     x1, x2, x3, x4, x5,
     >     t1, t2, t3, t4, t5,
     >     b1, b2, b3, b4, b5,
     >     p

      a1 = -10.91
      a2 = - 0.86
      a3 = - 8.60
      a4 = - 9.36
      a5 = -13.15

      x1 = 270.
      x2 = 302.
      x3 = 334.
      x4 = 367.
      x5 = 400.

      t1 = (w-x2)*(w-x3)*(w-x4)*(w-x5)
      t2 = (w-x1)*(w-x3)*(w-x4)*(w-x5)
      t3 = (w-x1)*(w-x2)*(w-x4)*(w-x5)
      t4 = (w-x1)*(w-x2)*(w-x3)*(w-x5)
      t5 = (w-x1)*(w-x2)*(w-x3)*(w-x4)

      b1 = (x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)
      b2 = (x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)
      b3 = (x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)
      b4 = (x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)
      b5 = (x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)

      p = a1*t1/b1 + a2*t2/b2 + a3*t3/b3 + a4*t4/b4 + a5*t5/b5

      futr  = EXP(p)

      RETURN
      END

