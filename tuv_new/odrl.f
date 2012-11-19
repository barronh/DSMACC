*=============================================================================*

      SUBROUTINE odrl(nz,z,nw,wl, c,
     $     dtrl)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Compute Rayleigh optical depths as a function of altitude and wavelength =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
*=            working wavelength grid                                        =*
*=  C       - REAL, number of air molecules per cm^2 at each specified    (O)=*
*=            altitude layer                                                 =*
*=  DTRL    - REAL, Rayleigh optical depth at each specified altitude     (O)=*
*=            and each specified wavelength                                  =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'


* input: 
      INTEGER nw, nz
      REAL wl(kw), z(kz)
      REAL c(kz)

* output:
* Rayleigh optical depths

      REAL dtrl(kz,kw)

* other:

      REAL srayl, wc, wmicrn, xx 
      INTEGER iz, iw
      
*_______________________________________________________________________

* compute Rayleigh cross sections and depths:

      DO 10, iw = 1, nw - 1
         wc = (wl(iw) + wl(iw+1))/2.

* Rayleigh scattering cross section from WMO 1985 (originally from
* Nicolet, M., On the molecular scattering in the terrestrial atmosphere:
* An empirical formula for its calculation in the homoshpere, Planet.
* Space Sci., 32, 1467-1468, 1984.

         wmicrn =  wc/1.E3
         IF( wmicrn .LE. 0.55) THEN
            xx = 3.6772 + 0.389*wmicrn + 0.09426/wmicrn
         ELSE
            xx = 4. + 0.04
         ENDIF
         srayl = 4.02e-28/(wmicrn)**xx

* alternate (older) expression from
* Frohlich and Shaw, Appl.Opt. v.11, p.1773 (1980).
C     xx = 3.916 + 0.074*wmicrn + 0.050/wmicrn
C     srayl(iw) = 3.90e-28/(wmicrn)**xx

         DO 20, iz = 1, nz - 1
            dtrl(iz,iw) = c(iz)*srayl
   20    CONTINUE

   10 CONTINUE
*_______________________________________________________________________

      RETURN
      END
