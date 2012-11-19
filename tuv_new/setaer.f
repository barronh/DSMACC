* This file contains the following subroutine, related to setting up the
* vertical profiles of atmospheric variables
*     setaer

*=============================================================================*

      SUBROUTINE setaer(ipbl, zpbl, aod330,
     $     tau550, ssaaer, alpha,
     $     nz, z, nw, wl, 
     $     dtaer, omaer, gaer)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of aerosols, and corresponding absorption     =*
*=  optical depths, single scattering albedo, and asymmetry factor.          =*
*=  Single scattering albedo and asymmetry factor can be selected for each   =*
*=  input aerosol layer (do not have to correspond to working altitude       =*
*=  grid).  See loop 27.                                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  DTAER   - REAL, optical depth due to absorption by aerosols at each   (O)=*
*=            altitude and wavelength                                        =*
*=  OMAER   - REAL, single scattering albedo due to aerosols at each      (O)=*
*=            defined altitude and wavelength                                =*
*=  GAER    - REAL, aerosol asymmetry factor at each defined altitude and (O)=*
*=            wavelength                                                     =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER kdata
      PARAMETER(kdata=51)

* input:

      REAL wl(kw)
      REAL z(kz)
      INTEGER nz
      INTEGER nw

      REAL tau550
      REAL ssaaer, alpha

* output: (on converted grid)
      REAL dtaer(kz,kw), omaer(kz,kw), gaer(kz,kw)

* local:
      REAL zd(kdata), aer(kdata)
      REAL cd(kdata), omd(kdata), gd(kdata)
      REAL womd(kdata), wgd(kdata)

      REAL cz(kz)
      REAL omz(kz)
      REAL gz(kz)

      REAL colold

      REAL wc, wscale
      INTEGER i, iw, nd

      REAL fsum
      EXTERNAL fsum

      REAL zpbl
      INTEGER ipbl

      REAL aod330

      REAL aodw(kw), ssaw(kw)
      REAL fract(kz)

*_______________________________________________________________________

* Aerosol data from Elterman (1968)
* These are vertical optical depths per km, in 1 km
* intervals from 0 km to 50 km, at 340 nm.
* This is one option.  User can specify different data set.

      DATA aer/
     1     2.40E-01,1.06E-01,4.56E-02,1.91E-02,1.01E-02,7.63E-03,
     2     5.38E-03,5.00E-03,5.15E-03,4.94E-03,4.82E-03,4.51E-03,
     3     4.74E-03,4.37E-03,4.28E-03,4.03E-03,3.83E-03,3.78E-03,
     4     3.88E-03,3.08E-03,2.26E-03,1.64E-03,1.23E-03,9.45E-04,
     5     7.49E-04,6.30E-04,5.50E-04,4.21E-04,3.22E-04,2.48E-04,
     6     1.90E-04,1.45E-04,1.11E-04,8.51E-05,6.52E-05,5.00E-05,
     7     3.83E-05,2.93E-05,2.25E-05,1.72E-05,1.32E-05,1.01E-05,
     8     7.72E-06,5.91E-06,4.53E-06,3.46E-06,2.66E-06,2.04E-06,
     9     1.56E-06,1.19E-06,9.14E-07/
*_______________________________________________________________________


* Altitudes corresponding to Elterman profile, from bottom to top:

      WRITE(kout,*)'aerosols:  Elterman (1968) continental profile'
      nd = 51
      DO 22, i = 1, nd
         zd(i) = FLOAT(i-1)
   22 CONTINUE

* assume these are point values (at each level), so find column
* increments

      DO 27, i = 1, nd - 1
         cd(i) = (aer(i+1) + aer(i)) / 2.
         omd(i) = ssaaer
         gd(i) = .61
   27 CONTINUE

*********** end data input.

* Compute integrals and averages over grid layers:
* for g and omega, use averages weighted by optical depth

      DO 29, i = 1, nd-1
         womd(i) = omd(i) * cd(i)
         wgd(i) = gd(i) * cd(i)
   29 CONTINUE
      CALL inter3(nz,z,cz, nd,zd,cd, 1)
      CALL inter3(nz,z,omz, nd, zd,womd, 1)
      CALL inter3(nz,z,gz , nd, zd,wgd, 1)
      DO 30, i = 1, nz-1
         IF (cz(i) .GT. 0.) THEN
            omz(i) = omz(i)/cz(i)
            gz(i)  = gz(i) /cz(i)
         ELSE
            omz(i) = 1.
            gz(i) = 0.
         ENDIF
   30 CONTINUE

* old column at 340 nm
*  (minimum value is pzero = 10./largest)

      colold = MAX(fsum(nz-1,cz),pzero)

* scale with new column tau at 550 nm

      IF(tau550 .GT. nzero) THEN
         DO i = 1, nz-1
            cz(i) = cz(i) * (tau550/colold) * (550./340.)**alpha 
         ENDDO
      ENDIF

* assign at all wavelengths
* (can move wavelength loop outside if want to vary with wavelength)

      DO 50, iw = 1, nw - 1
         wc = (wl(iw)+wl(iw+1))/2.

* Elterman's data are for 340 nm, so assume optical depth scales 
* inversely with first power of wavelength.

         wscale = (340./wc)**alpha

* optical depths:

         DO 40, i = 1, nz - 1
            dtaer(i,iw) = cz(i)  * wscale
            omaer(i,iw) = omz(i)
            gaer(i,iw) = gz(i)
   40    CONTINUE
   50 CONTINUE

*! overwrite for pbl:

      IF(ipbl .GT. 0) THEN	
         write (*,*) 'pbl aerosols, aod330 = ', aod330

* create wavelength-dependent optical depth and single scattering albedo:
	
         DO iw = 1, nw-1
            wc = (wl(iw)+wl(iw+1))/2.
            aodw(iw) = aod330*(wc/330.)**(-1.0)
            IF(wc .LT. 400.) THEN
               ssaw(iw) = 0.6
            ELSE
               ssaw(iw) = 0.9
            ENDIF	
         ENDDO

* divide aod among pbl layers, overwrite Elterman profile in pbl

         DO i = 1, ipbl
            fract(i) = (z(i+1) - z(i))/zpbl
         ENDDO
	
         DO iw = 1, nw-1
            DO i = 1, ipbl 
               dtaer(i, iw) = aodw(iw) * fract(i)
               omaer(i,iw) = ssaw(iw)
            ENDDO
         ENDDO

      ENDIF
*_______________________________________________________________________

      RETURN
      END
