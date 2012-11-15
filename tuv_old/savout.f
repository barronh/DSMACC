* This file contains the following subroutines, related to saving and writing
* some specific outputs:
*     saver1
*     saver2
*     outpt1
*     outpt2
*=============================================================================*

      SUBROUTINE saver1(it, itfix, iw, iwfix,  nz, izout,
     $           sirrad, saflux,
     $           svi_zw, svf_zw, svi_zt, svf_zt, svi_tw, svf_tw)

      IMPLICIT NONE

* Include parameter file

      INCLUDE 'params'

      INTEGER it, itfix, iw, iwfix, nz, izout
      REAL sirrad(kz,kw), saflux(kz,kw)


      REAL svf_zw(kz,kw), svf_tw(kt,kw), svf_zt(kz,kt)
      REAL svi_zw(kz,kw), svi_tw(kt,kw), svi_zt(kz,kt)

      INTEGER iz

* Save spectral irradiances and actinic fluxes into different arrays:
*  fn(iz,iw) at constant it = itfix
*  fn(iz,it) at constant iw = iwfix
*  fn(it,iw) at constant iz = izout

      IF(it .EQ. itfix) THEN
         DO iz = 1, nz
            svi_zw(iz,iw) = sirrad(iz,iw)
            svf_zw(iz,iw) = saflux(iz,iw)
         ENDDO
      ENDIF
      IF(iw .EQ. iwfix) THEN
         DO iz = 1, nz
            svi_zt(iz,it) = sirrad(iz,iw)
            svf_zt(iz,it) = saflux(iz,iw)
         ENDDO
      ENDIF
      svi_tw(it,iw) = sirrad(izout,iw)
      svf_tw(it,iw) = saflux(izout,iw)

      RETURN
      END

*=============================================================================*

      SUBROUTINE saver2(it,itfix, nz,izout, ns,isfix,ims, nj,ijfix,imj,
     $     rate, valj,
     $     svr_zs, svj_zj, svr_zt, svj_zt, svr_ts, svj_tj)


      IMPLICIT NONE

* Include parameter file

      INCLUDE 'params'

      INTEGER it, itfix, nz, izout, ns, isfix, nj, ijfix
      INTEGER ims(ks), imj(kj)

      REAL rate(ks,kz), valj(kj,kz)

      REAL svj_zj(kz,kj), svj_tj(kt,kj), svj_zt(kz,kt)
      REAL svr_zs(kz,ks), svr_ts(kt,ks), svr_zt(kz,kt)
      INTEGER iz, is, ij

* Save dose rates and j-values into arrays:
* fn(iz,is) and fn(iz,ij) at constant it = itfix
* fn(iz,it) at constant is = isfix and ij = ijfix
* fn(it,is) and fn(it,ij) at constant iz = izout

      IF(it .EQ. itfix) THEN
         DO iz = 1, nz
            DO is = 1, ns
               svr_zs(iz,is) = rate(is,iz)
            ENDDO
            DO ij = 1, nj
               svj_zj(iz,ij) = valj(ij,iz)
            ENDDO
         ENDDO
      ENDIF

      DO is = 1, ns
         IF(is .EQ. isfix) THEN
            DO iz = 1, nz
               svr_zt(iz,it) = rate(ims(is),iz)
            ENDDO
         ENDIF
      ENDDO
      
      DO ij = 1, nj
         IF(ij .EQ. ijfix) THEN
            DO iz = 1, nz
               svj_zt(iz,it) = valj(imj(ij),iz)
            ENDDO
         ENDIF
      ENDDO

      DO is = 1, ns
         svr_ts(it,is) = rate(is,izout)
      ENDDO
      DO ij = 1, nj
         svj_tj(it,ij) = valj(ij,izout)
      ENDDO

      RETURN
      END

*=============================================================================*

      SUBROUTINE outpt1( outfil, iout, 
     $     lirrad, laflux, lrates, ljvals, lmmech, lzenit,
     $     nms, ims, nmj, imj,
     $     nz, z, tlev, airden, izout,
     $     nw, wl, etf, iwfix,
     $     nt, t, sza, itfix,
     $     ns, slabel, isfix, nj, jlabel, ijfix,
     $     svj_zj, svj_tj, svj_zt,
     $     svr_zs, svr_ts, svr_zt,
     $     svf_zw, svf_tw, svf_zt,
     $     svi_zw, svi_tw, svi_zt )

      IMPLICIT NONE
      INCLUDE 'params'

* Wavelength grid:

      INTEGER nw, iw
      REAL wl(kw), wc(kw), wu(kw)
      REAL etf(kw)

* Altitude grid

      INTEGER nz, iz
      REAL z(kz)
      REAL tlev(kz), airden(kz)

* Time/sza grid

      INTEGER it, nt
      REAL t(kt), sza(kt)

* Radiation quantities

      REAL svi_zw(kz,kw), svi_tw(kt,kw), svi_zt(kz,kt)
      REAL svf_zw(kz,kw), svf_tw(kt,kw), svf_zt(kz,kt)

      INTEGER is, ns
      CHARACTER*50 slabel(ks)
      REAL svr_zs(kz,ks), svr_ts(kt,ks), svr_zt(kz,kt)

      INTEGER ij, nj
      CHARACTER*50 jlabel(kj)
      REAL svj_zj(kz,kj), svj_tj(kt,kj), svj_zt(kz,kt)

* output options

      INTEGER itfix, izout, iwfix, isfix, ijfix
      CHARACTER*20 outfil
      INTEGER iout

      LOGICAL lirrad, laflux, lrates, ljvals, lmmech, lzenit
      INTEGER i, nms, ims(ks), nmj, imj(kj)

      CHARACTER*6 finame
      INTEGER nlen

      DO iw = 1, nw - 1
         wu(iw) = wl(iw+1)
         wc(iw) = (wl(iw) + wu(iw))/2.
      ENDDO

      IF(iout .NE. 6) THEN
         CALL atrim(outfil,finame,nlen)
c         OPEN(UNIT=iout,FILE=finame(1:nlen),
c     $        STATUS='UNKNOWN')
         OPEN(UNIT=iout,
     $        FILE='../'//finame(1:nlen)//'.txt',
     $        STATUS='UNKNOWN')
      ENDIF

***** write out if looping over sza:

      IF(lzenit) THEN

* spectral irradiance:

      IF(lirrad) THEN
        WRITE(iout,100) 
 100    FORMAT('Spectral Irradiance, W m-2 nm-1')
         
         IF(itfix .GT. 0) THEN
            WRITE(iout,110) sza(itfix)
            WRITE(iout,120) 
            WRITE(iout,130) (wc(iw), iw = 1, nw - 1)
            DO iz = 1, nz
               WRITE(iout,140) z(iz), (svi_zw(iz,iw), iw = 1, nw - 1)
            ENDDO
         ENDIF
 110     FORMAT('values at sza = ',F10.3,' deg.')
 120     FORMAT('Columns: altitude (km), wavelengths (nm)')
 130     FORMAT(5X,'z, km',650(1X,F10.3))
 140     FORMAT(0PF10.4,650(1PE11.3))

         IF(izout .GT. 0) THEN
            WRITE(iout,150) z(izout)
            WRITE(iout,160) 
            WRITE(iout,170) (sza(it), it = 1, nt)
            DO iw = 1, nw - 1
               WRITE(iout,180) wc(iw), (svi_tw(it,iw), it = 1, nt)
            ENDDO
         ENDIF
 150     FORMAT('values at z = ',F10.3,' km') 
 160     FORMAT('Columns: wavelength (nm), solar zenith angles (deg.)')
 170     FORMAT(4X,'wc, nm',650(1X,F10.3))
 180     FORMAT(0PF10.4,650(1PE11.3))

         IF(iwfix .GT. 0) THEN
            WRITE(iout,190) wc(iwfix)
            WRITE(iout,200) 
            WRITE(iout,210)  (sza(it), it = 1, nt)
            DO iz = 1, nz
               WRITE(iout,220) z(iz), (svi_zt(iz,it), it = 1, nt)
            ENDDO
         ENDIF
 190     FORMAT('values at w = ',F10.4,' nm')
 200     FORMAT('Columns: altitude (km), solar zenith angles (deg.)')
 210     FORMAt(5X,'z, km',650(1X,F10.3))
 220     FORMAT(0PF10.4,650(1PE11.3))
      
         WRITE(iout,225)
 225     FORMAT(60('-'))

      ENDIF

* spectral actinic fluxes:

      IF(laflux) THEN
         WRITE(iout,230) 
 230     FORMAT('Spectral actinic flux, quanta cm-2 s-1 nm-1')

         IF(itfix .GT. 0) THEN
            WRITE(iout,240) sza(itfix)
            WRITE(iout,250) 
            WRITE(iout,260) (wc(iw), iw = 1, nw - 1)
            DO iz = 1, nz
               WRITE(iout,270) z(iz), (svf_zw(iz,iw), iw = 1, nw - 1)
            ENDDO
         ENDIF
 240     FORMAT('values at sza = ',F10.3,' deg.')
 250     FORMAT('Columns: altitude (km), wavelengths (nm)')
 260     FORMAT(5X,'z, km',650(1X,F10.3))
 270     FORMAT(0PF10.4,650(1PE11.3))

         IF(izout .GT. 0) THEN
            WRITE(iout,280) z(izout)
            WRITE(iout,290)
            WRITE(iout,300) (sza(it), it = 1, nt)
            DO iw = 1, nw - 1
               WRITE(iout,310) wc(iw), (svf_tw(it,iw), it = 1, nt)
            ENDDO
         ENDIF
 280     FORMAT('values at z = ',F10.3,' km') 
 290     FORMAT('Columns: wavelength (nm), solar zenith angles (deg.)')
 300     FORMAT(4X,'wc, nm',650(1X,F10.3))
 310     FORMAT(0PF10.4,650(1PE11.3))

         IF(iwfix .GT. 0) THEN
            WRITE(iout,320) wc(iwfix)
            WRITE(iout,330)
            WRITE(iout,340) (sza(it), it = 1, nt)
            DO iz = 1, nz
               WRITE(iout,350) z(iz), (svf_zt(iz,it), it = 1, nt)
            ENDDO
         ENDIF
 320     FORMAT('values at w = ',F10.4,' nm')
 330     FORMAT('Columns: altitude (km), solar zenith angles(deg.)')
 340     FORMAt(5X,'z, km',650(1X,F10.3))
 350     FORMAT(0PF10.4,650(1PE11.3))

         WRITE(iout,355)
 355     FORMAT(60('-'))
      ENDIF

* dose rates

      IF(lrates) THEN
         WRITE(iout,360) 
 360     FORMAT('Dose rates, W m-2')
         DO i = 1, nms
            WRITE(iout,370)  i, slabel(ims(i))
         ENDDO
 370     FORMAT(I4,' = ',A50)

         IF(itfix .GT. 0) THEN
            WRITE(iout,380) sza(itfix)
            WRITE(iout,390)
            WRITE(iout,400) (i, i = 1, nms)
            DO iz = 1, nz
               WRITE(iout,410) z(iz), (svr_zs(iz,ims(i)), i = 1, nms)
            ENDDO
         ENDIF
 380     FORMAT('values at sza = ',F10.3,' deg.')
 390     FORMAT('Columns: altitude, weighting spectra')
 400     FORMAT(5X,'z, km',650I11)
 410     FORMAT(0PF10.4,650(1PE11.3))

         IF(izout .GT. 0) THEN
            WRITE(iout,420) z(izout)
            write(iout,430)
            WRITE(iout,440) (i, i = 1, nms)
            DO it = 1, nt
               WRITE(iout,450) sza(it), (svr_ts(it,ims(i)), i = 1, nms)
            ENDDO
         ENDIF
 420     FORMAT('values at z = ',F10.3,' km') 
 430     FORMAT('Columns: sza, weighting spectra')
 440     FORMAT(1x,'sza, deg.',650I11)
 450     FORMAT(0PF10.4,650(1PE11.3))

         IF(isfix .GT. 0) THEN
            WRITE(iout,460) slabel(ims(isfix))
            WRITE(iout,470) 
            WRITE(iout,480) (sza(it), it = 1, nt)
            DO iz = 1, nz
               WRITE(iout,490) z(iz), (svr_zt(iz,it), it = 1, nt)
            ENDDO
         ENDIF
 460     FORMAT('values for ',A50)
 470     FORMAT('Columns: altitude (km), solar zenith angles (deg.)')
 480     FORMAT(5X,'z, km',1X,650(F10.3,1X))
 490     FORMAT(0PF10.4,650(1PE11.3))

         WRITE(iout,495)
 495     FORMAT(60('-'))
      ENDIF

* j-values

      IF(ljvals) THEN
         WRITE(iout,500) 
 500     FORMAT('Photolysis rate coefficients, s-1')
         DO i = 1, nmj
            WRITE(iout,510)  i, jlabel(imj(i))
         ENDDO
 510     FORMAT(I4,' = ',A50)
         
         IF(itfix .GT. 0) THEN
            WRITE(iout,520) sza(itfix)
            WRITE(iout,530) 
            WRITE(iout,540) (i, i = 1, nmj)
            DO iz = 1, nz
               WRITE(iout,550) z(iz), (svj_zj(iz,imj(i)), i = 1, nmj)
            ENDDO
         ENDIF
 520     FORMAT('values at sza = ',F10.3,' deg.')
 530     FORMAT('Columns: altitude, photo-reactions')
 540     FORMAT(5X,'z, km',650I11)
 550     FORMAT(0PF10.4,650(1PE11.3))

         IF(izout .GT. 0) THEN
            WRITE(iout,560) z(izout)
            write(iout,570) 
            WRITE(iout,580) (i, i = 1, nmj)
            DO it = 1, nt
               WRITE(iout,590) sza(it), (svj_tj(it,imj(i)), i = 1, nmj)
            ENDDO
         ENDIF
 560     FORMAT('values at z = ',F10.3,' km') 
 570     FORMAT('Columns: sza, photo-reactions')
 580     FORMAT(1x,'sza, deg.',650I11)
 590     FORMAT(0PF10.4,650(1PE11.3))

         IF(ijfix .GT. 0) THEN
            WRITE(iout,600) jlabel(imj(ijfix))
            WRITE(iout,610)
            WRITE(iout,620) (sza(it), it = 1, nt)
            DO iz = 1, nz
               WRITE(iout,630) z(iz), (svj_zt(iz,it), it = 1, nt)
            ENDDO
         ENDIF
 600     FORMAT('values for ',A50)
 610     FORMAT('Columns: altitude (km), solar zenith angles (deg.)')
 620     FORMAT(5X,'z, km',1X,650(F10.3,1X))
 630     FORMAT(0PF10.4,650(1PE11.3))

         WRITE(iout,635)
 635     FORMAT(60('-'))
      ENDIF

**** else write out if looping over time (lzenit = .false. )

      ELSE

* spectral irradiance:

      IF(lirrad) THEN
        WRITE(iout,1100) 
 1100   FORMAT('Spectral Irradiance, W m-2 nm-1')
         
         IF(itfix .GT. 0) THEN
            WRITE(iout,1110) t(itfix), sza(itfix)
            WRITE(iout,1120) 
            WRITE(iout,1130) (wc(iw), iw = 1, nw - 1)
            DO iz = 1, nz
               WRITE(iout,1140) z(iz), (svi_zw(iz,iw), iw = 1, nw - 1)
            ENDDO
         ENDIF
 1110     FORMAT('values at t = ',F10.3,' hrs.',
     $        2x, 'sza = ',F10.3)
 1120     FORMAT('Columns: altitude (km), wavelengths (nm)')
 1130     FORMAT(5X,'z, km',650(1X,F10.3))
 1140     FORMAT(0PF10.4,650(1PE11.3))

         IF(izout .GT. 0) THEN
            WRITE(iout,1150) z(izout)
            WRITE(iout,1160) 
            WRITE(iout,1170) (t(it), it = 1, nt)
            WRITE(iout,1175) (sza(it), it = 1, nt)
            DO iw = 1, nw - 1
               WRITE(iout,1180) wc(iw), (svi_tw(it,iw), it = 1, nt)
            ENDDO
         ENDIF
 1150     FORMAT('values at z = ',F10.3,' km') 
 1160     FORMAT(
     $         'Columns: wavelength (nm), times (hrs.) or sza (deg.)')
 1170     FORMAT(4X,'wc, nm',650(1X,F10.3))
 1175     FORMAT(4X,'sza = ',650(1X,F10.3))
 1180     FORMAT(0PF10.4,650(1PE11.3))

         IF(iwfix .GT. 0) THEN
            WRITE(iout,1190) wc(iwfix)
            WRITE(iout,1200) 
            WRITE(iout,1210)  (t(it), it = 1, nt)
            WRITE(iout,1215)  (sza(it), it = 1, nt)
            DO iz = 1, nz
               WRITE(iout,1220) z(iz), (svi_zt(iz,it), it = 1, nt)
            ENDDO
         ENDIF
 1190    FORMAT('values at w = ',F10.4,' nm')
 1200    FORMAT(
     $        'Columns: altitude (km), times (hrs.) or sza (deg.)')
 1210    FORMAt(5X,'z, km',650(1X,F10.3))
 1215    FORMAT(4X,'sza = ',650(1X,F10.3))
 1220    FORMAT(0PF10.4,650(1PE11.3))
      
         WRITE(iout,1225)
 1225    FORMAT(60('-'))
      ENDIF

* spectral actinic fluxes:

      IF(laflux) THEN
         WRITE(iout,1230) 
 1230     FORMAT('Spectral actinic flux, quanta cm-2 s-1 nm-1')

         IF(itfix .GT. 0) THEN
            WRITE(iout,1240) t(itfix), sza(itfix)
            WRITE(iout,1250) 
            WRITE(iout,1260) (wc(iw), iw = 1, nw - 1)
            DO iz = 1, nz
               WRITE(iout,1270) z(iz), (svf_zw(iz,iw), iw = 1, nw - 1)
            ENDDO
         ENDIF
 1240     FORMAT('values at t = ',F10.3,' hrs.',
     $        2x, 'sza = ',F10.3)
 1250     FORMAT('Columns: altitude (km), wavelengths (nm)')
 1260     FORMAT(5X,'z, km',650(1X,F10.3))
 1270     FORMAT(0PF10.4,650(1PE11.3))

         IF(izout .GT. 0) THEN
            WRITE(iout,1280) z(izout)
            WRITE(iout,1290)
            WRITE(iout,1300) (t(it), it = 1, nt)
            WRITE(iout,1305) (sza(it), it = 1, nt)
            DO iw = 1, nw - 1
               WRITE(iout,1310) wc(iw), (svf_tw(it,iw), it = 1, nt)
            ENDDO
         ENDIF
 1280     FORMAT('values at z = ',F10.3,' km') 
 1290     FORMAT(
     $         'Columns: wavelength (nm), times (hrs.) or sza (deg.)')
 1300     FORMAT(4X,'wc, nm',650(1X,F10.3))
 1305     FORMAT(4X,'sza = ',650(1x,F10.3))
 1310     FORMAT(0PF10.4,650(1PE11.3))

         IF(iwfix .GT. 0) THEN
            WRITE(iout,1320) wc(iwfix)
            WRITE(iout,1330)
            WRITE(iout,1340) (t(it), it = 1, nt)
            WRITE(iout,1345) (sza(it), it = 1, nt)
            DO iz = 1, nz
               WRITE(iout,1350) z(iz), (svf_zt(iz,it), it = 1, nt)
            ENDDO
         ENDIF
 1320    FORMAT('values at w = ',F10.4,' nm')
 1330    FORMAT(
     $        'Columns: altitude (km), times (hrs.) or sza (deg.)')
 1340    FORMAt(5X,'z, km',650(1X,F10.3))
 1345    FORMAT(4X,'sza = ',650(1X,F10.3))
 1350    FORMAT(0PF10.4,650(1PE11.3))

         WRITE(iout,1355)
 1355    FORMAT(60('-'))
      ENDIF

* dose rates

      IF(lrates) THEN
         WRITE(iout,1360) 
 1360     FORMAT('Dose rates, W m-2')
         DO i = 1, nms
            WRITE(iout,1370)  i, slabel(ims(i))
         ENDDO
 1370     FORMAT(I4,' = ',A50)

         IF(itfix .GT. 0) THEN
            WRITE(iout,1380) t(itfix), sza(itfix)
            WRITE(iout,1390)
            WRITE(iout,1400) (i, i = 1, nms)
            DO iz = 1, nz
               WRITE(iout,1410) z(iz), (svr_zs(iz,ims(i)), i = 1, nms)
            ENDDO
         ENDIF
 1380     FORMAT('values at t = ',F10.3,' hrs.',
     $        2x, 'sza = ',F10.3)
 1390     FORMAT('Columns: altitude, weighting spectra')
 1400     FORMAT(5X,'z, km',650I11)
 1410     FORMAT(0PF10.4,650(1PE11.3))

         IF(izout .GT. 0) THEN
            WRITE(iout,1420) z(izout)
            write(iout,1430)
            WRITE(iout,1440) (i, i = 1, nms)
            DO it = 1, nt
               WRITE(iout,1450) t(it), sza(it), 
     $              (svr_ts(it,ims(i)), i = 1, nms)
            ENDDO
         ENDIF
 1420     FORMAT('values at z = ',F10.3,' km') 
 1430     FORMAT('Columns: time, sza, weighting spectra')
 1440     FORMAT('time, hrs.',2x,'sza, deg.',650I11)
 1450     FORMAT(0PF10.4,1x,0PF10.3, 650(1PE11.3))

         IF(isfix .GT. 0) THEN
            WRITE(iout,1460) slabel(ims(isfix))
            WRITE(iout,1470) 
            WRITE(iout,1480) (t(it), it = 1, nt)
            write(iout,1485) (sza(it), it = 1, nt)
            DO iz = 1, nz
               WRITE(iout,1490) z(iz), (svr_zt(iz,it), it = 1, nt)
            ENDDO
         ENDIF
 1460    FORMAT('values for ',A50)
 1470    FORMAT('Columns: altitude (km), times (hrs.) or sza (deg.)')
 1480    FORMAT(5X,'z, km',1X,650(F10.3,1X))
 1485    FORMAT(4X,'sza = ',1X,650(F10.3,1x))
 1490    FORMAT(0PF10.4,650(1PE11.3))

         WRITE(iout,1495)
 1495    FORMAT(60('-'))
      ENDIF

* j-values

      IF(ljvals) THEN
         WRITE(iout,1500) 
 1500     FORMAT('Photolysis rate coefficients, s-1')
         DO i = 1, nmj
            WRITE(iout,1510)  i, jlabel(imj(i))
         ENDDO
 1510     FORMAT(I4,' = ',A50)
         
         IF(itfix .GT. 0) THEN
            WRITE(iout,1520) t(itfix), sza(itfix)
            WRITE(iout,1530) 
            WRITE(iout,1540) (i, i = 1, nmj)
            DO iz = 1, nz
               WRITE(iout,1550) z(iz), (svj_zj(iz,imj(i)), i = 1, nmj)
            ENDDO
         ENDIF
 1520     FORMAT('values at t = ',F10.3,' hrs.',
     $        2x, 'sza = ',F10.3)
 1530     FORMAT('Columns: altitude, photo-reactions')
 1540     FORMAT(5X,'z, km',650I11)
 1550     FORMAT(0PF10.4,650(1PE11.3))

         IF(izout .GT. 0) THEN
            WRITE(iout,1560) z(izout)
            write(iout,1570) 
            WRITE(iout,1580) (i, i = 1, nmj)
            DO it = 1, nt
               WRITE(iout,1590) t(it), sza(it),
     $              (svj_tj(it,imj(i)), i = 1, nmj)
            ENDDO
         ENDIF
 1560     FORMAT('values at z = ',F10.3,' km') 
 1570     FORMAT('Columns: time, sza, photo-reactions')
 1580     FORMAT('time, hrs.',2x,'sza, deg.',650I11)
 1590     FORMAT(0PF10.4,1X,0PF10.3, 650(1PE11.3))

         IF(ijfix .GT. 0) THEN
            WRITE(iout,1600) jlabel(imj(ijfix))
            WRITE(iout,1610)
            WRITE(iout,1620) (t(it), it = 1, nt)
            WRITE(iout,1625) (sza(it), it = 1, nt)
            DO iz = 1, nz
               WRITE(iout,1630) z(iz), (svj_zt(iz,it), it = 1, nt)
            ENDDO
         ENDIF
 1600     FORMAT('values for ',A50)
 1610     FORMAT('Columns: altitude (km), times (hrs.)')
 1620     FORMAT(5X,'z, km',1X,650(F10.3,1X))
 1625     FORMAT(4X,'sza = ',1X,650(F10.3,1X))
 1630     FORMAT(0PF10.4,650(1PE11.3))

         WRITE(iout,1635)
 1635    FORMAT(60('-'))
      ENDIF
      ENDIF

      IF(lmmech) THEN
         CALL outpt2(
     $        outfil, iout,
     $        izout, z, tlev, airden,
     $        nw, wl, wc, wu, svf_tw,
     $        nt, t, nj, jlabel, svj_tj)
      ENDIF
         
      RETURN
      END

*=============================================================================*

      SUBROUTINE outpt2(
     $     outfil, iout,
     $     izout, z, tlev, airden,
     $     nw, wl, wc, wu, savsaf,
     $     nt, t, nj, jlabel, savjvl)

* Output of actinic fluxes and J-values in format used by 
*   the NCAR Master Mechanism

      IMPLICIT NONE

      INCLUDE 'params'

* Wavelength grid:

      INTEGER nw, iw
      REAL wl(kw), wc(kw), wu(kw)

* Altitude grid

      REAL z(kz)

      INTEGER izout
      REAL tlev(kz), airden(kz)

* Time/sza grid

      INTEGER it, nt
      REAL t(kt)

      INTEGER ij, nj
      REAL savsaf(kt,kw), savjvl(kt,kj)
      CHARACTER*50 jlabel(kj)

* output options

      CHARACTER*20 outfil
      INTEGER iout

      EXTERNAL ftrim
      CHARACTER*6 ftrim

*  actinic flux

      WRITE(iout,100)
 100  FORMAT('Output for NCAR Master Mechanism')

      WRITE(iout,110) nw-1, nt
      WRITE(iout,120) izout, z(izout), tlev(izout), airden(izout)

      WRITE(iout,130) (3600.*t(it), it = 1, nt)
      DO iw = 1, nw - 1
         WRITE(iout,140) iw, wl(iw), wc(iw), wu(iw)
         WRITE(iout,130) (savsaf(it,iw), it = 1, nt)
      ENDDO

* j-values

      WRITE(iout,110) nt, nj
      DO ij = 1, nj
         WRITE(iout,150) jlabel(ij)
         WRITE(iout,130) (savjvl(it,ij), it = 1, nt)
      ENDDO

 110  FORMAT(I4,2x,I4)
 120  FORMAT(I4,2X,0PF10.3,2X,0PF10.3,2X,1PE11.4)
 130  FORMAT(7(1PE11.3))
 140  FORMAT(I4,3(1X,F10.3))
 150  FORMAT(A50)

      RETURN
      END

