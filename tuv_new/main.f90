      PROGRAM MAIN
        IMPLICIT None
        INCLUDE 'params'
        integer*4 ijlabel_local,ntime_local,itime_local
        REAL*4 svj_tj_local(kz,kj)
        CHARACTER*50 jlabel_local(kj)
        REAL*4 sza_local(kt)
        REAL*4, parameter :: Avogadro_local = 6.02214129e+23,
     &  R_local = 8.3144621
        character*50 jlabels_local(kj)
        INTEGER*4, parameter :: iyear_local=2012, imonth_local = 7,
     &                        iday_local = 1

        REAL*4, parameter  :: ro3col_local = 260., ralt_local = 1.,
     &                        rlat_local = 45., rlon_local = 0.,
     &                        rbox_temp_local = 298.,
     &                        rbox_pressure_local = 101325.,
     &  raird_local  = Avogadro_local * rbox_pressure_local /
     &                  R_local / rbox_temp_local * 1e-6
!        sza(1) = 0.
!        sza(2) = 15.
!        sza(3) = 30.
!        sza(4) = 45.
!        sza(5) = 60.
!        sza(6) = 75.
!        sza(7) = 90.
!        sza(8) = 105.
!        sza(9) = 120.
!        sza(10) = 135.
!        sza(11) = 150.
!        sza(12) = 165.
!        sza(13) = 180.
        ntime_local = 19
!        DO itime=1,ntime
!          sza(itime) = sza(itime) * pi / 180.
!        ENDDO

        call tuv(iyear_local, imonth_local, iday_local,
     &           ro3col_local, ralt_local, rlat_local,
     &           rlon_local, rbox_temp_local, raird_local,
     &           sza_local, svj_tj_local, jlabel_local)
        write(6,'(A50,100(X,ES8.2))') 
     &'jlabel', (sza_local(itime_local),itime_local=1,ntime_local)
        DO ijlabel_local=1,kj
          write(6,'(A50,100(X,ES8.2))')
     &    jlabel_local(ijlabel_local),
     &       (svj_tj_local(itime_local,ijlabel_local),
     &        itime_local=1,ntime_local)
        ENDDO
      END PROGRAM