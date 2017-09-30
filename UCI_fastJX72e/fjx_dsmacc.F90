!>>>>>>>>cloud-JX code (includes fractional cloud treatments) ver 7.2 (12/2013)<<<<<<<<<<<<
      subroutine DSMACC_FAST_JX(XLNG,YLAT,JDAY,PRESS_MB,TEMP_K,CLDFLAG,JLABELS,JVMAPS,JFACTORS,JVALUES,JINDS,NJ,SZA,DEBUG)
! USES:
      USE CMN_FJX_MOD
      USE FJX_SUB_MOD
      USE FJX_INIT_MOD
      USE CLD_SUB_MOD, ONLY : CLOUD_JX

! CALL sequence
!      call INIT_FJX (TITLJXX,JVN_,NJXX)
!      call ACLIM_FJX (YLAT, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1_)
!      call SOLAR_JX(PHOTAU,IDAY,YGRD,XGRD, SZA,U0,SOLF)
!      call CLOUD_JX (U0,SZA,REFLB,SOLF,FG0,LPRTJ,PPP,ZZZ,TTT, &
!             DDD,RRR,OOO,   LWP,IWP,REFFL,REFFI,     CLF,CWC,   &
!             AERSP,NDXAER,L1_,AN_,VALJXX,JVN_,   &
!             CLDFLAG,NRANDO,IRAN,L3RG,NICA,JCOUNT)

      implicit none
      integer, intent(in) :: DEBUG
      real*8, intent(in) :: XLNG, YLAT, JDAY, PRESS_MB, TEMP_K
      integer, intent(in) :: CLDFLAG
      real*8, dimension(JVN_), intent(OUT) :: JFACTORS, JVALUES
      character*6, dimension(JVN_), intent(OUT) :: JVMAPS
      character*50, dimension(JVN_), intent(OUT) :: JLABELS
      integer, dimension(JVN_), intent(out) :: JINDS
      integer, intent(out) :: NJ
      real*8, intent(inout) :: SZA
      real*8, dimension(L1_) :: ETAA,ETAB, ZOFL,RI,TI,CLDP,AER1,AER2
      real*8, dimension(L_) :: WLC,WIC
      real*8  GMTAU,PHOTAU,ALBEDO, XGRD,YGRD,PSURF, SCALEH
      real*8  CF,PMID,PDEL,ZDEL,ICWC,F1,ZKM
      integer, dimension(L1_):: NAA1,NAA2
      integer MYEAR, MONTH, IDAY, JLAT, ILON, NHRMET, IZOUT
      integer I,J,K,L,N,JP, NN
      integer NRAN, RANSEED, LTOP, NJXX
      real*8, dimension(LWEPAR) :: CLDFRW,CLDIWCW,CLDLWCW   ! WCW=Cloud Water Content (g/g)
      character*6, dimension(JVN_)   ::  TITLJXX
      real*8, dimension(21,4)    :: ERRJ, ERRJJ,ERRJ2
      integer JP04,JP09,JP11,JP15, ICLD
      integer, dimension(8)      :: JCNT
      character*11, dimension(4) ::  TITJX
      character*40, dimension(8), parameter :: TITCLD =  &
      ['clear sky - no clouds    ', &
       'avg cloud cover          ', &
       'avg cloud cover^3/2      ', &
       'ICAs - avg direct beam   ', &
       'ICAs - random N ICAs     ', &
       'QCAs - midpt of bins     ', &
       'QCAs - avg clouds in bins', &
       'ICAs - use all ICAs***   ']
      character*8, dimension(8), parameter :: TITERR =  &
      ['clearsky',&
       'avgcloud',&
       'cldf^3/2',&
       'ICA-beam',&
       'ICAs-ran',&
       'QCAs-mid',&
       'QCAs-avg',&
       'all ICAs']

!      integer, parameter                :: CLDFLDS = 128
      integer, parameter                :: CLDFLDS = 1

!--------------------key params sent to CLOUD_JX-------------------------
      real*8                     :: U0,REFLB,SOLF
      real*8,  parameter         :: FG0 = 1.10d0
      logical                    :: LPRTJ, FIRST=.true.
      real*8,  dimension(L1_+1)  :: PPP,ZZZ
      real*8,  dimension(L1_  )  :: TTT,DDD,RRR,OOO
      real*8,  dimension(L1_)    :: LWP,IWP,REFFL,REFFI
      real*8,  dimension(L1_)    :: CLF,CWC
      real*8,  dimension(L1_,AN_):: AERSP
      integer, dimension(L1_,AN_):: NDXAER
      real*8, dimension(L_,JVN_) :: VALJXX
      integer                    :: NRANDO,IRAN,L3RG
      integer                    :: NICA,JCOUNT
      SAVE FIRST
!---U0 = cos (SZA), SZA = solar zenith angle
!---REFLB = Lambertian reflectivity at the Lower Boundary
!---SOLF = solar flux factor for sun-earth distance
!---FG0 = scale for asymmetry factor to get equivalent isotropic (CLDFLAG=3 only)
!---LPRTJ = .true. = turn on internal print in both CLOUD_JX & PHOTO_JX
!--- P = edge press (hPa), Z = edge alt (m), T = layer temp (K)
!--- D = layer dens-path (# molec /cm2), O = layer O3 path (# O3 /cm2)
!--- R = layer rel.hum.(fraction)
!---LWP/IWP = Liquid/Ice water path (g/m2)
!---REFFL/REFFI = R-effective(microns) in liquid/ice cloud
!---CLF = cloud fraction (0.0 to 1.0)
!---CWC = cloud water content (g/m3) both liq & ice
!---AERSP = aerosol path (g/m2) & NDXAER = aerosol index type
!---  aerosols are dimensioned with up to AN_ different types in an ICA layer
!---L1_ = parameter, dim of profile variables, L_+1 for top (non CTM) layer
!---AN_ = parameter, dim of number of aerosols being passed
!---VALJXX = J-values from CLOUD_JX & PHOTO_JX
!---JVN_ = dim of max number of J-s reported out (in the order of fast-JX, not CTM)
!---CLDFLAG = integer index for type of cloud overlap
!---CLOUD_JX:   different cloud schemes
!---CLOUD_JX:   different cloud schemes (4:8 require max-ran overlap algorithm)
!       CLDFLAG = 1  :  Clear sky J's
!       CLDFLAG = 2  :  Averaged cloud cover
!       CLDFLAG = 3  :  cloud-fract**3/2, then average cloud cover
!       CLDFLAG = 4  :  Average direct solar beam over all ICAs, invert to get clouds
!       CLDFLAG = 5  :  Random select NRANDO ICA's from all(Independent Column Atmos.)
!       CLDFLAG = 6  :  Use all (up to 4) quadrature cloud cover QCAs (mid-pts of bin)
!       CLDFLAG = 7  :  Use all (up to 4) QCAs (average clouds within each Q-bin)
!       CLDFLAG = 8  :  Calculate J's for ALL ICAs (up to 20,000 per cell!)
!---NRANDO = number of random ICAs to do J's for (CLDFLAG=4)
!---IRAN = starting index for random number selection
!---L3RG = flag for setting max-ran overlap groups:
!---     =0   break max overlap groups at cloud fraction = 0
!---     else fixed layer (1:9, 9:last LWcloud, LWclud+1:LTOP)
!---NICA = total number of ICAs


!---fast-JX:  INIT_JX is called only once to read in & store all fast-JX data:
!              also sets up random sequence for cloud-JX
!-----------------------------------------------------------------------
      if (first) then
          call INIT_FJX (TITLJXX,JVN_,NJXX)
          first = .false.
      endif
!-----------------------------------------------------------------------
!--Set up atmosphere for a single column and time for J-values calculation
!--Nominally taken from CTM, but for standalone here is read in
      open (77,file='CTM_GrdCld.dat',status='old',err=91)
      read (77,*)
      ! BHH use MYEAR from args
      !read (77,'(i6)') MYEAR
      MYEAR = INT(JDAY/1000.)
      read (77,*)
      ! BHH use IDAY from args
      !read (77,'(i6)') IDAY
      IDAY = INT(MOD(JDAY,1000.))
      read (77,*)
      !read (77,'(i6)') MONTH
      SELECT CASE  (IDAY)
        CASE (:31)
          MONTH = 1
        CASE (32:59)
          MONTH = 2
        CASE (60:90)
          MONTH = 3
        CASE (91:120)
          MONTH = 4
        CASE (121:151)
          MONTH = 5
        CASE (152:181)
          MONTH = 6
        CASE (182:212)
          MONTH = 7
        CASE (213:243)
          MONTH = 8
        CASE (244:273)
          MONTH = 9
        CASE (274:304)
          MONTH = 10
        CASE (305:334)
          MONTH = 11
        CASE (335:)
          MONTH = 12
      END SELECT
      read (77,*)
      !read (77,'(i6)') NHRMET
      read (77,*)
      !read (77,'(2f6.1)') GMTAU,PHOTAU
      GMTAU = MOD(JDAY, 1.) * 24
      PHOTAU = GMTAU! + XLNG/15
      read (77,*)
      !read (77,'(f6.1,i4)') YLAT, JLAT
      read (77,*)
      !read (77,'(f6.1,i4)') XLNG, ILON
      read (77,*)
          YGRD = YLAT*CPI180
          XGRD = XLNG*CPI180
      read (77,'(f6.1)') PSURF
      read (77,'(f6.1)') ALBEDO
      !read(77,'(2i5)') CLDFLAG,NRANDO
      read (77,*)
      !if (CLDFLAG .le. 0) CLDFLAG = 9
      !  CLDFLAG = min(9, CLDFLAG)
      !  write (6,'(i5,1x,a,i5)') CLDFLAG, TITCLD(CLDFLAG),NRANDO
      !  write (6,'(2i5,5x,a,i5)') LPAR,LWEPAR, 'LPAR / LWEPAR', L1_
      read (77,*)
      ! BHH DSMACC LPAR include the calculated layer
      do L = 1,LPAR+1 - 1
        read (77,'(i3,1x,2f11.7,2x,f5.1,f5.2,f11.2,2(f7.3,i4))') &
                      J,ETAA(L),ETAB(L),TI(L),RI(L),ZOFL(L) &
                     ,AER1(L),NAA1(L),AER2(L),NAA2(L)
      enddo
      do L = 1,L1_
       PPP(L) = ETAA(L) + ETAB(L)*PSURF
      enddo
      ! BHH DSMACC finding output level based on pressure
      IZOUT=1
      DO L = 1,L1_
         !PRINT*,L,PPP(L),PRESS_MB
         IF (PRESS_MB > PPP(L)) THEN
            IZOUT = L
            IF (DEBUG>0) THEN
              print *,"BHH","IZ",IZOUT
              print *,"BHH","GMTAU",GMTAU
              print *,"BHH","PHOTAU",PHOTAU
            ENDIF
            EXIT
         ENDIF
      ENDDO
      ! BHH DSMACC moving all other levels up by one
      DO L=L1_,IZOUT+1,-1
         PPP(L) = PPP(L-1)
         TI(L) = TI(L-1)
         RI(L) = RI(L-1)
      ENDDO
      ! BHH DSMACC setting pressure at output point to explicit value
      PPP(IZOUT) = PRESS_MB
      ! BHH DSMACC setting temperature at output point to explicit value
      TI(IZOUT) = TEMP_K
      ! BHH DSMACC not specifying RH explicitly at this time
      RI(IZOUT) = 0.5*(RI(MAX(1, IZOUT-1))+RI(IZOUT+1))
      IF (DEBUG>1) THEN
        print*,'BHH',IZOUT
        DO L=1,L1_
           print*,'BHH',L,PPP(L),TI(L),RI(L)
        ENDDO
      ENDIF
!---ACLIM_FJX sets climatologies for O3, T, D & Z - overwrite with CTM data
      call ACLIM_FJX (YLAT, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1_)
      do L = 1,L_
       TTT(L) = TI(L)
       RRR(L) = RI(L)
      enddo
      ZZZ(1)  = 16.d5*log10(1013.25d0/PPP(1))        ! zzz in cm
      do L = 1,L_
       DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
       SCALEH      = 1.3806d-19*MASFAC*TTT(L)
       ZZZ(L+1) = ZZZ(L) -( log(PPP(L+1)/PPP(L)) * SCALEH )
      enddo
       ZZZ(L1_+1) = ZZZ(L1_) + ZZHT
      REFLB = ALBEDO
      IF (DEBUG>0) THEN
          LPRTJ = .true.
      ELSE
          LPRTJ = .false.
      ENDIF

!---following is readin for cloud data, currently has 160 atmospheres
!   from tropical T319 ECMWF atmosphere used in UCI CTM.
      JCNT(:) = 0
      do ICLD=1,CLDFLDS
        read (77,*)
       do L = LWEPAR,1,-1
        read (77,'(i3,1p,e14.5,28x,2e14.5)') &
                      J,CLDFRW(L),CLDLWCW(L),CLDIWCW(L)
        if ((CLDFRW(L) .gt. 0.001d0) .and. (CLDFLAG .ne. 1)) then
          write(6,'(A1,x,A6,x,A4,x,A4,x,i5,f10.4,x,2e10.2)')&
      'J','CLDFRW','CWCW','IWCW',J,CLDFRW(L),CLDLWCW(L),CLDIWCW(L)
        endif
       enddo !L
!---load CTM-based data on ozone and aerosols on top of climatology
!---to convert kg (STT) in grid cell w/AREA (m2) to # molecules/cm^2
!      D_O3(I,J,L) = 6.023d26*STT(I,J,L,1)/48.d0  *1.d-4 /AREAXY(I,J)
!---to convert kg (STT) in grid cell w/AREA (m2) to PATH (g/m^2)
!      P_AERSL(I,J,L) = STT(I,J,L,1)*1.d-3/AREAXY(I,J)
!---this data should be available from the CTM somewhere.
!---   ZZZ(1:L_+1) is geometric altitude (cm), approx is OK.
!---AERSP = aerosol path (g/m2) and NDXAER() must come from CTM
!    the index must match one in the std or UMich data sets (or add your own).
        AERSP(:,:)  = 0.d0
        NDXAER(:,:) = 0
        do L = 1,L_
          NDXAER(L,1) = NAA1(L)
          AERSP(L,1)  = AER1(L)
          NDXAER(L,2) = NAA2(L)
          AERSP(L,2)  = AER2(L)
        enddo !L
!---convert cloud data from our EC met fields into cloud data for cloud-JX
!---     init data = cloud fraction, and water content (g/m3) averaged over cell
!---     needs: cloud fraction, ice- and liq-water path (in cloud)
!---     and R-effective of ice and liquid clouds
        LTOP  = LWEPAR
      if (maxval(CLDFRW) .le. 0.005d0) then
        IWP(:) = 0.d0
        REFFI(:) = 0.d0
        LWP(:) = 0.d0
        REFFL(:) = 0.d0
      endif
      do L = 1,LTOP
          CF  = CLDFRW(L)
        if (CF .gt. 0.005d0) then
          CLF(L) = CF
          WLC(L) = CLDLWCW(L) / CF
          WIC(L) = CLDIWCW(L) / CF
        else
          CLF(L) = 0.d0
          WLC(L) = 0.d0
          WIC(L) = 0.d0
        endif
          CWC(L) = WIC(L)+WLC(L)
      enddo !L
!---derive R-effective for clouds:  the current UCI algorithm - use your own
      do L = 1,LTOP
!---ice clouds
        if (WIC(L) .gt. 1.d-12) then
            PDEL = PPP(L) - PPP(L+1)
            ZDEL = (ZZZ(L+1) - ZZZ(L))*0.01d0  ! m
          IWP(L) = 1000.d0*WIC(L)*PDEL*G100    ! g/m2
          ICWC =        IWP(L) / ZDEL          ! g/m3
          REFFI(L) = 164.d0 * (ICWC**0.23d0)
        else
          IWP(L) = 0.d0
          REFFI(L) = 0.d0
        endif
!---water clouds
        if (WLC(L) .gt. 1.d-12) then
            PMID = 0.5d0*(PPP(L)+PPP(L+1))
            PDEL = PPP(L) - PPP(L+1)
          F1   = 0.005d0 * (PMID - 610.d0)
          F1   = min(1.d0, max(0.d0, F1))
          LWP(L) = 1000.d0*WLC(L)*PDEL*G100     ! g/m2
          REFFL(L) = 9.6d0*F1 + 12.68d0*(1.d0-F1)
        else
          LWP(L) = 0.d0
          REFFL(L) = 0.d0
        endif
      enddo !L
!---cloud input as interpreted by fast_JX
!---extinction K(m2/g) = 3/4 * Q / [Reff(micron) * density(g/cm3)]
!           ODL = LWP(L) * 0.75d0 * 2.1d0 / REFFL(L)
!           ODI = IWP(L) * 0.75d0 * 2.0d0 / (REFFI(L) * 0.917d0)

!---fast-JX:  SOLAR_JX is called  once per grid-cell to set U0, SZA, SOLF
!--- your CTM code may have its own way of calculating and passing these quantities
!-----------------------------------------------------------------------
      call SOLAR_JX(PHOTAU,IDAY,YGRD,XGRD, SZA,U0,SOLF)
!-----------------------------------------------------------------------
      if (LPRTJ) then
        write(6,'(a,f8.3,3f8.5)')'solar zenith angle, solar-f' &
               ,SZA,SOLF,U0,REFLB
        write(6,'(a,f8.3,f8.3)') 'lat/lng',YLAT,XLNG
!        call JP_ATM0(PPP,TTT,DDD,OOO,ZZZ, L_)
      endif

!  locate the position of random number sequence based on year/day/hour
      IRAN = 13+ILON+3*JLAT+5*(MYEAR-1900)+7*IDAY + 11*nint(GMTAU)

!!!!!!!          if (ICLD .gt. 1) LPRTJ =.false.
       L3RG = PL3RG
      !BHH - control CLDFLAG from the input
      !do CLDFLAG = 1,8
      IF (DEBUG>0) write(6,*) ' call CLOUD_JX w/flag=',CLDFLAG,L3RG,NRANDO

!---special diag to find out what is wrong with QCA - mid ?????
!          if (ICLD .eq. 3) then
!             LPRTJ= .true.
!          else
!!!!             LPRTJ = .false.
!          endif
!=======================================================================

      call CLOUD_JX (U0,SZA,REFLB,SOLF,FG0,LPRTJ,PPP,ZZZ,TTT, &
             DDD,RRR,OOO,   LWP,IWP,REFFL,REFFI,     CLF,CWC,   &
             AERSP,NDXAER,L1_,AN_,VALJXX,JVN_,   &
             CLDFLAG,NRANDO,IRAN,L3RG,NICA,JCOUNT)

          JCNT(CLDFLAG) = JCNT(CLDFLAG) + JCOUNT
!=======================================================================
            if (CLDFLAG .eq. 8) then
              write(6,'(a,3i8)') ' cloud-JX v7.2: ICLD/L3RG/NICAs',ICLD,L3RG,NICA
            endif
!   4O3        PHOTON    O2        O(1D)     1.000 mapped to FJX:   3 O3(1D)
!   9NO2       PHOTON    N2        O         1.000 mapped to FJX:   9 NO2
!  11NO3       PHOTON    NO        O2        0.114 mapped to FJX:  10 NO3
!  15HNO3      PHOTON    NO2       OH        1.000 mapped to FJX:  13 HNO3
        IF (DEBUG>0) THEN
        JP04 = JIND(4)
        JP09 = JIND(9)
        JP11 = JIND(11)
        JP15 = JIND(15)
        TITJX(1) = 'J(O3>O1D)  '
        TITJX(2) = 'J(NO2)     '
        TITJX(3) = 'J(NO3>all) '
        TITJX(4) = 'J(HNO3)    '
        print *,'BHH',IZOUT,'PPP',PPP(IZOUT)
        print *,'BHH',IZOUT,'TTT',TTT(IZOUT)
        print *,'BHH',IZOUT,'RRR',RRR(IZOUT)
        print *,'BHH',IZOUT,'DDD',DDD(IZOUT)
        print *,'BHH',IZOUT,'OOO',OOO(IZOUT)
        PRINT *,'BHH',TITJX(1),VALJXX(IZOUT,JP04)
        PRINT *,'BHH',TITJX(2),VALJXX(IZOUT,JP09)
        PRINT *,'BHH',TITJX(3),VALJXX(IZOUT,JP11)
        PRINT *,'BHH',TITJX(4),VALJXX(IZOUT,JP15)
        ENDIF
      enddo ! ICLD
      JVMAPS(:) = JVMAP(:)
      JLABELS(:) = JLABEL(:)
      JFACTORS(:) = JFACTA(:)
      JVALUES(:) = VALJXX(IZOUT,:)
      JINDS(:) = JIND(:)
      NJ = NRATJ
      !BHH control CLDFLAG from INPUT menu
      !enddo   !CLDFLAG

      close (77)
      return

   91 stop 'error in opening CTM_GrdCld.dat file'
      end subroutine
