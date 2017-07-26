!<<<<<<<<<<<<<<<<<<fastJX codes<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!<<<<<<<<only access to external variable thru cmn_FJX.f and calls<<<<<<
!<<<<<<<<<<<<<<<<<<<<version 7.2  (6/2013, mjp) w/ RANSET for clouds<<<<

!
! !MODULE: FJX_INIT
!
! !DESCRIPTION: FJX_INIT contains variables and routines to input fast-JX data
!
!
! !INTERFACE:
!
      MODULE FJX_INIT_MOD
!
! !USES:
!
      USE CMN_FJX_MOD

      USE FJX_SUB_MOD, ONLY : EXITC

      IMPLICIT NONE
      PRIVATE

!
! !PUBLIC SUBROUTINES:
!
      PUBLIC  :: INIT_FJX

      CONTAINS

!-----------------------------------------------------------------------
      subroutine INIT_FJX (TITLEJXX,NJXU,NJXX)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in)  ::NJXU
      integer, intent(out) ::NJXX
      character*6, intent(out), dimension(NJXU) :: TITLEJXX

      integer  JXUNIT,J, RANSEED

      write(6,*) ' fast-JX ver-7.2 standalone CTM code'

      if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
        call EXITC(' INIT_JX: invalid no. wavelengths')
      endif

! Use channel 8 to read fastJX data files:
      JXUNIT  = 8

! Read in fast-J X-sections (spectral data)
      call RD_XXX(JXUNIT,'FJX_spec.dat')

! Read in cloud scattering data
      call RD_CLD(JXUNIT,'FJX_scat-cld.dat')

! Read in aerosols scattering data
      call RD_MIE(JXUNIT,'FJX_scat-aer.dat')

! Read in UMich aerosol scattering data
      call RD_UM (JXUNIT,'FJX_scat-UMa.dat')

! Read in T & O3 climatology used to fill e.g. upper layers or if O3 not calc.
      call RD_PROF(JXUNIT,'atmos_std.dat')

        NJXX = NJX
      do J = 1,NJX
        TITLEJXX(J) = TITLEJX(J)
      enddo

! Read in photolysis rates used in chemistry code and mapping onto FJX J's
!---CTM call:  read in J-values names and link to fast-JX names
      call RD_JS_JX(JXUNIT,'FJX_j2j.dat', TITLEJXX,NJXX)

!---setup the random number sequence RAN4
         RANSEED = 66
      call RANSET (NRAN_,RAN4,RANSEED)


      END SUBROUTINE INIT_FJX


!-----------------------------------------------------------------------
      subroutine RD_XXX(NUN,NAMFIL)
!-----------------------------------------------------------------------
!  Read in wavelength bins, solar fluxes, Rayleigh, T-dep X-sections.
!
!>>>>NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!>>>>NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-only:
!           WX_ = 18 should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (cmn_FXJ.f).
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) scaled
!
!-----------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (JX_spec.dat) >> j2 for fast-J2
!     NUN      Channel number for reading data file
!
!     NJX    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, JJ, K, IW, NQRD, NWWW,   LQ
      real*8  QQ2(199), TQQ2

      character*78 TITLE0
      character*6  TITLEJ2,TITLEJ3
      character*1  TSTRAT

      TQQ(:,:) = 0.d0

!----------spectral data----set for new format data------------------
!         note that X_ = max # Xsects read in
!                   NJX = # fast-JX J-values derived from this (.le. X_)

! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects

      open (NUN,FILE=NAMFIL,status='old',form='formatted')
      read (NUN,100) TITLE0

!----note that NQRD is not used any more, a read until 'endofJ' is performed
      read (NUN,101) NQRD,NWWW
         NW1 = 1
         NW2 = NWWW
      write(6,'(1x,a)') TITLE0
      write(6,'(i8)') NWWW
!----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
      read (NUN,102) (WL(IW),IW=1,NWWW)
      read (NUN,102) (FL(IW),IW=1,NWWW)
!      read (NUN,102) (FW(IW),IW=1,NWWW)
!      read (NUN,102) (FP(IW),IW=1,NWWW)
      read (NUN,102) (QRAYL(IW),IW=1,NWWW)

!---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
!---NB the O3 and q-O3-O1D are at different temperatures and cannot be combined
      read (NUN,103) TITLEJX(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,1), (QO2(IW,3),IW=1,NWWW)

      read (NUN,103) TITLEJX(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

      read (NUN,103) TITLEJX(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

      SQQ(1) = ' '
      SQQ(2) = ' '
      SQQ(3) = ' '

      LQQ(1) = 3
      LQQ(2) = 3
      LQQ(3) = 3

!---Read remaining species:  X-sections at 1-2-3 T_s
!---read in 1 to 3 X-sects per J-value (JJ)
        JJ = 3
!-- read new Xsection block
    3 continue
       read (NUN,104) TITLEJ2,TSTRAT,TQQ2,(QQ2(IW),IW=1,NWWW)
        if (TITLEJ2 .eq. 'endofJ') goto 1
!---try to add a new Xsect
    2 continue
        JJ = JJ+1
        if (JJ .gt. X_) then
        call EXITC(' RD_XXX:  X_ not large enough for Xsects read in')
        endif
         TITLEJX(JJ) = TITLEJ2
          LQQ(JJ) = 1
          SQQ(JJ) = TSTRAT
           LQ = 1
          TQQ(LQ,JJ) = TQQ2
         do IW = 1,NWWW
          QQQ(IW,LQ,JJ) = QQ2(IW)
         enddo
!---read in a maximum (at present) of 3 Xsects for TITLEJX(JJ)
       do LQ = 2,3
        read (NUN,104) TITLEJ2,TSTRAT,TQQ2,(QQ2(IW),IW=1,NWWW)
          if (TITLEJ2 .eq. 'endofJ') goto 1
          if (TITLEJ2 .ne. TITLEJX(JJ)) goto 2
!--- OK have found another Temperature or Pressure for this Xsect
          LQQ(JJ) = LQ
          TQQ(LQ,JJ) = TQQ2
         do IW = 1,NWWW
          QQQ(IW,LQ,JJ) = QQ2(IW)
         enddo
       enddo
       goto 3
    1 continue
       NJX = JJ

!---TROP-ONLY (W_ = 12 or 8) then drop the strat Xsects (labeled 'x')
      if (W_ .eq. 12 .or. W_ .eq. 8) then
        write(6,'(a)')  &
         ' >>>TROP-ONLY reduced wavelengths, drop strat X-sects'
        JJ = 3
        do J = 4,NJX
         if (SQQ(J) .ne. 'x') then
!---collapse Xsects
          JJ = JJ+1
          if (JJ .lt. J) then
             TITLEJX(JJ) = TITLEJX(J)
             LQQ(JJ) = LQQ(J)
             SQQ(JJ) = SQQ(J)
           do LQ = 1,LQQ(J)
             TQQ(LQ,JJ) = TQQ(LQ,J)
            do IW = 1,NWWW
             QQQ(IW,LQ,JJ) = QQQ(IW,LQ,J)
            enddo
           enddo
          endif
         endif
        enddo
         NJX = JJ
      endif

      do J = 1,NJX
        write(6,200) J,TITLEJX(J),SQQ(J),LQQ(J),(TQQ(I,J),I=1,LQQ(J))
      enddo

!---need to check that TQQ (= T(K) or p(hPa)) is monotonically increasing:
      do J = 1,NJX
        if (LQQ(J) .eq. 3) then
          if (TQQ(2,J) .ge. TQQ(3,J)) then
            call EXITC ('TQQ out of order')
          endif
          if (TQQ(1,J) .ge. TQQ(2,J)) then
            call EXITC ('TQQ out of order')
          endif
        endif
        if (LQQ(J) .eq. 2) then
          if (TQQ(1,J) .ge. TQQ(2,J)) then
            call EXITC ('TQQ out of order')
          endif
        endif
      enddo

!---now collapse all the wavelengths for TROP-ONLY (W_ = 12 or 8)
      if (W_ .ne. WX_) then
!---TROP-ONLY
       if (W_ .eq. 12) then
        write(6,'(a)') &
         ' >>>TROP-ONLY reduce wavelengths to 12, drop strat X-sects'
        NW2 = 12
        do IW = 1,4
          WL(IW) = WL(IW+4)
          FL(IW) = FL(IW+4)
          QRAYL(IW) = QRAYL(IW+4)
         do K = 1,3
          QO2(IW,K) = QO2(IW+4,K)
          QO3(IW,K) = QO3(IW+4,K)
          Q1D(IW,K) = Q1D(IW+4,K)
         enddo
         do J = 4,NJX
           do K = 1,LQQ(J)
            QQQ(IW,K,J) = QQQ(IW+4,K,J)
           enddo
         enddo
        enddo
        do IW = 5,12
          WL(IW) = WL(IW+6)
          FL(IW) = FL(IW+6)
          QRAYL(IW) = QRAYL(IW+6)
         do K = 1,3
          QO2(IW,K) = QO2(IW+6,K)
          QO3(IW,K) = QO3(IW+6,K)
          Q1D(IW,K) = Q1D(IW+6,K)
         enddo
         do J = 4,NJX
           do K = 1,LQQ(J)
            QQQ(IW,K,J) = QQQ(IW+6,K,J)
           enddo
         enddo
        enddo
!---TROP-QUICK  (must scale solar flux for W=5)
       elseif (W_ .eq. 8) then
        write(6,'(a)') &
         ' >>>TROP-QUICK reduce wavelengths to 8, drop strat X-sects'
        NW2 = 8
        do IW = 1,1
          WL(IW) = WL(IW+4)
          FL(IW) = FL(IW+4)  * 2.d0
          QRAYL(IW) = QRAYL(IW+4)
         do K = 1,3
          QO2(IW,K) = QO2(IW+4,K)
          QO3(IW,K) = QO3(IW+4,K)
          Q1D(IW,K) = Q1D(IW+4,K)
         enddo
         do J = 4,NJX
           do K = 1,LQQ(J)
            QQQ(IW,K,J) = QQQ(IW+4,K,J)
           enddo
         enddo
        enddo
        do IW = 2,8
          WL(IW) = WL(IW+10)
          FL(IW) = FL(IW+10)
          QRAYL(IW) = QRAYL(IW+10)
         do K = 1,3
          QO2(IW,K) = QO2(IW+10,K)
          QO3(IW,K) = QO3(IW+10,K)
          Q1D(IW,K) = Q1D(IW+10,K)
         enddo
         do J = 4,NJX
           do K = 1,LQQ(J)
            QQQ(IW,K,J) = QQQ(IW+10,K,J)
           enddo
         enddo
        enddo

       else
         call EXITC(' no. wavelengths wrong: W_ .ne. 8,12,18')
       endif
      endif

      close(NUN)

  100 format(a)
  101 format(10x,5i5)
  102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))
  103 format(a6,1x,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  104 format(a6,a1,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  200 format(1x,' x-sect:',i3,a10,a4,i5,3(3x,f6.2))
  201 format(' Number of x-sections supplied to Fast-J2: ',i3,/, &
             ' Maximum number allowed (X_) only set to: ',i3, &
             ' - increase in cmn_FJX.f')

      END SUBROUTINE RD_XXX


!-----------------------------------------------------------------------
      subroutine RD_CLD(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols/cloud scattering data set for fast-JX (ver 5.3+)
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NCC      Number of categories for cloud scattering phase functions
!     QCC      Cloud scattering phase functions
!     WCC      5 Wavelengths for supplied phase functions
!     PCC      Phase function: first 8 terms of expansion
!     RCC      Effective radius associated with cloud type
!     SCC      Single scattering albedo
!     DCC      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K
      character*78 TITLE0
      character*20 TITLAA(A_)   ! TITLAA: Title for scatering data

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(i2,a78)') NCC,TITLE0
        if (NCC .gt. C_) then
          call EXITC(' too many cld-data sets: NCC > C_')
        endif

      write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX

      read (NUN,*)
      read (NUN,*)
      do J = 1,NCC
          read (NUN,'(3x,a8,1x,2f6.3)') TITLAA(J),RCC(J),DCC(J)
        do K = 1,5
          read (NUN,'(f4.0,f7.4,f7.4,7f6.3)') &
        WCC(K,J),QCC(K,J),SCC(K,J),(PCC(I,K,J),I=2,8)
          PCC(1,K,J) = 1.d0
        enddo
      enddo

      close(NUN)

      write(6,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):' &
                   ,(WCC(K,1),K=1,5)
      write(6,*) TITLE0
      do J=1,NCC
      write(6,'(i3,1x,a8,7f8.3)') &
                    J,TITLAA(J),RCC(J),DCC(J),(QCC(K,J),K=1,5)
      enddo

      END SUBROUTINE RD_CLD


!-----------------------------------------------------------------------
      subroutine RD_MIE(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------aerosols scattering data set for fast-JX (ver 5.3+)
!  >>>>>>>>>>>>>>>>spectral data rev to J-ref ver8.5 (5/05)<<<<<<<<<<<<
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     WAA      5 Wavelengths for the supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!     DAA      density (g/cm^3)
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K
      character*78 TITLE0
      character*20 TITLAA(A_)   ! TITLAA: Title for scatering data

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(i2,a78)') NAA,TITLE0
        if (NAA .gt. A_) then
          call EXITC(' too many aerosol-data sets: NAA > A_')
        endif

      write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX

      read (NUN,*)
      read (NUN,*)
      do J = 1,NAA
          read (NUN,'(3x,a8,1x,2f6.3)') TITLAA(J),RAA(J),DAA(J)
        do K = 1,5
          read (NUN,'(f4.0,f7.4,f7.4,7f6.3)') &
        WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
          PAA(1,K,J) = 1.d0
        enddo
      enddo

      close(NUN)

      write(6,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):' &
                   ,(WAA(K,1),K=1,5)
      write(6,*) TITLE0
      do J=1,NAA
      write(6,'(i3,1x,a8,7f8.3)') &
                    J,TITLAA(J),RAA(J),DAA(J),(QAA(K,J),K=1,5)
      enddo

      END SUBROUTINE RD_MIE


!-----------------------------------------------------------------------
      subroutine RD_UM(NUN,NAMFIL)
!-----------------------------------------------------------------------
!-------UMich aerosol optical data for fast-JX (ver 6.1+)
!-----------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) :: NUN
      character(*), intent(in) ::  NAMFIL

      integer  I, J, K, L
      character*78 TITLE0
      character*20 TITLUM(33)   ! TITLUM: Title for U Michigan aerosol data set

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(a78)') TITLE0
        write(6,*) 'UMichigan Aerosols', TITLE0
      read(NUN,'(5x,10f5.0)') WMM
        write(6,'(a,10f7.1)') ' UMIchigan aerosol wavelengths:',WMM

!---33 Different UM Aerosol Types:  SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4,
!---      FF00(0%BC), FF02, ...FF14(14%BC),  BB00, BB02, ...BB30(30%BC)
      do L=1,33
          read(NUN,'(a4)') TITLUM(L)
!---21 Rel Hum:    K=1=0%, =2=5%, ... =20=95%, =21=99%
        do K=1,21
!---6 wavelengths: J=1=200nm, 2=300nm, 3=400nm, (4'=550nm) 5=600nm, 6=1000nm
!---3 optic vars:  I=1=SSAlbedo,  =2=g,  =3=k-ext
          read(NUN,'(18f9.5)')  ((UMAER(I,J,K,L),I=1,3),J=1,6)
        enddo
      enddo

      close(NUN)

        write(6,'(a)') 'collapse UM wavelengths, drop 550 nm'
          WMM(4) = WMM(5)
          WMM(5) = WMM(6)
       do L=1,33
       do K=1,21
       do I=1,3
          UMAER(I,4,K,L) = UMAER(I,5,K,L)
          UMAER(I,5,K,L) = UMAER(I,6,K,L)
       enddo
       enddo
       enddo

        write(6,'(7(i5,1x,a4))') (L,TITLUM(L), L=1,33)

      END SUBROUTINE RD_UM


!-----------------------------------------------------------------------
      subroutine RD_PROF(NJ2,NAMFIL)
!-----------------------------------------------------------------------
!  Routine to input T and O3 reference profiles 'atmos_std.dat'
!-----------------------------------------------------------------------
      implicit none

      integer, intent(in) ::  NJ2
      character(*), intent(in) ::  NAMFIL
!
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
      real*8  OFAC, OFAK

      character*78 TITLE0
!
      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
!      write(6,'(1X,A)') TITLE0
      write(6,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
        read (NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
      enddo
      close (NJ2)

!  Extend climatology to 100 km
      OFAC = exp(-2.d5/5.d5)
      do I = 32,LREF
        OFAK = OFAC**(I-31)
        do M = 1,NTMONS
        do L = 1,NTLATS
          OREF(I,L,M) = OREF(31,L,M)*OFAK
        enddo
        enddo
      enddo
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 42,LREF
        TREF(I,L,M) = TREF(41,L,M)
      enddo
      enddo
      enddo

 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')

      END SUBROUTINE RD_PROF


!-----------------------------------------------------------------------
      subroutine RD_JS_JX(NUNIT,NAMFIL,TITLEJX,NJX)
!-----------------------------------------------------------------------
!  Read 'FJX_j2j.dat' that defines mapping of fast-JX J's (TITLEJX(1:NJX))
!    onto the CTM reactions:  react# JJ, named T_REACT, uses fast-JX's JVMAP
!    including scaling factor JFACTA
!-----------------------------------------------------------------------
!---mapping variables stored in  block /jvchem/JFACTA,JIND,NRATJ,JLABEL,JVMAP
!           real*8  JFACTA(JVN_)          integer JIND(JVN_), NRATJ
!           character*50 JLABEL(JVN_)     character*6  JVMAP(JVN_)
!     JFACTA    multiplication factor for fast-JX calculated J
!     JLABEL    label(*50) of J-value used in the main chem model
!     JVMAP     label(*6) of J-value used to match with fast-JX J's
!     NRATJ     number of Photolysis reactions in CTM chemistry, derived here
!                   NRATJ must be .le. JVN_
!-----------------------------------------------------------------------
      implicit none
!
      integer, intent(in)                    ::  NUNIT, NJX
      character(*), intent(in)               ::  NAMFIL
      character*6, intent(in),dimension(NJX) :: TITLEJX
      integer   J,JJ,K
      character*120 CLINE
      character*50 T_REACT
      character*6  T_FJX
      real*8 F_FJX

! Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
! The chemistry code title describes fully the reaction (a50)
! Blank (unfilled) chemistry J's are unmapped
! The number NRATJ is the last JJ readin that is .le. JVN
!   include fractional quantum yield for the fast-JX J's

      JLABEL(:) = '------'
      JVMAP(:) = '------'
      JFACTA(:) = 0.d0

      open (NUNIT,file=NAMFIL,status='old',form='formatted')

       read (NUNIT,'(a)') CLINE
         write(6,'(a)') CLINE
      do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       if (JJ .gt. JVN_) exit
        JLABEL(JJ) = T_REACT
        JFACTA(JJ) = F_FJX
        JVMAP(JJ) = T_FJX
        NRATJ = JJ
      enddo

      close(NUNIT)

!---Zero / Set index arrays that map Jvalue(j) onto rates
      do K = 1,NRATJ
         JIND(K) = 0
       do J = 1,NJX
        if (JVMAP(K) .eq. TITLEJX(J)) then
         JIND(K) = J
        endif
       enddo
      enddo

      write(6,'(a,i4,a)')' Photochemistry Scheme with',NRATJ,' J-values'
      do K=1,NRATJ
       if (JVMAP(K) .ne. '------' ) then
        J = JIND(K)
        if (J.eq.0) then
         write(6,'(i5,a50,f6.3,a,1x,a6)') K,JLABEL(K),JFACTA(K), &
               ' no mapping onto fast-JX',JVMAP(K)
        else
         write(6,'(i5,a50,f6.3,a,i4,1x,a6)') K,JLABEL(K),JFACTA(K), &
               ' mapped to FJX:',J,TITLEJX(J)
        endif
       endif
      enddo

      END SUBROUTINE RD_JS_JX


!-----------------------------------------------------------------------
      SUBROUTINE RANSET (ND,RAN4L,ISTART)
!-----------------------------------------------------------------------
!  generates a sequence of real*4 pseudo-random numbers RAN4L(1:ND)
!     program RAN3 from Press, based on Knuth
      implicit none
      integer, parameter ::  MBIG=1000000000
      integer, parameter ::  MSEED=161803398
      integer, parameter ::  MZ=0
      real*4 , parameter ::  FAC=1.e-9
      integer,intent(in)    :: ND
      real*4, intent(out)   :: RAN4L(ND)
      integer,intent(inout) :: ISTART
      integer :: MA(55),MJ,MK,I,II,J,K,INEXT,INEXTP
!---initialization and/or fix of ISEED < 0
        MJ = MSEED - abs(ISTART)
        MJ = mod(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        do I=1,54
          II = mod(21*I,55)
          MA(II) = MK
          MK = MJ-MK
          if (MK.lt.MZ) then
            MK=MK+MBIG
          endif
          MJ = MA(II)
        enddo
        do K=1,4
         do I=1,55
           MA(I)=MA(I)-MA(1+MOD(I+30,55))
           if (MA(I) .lt. MZ) then
             MA(I) = MA(I)+MBIG
           endif
         enddo
        enddo
        INEXT = 0
        INEXTP = 31
        ISTART = 1
!---generate next ND pseudo-random numbers
      do J=1,ND
         INEXT = mod(INEXT,55) +1
         INEXTP = mod(INEXTP,55) +1
         MJ = MA(INEXT) - MA(INEXTP)
        if (MJ .lt. MZ) then
          MJ=MJ+MBIG
        endif
         MA(INEXT) = MJ
         RAN4L(J) = MJ*FAC
      enddo

      END SUBROUTINE RANSET


      END MODULE FJX_INIT_MOD
