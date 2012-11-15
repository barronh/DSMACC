* This file contains the following subroutines, related to reading
* simple input parameters from an input file, and interactive control. 
*     rdinp
*     write1
*     readin
*     chkval
*     newval
*     gethlp
*     select
*     atrim
*=============================================================================*

      SUBROUTINE rdinp(intrct, 
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ims,    slabel, imj,    jlabel)

*-----------------------------------------------------------------------------*
*= This subroutine allows user input of some simple variables from a         =*
*= pre-exisiting file, or interactively.  It also can be used to choose      =*
*= some simple outputs.                                                      =*
*= Some help and checking of the user input values is performed, but is not  =*
*= comprehensive.                                                            =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* Input/output variables

      LOGICAL intrct
      CHARACTER*6 inpfil, outfil
      INTEGER nstr
      REAL lat, lon, tmzone
      INTEGER iyear, imonth, iday
      REAL zstart, zstop, wstart, wstop, tstart, tstop
      INTEGER nz, nwint, nt
      LOGICAL lzenit
      REAL alsurf, psurf, o3col, so2col, no2col
      REAL taucld, zbase, ztop, tauaer, ssaaer, alpha
      REAL dirsun, difdn, difup, zout, zaird, ztemp
      LOGICAL lirrad, laflux, lrates, ljvals, lmmech

      INTEGER isfix, ijfix, itfix, izfix, iwfix
      INTEGER nms, ims(ks), nmj, imj(kj)
      CHARACTER*50 slabel(ks), jlabel(kj)

* Internal

      LOGICAL ls(ks), lj(kj)
      INTEGER ns, nj
      CHARACTER*7 avar, aerror
      INTEGER i
      LOGICAL okvalu, oktype, okfind, oknew
      LOGICAL lexist

      CHARACTER*6 savfil, tmpfil
      CHARACTER*20 finame
      INTEGER nlen

* Headers and help lines

      INTEGER nhead, nhelp
      CHARACTER*70 ahline(400)

*----------------------------------------------------------

** initalize

      okvalu = .TRUE.
      oktype = .TRUE.
      okfind = .TRUE.
      oknew  = .TRUE.

      IF(.NOT. intrct) GO TO 12

      WRITE(*,100) 
 100  FORMAT(66('='))

* read headers and help file
         
      OPEN(UNIT=kin,FILE='helpin', STATUS='OLD')
      READ(kin,*) nhead, nhelp
      DO i = 1, nhelp-1
         READ(kin,110) ahline(i)
      ENDDO
 110  FORMAT(a70)
      CLOSE(kin)

      DO i = 1, nhead
         WRITE(*,*) ahline(i)
      ENDDO
      WRITE(*,100) 

      WRITE(*,*) ' Type ?? for general information'
      WRITE(*,*) ' or <enter> to continue'
      READ(*,115) avar
 115  FORMAT(A7)
      IF(avar(1:2) .EQ. '??' ) THEN
         CALL gethlp(avar, nhelp, ahline)
         PAUSE
      ENDIF

      WRITE(*,100) 

* choose input file

 10   CONTINUE
      WRITE(*,*) 'select input file'
      WRITE(*,*) ' <enter>: usrinp (if created before)'
      WRITE(*,*)' 1: defin1 (default No. 1, optimized for surface UV)' 
      WRITE(*,*)' 2: defin2 (default No. 2, optimized for photochem )' 
      WRITE(*,*)' 3: defin3 (default No. 3, optimized for master mech)' 
      WRITE(*,*)' 4: defin4 (default No. 4, sample of all outputs)'
      
      WRITE(*,*) ' file-name for others '
      READ(*,120) inpfil
 120  FORMAT(a6)
      
      IF(inpfil .EQ. ' ') inpfil = 'usrinp'
      IF(inpfil(1:1) .EQ. '1') inpfil = 'defin1'
      IF(inpfil(1:1) .EQ. '2') inpfil = 'defin2'
      IF(inpfil(1:1) .EQ. '3') inpfil = 'defin3'
      IF(inpfil(1:1) .EQ. '4') inpfil = 'defin4'
      
 12   CONTINUE
      INQUIRE(file=inpfil,exist=lexist)
      IF(.NOT. lexist) THEN
         WRITE(*,*) '****** file does not exist: ', inpfil
         GO TO 10
      ENDIF

      write(*,*)"input file = ", inpfil
      IF(inpfil .EQ. 'defin1' .OR. inpfil .EQ. 'defin2' .OR. 
     $     inpfil .EQ. 'defin3') THEN
         OPEN(UNIT=kin,FILE=inpfil,STATUS='OLD')
         write(*,*)"old input file ?"
      ELSE
         write(*,*)"unknown input file ?"
         OPEN(UNIT=kin,FILE=inpfil,STATUS='UNKNOWN')
      ENDIF   
      READ(kin,130) avar
 130  FORMAT(A7)
      IF(avar .NE. 'TUV inp') THEN
         WRITE(*,*) 'This is not a legal TUV input file'
         GO TO 10
      ENDIF

      WRITE(*,*)
      WRITE(*,*)

* read input file:
      write(*,*)"rdinp -1 = ",nstr
      CALL readin(
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ns, ls, ims, slabel, nj, lj, imj, jlabel)
      write(*,*)"rdinp -2 = ",nstr
      CLOSE (kin)

************* check and report values:

 20   CONTINUE

      CALL chkval(aerror, okvalu, oknew,
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ns, ls, ims, slabel, nj, lj, imj, jlabel)

      IF (.NOT. intrct .AND. okvalu) GO TO 99

      CALL write1(6,
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ns, ls, ims, slabel, nj, lj, imj, jlabel)

      IF(.NOT. okvalu) THEN
         WRITE(*,*) 
     $        '****** incorrect value for variable: ', aerror
         IF(.NOT. oknew) THEN
            WRITE(*,*) 'New value was put in table, re-type to confirm'
         ENDIF
         avar = aerror
         GO TO 40
      ENDIF

* confirm or change or get help

 30   CONTINUE
      WRITE(*,*) 'Type ?variable for help on a variable, or'
      WRITE(*,*) '<enter> = keep these settings, or'
      WRITE(*,*) 'Type variable name to change (lower case):'
      READ(*,300) avar
 300  FORMAT(A7)

* confirm

      IF(oktype .AND. okfind .AND. oktype .AND. avar .EQ. ' ') GO TO 50
         
* get  help

      IF(avar(1:1) .EQ. '?') THEN
         CALL gethlp(avar, nhelp, ahline)
         PAUSE
         GO TO 20
      ENDIF

* change values

 40   CONTINUE
      CALL newval(avar,okfind,oktype,
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ns, ls, ims, slabel, nj, lj, imj, jlabel)

      IF(.NOT. okfind) THEN
         WRITE(*,*) 
     $        '****** no such variable - try again (r = refresh screen)'
         GO TO 30
      ENDIF
      IF(.NOT. oktype) THEN
         WRITE(*,*) '****** improper type for variable ', avar
         GO TO 40
      ENDIF

      IF(avar .EQ. 'inpfil') GO TO 12
      GO TO 20

** check output file name

      CALL atrim(outfil, tmpfil, nlen)
c      finame = tmpfil(1:nlen)
      finame = '../'//tmpfil(1:nlen)//'.txt'
      OPEN(UNIT=20,FILE=finame,STATUS='UNKNOWN',ERR=60)
      CLOSE(20)
      GO TO 65
 60   CONTINUE
      CLOSE(20)
      WRITE(*,*) '****** improper file name', outfil
      avar = 'outfil'
      GO TO 40
 65   CONTINUE

**************** write inputs to new file:

 50   CONTINUE

      WRITE(*,*) 'save inputs to file?'
      WRITE(*,*) ' <enter> = do not save'
      WRITE(*,*) ' 1 = save to file: usrinp '
      WRITE(*,*) ' or write new file name'
      READ(*,500) savfil
 500  FORMAT(A6) 
      IF (savfil .EQ. ' ') GO TO 99
      IF (savfil(1:1) .EQ. '1') savfil = 'usrinp'

* check new file name:
* (cannot overwrite default input files: defin1, defin2, defin3 )

      CALL atrim(savfil, tmpfil, nlen)
      finame = tmpfil(1:nlen)
      IF(finame .EQ. 'defin1' .OR. finame .EQ. 'defin2' .OR.
     $     finame .EQ. 'defin3' .OR. finame .EQ. 'defin4') THEN
         WRITE(*,*) '****** Cannot overwrite default input files'
         GO TO 70
      ENDIF
      OPEN(UNIT=20,FILE=finame,STATUS='UNKNOWN',ERR=70)
      GO TO 75
 70   CONTINUE
      CLOSE(20)
      WRITE(*,*) '****** improper file name', savfil
      GO TO 50
 75   CONTINUE

      CALL write1(20,
     $     savfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ns, ls, ims, slabel, nj, lj, imj, jlabel)
      CLOSE (20)

***** done

 99   CONTINUE

C      WRITE(*,*) 'done: loading inputs'

      CALL write1(kout,
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ns, ls, ims, slabel, nj, lj, imj, jlabel)

      END
      
*=============================================================================*

      SUBROUTINE write1(iout,
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ns, ls, ims, slabel, nj, lj, imj, jlabel)

      IMPLICIT NONE
      INCLUDE 'params'

* Input/output variables

      INTEGER iout

      CHARACTER*6 inpfil, outfil
      INTEGER nstr
      REAL lat, lon, tmzone
      INTEGER iyear, imonth, iday
      REAL zstart, zstop, wstart, wstop, tstart, tstop
      INTEGER nz, nwint, nt
      LOGICAL lzenit
      REAL alsurf, psurf, o3col, so2col, no2col
      REAL taucld, zbase, ztop, tauaer, ssaaer, alpha
      REAL dirsun, difdn, difup, zout, zaird, ztemp
      LOGICAL lirrad, laflux, lrates, ljvals, lmmech

      INTEGER isfix, ijfix, itfix, izfix, iwfix
      INTEGER nms, ims(ks), nmj, imj(kj)
      CHARACTER*50 slabel(ks), jlabel(kj)
      LOGICAL ls(ks), lj(kj)
      INTEGER ns, nj

      INTEGER i

      INTEGER ninp, nout
      CHARACTER*6 inpf, outf

* right-justify 6-character file names

      CALL atrim(inpfil,inpf,ninp)
      CALL atrim(outfil,outf,nout)
      DO i = 1, 6 - ninp
         inpf = ' '//inpf
      ENDDO

      DO i = 1, 6 - nout
         outf = ' '//outf
      ENDDO

      IF(iout .EQ. 6) then
         WRITE(*,100)
         WRITE(*,1000)
         WRITE(*,105) inpf, outf, nstr
         WRITE(*,110) lat, lon, tmzone
         WRITE(*,115) iyear, imonth, iday
         WRITE(*,120) zstart, zstop, nz
         WRITE(*,125) wstart, wstop, nwint
         WRITE(*,130) tstart, tstop, nt
         WRITE(*,135) lzenit, alsurf, psurf
         WRITE(*,140) o3col, so2col, no2col
         WRITE(*,145) taucld, zbase, ztop
         WRITE(*,150) tauaer, ssaaer, alpha
         WRITE(*,155) dirsun, difdn, difup
         WRITE(*,160) zout, zaird, ztemp
         WRITE(*,165) lirrad, laflux, lmmech
         WRITE(*,170) lrates, isfix, nms
         WRITE(*,175) ljvals, ijfix, nmj
         WRITE(*,180) iwfix, itfix, izfix
         WRITE(*,1000)
      ELSE
         WRITE(iout,100) 
         WRITE(iout,1000)
         WRITE(iout,105) inpf, outf, nstr
         WRITE(iout,110) lat, lon, tmzone
         WRITE(iout,115) iyear, imonth, iday
         WRITE(iout,120) zstart, zstop, nz
         WRITE(iout,125) wstart, wstop, nwint
         WRITE(iout,130) tstart, tstop, nt
         WRITE(iout,135) lzenit, alsurf, psurf
         WRITE(iout,140) o3col, so2col, no2col
         WRITE(iout,145) taucld, zbase, ztop
         WRITE(iout,150) tauaer, ssaaer, alpha
         WRITE(iout,155) dirsun, difdn, difup
         WRITE(iout,160) zout, zaird, ztemp
         WRITE(iout,165) lirrad, laflux, lmmech
         WRITE(iout,170) lrates, isfix, nms
         WRITE(iout,175) ljvals, ijfix, nmj
         WRITE(iout,180) iwfix, itfix, izfix
         WRITE(iout,1000)
 1000    FORMAT(66('='))

         IF(iout .EQ. kout ) then
            WRITE(iout,1012)
 1012       FORMAT('==== Spectral weighting functions used:')
            DO i = 1, nms
               WRITE(iout,200) ls(ims(i)), i, slabel(ims(i))
            ENDDO
            WRITE(iout,1022)
 1022       FORMAT('==== Photolysis reactions used:')
            DO i = 1, nmj
               WRITE(iout,200) lj(imj(i)), i, jlabel(imj(i))
            ENDDO
            WRITE(iout,1000)

         ELSE

            WRITE(iout,1010)
 1010       FORMAT('===== Avaliable spectral weighting functions:')
            DO i = 1, ns
               WRITE(iout,200) ls(i), i, slabel(i)
            ENDDO
            WRITE(iout,1020)
 1020       FORMAT('===== Available photolysis reactions:')
            DO i = 1, nj
               WRITE(iout,200) lj(i), i, jlabel(i)
            ENDDO
 200        FORMAT(L1,I3,1X,A50)
            WRITE(iout,1000)
         
         ENDIF
      ENDIF

 100  FORMAT('TUV inputs:')
 105  FORMAT('inpfil = ',5x,A6,  3X,'outfil = ',5x,A6,  
     $     3X,'nstr = ',2X,I11)
 110  FORMAT('lat = ',3X,F11.3,3X,'lon = ',3X,F11.3,
     $     3X,'tmzone = ',F11.1)
 115  FORMAT('iyear = ',1X,I11,  3X,'imonth = ',I11,
     $     3X,'iday = ',2X,I11)
 120  FORMAT('zstart = ',F11.3,3X,'zstop = ',1X,F11.3,
     $     3X,'nz = ',4X,I11)
 125  FORMAT('wstart = ',F11.3,3X,'wstop = ',1X,F11.3,
     $     3X,'nwint = ',1X,I11)
 130  FORMAT('tstart = ',F11.3,3X,'tstop = ',1X,F11.3,
     $     3X,'nt = ',4X,I11)
 135  FORMAT('lzenit = ',10X,L1,3X,'alsurf = ',F11.3,
     $     3X,'psurf = ',1X,F11.1)
 140  FORMAT('o3col = ',1X,F11.3,3X,'so2col = ',F11.3,
     $     3X,'no2col = ',F11.3)
 145  FORMAT('taucld = ',F11.3,3X,'zbase = ',1X,F11.3,
     $     3X,'ztop = ',2X,F11.3)
 150  FORMAT('tauaer = ',F11.3,3X,'ssaaer = ',F11.3,
     $     3X,'alpha = ',1X,F11.3)
 155  FORMAT('dirsun = ',F11.3,3X,'difdn = ',1X,F11.3,
     $     3X,'difup = ',1X,F11.3)
 160  FORMAT('zout = ',2X,0pF11.3,3X,'zaird = ',1X,1pE11.3,
     $     3X,'ztemp = ',1X,0pF11.3)
 165  FORMAT('lirrad = ',10X,L1,  3X,'laflux = ',10X,L1,  
     $     3X,'lmmech = ',10X,L1)
 170  FORMAT('lrates = ',10X,L1,  3X,'isfix = ',1X,I11,  
     $     3X,'nms = ',3X,I11)
 175  FORMAT('ljvals = ',10X,L1,  3X,'ijfix = ',1X,I11,  
     $     3X,'nmj = ',3X,I11)
 180  FORMAT('iwfix = ',1X,I11,  3X,'itfix = ',1X,I11,
     $     3X,'izfix = ',1X,I11)

      RETURN
      END

*=============================================================================*

      SUBROUTINE readin(
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ns, ls, ims, slabel, nj, lj, imj, jlabel)

      IMPLICIT NONE
      INCLUDE 'params'

* Input/output variables

      CHARACTER*6 inpfil, outfil
      INTEGER nstr
      REAL lat, lon, tmzone
      INTEGER iyear, imonth, iday
      REAL zstart, zstop, wstart, wstop, tstart, tstop
      INTEGER nz, nwint, nt
      LOGICAL lzenit
      REAL alsurf, psurf, o3col, so2col, no2col
      REAL taucld, zbase, ztop, tauaer, ssaaer, alpha
      REAL dirsun, difdn, difup, zout, zaird, ztemp
      LOGICAL lirrad, laflux, lrates, ljvals, lmmech
      INTEGER isfix, ijfix, itfix, izfix, iwfix
      INTEGER nms, ims(ks), nmj, imj(kj)
      CHARACTER*50 slabel(ks), jlabel(kj)
      LOGICAL ls(ks), lj(kj)
      INTEGER ns, nj

* internal

      INTEGER i, idum

      READ(kin,*)
      READ(kin,105) inpfil, outfil, nstr
      write(*,*)"nstr from read = ",nstr
      READ(kin,110) lat, lon, tmzone
      READ(kin,115) iyear, imonth, iday
      READ(kin,120) zstart, zstop, nz
      READ(kin,125) wstart, wstop, nwint
      READ(kin,130) tstart, tstop, nt
      READ(kin,135) lzenit, alsurf, psurf
      READ(kin,140) o3col, so2col, no2col
      READ(kin,145) taucld, zbase, ztop
      READ(kin,150) tauaer, ssaaer, alpha
      READ(kin,155) dirsun, difdn, difup
      READ(kin,160) zout, zaird, ztemp
      READ(kin,165) lirrad, laflux, lmmech
      READ(kin,170) lrates, isfix, nms
      READ(kin,175) ljvals, ijfix, nmj
      READ(kin,180) iwfix, itfix, izfix
      READ(kin,*)

 105  FORMAT(14X,A6,  3X,14X,A6,  3X,9X,I11)
 110  FORMAT(9X,F11.3,3X,9X,F11.3,3X,9X,F11.1)
 115  FORMAT(9X,I11,  3X,9X,I11  ,3X,9X,I11)
 120  FORMAT(9X,F11.3,3X,9X,F11.3,3X,9X,I11)
 125  FORMAT(9X,F11.3,3X,9X,F11.3,3X,9X,I11)
 130  FORMAT(9X,F11.3,3X,9X,F11.3,3X,9X,I11)
 135  FORMAT(19X,L1,  3X,9X,F11.3,3X,9X,F11.1)
 140  FORMAT(9X,F11.3,3X,9X,F11.3,3X,9X,F11.3)
 145  FORMAT(9X,F11.3,3X,9X,F11.3,3X,9X,F11.3)
 150  FORMAT(9X,F11.3,3X,9X,F11.3,3X,9X,F11.3)
 155  FORMAT(9X,F11.3,3X,9X,F11.3,3X,9X,F11.3)
 160  FORMAT(9X,F11.3,3X,9X,E11.3,3X,9X,F11.3)
 165  FORMAT(19X,L1,  3X,19X,L1,  3X,19X,L1)
 170  FORMAT(19X,L1,  3X,9X,I11,  3X,9X,I11)
 175  FORMAT(19X,L1,  3X,9X,I11,  3X,9X,I11)
 180  FORMAT(9X,I11,  3X,9X,I11  ,3X,9X,I11)

* read spectral weighting function list,
* and identify those to be used:

      READ(kin,*)
      ns = 0
      nms = 0
      DO i = 1, ks
         READ(kin,200,err=10) ls(i), idum, slabel(i)
         IF (idum .NE. i) STOP 'error in spectral function list'
         ns = ns + 1
         IF(ls(i)) THEN
            nms = nms + 1
            ims(nms) = idum
         ENDIF
      ENDDO

* read photolysis reaction list,
* and identify those to be used:

 10   CONTINUE
      nj = 0
      nmj = 0
      DO i = 1, kj
         READ(kin,200,err=20) lj(i), idum, jlabel(i)
         IF (idum .NE. i) STOP 'error in photolysis reaction list'
         nj = nj  + 1
         IF(lj(i)) THEN
            nmj = nmj + 1
            imj(nmj) = idum
         ENDIF
      ENDDO

 200  FORMAT(L1,I3,1X,A50)

 20   CONTINUE
      RETURN
      END

*=============================================================================*

      SUBROUTINE chkval(aerror, okvalu, oknew,
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ns, ls, ims, slabel, nj, lj, imj, jlabel)

      IMPLICIT NONE
      INCLUDE 'params'

* Input/output variables

      LOGICAL okvalu, oknew
      CHARACTER*7 aerror

      CHARACTER*6 inpfil, outfil
      INTEGER nstr
      REAL lat, lon, tmzone
      INTEGER iyear, imonth, iday
      REAL zstart, zstop, wstart, wstop, tstart, tstop
      INTEGER nz, nwint, nt
      LOGICAL lzenit
      REAL alsurf, psurf, o3col, so2col, no2col
      REAL taucld, zbase, ztop, tauaer, ssaaer, alpha
      REAL dirsun, difdn, difup, zout, zaird, ztemp
      LOGICAL lirrad, laflux, lrates, ljvals, lmmech

      INTEGER isfix, ijfix, itfix, izfix, iwfix
      INTEGER nms, ims(ks), nmj, imj(kj)
      CHARACTER*50 slabel(ks), jlabel(kj)
      LOGICAL ls(ks), lj(kj)
      INTEGER ns, nj

* Internal

      INTEGER jday, nday
      LOGICAL oky, okm, okd
      CHARACTER*6 finame
      INTEGER nlen

* Initialize

      okvalu = .TRUE.
      oknew = .TRUE.

**** Check All Values:

* skip:
*     inpfil
*     outfil

* number of streams must be even if larger than 2:
* (if less than 2, uses 2-stream delta-Eddington)

      IF(nstr .GT. 2) THEN
         IF (2*(nstr/2) .NE. nstr .OR. nstr .GT. 32) THEN
            aerror = 'nstr'
            okvalu = .FALSE.
            RETURN
         ENDIF
      ENDIF
            
* latitude must be between -90 and +90 degrees:

      IF(ABS(lat) .GT. 90.) THEN
         aerror = 'lat'
         okvalu = .FALSE.
         RETURN
      ENDIF

* longitude must be between -180 and +180 degrees:

      IF(ABS(lon) .GT. 180.) THEN
         aerror = 'lon'
         okvalu = .FALSE.
         RETURN
      ENDIF

* time zone must be between -12 and +12 hours:

      IF(ABS(tmzone) .GT. 12.) THEN
         aerror = 'tmzone'
         okvalu = .FALSE.
         RETURN
      ENDIF

* check for legal dates:

      CALL calend(iyear, imonth, iday,
     $     jday, nday, oky, okm, okd)
      IF( .NOT. oky) THEN
         aerror = 'iyear'
         okvalu = .FALSE.
         RETURN
      ENDIF
      IF( .NOT. okm) THEN
         aerror = 'imonth'
         okvalu = .FALSE.
         RETURN
      ENDIF
      IF( .NOT. okd) THEN
         aerror = 'iday'
         okvalu = .FALSE.
         RETURN
      ENDIF

* surface elevation must be equal to or larger than zero

      IF(zstart .LT. 0.) THEN
         aerror = 'zstart'
         okvalu = .FALSE.
         RETURN
      ENDIF

* top of atmosphere should not exceed 120 km.

      IF(zstop .GT. 120.) THEN
         aerror = 'zstop'
         okvalu = .FALSE.
         RETURN
      ENDIF

* number of altitude levels must not exceed kz-1 (set in include file)

      IF(nz .gt. kz - 1) THEN
         aerror = 'nz'
         okvalu = .FALSE.
         RETURN
      ENDIF

* starting wavelength

      IF(wstart .lt. 100. .OR. wstart .GE. wstop) THEN
         aerror = 'wstart'
         okvalu = .FALSE.
         RETURN
      ENDIF
      
* ending wavelength

      IF(wstop .gt. 1000.) THEN
         aerror = 'wstop'
         okvalu = .FALSE.
         RETURN
      ENDIF

* number of wavelength intervals

      IF(nwint .EQ. 0 .OR. ABS(nwint) .GT. kw-1) THEN
         aerror = 'nwint'
         okvalu = .FALSE.
         RETURN
      ENDIF

* starting and stop time (or solar zenith angle).
* if zenith angle, must be between 0 and 180 deg.
* if time, must be between -13 and +35 hrs.

      IF(lzenit) THEN
         IF(tstart .LT. 0. .OR. tstart .GT. 180.) THEN
            aerror = 'tstart'
            okvalu = .FALSE.
            RETURN
         ENDIF
         IF(tstop .LT. tstart .OR. tstop .GT. 180.) THEN
            aerror = 'tstop'
            okvalu = .FALSE.
            RETURN
         ENDIF
      ELSE
         IF(tstart .LT. -13. .OR. tstart .GT. 36.) THEN
            aerror = 'tstart'
            okvalu = .FALSE.
            RETURN
         ENDIF
         IF(tstop .LT. tstart .OR. tstop .GT. 36.) THEN
            aerror = 'tstop'
            okvalu = .FALSE.
            RETURN
         ENDIF
      ENDIF

* number of time steps

      IF(nt .LT. 1 .OR. nt .GT. kt) THEN
         aerror = 'ntz'
         okvalu = .FALSE.
         RETURN
      ENDIF

* skip:
*     lzenith

* surface albedo

      IF(alsurf .LT. 0. .OR. alsurf .GT. 1.) THEN
         aerror = 'albnew'
         okvalu = .FALSE.
         RETURN
      ENDIF

* skip:
*     surface pressure

* ozone, SO2, and NO2 columns, cannot be negative

      IF(o3col .LT. 0.) THEN
         aerror = 'o3col'
         okvalu = .FALSE.
         RETURN
      ENDIF

      IF(so2col .LT. 0.) THEN
         aerror = 'so2col'
         okvalu = .FALSE.
         RETURN
      ENDIF

      IF(no2col .LT. 0.) THEN
         aerror = 'no2col'
         okvalu = .FALSE.
         RETURN
      ENDIF

* cloud optical depth, base, and top

      IF(taucld .LT. 0.) THEN
         aerror = 'taucld'
         okvalu = .FALSE.
         RETURN
      ENDIF

      IF(taucld .GT. 0.) THEN

         IF(zbase .LT. zstart) THEN
            aerror = 'zbase'
            okvalu = .FALSE.
            write(*,*) '****** zbase cannot be less than zstart'
            write(*,*) '****** changed to zbase = zstart'
            zbase = zstart
            oknew = .FALSE.
            RETURN
         ENDIF

         IF(ztop .LE. zbase) THEN
            aerror = 'ztop'
            okvalu = .FALSE.
            write(*,*) '****** ztop cannot be less than zbase'
            write(*,*) '****** changed to ztop = zbase + 0.01 km'
            ztop = zbase + 0.01
            oknew = .FALSE.
            RETURN
         ENDIF

      ENDIF

* aerosol optical depth, single scattering albedo, Angstrom coefficient

      IF(tauaer .LT. 0.) THEN
         aerror = 'tauaer'
         okvalu = .FALSE.
         RETURN
      ENDIF

      IF(ABS(ssaaer) .GT. 1.) THEN
         aerror = 'tssaaer'
         okvalu = .FALSE.
         RETURN
      ENDIF

      IF(alpha .LT. 0.) THEN
         aerror = 'alpha'
         okvalu = .FALSE.
         RETURN
      ENDIF

* weighting factors for direct, diffuse radiation components

      IF(ABS(dirsun) .GT. 1.) THEN
         aerror = 'dirsun'
         okvalu = .FALSE.
         RETURN
      ENDIF

      IF(ABS(difdn) .GT. 1.) THEN
         aerror = 'difdn'
         okvalu = .FALSE.
         RETURN
      ENDIF

      IF(ABS(difup) .GT. 1.) THEN
         aerror = 'difup'
         okvalu = .FALSE.
         RETURN
      ENDIF

* altitude for specific outputs.  Cannot be less than zstart or more than
*  zstop

      IF(zout .LT. zstart) THEN
         aerror = 'zout'
         okvalu = .FALSE.
         write(*,*) '****** zout cannot be less than zstart'
         write(*,*) '****** changed to zout = zstart'
         oknew = .FALSE.
         zout = zstart
         RETURN
      ENDIF
      IF(zout .GT. zstop) THEN
         aerror = 'zout'
         okvalu = .FALSE.
         write(*,*) '****** zout cannot be more than zstop'
         write(*,*) '****** changed to zout = zstop'
         zout = zstop
         RETURN
      ENDIF

* skip:
*     lirrad
*     laflux
*     lmmech
*     lrates

* index of a spectral weighting function selected for detailed output

      IF(isfix .LT. 0 .OR. isfix .GT. nms) THEN
         aerror = 'isfix'
         okvalu = .FALSE.
         RETURN
      ENDIF

* skip:
*     nms
*     ns
*     lrates
*     ls
*     slabel

* Index of a spectral weighting function selected for detailed output

      IF(ijfix .LT. 0 .OR. ijfix .GT. nmj) THEN
         aerror = 'ijfix'
         okvalu = .FALSE.
         RETURN
      ENDIF

* skip:
*     nmj
*     nj
*     ljvals
*     lj
*     jlabel

* Index for output at fixed wavelength

      IF(iwfix .lt. 0 .OR. iwfix .GT. ABS(nwint)) THEN
         aerror = 'iwfix'
         okvalu = .FALSE.
         RETURN
      ENDIF

* Index for output at fixed time

      IF(itfix .lt. 0 .OR. itfix .GT. nt) THEN
         aerror = 'itfix'
         okvalu = .FALSE.
         RETURN
      ENDIF
         
* Index for output at fixed altitude

      IF(izfix .LT. 0 .OR. izfix .GT. nz) THEN
         aerror = 'izfix'
         okvalu = .FALSE.
         RETURN
      ENDIF
      IF(izfix .GT. 0) THEN
         zout = zstart + (zstop - zstart) * 
     $        float(izfix-1)/float(nz-1)
         WRITE(*,*) '****** changed to zout = ', zout
      ENDIF

      RETURN
      END

*=============================================================================*

      SUBROUTINE  newval(avar,okfind,oktype,
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3col,  so2col, no2col,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ns, ls, ims, slabel, nj, lj, imj, jlabel)

      IMPLICIT NONE
      INCLUDE 'params'

* Input/output variables

      CHARACTER*7 avar
      LOGICAL oktype, okfind

      CHARACTER*6 inpfil, outfil
      INTEGER nstr
      REAL lat, lon, tmzone
      INTEGER iyear, imonth, iday
      REAL zstart, zstop, wstart, wstop, tstart, tstop
      INTEGER nz, nwint, nt
      LOGICAL lzenit
      REAL alsurf, psurf, o3col, so2col, no2col
      REAL taucld, zbase, ztop, tauaer, ssaaer, alpha
      REAL dirsun, difdn, difup, zout, zaird, ztemp
      LOGICAL lirrad, laflux, lrates, ljvals, lmmech

      INTEGER isfix, ijfix, itfix, izfix, iwfix
      INTEGER nms, ims(ks), nmj, imj(kj)
      CHARACTER*50 slabel(ks), jlabel(kj)
      LOGICAL ls(ks), lj(kj)
      INTEGER ns, nj

*  Internal

      INTEGER n

      oktype = .TRUE.
      okfind = .TRUE.
      n = INDEX(avar,' ')

      IF (avar(1:n-1) .EQ. 'r') RETURN

      IF (avar(1:n-1) .EQ. 'inpfil') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,100,ERR=32)inpfil
         RETURN
      ENDIF
 100  FORMAT(A6)

      IF (avar(1:n-1) .EQ. 'outfil') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,100,ERR=32) outfil
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'nstr') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) nstr
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'lat') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) lat
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'lon') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) lon
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'tmzone') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) tmzone
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'iyear') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) iyear
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'imonth') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) imonth
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'iday') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) iday
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'zstart') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) zstart
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'zstop') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) zstop
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'nz') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) nz
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'wstart') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) wstart
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'wstop') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) wstop
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'nwint') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) nwint
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'tstart') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) tstart
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'tstop') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) tstop
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'nt') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) nt
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'lzenit') THEN
         IF(lzenit) THEN
            lzenit = .FALSE.
         ELSE
            lzenit = .TRUE.
         ENDIF
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'alsurf') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) alsurf
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'psurf') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) psurf
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'o3col') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) o3col
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'so2col') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) so2col
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'no2col') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) no2col
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'taucld') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) taucld
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'zbase') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) zbase
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'ztop') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) ztop
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'tauaer') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) tauaer
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'ssaaer') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) ssaaer
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'alpha') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) alpha
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'dirsun') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) dirsun
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'difdn') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) difdn
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'difup') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) difup
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'zout') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) zout
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'zaird') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) zaird
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'ztemp') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) ztemp
         RETURN
      ENDIF


      IF (avar(1:n-1) .EQ. 'lirrad') THEN
         IF(lirrad) THEN
            lirrad = .FALSE.
         ELSE
            lirrad = .TRUE.
         ENDIF
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'laflux') THEN
         IF(laflux) THEN
            laflux = .FALSE.
         ELSE
            laflux = .TRUE.
         ENDIF
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'lmmech') THEN
         IF(lmmech) THEN
            lmmech = .FALSE.
         ELSE
            lmmech = .TRUE.
         ENDIF
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'lrates') THEN
         IF(lrates) THEN
            lrates = .FALSE.
         ELSE
            lrates = .TRUE.
         ENDIF
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'isfix') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) isfix
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'nms') THEN
         CALL select(ns,nms,ims,ls,slabel)
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'ljvals') THEN
         IF(ljvals) THEN
            ljvals = .FALSE.
         ELSE
            ljvals = .TRUE.
         ENDIF
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'ijfix') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) ijfix
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'nmj') THEN
         CALL select(nj,nmj,imj,lj,jlabel)
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'iwfix') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) iwfix
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'itfix') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) itfix
         RETURN
      ENDIF

      IF (avar(1:n-1) .EQ. 'izfix') THEN
         WRITE(*,*) 'write new value for ', avar
         READ(*,*,ERR=32) izfix
         RETURN
      ENDIF

      okfind = .FALSE.
      RETURN

 32   oktype = .FALSE.
      RETURN

      END

*=============================================================================*

      SUBROUTINE gethlp(avar, nhelp, ahline)

      IMPLICIT NONE

      INTEGER nhelp
      CHARACTER*70 ahline(400)
      CHARACTER*7 avar
      INTEGER i, j, k

      IF(avar(1:2) .EQ. '??') THEN      

         DO i = 1, nhelp
            IF(ahline(i)(1:2) .EQ. '??') THEN
               READ(ahline(i+1),*) j
               WRITE(*,2000)
               DO k = 1, j
                  WRITE(*,*) ahline(i+1+k)
               ENDDO
               WRITE(*,2000)
               RETURN
            ENDIF
         ENDDO

      ELSEIF (avar(1:1) .EQ. '?') THEN

         DO i = 1, nhelp
            IF(ahline(i)(1:1) .EQ. '?') THEN
               IF(avar(2:7) .EQ. ahline(i)(2:7)) THEN
                  READ(ahline(i+1),*) j
                  WRITE(*,2000)
                  DO k = 1, j
                     WRITE(*,*) ahline(i+1+k)
                  ENDDO
                  WRITE(*,2000)
                  RETURN
               ENDIF
            ENDIF
         ENDDO

         WRITE(*,2000)
         WRITE(*,*) '****** no help found for variable ', avar(2:7)

      ENDIF

 2000 FORMAT(70('-'))

      RETURN
      END

*=============================================================================*

      SUBROUTINE select(n,nm,im,l,label)

      IMPLICIT NONE

* Input/output

      INTEGER n, nm, im(n)
      LOGICAL l(n)
      CHARACTER*50 label(n)

*  Internal

      INTEGER i, j, ii, jj, m

      CHARACTER*6 ainp

      m = 1 + n/10
 10   CONTINUE
      DO 20 jj = 1, m

 30      WRITE(*,2000)
         DO j = 1, 10
            i = 10*(jj-1) + j
            IF (i .GT. n) GO TO 40
            WRITE(*,200) l(i), i, label(i)
         ENDDO
 40      CONTINUE

         WRITE(*,2000)
         WRITE(*,*) 
     $        'Number = toggle T/F'
         WRITE(*,*) 
     $        'n = next page'
         WRITE(*,*) 
     $        'r = return to first page'
         WRITE(*,*) 
     $        '<Enter> = done'
         
         READ(*,100) ainp
         IF(ainp .EQ. ' ') GO TO 70
         IF(INDEX(ainp,'n') .NE. 0) GO TO 20
         IF(INDEX(ainp,'r') .NE. 0) GO TO 10

         READ(ainp,110,ERR=50) ii
         IF(ii .LT. 1 .OR. ii .GT. n) GO TO 50
         GO TO 60

 50      CONTINUE
         WRITE(*,*) 'input error, try again'
         GO TO 30

 60      CONTINUE
         IF(l(ii)) THEN
            l(ii) = .FALSE.
         ELSE
            l(ii) = .TRUE.
         ENDIF
         GO TO 30

 20   CONTINUE
 70   CONTINUE

      nm = 0
      DO i = 1, n
         IF(l(i)) THEN
            nm = nm + 1
            im(nm) = i
         ENDIF
      ENDDO

 100  FORMAT(A6)
 110  FORMAT(I6)
 200  FORMAT(L1,I3,1X,A50)
 2000 FORMAT(60('-'))

      RETURN 
      END

*=============================================================================*

      SUBROUTINE atrim(a,aa,n)

* Trim blanks from character string
* Input: a
* Output:  aa, n
* Internal: i

      IMPLICIT NONE
      CHARACTER*6 a, aa
      INTEGER n, i

      aa = ' '
      n = 0
      DO i = 1, 6
         IF(a(i:i) .NE. ' ') THEN
            n = n + 1
            aa(n:n) = a(i:i)
         ENDIF
      ENDDO

      RETURN
      END

