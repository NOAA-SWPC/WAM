PROGRAM extract_atcf

! called by plot_ellipses.pro, this should extract the desired GFS ensemble forecast 
! information and write it to idl_fcst.dat, which idl will read back in.

! NOTE: will need to change directory name of output file, possibly.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HARDWIRED CHANGES TO MAKE !!!!!
! PARAMETER (nmembers = 4)
! DATA cmem /'01','02','03','04'/
! CHARACTER*2, DIMENSION(4) :: cmem
! build_filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PARAMETER (nmembers = 4)
PARAMETER (nleads   = 29)  ! 28

TYPE ens_fcst
  CHARACTER(LEN=2)  :: basin
  INTEGER           :: id
  INTEGER           :: ibasetime
  INTEGER           :: iflead
  REAL              :: ctr_lat(nmembers)
  REAL              :: ctr_lon(nmembers)
  REAL              :: centralpressure(nmembers)
  REAL              :: windspeed(nmembers)
END TYPE

TYPE (ens_fcst), DIMENSION(nleads) :: hurricane_ens_fcst
TYPE (ens_fcst), DIMENSION(nmembers, nleads) :: GEFS

CHARACTER*120 :: infile, outfile, outfile_ens
CHARACTER*2, DIMENSION(nmembers)  :: cmem
CHARACTER*2, DIMENSION(nmembers)  :: cmem_avail
CHARACTER*2   :: stormnumber

CHARACTER*10  :: cyyyymmddhh_in
CHARACTER*2   :: cbasin_in
CHARACTER*2   :: cstormno_in
CHARACTER*80  :: rundir_in
CHARACTER*2   :: glvl_in
CHARACTER*3   :: fcst_len_in
CHARACTER*3   :: fcst_output_int_in
CHARACTER*4, DIMENSION(183) :: cmmdd
CHARACTER*4   :: cmmddt

INTEGER       :: fcst_len
INTEGER       :: fcst_output_int
INTEGER       :: ios
INTEGER       :: imemktr
INTEGER       :: nmember_out

LOGICAL iex

DATA cmem /'01','02','03','04'/

DATA cmmdd /&
 '0601','0602','0603','0604','0605','0606','0607','0608','0609','0610',&
 '0611','0612','0613','0614','0615','0616','0617','0618','0619','0620',&
 '0621','0622','0623','0624','0625','0626','0627','0628','0629','0630',&
 '0701','0702','0703','0704','0705','0706','0707','0708','0709','0710',&
 '0711','0712','0713','0714','0715','0716','0717','0718','0719','0720',&
 '0721','0722','0723','0724','0725','0726','0727','0728','0729','0730',&
 '0731',&
 '0801','0802','0803','0804','0805','0806','0807','0808','0809','0810',&
 '0811','0812','0813','0814','0815','0816','0817','0818','0819','0820',&
 '0821','0822','0823','0824','0825','0826','0827','0828','0829','0830',&
 '0831',&
 '0901','0902','0903','0904','0905','0906','0907','0908','0909','0910',&
 '0911','0912','0913','0914','0915','0916','0917','0918','0919','0920',&
 '0921','0922','0923','0924','0925','0926','0927','0928','0929','0930',&
 '1001','1002','1003','1004','1005','1006','1007','1008','1009','1010',&
 '1011','1012','1013','1014','1015','1016','1017','1018','1019','1020',&
 '1021','1022','1023','1024','1025','1026','1027','1028','1029','1030',&
 '1031',&
 '1101','1102','1103','1104','1105','1106','1107','1108','1109','1110',&
 '1111','1112','1113','1114','1115','1116','1117','1118','1119','1120',&
 '1121','1122','1123','1124','1125','1126','1127','1128','1129','1130'/

! ----------------------------------------------------------------
! input beginning and end date to process, and then find indices of
! dates to process in cmmdd array
! ----------------------------------------------------------------

CALL getarg (1, cyyyymmddhh_in)
CALL getarg (2, cstormno_in)
CALL getarg (3, cbasin_in)   
CALL getarg (4, rundir_in)   
CALL getarg (5, glvl_in)   
CALL getarg (6, fcst_len_in)   
CALL getarg (7, fcst_output_int_in)   
read (fcst_len_in,*)fcst_len  
read (fcst_output_int_in,*)fcst_output_int

! ----------------------------------------------------------------
! set up forecast structure, and write a warning to the output file
! that will be overwritten if the program later runs successfully
! to conclusion
! ----------------------------------------------------------------

! print *, 'in extract: fcst_len: ',fcst_len,' fcst_output_int: ',fcst_output_int
DO i = 1, nleads
  CALL init_ens_fcst_structure (hurricane_ens_fcst(i))
END DO

! print *, 'in extract: before outfile '
outfile = trim(rundir_in) // '/fimens/track/idl_fcst.dat'
! print *, 'in extract: outfile: ',outfile
OPEN (UNIT=1, FILE=outfile, STATUS='replace',FORM='formatted')
WRITE (1, 460) 666   ! this will indicate an error unless replaced later
CLOSE (1)

! ----------------------------------------------------------------
! read in the ensemble forecast data for this lead
! ----------------------------------------------------------------
DO ilead = 0,fcst_len,fcst_output_int   
  !PRINT *,'calling get_ens_fcstinfo, lead = ',ilead,cyyyymmddhh_in
  CALL get_ens_fcstinfo(ilead, cyyyymmddhh_in, cbasin_in, &
     cstormno_in,rundir_in, glvl_in, fcst_len_in ,hurricane_ens_fcst,iex, &
     imemktr)
  !PRINT *,'returning from get_ens_fcstinfo'
END DO    

! --------------------------------------------------------------
! write to file
! --------------------------------------------------------------

! outfile = '/lfs1/projects/rtfim/FIMYENS/FIMwfm/ensplots/idl_fcst.dat'
! PRINT *,'Program extract_atcf.x writing to ',TRIM(outfile)
OPEN (UNIT=1, FILE=outfile, STATUS='replace',FORM='formatted', IOSTAT=ios)
if (ios .ne. 0) then
   print *, 'error: ',ios,' opening ',outfile
endif

DO ilead = 0, fcst_len, fcst_output_int
  ileadidx = 1 + ilead/fcst_output_int
  WRITE (1, 460, iostat=ios) ileadidx
  if (ios .ne. 0) then
   print *, 'error: ',ios,' writing ',ileadidx
  endif
  460 FORMAT(i3)
  DO imem = 1, imemktr
    WRITE (1,461,iostat=ios) imem, hurricane_ens_fcst(ileadidx)%iflead, &
      hurricane_ens_fcst(ileadidx)%ctr_lat(imem),&
      hurricane_ens_fcst(ileadidx)%ctr_lon(imem),&
      hurricane_ens_fcst(ileadidx)%centralpressure(imem),&
      hurricane_ens_fcst(ileadidx)%windspeed(imem)
    461 FORMAT (i2,1x,i3,1x,4(f8.2,1x))
    if (ios .ne. 0) then
     print *, 'error: ',ios,' writing values'
    endif
!    print *,  'writing to :',trim(outfile),' ',ileadidx,' imem: ',imem, hurricane_ens_fcst(ileadidx)%iflead, &
!      hurricane_ens_fcst(ileadidx)%ctr_lat(imem),&
!      hurricane_ens_fcst(ileadidx)%ctr_lon(imem),&
!      hurricane_ens_fcst(ileadidx)%centralpressure(imem),&
!      hurricane_ens_fcst(ileadidx)%windspeed(imem)
    462 FORMAT (i2,1x,i2,1x,i3,1x,4(f8.2,1x))
  END DO
END DO

! PRINT *,'end program extract_atcf.x'
! PRINT *,'BEFORE END OF PROGRAM - IMEMKTR: ',imemktr
 
  WRITE (*, '(i2)') imemktr
  ! PRINT (*, '(i2)') imemktr

CONTAINS

! ===================================================================

SUBROUTINE get_ens_fcstinfo(ilead_in, cyyyymmddhh_in, cbasin_in, &
  cstormno_in, rundir_in, glvl_in, fcst_len_in, hurricane_ens_fcst, iex, &
  imemktr)

! --------------------------------------------------------------------
! read ATCF hurricane ens forecast files.
! the chosen lead in hours (ilead).  Return the forecast information 
! for all in the storms that have been tracked for this lead time.
! --------------------------------------------------------------------

INTEGER, INTENT(IN)      :: ilead_in ! forecast lead, h
CHARACTER*10, INTENT(IN) :: cyyyymmddhh_in
CHARACTER*2, INTENT(IN)  :: cbasin_in
CHARACTER*2, INTENT(IN)  :: cstormno_in
CHARACTER*80, INTENT(IN) :: rundir_in
CHARACTER*2, INTENT(IN)  :: glvl_in
CHARACTER*3, INTENT(IN)  :: fcst_len_in

TYPE (ens_fcst), INTENT(OUT), DIMENSION(29)  :: hurricane_ens_fcst
LOGICAL, INTENT(OUT)     :: iex ! did the forecast files exist?

CHARACTER*120 :: infile
CHARACTER*112 cline
CHARACTER*2 cbasin
CHARACTER*2 cstormno
CHARACTER*2, DIMENSION(4) :: cmem
CHARACTER*2, DIMENSION(4) :: cmem_avail

INTEGER imemktr, ileady

DATA cmem /'01','02','03','04'/
           

! ---------------------------------------------------------------------
! check to make sure that all ensemble forecast members were computed;
! if not all 10, don't make plot...
! ---------------------------------------------------------------------

! PRINT *,' ------ ilead = ',ilead_in  
iex = .TRUE.
imemktr = 0
DO imem = 1,nmembers
  CALL build_filename(cyyyymmddhh_in, cmem(imem), rundir_in, glvl_in, fcst_len_in, infile)
  INQUIRE (file=infile,exist=iex)
  ! PRINT *,TRIM(infile),' iex: ',iex
  IF (iex) THEN
     imemktr = imemktr + 1
     cmem_avail(imemktr) = cmem(imem)
  ENDIF
END DO
! PRINT *,imemktr, ' member forecast files found'

ileadidx = 1 + ilead/fcst_output_int
!print *,ileadidx,ilead,trim(infile)

! IF (imemktr .eq. nmembers) THEN   ! only process if all members available

  ifound = 0
    ! DO imem = 1, nmembers
    DO imem = 1, imemktr
    ! CALL build_filename(cyyyymmddhh_in, cmem(imem), infile)
    CALL build_filename(cyyyymmddhh_in, cmem_avail(imem), rundir_in, glvl_in, fcst_len_in, infile)
    ! PRINT *, '************* infile: ',TRIM(infile)
    OPEN (UNIT=1, FILE=infile, STATUS='old', FORM='formatted')
    DO
      READ (1,'(a112)', END=2000) cline
      ! PRINT *,'cline: ',cline
      READ (cline(31:33), '(i3)', ERR=234) ileady
      cstormno = cline(5:6)
      READ (cline(5:6), '(i2)') id
      cbasin = cline(1:2)

      ifoundtc = 0
      IF (cbasin .eq. cbasin_in .and. cstormno .eq. cstormno_in .and. &
      ileady .eq. ilead_in) THEN

        ! --------------------------------------------------------
        ! This member is tracking the storm.  Get the central pressure, 
        ! max wind speed, center's lat/lon for this member
        ! --------------------------------------------------------

        ifoundtc = 1
        READ (cline(36:38), '(i3)') ilat
! print *, " *******ifoundtc: ",ifoundtc, 'ilat: ',ilat
        IF (ilat .ne. 0) THEN
          IF (cline(39:39) .eq. 'N' .or. cline(39:39) .eq. 'n') THEN
            hurricane_ens_fcst(ileadidx)%ctr_lat(imem) = REAL(ilat)/10.
          ELSE
           hurricane_ens_fcst(ileadidx)%ctr_lat(imem) = - REAL(ilat)/10.
          END IF

          READ (cline(42:45), '(i4)') ilon
          IF (cline(46:46) .eq. 'E' .or. cline(46:46) .eq. 'e') THEN
            hurricane_ens_fcst(ileadidx)%ctr_lon(imem) = REAL(ilon)/10.
          ELSE
            hurricane_ens_fcst(ileadidx)%ctr_lon(imem) = 360. - REAL(ilon)/10.
          END IF

          IF (cline(54:57) .NE. ' -99') THEN
            READ (cline(54:57), '(i4)') imslp
            hurricane_ens_fcst(ileadidx)%centralpressure(imem) = REAL(imslp)
          ENDIF

          IF (cline(49:51) .NE. '***') THEN
            READ (cline(49:51), '(i3)') iwindkt  ! in knots
            hurricane_ens_fcst(ileadidx)%windspeed(imem) = REAL(iwindkt)*.514444
          ENDIF
          GOTO 2000
        END IF
      END IF
      234 CONTINUE
    END DO
    2000 CLOSE (1)
! print *, " mem: ",imem," *******windspeed: ",hurricane_ens_fcst(ileadidx)%windspeed(imem)
  END DO   ! imem = 1, 10

  hurricane_ens_fcst (ileadidx)%basin      = cline(1:2)
  hurricane_ens_fcst (ileadidx)%id         = idsv
  READ (cline(9:18), '(i10)') ibasetime
  hurricane_ens_fcst (ileadidx)%ibasetime  = ibasetime
  hurricane_ens_fcst (ileadidx)%iflead     = ilead

! ENDIF ! imemktr = 10

RETURN
END SUBROUTINE get_ens_fcstinfo
  
! ======================================================================

! SUBROUTINE build_filename(cyyyymmddhh, cmem, infile)
SUBROUTINE build_filename(cyyyymmddhh, cmem, rundir_in, glvl_in, fcst_len_in, infile)

CHARACTER*10, INTENT(IN) :: cyyyymmddhh
CHARACTER*2, INTENT(IN) :: cmem     ! member number
CHARACTER*80, INTENT(IN) :: rundir_in
CHARACTER*2, INTENT(IN)  :: glvl_in
CHARACTER*3, INTENT(IN)  :: fcst_len_in
CHARACTER*120, INTENT(OUT) :: infile

CHARACTER*90 :: cdir
LOGICAL iex

!cdir = '/Users/thamill/hfip/2010/'
!cdir = '/lfs1/projects/fim/fiorino/w21/dat/tc/adeck/esrl/2010/gfsenkf/'

! cdir = '/lfs1/projects/gfsenkf/gfsenkf_t574/' // cyyyymmddhh // &
!   '/fimens/mem0' //cmem//'/fim_C/'
! infile = TRIM(cdir)//'track.'//cyyyymmddhh//'.FIM'//cmem

cdir = TRIM(rundir_in) // '/tracker_0' // cmem // '/' // fcst_len_in
infile = TRIM(cdir)//'/track.'//cyyyymmddhh//'00.F'// TRIM(glvl_in) // '0' // cmem

!print*,'in build_filename: reading ',infile
INQUIRE (file=infile,exist=iex)
!print*,'reading iex: ',iex

RETURN
END SUBROUTINE build_filename

! =============================================================

SUBROUTINE init_ens_fcst_structure(eforecast)
TYPE (ens_fcst), INTENT(OUT) :: eforecast

eforecast%basin = 'ZZ'
eforecast%id = -999

eforecast%ibasetime            = -999
eforecast%iflead               = -999

eforecast%ctr_lat(:)           = -999.99
eforecast%ctr_lon(:)           = -999.99
eforecast%centralpressure(:)   = -999.99
eforecast%windspeed(:)         = -999.99

RETURN
END SUBROUTINE init_ens_fcst_structure

! ===============================================================

END PROGRAM extract_atcf

