!23456
      module idea_das_saber
!=======================================================
!ncdump -h  /scratch3/NCEPDEV/swpc/noscrub/Svetlana.Karol/SABER_NC/DAILY/SABERT_2016d025.nc
!netcdf SABERT_2016d025 {
!dimensions:
!        Npixels = 1355 ;
!        Nlev = 72 ;
!        Nlevp = 73 ;
!variables:
!        float T_obs(Npixels, Nlev) ;
!        float eT_obs(Npixels, Nlev) ;
!        float lon_obs(Npixels) ;
!        float lat_obs(Npixels) ;
!        float sza_obs(Npixels) ;
!        double time_obs(Npixels) ;
!        double time_lst(Npixels) ;
!        int day_obs ;
!        float pmb_obs(Nlev) ;
!        float pint_obs(Nlevp) ;
!=========================================================
      use das_datatypes, only : str_conlimb
      implicit NONE
!  
!
      TYPE(str_conlimb),  pointer    :: Stsab_T, Stsab_O3, Stsab_OP
!
      integer             :: day_obs, npix_obsT, npix_obsO3, npix_obsOP
      integer             ::          nz_obsT, nz_obsO3, nz_obsOP
      integer             ::   file_status

      real, dimension(:,:), allocatable :: OP_obs, eOP_obs, dens_obs
      real, dimension(:), allocatable   :: pmb_obsOP, pint_obsOP, lon_obsOP, lat_obsOP
      real, dimension(:), allocatable   :: sza_obsOP, time_obsOP, time_lstOP
!
      real, dimension(:,:), allocatable :: O3_obs, eO3_obs
      real, dimension(:), allocatable   :: pmb_obsO3, pint_obsO3, lon_obsO3, lat_obsO3
      real, dimension(:), allocatable   :: sza_obsO3, time_obsO3, time_lstO3
!
      real, dimension(:,:), allocatable :: T_obs,  eT_obs
      real, dimension(:), allocatable   :: pmb_obsT, pint_obsT, lon_obsT, lat_obsT
      real, dimension(:), allocatable   :: sza_obsT, time_obsT, time_lstT
!
      character(len=132)  :: DirTK_list='/scratch3/NCEPDEV/swpc/noscrub/Svetlana.Karol/SABER_NC/DAILY/SABERT_'
      character(len=132)  :: DirO3_list='/scratch3/NCEPDEV/swpc/noscrub/Svetlana.Karol/SABER_NC/DAILY/SABERO3_'
      character(len=132)  :: DirO1_list='/scratch3/NCEPDEV/swpc/noscrub/Svetlana.Karol/SABER_NC/DAILY/SABERO3P_'
!
      character(len=11)   :: File_day='2016d025' 
      character(len=3)    :: Fend ='.nc'
      character(len=132)  :: File_saberT   !=trim(Dir_list{TK, O3, O3p})//trim(File_day)//Fend
      character(len=132)  :: File_saberO3
      character(len=132)  :: File_saberOP

      contains
!
      subroutine mdd_2_yddd(yyyy, mm, dd, ddd)   !strymd)
!
! Make (mm, dd) => ddd of year
!
       integer :: yyyy, mm, dd
       character(len=8) strymd
!
       integer :: ddd
       character(len=4) syy 
       character(len=3) sddd
!
       integer :: mmdd, Leap
       real :: dinmm(12), Ldinmm(12), daysmm(12)
       dinmm = (/31, 28, 31,  30, 31, 30, 31, 31, 30, 31, 30, 31/)
       Ldinmm = (/31, 28, 31,  30, 31, 30, 31, 31, 30, 31, 30, 31/)

        Leap =0
        if((yyyy/4)*4.eq.yyyy) Leap=1
        if (leap.eq.0) daysmm = dinmm
        if (leap.eq.1) daysmm = Ldinmm

        if (mm.eq.1) ddd = dd
        if (mm.gt.1) then
         ddd =sum( daysmm(1:mm-1)) + dd
        endif 


!write(syy, fmt='(i4.4)') yyyy
!write(sddd, fmt='(i3.3)') ddd
!strymd = trim(syy)//'d'//trim(sddd)


        end subroutine  mdd_2_yddd
!===============================================================
      subroutine datafile_status(file, file_status)
      use netcdf
     
      implicit none
!
! check existance of datafile nc-file only
!
      integer :: file_status
      character(len=*) file
      integer :: ncid,  ierNC
      file_status = -1
!
        ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)     
        if(iernc.ne.0) then
        print * 
         print *, ' cannot find  data-file',  file
        print * 
        else 
             file_status =  IerNC
             iernc=nf90_close(ncid)             
        endif

       end subroutine datafile_status
!===============================================================
      subroutine handle_err_das( status, message )
      implicit none

      integer,   intent(in)        :: status
      character(len=*), intent(in) :: message
      print*,'** DAS handle_err status = ',status, trim(message)
      
!      call mpi_wam_quit(23999, message)

      end subroutine handle_err_das
!
      SUBROUTINE IDEA_INIT_SABER_DAY(DATE_8, Irec)
!==================================================================
!
! read entire SABER_DAY data: three types SABERO3-SABERO3P-SABERT
! dimensions:
!	Npixels = 1423 ;
!	Nlev =    72 ;
!	Nlevp =   73 
!==================================================================
!                    
      use netcdf
      use das_datatypes, only : get_conv_size, READ_CONV_LIMB
      use das_datatypes, only : str_conlimb,   alloc_conv_data
!
! input: file, 
!        namevar,
!        kz1, kz2
!
! out
      implicit NONE   
      integer, intent(in) ::   DATE_8
      integer, intent(inout) ::  Irec 

!      TYPE(str_conlimb), pointer :: Stsab_T, Stsab_O3, Stsab_OP
     
      
      integer             :: kz1, kz2
     
      character(len=8)    :: s8_date

      CALL get_saber_y4ddd(DATE_8, S8_DATE)
      print *, S8_DATE, ' STR-date '

      File_SABERT =trim(DirTK_list)//trim(S8_DATE)//trim(Fend)
      File_SABERO3=trim(DirO3_list)//trim(S8_DATE)//trim(Fend)
      File_SABEROP=trim(DirO1_list)//trim(S8_DATE)//trim(Fend)

      print *, File_SABERT, ' VAY File_SABERT '

      CALL  datafile_status(file_sabert, file_status)
      CALL  get_conv_size(file_sabert,  nz_obsT, npix_obsT)
       kz1 = 1
       kz2 = nz_obsT
      CALL  alloc_conv_data(stsab_T, nz_obsT, npix_obsT)
      call  READ_CONV_LIMB(file_sabert,  'T',    stsab_T,    kz1, kz2)
!
      print *, File_SABERO3, ' VAY File_SABERO3 '
      CALL  datafile_status(file_saberO3, file_status)
      CALL  get_conv_size(file_saberO3,  nz_obsO3, npix_obsO3)
       kz1 = 1
       kz2 = nz_obsO3
      CALL  alloc_conv_data(stsab_O3, nz_obsO3, npix_obsO3)
      call  READ_CONV_LIMB(file_saberO3,  'O3',    stsab_O3,    kz1, kz2)
!
      print *, File_SABEROP,' VAY File_SABERO3P '
      CALL  datafile_status(file_saberOP, file_status)
      CALL  get_conv_size(file_saberOP,  nz_obsOP, npix_obsOP)
       kz1 = 1
       kz2 = nz_obsOP
      CALL  alloc_conv_data(stsab_OP, nz_obsOP, npix_obsOP)
      call  READ_CONV_LIMB(file_saberOP,  'OP',    stsab_OP,    kz1, kz2)!
! 
! Read in Daily SABERT, SABER03, SABERO3P
!  Check ingested values:  
!
         print *, maxval(Stsab_T%val), minval(Stsab_T%val), ' VAY-Stsab_T%val'
         print *, maxval(Stsab_O3%val), minval(Stsab_O3%val), ' VAY-Stsab_O3%val'
         print *, maxval(Stsab_OP%val), minval(Stsab_OP%val), ' VAY-Stsab_OP%val'
!     
!
!
         irec = 1    !  flag for allocated SABER-day data
                     !  to read new day Data Structures Should be Released
      END SUBROUTINE IDEA_INIT_SABER_DAY
!
      SUBROUTINE get_saber_y4ddd(DATE_8, S8_DATE)
!
      integer, intent(in)              :: DATE_8
!
      character(len=8), intent(out)    :: s8_date
!
      character(len=4)    ::   syy
      character(len=3)    ::   sddd
      integer             :: mmdd, yyyy,mm, dd, ddd
      yyyy = int(Date_8/10000)
      mmdd =Date_8-yyyy*10000
      mm = int(mmdd)/100
      dd = mmdd -mm*100
      call mdd_2_yddd(yyyy, mm, dd, ddd)
      write(syy, fmt='(i4.4)') yyyy
      write(sddd, fmt='(i3.3)') ddd
      S8_DATE = trim(syy)//'d'//trim(sddd)
!
      END SUBROUTINE get_saber_y4ddd  
!
       end module idea_das_saber
!
       SUBROUTINE read_saber_day(CDATE8, IREC)
       use netcdf
!
      use das_datatypes, only : get_conv_size, READ_CONV_LIMB
      use das_datatypes, only : str_conlimb,   alloc_conv_data, dealloc_saber
      use idea_das_saber, only:  DirTK_list, DirO3_list,DirO1_list
      use idea_das_saber, only:  Fend, datafile_status
      use idea_das_saber, only:  get_saber_y4ddd 
      use idea_das_saber, only: day_obs, npix_obsT, npix_obsO3, npix_obsOP
      use idea_das_saber, only:      nz_obsT, nz_obsO3, nz_obsOP, file_status
      use idea_das_saber, only:   Stsab_T, Stsab_O3, Stsab_OP
       implicit none
       integer, intent(in) :: CDATE8
       integer, intent(inout) :: IREC
!
!      TYPE(str_conlimb), pointer :: Stsab_T, Stsab_O3, Stsab_OP
!      
      integer             :: kz1, kz2

      character(len=132)  :: File_saberT   !=trim(Dir_list{TK, O3, O3p})//trim(File_day)//Fend
      character(len=132)  :: File_saberO3
      character(len=132)  :: File_saberOP 
! 
      character(len=8)    :: s8_date
     
      call get_saber_y4ddd(CDATE8, S8_DATE)
!
      print *, S8_DATE, ' STR-date '
      File_SABERT =trim(DirTK_list)//trim(S8_DATE)//trim(Fend)
      File_SABERO3=trim(DirO3_list)//trim(S8_DATE)//trim(Fend)
      File_SABEROP=trim(DirO1_list)//trim(S8_DATE)//trim(Fend)

      print *, File_SABERT, ' VAY File_SABERT '

      CALL  datafile_status(file_sabert, file_status)
!      print *, file_status
      CALL  get_conv_size(file_sabert,  nz_obsT, npix_obsT)
       kz1 = 1
       kz2 = nz_obsT
      if (irec.ne.0) call dealloc_saber(stsab_T)
!      print *, 'after call dealloc_saber T '
      CALL  alloc_conv_data(stsab_T, nz_obsT, npix_obsT)
!      print *, 'after call alloc_saber T '
      call  READ_CONV_LIMB(file_sabert,  'T',    stsab_T,    kz1, kz2)
!
      print *, File_SABERO3, ' VAY File_SABERO3 '
      CALL  datafile_status(file_saberO3, file_status)
      CALL  get_conv_size(file_saberO3,  nz_obsO3, npix_obsO3)
       kz1 = 1
       kz2 = nz_obsO3
      if (irec.ne.0) call dealloc_saber(stsab_O3)
      CALL  alloc_conv_data(stsab_O3, nz_obsO3, npix_obsO3)
      call  READ_CONV_LIMB(file_saberO3,  'O3',    stsab_O3,    kz1, kz2)
!
      print *, File_SABEROP,' VAY File_SABERO3P '
      CALL  datafile_status(file_saberOP, file_status)
      CALL  get_conv_size(file_saberOP,  nz_obsOP, npix_obsOP)
       kz1 = 1
       kz2 = nz_obsOP
      if (irec.ne.0) call dealloc_saber(stsab_OP)
      CALL  alloc_conv_data(stsab_OP, nz_obsOP, npix_obsOP)
      call  READ_CONV_LIMB(file_saberOP,  'OP',    stsab_OP,    kz1, kz2)!
!     
! Read in Daily SABERT, SABER03, SABERO3P
!  Check ingested values:  
!
         print *, maxval(Stsab_T%val),  minval(Stsab_T%val), ' VAY-Stsab_T%val'
!         print *, maxval(Stsab_O3%val), minval(Stsab_O3%val), ' VAY-Stsab_O3%val'
!         print *, maxval(Stsab_OP%val), minval(Stsab_OP%val), ' VAY-Stsab_OP%val'
!     
!
        IREC = 1
!
         print *, ' VAY-done :  read_saber_day', CDATE8
       END  SUBROUTINE read_saber_day
