      module  das_datatypes
!
      use IDEA_MPI_def,      only : mpi_WAM_quit, mpi_id
      use netcdf
!
      implicit NONE

       integer, parameter :: nzd_max=150
       real, parameter    :: nodata=-9999.
!
       type Hsxy
        integer :: nph
        real  ::  Wxy(4)
        integer :: Ix(4), Iy(4)
       end type Hsxy

       type Hzdir
        integer :: npz  ! 2-point inter-n
        integer :: nzd  ! actual # of vert. levels
                       ! 100 - max # of vert. levels
         real, dimension(2,nzd_max)      :: Wz !  Wz(npz, nz)
         integer,  dimension(2,nzd_max) ::  Iz !  Iz(npz, nz)
       end type Hzdir
!
       TYPE str_conlimb
            character(len=12)  :: nameSens
            character(len=6)   :: nameVar
            integer :: nz
            integer :: npix
        real, pointer ::  zkm(:)
        real, pointer ::  pmb(:)
        real, pointer :: dzkm(:), dpmb(:)
        real, pointer ::  time(:)
        real, pointer ::  lon(:)
        real, pointer ::  lat(:)
        real, pointer ::  val(:,:)
        real, pointer ::  err(:,:)
        real, pointer ::  Lcor(:)
        real, pointer ::  Ak(:,:)       ! single kernel nlev x nlev
        real, pointer ::  vapr(:)      ! a priori  ~ 5% deviation from the VAL or last iteration
!
! extra HOBS ..... errF ....OmTF
!
        real, pointer ::  omtf(:,:)                 ! data space (nzd, nhd)
        real, pointer ::  varf(:,:)                 ! (nzd, nhd)
!
! obs oper-rs
!    
       type(Hsxy),  allocatable, dimension(:) :: Hs
       type(Hzdir), allocatable, dimension(:) :: Hz 
!
! output-back (it1:it2)
! and marking "bad" data
!
        integer          :: it1, it2
        integer, pointer ::  gmark(:)  ! 1-for accepted data 0-mark for "non-analyzed" 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   48 records per day
!        damping only "obsdt"
!               or 
!      collecting them in 'obsday'
!   extra-diagnsotics at (nzfd, npix)
!
! integer           ::   nzfd  !.... or ... plev
!        real, pointer ::   Fval(:,:)
!        real, pointer ::   Aval(:,:)
!        real, pointer ::   Aerr(:,:)  
!        real, pointer ::   FPs(:)
!        real, pointer ::   Chi2(:)       
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      END TYPE str_conlimb
!
      TYPE var_data_str
        integer   :: maxd = 10      ! maximum of data-sources
        character(len=12)  :: varName
        integer  :: nobs_src
        integer  :: ymd 
! character(len=256), pointer :: dins(:)
! character(len=256), pointer :: dpaths(:)
! character(len=256), pointer :: dprefix(:)
! logical,  pointer :: file_exist(:)
!
! 4-possible types
       integer  ::  nst1  ! nst_conlimb
       integer  ::  nst2  ! nst_nadcol
       integer  ::  nst3  ! nst_nadpak
       integer  ::  nst4  ! nst_nadrad
       integer  ::  nstot ! sum(1-4)
!
       type(str_conlimb) , allocatable, dimension(:)   ::  conlimb
! type(str_nadircol), allocatable, dimension(:)   ::  nad_col
! type(str_nadirpak), allocatable, dimension(:)   ::  nad_pak
! type(str_nadrad),   allocatable, dimension(:)   ::  nad_rad

       END TYPE    var_data_str


       type obs_day
          integer :: ymd
          integer :: nvar
         type(var_data_str), allocatable  :: D(:)
       end type obs_day
!
!
!

      contains
!
      subroutine  alloc_conv_data(C, Nz, Npix)
!
      integer Nz, Npix
     
      type(str_conlimb),pointer :: C
      integer K
      write(6,*)' Allocating  npix, nz for limb-sensor '
    
     
!
! Allocate data-arrays
!  
       allocate(C)
       write(6,*) ' npix, nz', npix, nz, 'alloc_conlimb '

       C%npix= npix
       C%nz= nz

       if (npix.ge.1) then 

        allocate(C%lon(npix), C%lat(npix), C%time(npix))
        allocate(C%zkm(nz), C%pmb(nz), C%Lcor(nz))
        allocate(C%dzkm(nz), C%dpmb(nz) )
        allocate(C%val(nz, npix),  C%ERR(nz, npix))
        allocate(C%OmTF(nz, npix), C%VarF(nz, npix))
        allocate(C%AK(nz,nz), C%VAPR(nz) )

        allocate(C%gmark(npix))

!       allocate(C%Hs(npix), C%Hz(npix))
       endif


!
       write(6,*)' Allocate  alloc_conlimb '
       end subroutine  alloc_conv_data

!
       SUBROUTINE dealloc_saber(C)
!
       type(str_conlimb),   pointer ::  C
!       integer N
!       integer K, istat
!       istat=0   
!        deallocate(C)                  
        deallocate(C%lon, C%lat, C%time)
        deallocate(C%zkm, C%pmb, C%Lcor)
        deallocate(C%dzkm, C%dpmb )
        deallocate(C%val,  C%ERR)
        deallocate(C%OmTF, C%VarF)
        deallocate(C%AK, C%VAPR )
        deallocate(C%gmark)
         C%npix= 0
         C%nz= 0      

!       do k=1, N
!        if(associated(C%Val) deallocate(C%Val)
!        if(allocated(C%ERR)) deallocate(C%ERR)
!        if(allocated(C%LON)) deallocate(C%LON)
!        if(allocated(C%LAT)) deallocate(C%LAT)
!        if(allocated(C%Time)) deallocate(C%Time)

!        if(allocated(C%gmark)) deallocate(C%gmark)

!        if(allocated(C%PMB)) deallocate(C%PMB)
!        if(allocated(C%ZKM)) deallocate(C%ZKM)
!        if(allocated(C%dPMB)) deallocate(C%dPMB)
!        if(allocated(C%dZKM)) deallocate(C%dZKM)
!        if(allocated(C%Lcor)) deallocate(C%Lcor)
!
!        if(allocated(C%OmTF)) deallocate(C%OmTF)
!        if(allocated(C%VarF)) deallocate(C%VarF,  stat=istat)
!
!        if(allocated(C%Hs)) deallocate(C%Hs, stat=istat)
!        if (istat.ne.0) print*,   ' deallocate(C%Hs  ....',  istat
!        if(allocated(C%Hz)) deallocate(C%Hz, stat=istat)
!        if (istat.ne.0) print *,   ' deallocate(C%Hs  ....',  istat
!        enddo

        END  SUBROUTINE dealloc_saber 
!
      subroutine get_conv_size(file, nz, npix)
      use netcdf     
      implicit none
     
      character(len=*) file
      integer :: nz, npix
!
      integer :: ncid, vid, ierNC, ier_nc
      integer :: nlev
      integer, dimension(nf90_max_var_dims) :: dimidT
!
        ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)

        if(iernc.ne.0) then
          npix = 0
          nz=0
          print * 
          print *, ' cannot open data-file', file  
          print *, ' get_conv_size(file, npix, nz) '
          print * 
          return
        endif
!
        ier_nc=nf90_inq_varid( ncid, 'lon_obs', vid )
        ierNC=nf90_inquire_variable(ncid, vid, dimids=dimidT)
        iernc = nf90_inquire_dimension(ncid, dimidT(1), len=npix)
        ier_nc=nf90_inq_varid( ncid, 'pmb_obs', vid )
        ierNC=nf90_inquire_variable(ncid, vid, dimids=dimidT)
        iernc = nf90_inquire_dimension(ncid, dimidT(1), len=nz)
        iernc=nf90_close(ncid)
        print *, nz, npix, ' nz-npix file SABER', file
!     
       end   subroutine get_conv_size
!
!
       subroutine READ_CONV_LIMB(file,  name_var,    st,    kz1, kz2)
!
! input: file, 
!        namevar,
!        kz1, kz2
!
! output: ST - filled structure of data for the limb-sensor
!     Read ALL DAY
!  
      use netcdf

      implicit none

      character(len=*) file, name_var
      type(str_conlimb) ::  St

      integer :: kz1, kz2, nz, ier_nc
      integer ::  ierror
      character(len=12) :: str_var,str_err, str_e, str_obs
! Read-nc files and fill-in "whole-day struct"
!
      integer :: ncid, vid, ierNC
      integer :: nlev, npix
      integer, dimension(nf90_max_var_dims) :: dimidT
      real, allocatable :: wrk1(:), wrk2(:,:)
      integer :: i, k
        str_e='e'
        str_obs='_obs'
       str_var=trim(name_var)//trim(str_obs)
       str_err=trim(str_e)//str_var
        print *,  str_var, ' .....str_var **_obs '
        print *,  str_err, ' .....str_err  e**_obs '
!       pause  ' conv_limb data '
!
        ierNC=NF90_OPEN(trim(File), nf90_nowrite, ncid)
!
!        print *, ' ncid ', file
        iernc=nf90_inq_varid( ncid, Str_var, vid )

        ierNC=nf90_inquire_variable(ncid, vid, dimids=dimidT)
        iernc = nf90_inquire_dimension(ncid, dimidT(1), len=nlev)
        iernc = nf90_inquire_dimension(ncid, dimidT(2), len=npix)
        print *,  St%nz , St%npix , 'St%nz , St%npix - allocated'
        print *,  nlev, npix, ' data-sizes from nc-file '

       if (St%npix.ne.npix) then

        print *, 'das_datatypes.F90/ read_conv_limb   ..... ERROR '
               print *, ' Check OBSD-allocation for data ', trim(File)
               print *,  St%nz , St%npix , 'St%nz , St%npix - allocated'
        print *, nlev, npix, ' data-sizes from nc-file ', File
!        stop  ' in      das_datatypes.F90/ read_conv_limb '             
!        
       endif
       St%npix = npix
!
!      in das_observer   of das_comp.F90
!      St%nameSens=trim(name_var)
!
        St%nameVar=trim(name_var)
!        print *, St%nameVar, trim(name_var)
!      
      if (St%nz.ne.nlev) then
            print *, ' Restricted use of data in Vertical '
            print *,  ' between data layers: ', kz1, kz2   
            print *,      nlev , ' nlev '
      endif
      nz = St%nz
      allocate(wrk1(nlev),  stat=ierror )
!  if ( ierror /= 0 ) call mpi_wam_quit(239991,' read_conv_limb wrk1(nlev) ')
      allocate( wrk2(nlev,npix), stat=ierror )
! if ( ierror /= 0 )call mpi_wam_quit(239991, ' read_conv_limb wrk2 ')
       print *,  ' name_var   ', name_var
       print *,  '  file ', file
       print *,  '  data-size ', St%nz, St%npix
       wrk1 = 0.
       wrk2 = 0.
!     
!==============================================
!          allocate(St%Val(nlev, npix), St%ERR(nlev, npix))
!          allocate( St%lon(npix), St%lat(npix), St%time(npix))
!          allocate(St%pmb(nlev), St%zkm(nlev), St%dpmb(nlev), St%dzkm(nlev))
!          allocate(St%Lcor(nlev))
!          allocate(St%omtf(nlev, npix), St%varF(nlev, npix))
!
! Obs. operators : Hs and Hz
!!
!          allocate(St%Hs(npix), St%Hz(npix))
!          allocate(St%Hs(npix)%Ix(nphT), St%Hs(npix)%Iy(nphT), St%Hs(npix)%wxy(nphT))
!
! Allocate Forecast-arrays at NPIX-locations
!       
!          allocate(St%XF(plev, npix), St%VF(plev, npix))
!          allocate(St%ZF(plev, npix), St%PF(plev, npix))
!          allocate(St%dZF(plev, npix), St%dPF(plev, npix))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        iernc=nf90_inq_varid( ncid, Str_var, vid ) 

            if (iernc.ne.0) 
     &  print *, Str_var, ' vid is not found in read_conv_limb ', iernc
    
        iernc= nf90_get_var( ncid, vid, wrk2)
!   
        if (iernc.ne.0) 
     &  print *,  Str_var, ' cannot be read in  read_conv_limb ', iernc
        print *,   maxval(wrk2), minval(wrk2) ,  Str_var    

        St%val(1:nz,1:npix) =wrk2(kz1:kz2,1:npix)

!
        iernc=nf90_inq_varid( ncid, Str_err, vid )   
        iernc= nf90_get_var( ncid, vid, wrk2)
        print *,  maxval(wrk2), minval(wrk2) , Str_err
 
       St%err(1:nz,1:npix) =wrk2(kz1:kz2,1:npix)
       print *,   maxval(St%val), maxval(sT%Err), ' St%Saber  ' 


        iernc=nf90_inq_varid( ncid, 'lon_obs', vid )
        iernc= nf90_get_var( ncid, vid, St%lon)

        where(st%lon.lt.0.)
           st%lon = st%lon +360.0
        endwhere

        iernc=nf90_inq_varid( ncid, 'lat_obs', vid )
        iernc= nf90_get_var( ncid, vid, St%lat)
        iernc=nf90_inq_varid( ncid, 'time_obs', vid )
        iernc= nf90_get_var( ncid, vid, St%time)
        print *, maxval(St%lat), maxval(St%lon), maxval(St%time)/86400.
       
        iernc=nf90_inq_varid( ncid, 'pmb_obs', vid )
        iernc= nf90_get_var( ncid, vid, wrk1)
        St%pmb(1:nz) = wrk1(kz1:kz2)
        St%zkm = -7.0*alog(St%pmb/1000.0)

         St%Lcor(1:nz) =5.
         St%omtf(1:nz, 1:npix) = nodata
         St%varF(1:nz, 1:npix) =0.


         St%it1 =1
         St%it2 =npix
!
!         call check_orbits(St%lat, St%lon, St%time, St%gmark, St%val, nz, npix)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! compute dpmb ....... dzkm .....
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        call compute_dthick_obs(nZ, St%pmb, St%zkm, St%dpmb, St%dzkm)
!
!        print * , St%dpmb(1:nlev),  '    St%dPthick-TMLS '

        iernc=nf90_close(ncid)
!        if (str_var.eq.'T_obs') pause

        deallocate(wrk1, wrk2)
!
       end subroutine read_conv_limb
!23456
       subroutine compute_dthick_obs(nz, pmb, zkm, dpmb, dzkm)
       implicit none
       integer :: nz
       real, dimension(nz) :: pmb, zkm,  dpmb, dzkm
       integer :: k
       real :: p1
       if(pmb(1).ge.pmb(2)) then
          do k=2, nz-1
          dpmB = .5*(pmb(k-1)-pmb(k+1))
          dzkm(k) = .5*(zkm(k+1)-zkm(k-1))
          enddo
         P1 =pmb(1)  !(1.+exp(dzkm(2)/7.))
         dzkm(1)  = dzkm(2)
         dpmb(1) = (P1-pmb(2))
         dzkm(nz)  = dzkm(nz-1)
         dpmb(nz) = .5*pmb(nz-1)
       else
         do k=2, nz-1
         dpmB = .5*(pmb(k+1)-pmb(k-1))
         dzkm(k) = .5*(zkm(k-1)-zkm(k+1))
         enddo
! reverse-WACCM-GCM style
         P1 =pmb(nz)  !(1.+exp(dzkm(2)/7.))
         dzkm(nz)  = dzkm(nz-1)
         dpmb(nz) = .5*pmb(nz-1)
         dzkm(1)  = dzkm(2)
         dpmb(1) = .5*pmb(1)

      endif
!     print *
!     do k=1, nz
!     print *, pmB, dpmB, zkm(k), dzkm(k), k
!     enddo
!     print *
      end subroutine compute_dthick_obs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Very simple DQC => Avoid errors in "lon-lat-time"
!                     gmark(k) = 0
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       subroutine  check_orbits(lat, lon, time, gmark, val, nz, npix)
!

      implicit NONE
!23456
       integer ::  nz, npix
       real, dimension(npix) :: lat, lon, time
       real, dimension(nz, npix) :: val
       integer, dimension(npix) :: gmark

       integer, parameter :: lats=-86., latn=+86.
       integer, parameter :: lonw= 0., lone=360.
       integer, parameter ::  ts0=-1., ts24=86401.
       integer :: k, nbad
       gmark(1:npix) =1
       nbad = 0
       do k=1, npix
       if(lat(k).lt.lats.or.lat(k).gt.latn.or.lon(k).lt.lone.
     & or.lon(k).gt.lonw.or. time(k).lt.ts0.or. time(k).gt.ts24) then
         gmark(k) = 0
         val(1:nz,k) = -9999.
!         print *, k, lat(k), lon(k), time(k), ' genobs_check orbits lat=lon=time '
!
! put "extreme" values and gmark(k)=0
!
!         lat(k)    =  90.
!         lon(k)   = 360.
!         time(k) = 86400.
         nbad = nbad+1
        ENDIF
       enddo
       if (nbad.ge.1) then
  
       print *,nbad,' @ of bad data points in subroutine check_orbits '
                 print *, ' time max/min ', maxval(time), minval(time)
                 print *, ' lat max/min ', maxval(lat), minval(lat)
                 print *, ' lon max/min ', maxval(lon), minval(lon) 
        
        endif
          where (gmark.eq.0) 
!
! data points with such "coord" [90. 360, 86400] is a "bad" data-point
! it still allows don't copy All data => Good data and use
!
! "Local grid/time step" sub-types in filters 
!  Obs with gmark(iobs) = 0 is not used in solution
!
                lat = latn
                lon = lonw
                time =86400.
          endwhere
         if (nbad.eq.npix) then 
          print *,  ' Flag data-sets as a bad-source of data '
          call mpi_wam_quit( 23999, ' stop in check_orbits of das_data ')
         endif
       end subroutine check_orbits
!
       SUBROUTINE dt_conlimb(B, C, tw1, tw2)
       use idea_das_utils, only : WHERE_1D  
       implicit NONE
!
!      integer ::   N                   ! N - number of limb-sensors: MLS/SABER/TIDI
!
!       integer ::   t1,   t2            ! time of day in seconds
       real :: tw1, tw2
       type(str_conlimb), pointer     :: B
       type(str_conlimb), pointer     :: C
       integer  :: K, i1, i2
       integer  :: nz, npix, npixB
       integer, allocatable, dimension(:) :: IP 
       real,   allocatable :: Xt(:)

     
       real  :: t1gap, t2gap
!
!        tw1 =float(t1)
!        tw2 =float(t2)
!
!         do k=1, N
            nz   = B%nz
            npixb =B%npix
!             write(6, *)  ' nz, npixb ', nz, npixb
             allocate(Xt(npixb), IP(npixb))
             Xt(1:npixB) =B%time
!
!      write(6, *)  ' tw1 - tw2 ', tw1, tw2, minval(xt), maxval(xt)
            t2gap= minval(xt, mask=XT>tw1)
            t1gap = maxval(xt, mask=XT<tw2)
            if(t2gap.gt.tw2.and.t1gap.lt.tw1) then
!
! data gap
!
              write(6,*) '  data gap  for tw1 - tw2 ', tw1, tw2

              C%npix= 0
              C%nz= 0

             return
             endif
!23456
         call where_1d(tw1, Xt, tw2, npixb, IP, npix, i1, i2)
       
       write(6, *)  ' i1, i2, npix, nz', i1, i2, npix, nz,
     &    Xt(i1), xt(i2)
!
         allocate(C)
!
         C%npix=  npix
         C%nz   = nz
!
! needs extra wor for MLS(1) + SABER(2)  + => C(k)%
!
!         C%namesens =B%namesens
!
!         C%nameVAR =B%nameVAR


         allocate(C%lon(npix), C%lat(npix), C%time(npix))
         allocate(C%gmark(npix))

         allocate(C%zkm(nz), C%pmb(nz), C%Lcor(nz))
         allocate(C%dzkm(nz), C%dpmb(nz) )
         allocate(C%val(nz, npix), C%ERR(nz, npix))
         allocate(C%OmTF(nz, npix), C%varF(nz, npix))
         allocate(C%AK(nz,nz), C%VAPR(nz) )
!   allocate(C%Hs(npix), C%Hz(npix))
!
          C%zkm = B%zkm
          C%pmb = B%pmb
          C%dzkm = B%dzkm
          C%dpmb = B%dpmb
          C%Lcor   = B%Lcor
!          C%AK =  B%AK 
!          C%VAPR= B%Vapr
!
! copy
!
           C%lon(1:npix) =B%lon(i1:i2)
           C%lat(1:npix) =B%lat(i1:i2)
           C%time(1:npix) =B%time(i1:i2)

!           C%gmark(1:npix) =B%gmark(i1:i2)

           C%val(:, 1:npix) =B%val(:, i1:i2) 
           C%ERR(:, 1:npix) =B%ERR(:, i1:i2) 

             C%OmTF =-999.   
             C%VarF =   0.   

             C%it1 =  i1
             C%it2 =  i2

             deallocate(Xt, IP)
!   C(k)%         enddo
!
           END  SUBROUTINE dt_conlimb
!
! next:    SUBROUTINE compute_omtf(Xf, Ef, lonf, latf, Cdt)
!
      end module  das_datatypes
