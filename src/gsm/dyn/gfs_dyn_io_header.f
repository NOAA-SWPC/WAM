      module gfs_dyn_io_header
      use gfs_dyn_machine
      implicit none
      save
      integer              ifin
      integer              icen
!jw      integer              icen2
!jw      integer              ienst
!jw      integer              iensi
!jw      integer              itrun
       integer,target  ::   icen2
       integer,target  ::   ienst
       integer,target  ::   iensi
       integer,target  ::   itrun
!
      integer lonb, latb, iens(5), idpp, idvt, idrun
     &,       idusr, ncldt, irealf, iorder, lonr, latr
!
!jw--- save iniital oro data
      REAL(KIND=KIND_IOr), allocatable :: Z_R(:)
      REAL(KIND=KIND_IO4), allocatable :: Z(:)
!     real(kind=kind_io4), allocatable :: gz_grid(:,:)      !  Moorthi
!!
      character*8   lab(4)
      end module gfs_dyn_io_header
