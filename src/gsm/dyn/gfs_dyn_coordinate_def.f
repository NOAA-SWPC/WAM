      module gfs_dyn_coordinate_def
!
! program log
! 2011 02 20 : henry juang, add some matrixes for mass_dp and ndsl options
!
      use gfs_dyn_machine
      implicit none
      
       real(kind=kind_evod) , allocatable ::
!jw     . AK5(:),BK5(:),CK5(:),CK(:),DBK(:),bkl(:),   		! hmhj
     .                      CK(:),DBK(:),bkl(:),   		! hmhj
     . AMHYB(:,:),BMHYB(:,:),SVHYB(:),tor_hyb(:),
     . SMHYB(:,:),HMHYB(:,:),
     . D_HYB_m(:,:,:),THREF(:),dm205_hyb(:,:,:)			! hmhj
!jws
       real(kind=kind_evod) eps_si                              ! hmhj
       integer(kind=kind_io4),target :: vertcoord_id            ! hmhj
       real(kind=kind_evod),allocatable,target :: AK5(:),BK5(:),CK5(:)
!jwe

!
      real(kind=kind_evod), allocatable :: vcoord(:,:)
      real(kind=kind_evod), allocatable :: y_ecm(:,:),
     &                                     t_ecm(:,:),sv_ecm(:),
     &                                     AM_slg(:,:),BM_slg(:,:),
     &                                     D_slg_m(:,:,:),
     &                                     SV_slg(:),tor_slg(:)
      integer nvcoord, idsl, idvc, idvm
      end module gfs_dyn_coordinate_def
