       module akbk_hyb_def
!     use resol_def
      use machine
      implicit none
      save
       real(kind=kind_evod) , allocatable ::
     . AK5(:),BK5(:),CK(:),DBK(:),bkl(:),   
     . AMHYB(:,:),BMHYB(:,:),SVHYB(:),tor_hyb(:),
     . D_HYB_m(:,:,:)

       end module akbk_hyb_def
