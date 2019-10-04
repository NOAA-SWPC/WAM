      module deldifs_def
      use gfs_dyn_machine
      implicit none
      save
      real(kind=kind_evod),allocatable :: dne(:),dno(:),
     &                                    sf(:),rtrd(:),rthk(:),
     &                                    bkly(:),ckly(:)
     &,                                   dneh(:), dnoh(:),cthk(:)
      real(kind=kind_evod) rtnp
      integer jdel
      end module deldifs_def
