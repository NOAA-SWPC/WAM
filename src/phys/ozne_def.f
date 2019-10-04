      module ozne_def
      use machine , only : kind_phys
      implicit none
      
      integer, parameter :: kozpl=28, kozc=48
      integer latsozp, levozp, timeoz, latsozc, levozc, timeozc
     &,       PL_Coeff
      real (kind=kind_phys) blatc, dphiozc
      real (kind=kind_phys), allocatable :: PL_LAT(:), PL_Pres(:)
     &,                                     PL_TIME(:)
      end module ozne_def
