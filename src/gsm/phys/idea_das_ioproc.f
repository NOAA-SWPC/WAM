!
! DAS on single "IO"
!
      SUBROUTINE make_das(Tf, Pf, X3f, Xaf, nx, ny, nz, ioproc, 
     & tw1, tw2, fhour,  Curr_NC_WAMDAY)
       use idea_das_saber, only:   Stsab_T, Stsab_O3, Stsab_OP
       use das_datatypes, only : str_conlimb
       use das_datatypes, only : dt_conlimb, dealloc_saber
!       use idea_das_utils, only:   WHERE_1D  
!
! Need to add lats/lons in addition to Pf
!      use saber_data_1hr, only : Pd, Td, X3d, Xpf
!      use saber_data_1hr, lond, latd, pdata, nxd, nzs
!      vertical grid of WAM
!
      IMPLICIT NONE
      integer , intent(in)        :: ioproc
      integer , intent(in)        :: nx, ny, nz
      real, dimension(nx, ny,nz)  :: Tf, Pf, X3f, Xaf 
      real, intent(in)            :: tw1, tw2, fhour
      integer, intent(in)         :: Curr_NC_WAMDAY
!
! locals
!
      integer :: i, j, k, ipix
      integer :: Nobs, Nzd
      real    :: XL, YL
      real, allocatable, dimension(:,:) ::  oTf, o3f, oAf 
      TYPE(str_conlimb), pointer :: Stsab_dT,  Stsab_dO3, Stsab_dOP
      real :: ts1, ts2,rhr
!----------------------------------------------------------
! take from Stsab_T, Stsab_O3, Stsab_OP and allocate
!         Pd, Td, X3d, Xpf,    latd, pdata, nxd, nzs
!-----------------------------------------------------------
       ts1 = 3600.*tw1
       ts2 = 3600.*tw2
       rhr = 1./3600.
!
       print *, 'VAY:make_das on ioproc=', ioproc
       print *, 'VAY-DAS-day', Curr_NC_WAMDAY
       print *, 'VAY-DAS-wind', tw1, fhour, tw2

       call dt_conlimb(Stsab_T,  Stsab_dT,  ts1, ts2)
       print *, 'after dt_conlimb VAY:das_dT'

       call dt_conlimb(Stsab_O3, Stsab_dO3, ts1, ts2)
       call dt_conlimb(Stsab_OP, Stsab_dOP, ts1, ts2)
!
       print *,'VAY:das_dT', 
     & maxval(Stsab_dT%Val),minval(Stsab_dT%Val,mask=Stsab_dT%Val>0.)
       print *,'VAY:das_Npix_Nz', Stsab_dT%npix, Stsab_dT%nz
       print *,'VAY:das_Lon',maxval(Stsab_dT%lon), minval(Stsab_dT%lon)
       print *,'VAY:das_Lat',maxval(Stsab_dT%lat), minval(Stsab_dT%lat)
!
       print *, 'VAY:das_THR',
     & maxval(Stsab_dT%time)*rhr, minval(Stsab_dT%time)*rhr

!      integer, allocatable, :: NobsG1(:), NobsG2(:), NobsG3(:)    
!
! simple 3dVAR for SABER T, O3, O3P  lond, latd, pdata
!  FWD-Step
!      do k=1, nzd
!        do ipix = 1, Nobs
!      XL = lond(ipix)
!      YL = latd(ipix)     
! 
!       oTf(k, i) = linterp2(XL, YL, Tf(:,:,k), lonf, latf)- dT(k,i)
!       o3f(k, i) = linterp2(XL, YL, x3f(:,:,k), lonf, latf)- dO3(k,i)
!       oAf(k, i) = linterp2(XL, YL, xAf(:,:,k), lonf, latf)- dOP(k,i)
!
!     enddo  
!     enddo
!  INV-step
!      do k=kf1, Kf1+nzd-1
!       do j=1,ny
!       do i=1,nx
!        if (wg1(i,j,k).gt.1 ) then
!          do ipix = 1, NobsG1(k)
!           enddo
!         endif
!          do ipix = 1, NobsG2(k)
!           enddo
!          do ipix = 1, NobsG3(k)
!           enddo
!       enddo
!       enddo
!      enndo 

       call dealloc_saber(Stsab_dT)
       call dealloc_saber(Stsab_dO3)
       call dealloc_saber(Stsab_dOP)

       print *, 'VAY:make_das on ioproc=', ioproc
       print *, 'VAY-DAS-day', Curr_NC_WAMDAY
       print *, 'VAY-DAS-wind', tw1, fhour, tw2
      END SUBROUTINE MAKE_DAS
