      module resol_def
      use machine
      implicit none
      save
      integer              jcap,jcap1,jcap2,latg,latg2,latr,latr2
      integer              levh,levm1,levp1,levs,lnt,lnt2,lnt22
      integer              lnte,lnted,lnto,lntod,lnuv
      integer              lonf,lonfx,lonr,lonrx
      integer              ntrac
      integer              nxpt,nypt,jintmx,latgd
      integer ntoz,ntcw
      integer lsoil,nmtvr,num_p3d,num_p2d
      integer ngrids_sfcc, ngrids_flx
! jbao new gfs phys
      integer ngrids_nsst,nr_nsst,nf_nsst
!     real(kind=kind_evod) rerth

      INTEGER   P_GZ,P_ZEM,P_DIM,P_TEM,P_RM,P_QM
      INTEGER   P_ZE,P_DI,P_TE,P_RQ,P_Q,P_DLAM,P_DPHI,P_ULN,P_VLN
      INTEGER   P_W,P_X,P_Y,P_RT,P_ZQ
      INTEGER                LOTS,LOTD,LOTA

      integer kwq,kwte,kwdz,kwrq
!jbao
!JFM  parameter (levs = NVL_VALUE, levp1 = levs+1, latr = 1, lonr = 1)
!JFM  levs and levp1 are set in physics.F90
      parameter (latr = 1, lonr = 1)
      integer   thermodyn_id, sfcpress_id                       ! hmhj
      end module resol_def

! jbao new gfs physics
      module ozne_def
      use machine , only : kind_phys
      implicit none
      save
      integer, parameter :: kozpl=28, kozc=48
      integer latsozp, levozp, timeoz, latsozc, levozc, timeozc &
     &,       PL_Coeff
      real (kind=kind_phys) blatc, dphiozc
      real (kind=kind_phys), allocatable :: PL_LAT(:), PL_Pres(:) &
     &,                                     PL_TIME(:)
      end module ozne_def
! end jbao new gfs physics
