!#include "../../../ESMFVersionDefine.h"

!-----------------------------------------------------------------------
!
      module module_export_import_data
!
!-----------------------------------------------------------------------
!
!***  list the roots of the field names of the arrays that will be
!***  transferred between export and import states during the
!***  integration.  these lists can then be used in simple do loops
!***  to redirect the data's pointers.
!
! Oct 17 2009       Sarah Lu, remove shum, oz, cld; modify ndata count
! Feb 20 2010       Sarah Lu, p, dp added to phy export state
! Feb 29 2011       Henry Juang, add dp as prognostic variable for NDSL
!-----------------------------------------------------------------------
!
      USE ESMF
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  set the number of 2- and 3-dimensional field names
!***  and how many to be transferred between components.
!-----------------------------------------------------------------------
!
!---------------------------------
!***  total number of field names
!---------------------------------
!
      integer,parameter :: ndata_1d=1
      integer,parameter :: ndata_2d=3
!     integer,parameter :: ndata_2d=2
!*    integer,parameter :: ndata_3d=9
      integer,parameter :: ndata_3d=6 ! For WAM, no stochastic_wts
!      integer,parameter :: ndata_3d=12
!
!-----------------------------------------------------
!***  number of fields moved from dynamics 
!-----------------------------------------------------
!
      integer,parameter :: ndata_1d_dyn_imp=1
      integer,parameter :: ndata_2d_dyn_imp=3
!*    integer,parameter :: ndata_3d_dyn_imp=6
      integer,parameter :: ndata_3d_dyn_imp=6
      integer,parameter :: ndata_1d_dyn_exp=1
      integer,parameter :: ndata_2d_dyn_exp=2
!*    integer,parameter :: ndata_3d_dyn_exp=9
      integer,parameter :: ndata_3d_dyn_exp=6  ! For WAM, no stochastic_wts
!      integer,parameter :: ndata_3d_dyn_exp=12  ! added stochastic_wts
!
!-----------------------------------------------------
!***  number of fields moved from physics 
!-----------------------------------------------------
!
      integer,parameter :: ndata_1d_phy_imp=1
      integer,parameter :: ndata_2d_phy_imp=2
!*    integer,parameter :: ndata_3d_phy_imp=9
      integer,parameter :: ndata_3d_phy_imp=6   ! For WAM, no stochastic_wts
!      integer,parameter :: ndata_3d_phy_imp=12 ! added stochastic wts
      integer,parameter :: ndata_1d_phy_exp=1
      integer,parameter :: ndata_2d_phy_exp=3
!*    integer,parameter :: ndata_3d_phy_exp=6
!*    integer,parameter :: ndata_3d_phy_exp=3     ! q,cld,o3 removed
!     integer,parameter :: ndata_3d_phy_exp=5     ! p, dp added
      integer,parameter :: ndata_3d_phy_exp=6     ! p, dp, dpdt  added Moorthi - Nov23, 2015
!
!----------------------------------------------------------------
!***  the names of the fields that will move through the coupler
!----------------------------------------------------------------
!
      character(esmf_maxstr),dimension(ndata_1d) :: datanames_1d        &
                                                     =(/'date'/)
!
      character(esmf_maxstr),dimension(ndata_2d) :: datanames_2d        &
                                                     =(/'hs  '          &
                                                       ,'ps  '          &
                                                       ,'rqtk'          &
                                                            /)
!
      character(esmf_maxstr),dimension(ndata_3d) :: datanames_3d        &
                                                     =(/'t     '        &
                                                       ,'u     '        &
                                                       ,'v     '        &
!*                                                     ,'shum  '        &
!*                                                     ,'oz    '        &
!*                                                     ,'cld   '        &
                                                       ,'p     '        &
                                                       ,'dp    '        &
                                                       ,'dpdt  '        &
!                                                       ,'shum_wts'      &
!                                                       ,'sppt_wts'      &
!                                                       ,'skebu_wts'     &
!                                                       ,'skebv_wts'     &
!                                                       ,'vcu_wts'       &
!                                                       ,'vcv_wts'       &
                                                            /)
!
!-----------------------------------------------------------------------
!
      end module module_export_import_data
!
!-----------------------------------------------------------------------
