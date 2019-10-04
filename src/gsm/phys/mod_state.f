      module mod_state
!
!    New module to supply domain information to the GFS output routines
!    called by wrtout.
!
!! REVISION LOG:
!  Jul 2010 Sarah Lu, add buff_mult_pieceg for 2d aerosol diag fields
!  Aug 2010 Jun Wang  add buff_mult_piecenst for 2d nst fields
!
      use machine,    ONLY: kind_io4
      implicit none
!
!
      real(kind=kind_io4), allocatable :: buff_mult_piece(:),
     1                                    buff_mult_pieces(:)
      real(kind=kind_io4),allocatable,target :: buff_mult_piecef(:,:,:),
     1                                    buff_mult_piecesf(:,:,:,:)
      real(kind=kind_io4),allocatable,target :: buff_mult_pieceg(:,:,:)
      real(kind=kind_io4),allocatable,target :: 
     1                                    buff_mult_piecenst(:,:,:)
!
      real(kind=kind_io4), allocatable,target ::
     1                                    buff_mult_piecea2d(:,:,:),
     1                                    buff_mult_piecea3d(:,:,:)
!
      integer , allocatable :: ivar_global(:),ivar_global_a(:,:)
     &,                        ivarg_global(:),ivarg_global_a(:,:)
      integer , allocatable :: maskss(:,:,:)
!
      integer ngrid ,ngrida,ngridg
!jw
      integer ngrid2d,ngrid3d,ngridnst

      end module mod_state

