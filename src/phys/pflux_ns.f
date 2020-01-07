       Module Pflux_NS_Bx_mod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%   Author: Russell Cosgrove, Sept. 26, 2019
!%
!%   Translated from January 20th, 2014 Matlab version.
!%   Changes from January 20th, 2014 Matlab version: rounding was removed from totalEnergyFlux_GWatts_1 and totalEnergyFlux_GWatts_2, 
!%		which may create changes in the third significant figure for low activity.  This version is deemed slightly superior.
!%
!%   DESCRIPTION:
!%
!%       This function takes user input of model parameters to produce
!%       polar plots of Poytning flux over the Northern Hemisphere.
!%
!%   INPUTS:
!%       clka : (IMF clock angle) is in (degrees).
!%              (0 is north, +-180 is south.)
!%
!%       BT : (IMF magnitude in GSM plane) is in (nT).
!%            (median is around 4 nT, standard deviation is also around 4 nT)
!%
!%       Vsw : (Solar wind speed) is in (km/s).
!%             (median is around 400 km/s. Standard deviation is around 100 km/s.)
!%
!%       Np : (Solar wind number density) is in (cm^{-3}).
!%            (median is around 4 cm^-3, Standard deviation is also around 4 cm^-3.)
!%
!%       sin_dipt : (Sine of the dipole tilt angle).
!%                  (median is 0.  Limits are +-0.55.)
!%
!%       Bx : (GSM Bx) nT.
!%            (median is 0. Range is -10 to 10)
!%
!%       AL : (Al index).
!%            (median is about -51.  Range is 0 to -500.)
!%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      implicit none
! File sizes for reading
      integer, parameter :: nN=7790
      real (kind=8)    :: altitude_N(nN), geoParam_N(nN,8),stdNorth(6),
     & a(34,16), a_AL(34,16)
      integer            :: NumOrbitsInSupport
      real (kind=8)      :: totalEnergyFlux_GWatts
      real (kind=4)      :: BV(12000,16), mu(12000), fac(12000),
     & pp_out(60,200), ILAT(60,200), MLT(60,200)
      real (kind=8), allocatable :: omega(:,:)

!---------------SUBROUTINES and FUNCTIONS DEFINED BELOW-------------------------
      contains
!-------------------------------------------------------------------------------
!---Main Program with Bx--------------------------
      Subroutine Pflux_NS_Bx(clka,BT,Vsw,Np,sin_dipt,Bx)
       implicit none
      real (kind=8) :: clka,BT,Vsw,Np,sin_dipt,Bx
      real (kind=8) :: Alt,efld,momentum,test(2,2), pi = 4.D0*atan(1.D0)
      real (kind=8) :: ff, r1, r2, r3, r4, r5, n1, n2, n3, n4, n5, n6,
     &                 n7,  bv_var(34), bv_y(16)
      real (kind=8) ::totalEnergyFlux_GWatts_1,totalEnergyFlux_GWatts_2
     &, cellArea
      real (kind=4) :: ppp(12000), pp(12000)
      logical :: mask1(nN), mask2(nN), mask3(nN), mask4(nN), mask5(nN),
     & mask6(nN), mask(nN)
      integer :: ones(nN), ind

! Calculating efield and momentum flux from inputs
      efld = BT**(2.0D0/3.0D0)*Vsw
      momentum = Np*(Vsw**2.0D0)
      Alt = (2650.0D0+3300.0D0)/2.0D0

! Fractional range for search
      ff = 1.0D0/8.0D0
! Typical ranges for parameters
      r1 = 2.0D0
      r2 = 2905.0D0
      r3 = (4.0e6-7.5e4)
      r4 = 1.1D0
      r5 = 20.0D0

      mask1 = abs(cos(geoParam_N(:,1))-cos(clka*pi/180.0D0)) < r1*ff
      mask2 = abs(sin(geoParam_N(:,1))-sin(clka*pi/180.0D0)) < r1*ff
      mask3 = abs(geoParam_N(:,2)-efld) < r2*ff
      mask4 = abs(geoParam_N(:,3)-momentum) < r3*ff
      mask5 = abs(geoParam_N(:,4)-sin_dipt) < r4*ff
      mask6 =  abs(geoParam_N(:,8)-Bx) < r5*ff
      mask = mask1 .and. mask2 .and. mask3 .and. mask4 .and. mask5 .and.
     & mask6
      ones = 1
      NumOrbitsInSupport = sum(pack(ones,mask))

      n1 = cos(clka*pi/180.0D0)
      n2 = sin(clka*pi/180.0D0)
      n3 = efld
      n4 = momentum
      n5 = sin_dipt
      n6=Bx
      n7 = Alt
      n3 = n3/stdNorth(1)
      n4 = n4/stdNorth(2)
      n5 = n5/stdNorth(3)
      n6 = n6/stdNorth(5)
      n7 = n7/stdNorth(6)

      bv_var = (/ n1, n2, n3, n4, n5, n6, n7, cos(2*clka*pi/180.0D0),
     &  sin(2*clka*pi/180.0D0),
     &  n3**2.0D0, n4**2.0D0, n5**2.0D0, n6**2.0D0, n7**2.0D0,
     &  n1*n3, n1*n4, n1*n5, n1*n6, n1*n7, n2*n3, n2*n4, n2*n5, n2*n6,
     &  n2*n7, n3*n4,
     &  n3*n5, n3*n6, n3*n7, n4*n5, n4*n6, n4*n7, n5*n6, n5*n7, 1.0D0 /)

      do ind=1,16
       bv_y(ind) = sum(bv_var*a(:,ind))
      end do

      do ind=1,12000
       pp(ind) = sum(BV(ind,:)*bv_y) + mu(ind)
      end do

!Estimate of the missing component
      ppp = abs(pp*fac)

!Calculate area of a cell in the polar mesh (all cells have the same area)
      cellArea = (2.0D0*pi*(1.0D0-cos(30.0D0*pi/180.0D0))
     &           *6471.0e3**2.0D0)/(60.0D0*200.0D0)

!Compute total Poynting flux in the uncorrected Poynting flux map and
!in the estimate of the missing component.
      totalEnergyFlux_GWatts_1 = sum(cellArea*ppp)/1.0e9
      totalEnergyFlux_GWatts_2 = sum(cellArea*abs(pp))/1.0e9

!Form final Poynting flux estimate as explained in Cosgrove et al. (2013)
      pp = pp +abs(pp)*totalEnergyFlux_GWatts_1/totalEnergyFlux_GWatts_2
      totalEnergyFlux_GWatts = sum(cellArea*pp)/1.0e9

      pp_out = reshape(pp, (/ 60, 200 /) )



      End Subroutine Pflux_NS_Bx

!-------------------------------------------------------------------------------
!---Main Program with AL--------------------------
      Subroutine Pflux_NS_AL(clka,BT,Vsw,Np,sin_dipt,AL)
        implicit none
        real (kind=4) :: clka,BT,Vsw,Np,sin_dipt,AL
      real (kind=8) :: Alt,efld,momentum,test(2,2), pi = 4.D0*atan(1.D0)
        real (kind=8) :: ff, r1, r2, r3, r4, r5, n1, n2, n3, n4, n5, n6,
     &                   n7,  bv_var(34), bv_y(16)
        real (kind=8) :: totalEnergyFlux_GWatts_1,
     &                   totalEnergyFlux_GWatts_2, cellArea
        real (kind=4) :: ppp(12000), pp(12000)
        logical       :: mask1(nN), mask2(nN), mask3(nN), mask4(nN),
     &                   mask5(nN), mask6(nN), mask(nN)
        integer       :: ones(nN), ind

!Calculating efield and momentum flux from inputs
        efld = BT**(2.0D0/3.0D0)*Vsw
        momentum = Np*(Vsw**2.0D0)
        Alt = (2650.0D0+3300.0D0)/2.0D0

! Fractional range for search
        ff = 1.0D0/8.0D0
! Typical ranges for parameters
        r1 = 2.0D0
        r2 = 2905.0D0
        r3 = (4.0e6-7.5e4)
        r4 = 1.1D0
        r5 = 500.0D0

       mask1 = abs(cos(geoParam_N(:,1))-cos(clka*pi/180.0D0)) < r1*ff
       mask2 = abs(sin(geoParam_N(:,1))-sin(clka*pi/180.0D0)) < r1*ff
       mask3 = abs(geoParam_N(:,2)-efld) < r2*ff
       mask4 = abs(geoParam_N(:,3)-momentum) < r3*ff
       mask5 = abs(geoParam_N(:,4)-sin_dipt) < r4*ff
       mask6 =  abs(geoParam_N(:,5)-AL) < r5*ff
       mask = mask1 .and. mask2 .and. mask3 .and. mask4 .and. mask5.and.
     & mask6
       ones = 1
       NumOrbitsInSupport = sum(pack(ones,mask))

       n1 = cos(clka*pi/180.0D0)
       n2 = sin(clka*pi/180.0D0)
       n3 = efld
       n4 = momentum
       n5 = sin_dipt
       n6 = AL
       n7 = Alt
       n3 = n3/stdNorth(1)
       n4 = n4/stdNorth(2)
       n5 = n5/stdNorth(3)
       n6 = n6/stdNorth(4)
       n7 = n7/stdNorth(6)

       bv_var = (/ n1, n2, n3, n4, n5, n6, n7, cos(2*clka*pi/180.0D0),
     &  sin(2*clka*pi/180.0D0),
     &  n3**2.0D0, n4**2.0D0, n5**2.0D0, n6**2.0D0, n7**2.0D0,
     &  n1*n3, n1*n4, n1*n5, n1*n6, n1*n7, n2*n3, n2*n4, n2*n5, n2*n6,
     &  n2*n7, n3*n4,
     &  n3*n5, n3*n6, n3*n7, n4*n5, n4*n6, n4*n7, n5*n6, n5*n7, 1.0D0 /)

       do ind=1,16
        bv_y(ind) = sum(bv_var*a_AL(:,ind))
       end do

       do ind=1,12000
       pp(ind) = sum(BV(ind,:)*bv_y) + mu(ind)
       end do

!Estimate of the missing component
       ppp = abs(pp*fac)

!Calculate area of a cell in the polar mesh (all cells have the same area)
      cellArea = (2.0D0*pi*(1.0D0-cos(30.0D0*pi/180.0D0))*
     &6471.0e3**2.0D0)/(60.0D0*200.0D0)

!Compute total Poynting flux in the uncorrected Poynting flux map and
!in the estimate of the missing component.
      totalEnergyFlux_GWatts_1 = sum(cellArea*ppp)/1.0e9
      totalEnergyFlux_GWatts_2 = sum(cellArea*abs(pp))/1.0e9

!Form final Poynting flux estimate as explained in Cosgrove et al. (2013)
      pp = pp +abs(pp)*totalEnergyFlux_GWatts_1/totalEnergyFlux_GWatts_2

      totalEnergyFlux_GWatts = sum(cellArea*pp)/1.0e9

      pp_out = reshape(pp, (/ 60, 200 /) )

      End Subroutine Pflux_NS_AL

!-------------------------------------------------------------------------------
!---Read in the model parameters--------------------------
      Subroutine Read_model_parameters

      call Read_mesh_txt('epf_geoParam_N.txt', nN, 8)
      geoParam_N = omega
      deallocate(omega)

      call Read_mesh_txt('epf_altitude_N.txt', nN, 1)
      altitude_N = omega(:,1)
      deallocate(omega)

      call Read_mesh_txt('epf_stdNorth.txt', 1, 6)
      stdNorth = omega(1,:)
      deallocate(omega)

      call Read_mesh_txt('epf_NS_NS_noAl_BxandAlt.txt', 34, 16)
      a = omega
      deallocate(omega)

      call Read_mesh_txt('epf_NS_NS_noBx_AlandAlt.txt', 34, 16)
      a_AL = omega
      deallocate(omega)

      call Read_mesh_txt('epf_Basis_NS_AlandBx_BV.txt', 12000, 16)
      BV = omega
      deallocate(omega)

      call Read_mesh_txt('epf_Basis_NS_AlandBx_mu.txt', 12000, 1)
      mu = omega(:,1)
      deallocate(omega)

      call Read_mesh_txt('epf_ffmc.txt', 12000, 1)
      fac = omega(:,1)
      deallocate(omega)

      call Read_mesh_txt('epf_ILAT.txt', 60, 200)
      ILAT = omega
      deallocate(omega)

      call Read_mesh_txt('epf_MLT.txt', 60, 200)
      MLT = omega
      deallocate(omega)


      End Subroutine Read_model_parameters

!-------------------------------------------------------------------------------
!---Read a single 2D matrix from text file--------------------------
      Subroutine Read_mesh_txt(filename,nrow,ncolumn)
      implicit none
      integer :: unitno=10
      integer :: kx,ncolumn,nrow
      character (len=*) :: filename
!call getarg(1,filename)
!call get_command_argument(1,filename)
!call get_command_argument(2,nrow_t)
!call get_command_argument(3,ncolumn_t)
!read(nrow_t,*)nrow
!read(ncolumn_t,*)ncolumn
      open(unit=unitno,file=trim(filename),form="formatted",
     &status="old",action="read")
      allocate (omega(nrow,ncolumn))
       do kx=1,nrow
       read(unitno,*) omega(kx,:)
       end do
      close(unitno)
      End Subroutine Read_mesh_txt
!-------------------------------------------------------------------------------

      END Module Pflux_NS_Bx_mod

!Program test
!	use Pflux_NS_Bx_mod
!	call Read_model_parameters
!	call Pflux_NS_Bx(-67.3,5.2,522.3,3.6,-0.14,4.7)
!	!call Pflux_NS_AL(-67.3,5.2,522.3,3.6,-0.14,-61.7)
!	print*,'NumOrbitsInSupport',NumOrbitsInSupport
!	print*,'totalEnergyFlux_GWatts',totalEnergyFlux_GWatts
!	print*,'pp_out(16,50)',pp_out(16,50)
!	print*,'ILAT(16,50)',ILAT(16,50)
!	print*,'MLT(16,50)',MLT(16,50)
!End Program test
