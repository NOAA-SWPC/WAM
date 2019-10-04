module module_quasipost
  use module_RELAX4E
  use module_REDUCTION
  implicit none
  private

  public :: quasipost

  integer, parameter :: npres = 33
  real, parameter :: badheight=-9e9

  ! NCEP Unified Post standard pressure levels (SLPDEF) used for this
  ! membrane MSLP calculation.  These are ALL of the post pressure
  ! levels up to 200mbar:
  real, parameter :: post_stdpres(npres) = (/ 20000.,          &
       22500., 25000., 27500., 30000., 32500., 35000., 37500., 40000., &
       42500., 45000., 47500., 50000., 52500., 55000., 57500., 60000., &
       62500., 65000., 67500., 70000., 72500., 75000., 77500., 80000., &
       82500., 85000., 87500., 90000., 92500., 95000., 97500.,100000./)

  ! index within post_stdpres of the 850mbar, 700mbar and 500mbar
  ! levels, respectively:
  integer, parameter :: k850 = 27, k700=21, k500=13

  ! Pressure "interface" levels, used only for interpolation.  These
  ! are half-way between pressure levels (post_stdpres) in pressure
  ! space (instead of z, Z or density), to match assumptions made in
  ! the post's Memberane MSLP calculation:
  real, parameter :: post_istdpres(npres+1) = (/ 18750., &
       21250., 23750., 26250., 28750., 31250., 33750., 36250., 38750., &
       41250., 43750., 46250., 48750., 51250., 53750., 56250., 58750., &
       61250., 63750., 66250., 68750., 71250., 73750., 76250., 78750., &
       81250., 83750., 86250., 88750., 91250., 93750., 96250., 98750., &
       101250./)

  ! Constants from the NCEP Unified Post used for interpolation and
  ! extrapolation:
  real, parameter :: post_H1=1.0
  real, parameter :: post_PQ0=379.90516
  real, parameter :: post_A2=17.2693882
  real, parameter :: post_A3=273.16
  real, parameter :: post_A4=35.86
  real, parameter :: post_D608=0.608
  real, parameter :: post_RD=287.04
  real, parameter :: post_G=9.81
  real, parameter :: post_GAMMA=6.5E-3
  real, parameter :: post_RGAMOG=post_RD*post_GAMMA/post_G
  real, parameter :: post_RHmin=1.0E-6     ! minimal RH bound
  real, parameter :: post_smallQ=1.E-12

  real, parameter :: post_slope=-6.6e-4 ! K/km

  REAL, PARAMETER :: old_COEF3=post_RD*post_SLOPE
  REAL, PARAMETER :: old_COEF2=-1./old_COEF3

  type TRACKER_DATA
     integer :: ids,ide,jds,jde,kds,kde
     integer :: ims,ime,jms,jme,kms,kme
  end type TRACKER_DATA


contains
  subroutine quasipost_message(what)
    character*(*), intent(in) :: what
    print "('Quasipost: ',A)",trim(what)
  end subroutine quasipost_message

  subroutine quasipost(grid)
    use MODULE_SOLVER_INTERNAL_STATE, only : SOLVER_INTERNAL_STATE
    type(solver_internal_state), intent(inout) :: grid
    integer :: ids,ide,jds,jde,kds,kde, &
               ims,ime,jms,jme,kms,kme, &
               ips,ipe,jps,jpe,kps,kpe
    ids=grid%ids ; jds=grid%jds ; kds=1
    ide=grid%ide ; jde=grid%jde ; kde=grid%LM
    ims=grid%ims ; jms=grid%jms ; kms=1
    ime=grid%ime ; jme=grid%jme ; kme=grid%LM
    ips=grid%its ; jps=grid%jts ; kps=1
    ipe=grid%ite ; jpe=grid%jte ; kpe=grid%LM
    call quasipost_impl(grid,ids,ide,jds,jde,kds,kde,&
                            ims,ime,jms,jme,kms,kme,&
                            ips,ipe,jps,jpe,kps,kpe)
  end subroutine quasipost

  subroutine quasipost_impl(grid,ids,ide,jds,jde,kds,kde,&
                            ims,ime,jms,jme,kms,kme,&
                            ips,ipe,jps,jpe,kps,kpe)
    use MODULE_SOLVER_INTERNAL_STATE, only : SOLVER_INTERNAL_STATE
    use module_exchange, only: halo_exch
    use MPI
    implicit none

    ! ------------------------------
    
    type(SOLVER_INTERNAL_STATE), intent(inout) :: grid
    integer :: ids,ide,jds,jde,kds,kde, &
               ims,ime,jms,jme,kms,kme, &
               ips,ipe,jps,jpe,kps,kpe
    real :: presTv(ims:ime,jms:jme,npres), Pmsl(ims:ime,jms:jme)
    real :: presZ(ims:ime,jms:jme,npres)
    integer :: relaxmask(ims:ime,jms:jme)
    real :: relaxwork(ims:ime,jms:jme)
    real :: interP(ims:ime,jms:jme,npres+1)
    real(kind=8) :: dum8
    integer :: ground_mask(ims:ime,jms:jme,npres)
    integer :: ground_level(ims:ime,jms:jme)
    integer :: ipres,i,j,mpres,imin,jmin,k,need_to_relax,imax,jmax,ierr
    real :: pmin

    ! Make sure the three constant pressure level values are right:
100 format('In module_membrane_mslp, post_stdpres(',A,')=',F0.3,&
           ' but should be ',F0.3)
    if(abs(post_stdpres(k850)-85000.)>1) then
       write(0,100) 'k850',post_stdpres(k850),85000.
       call mpi_abort(MPI_COMM_WORLD,1,ierr)
    endif
    if(abs(post_stdpres(k700)-70000.)>1) then
       write(0,100) 'k700',post_stdpres(k700),70000.
       call mpi_abort(MPI_COMM_WORLD,1,ierr)
    endif
    if(abs(post_stdpres(k500)-50000.)>1) then
       write(0,100) 'k500',post_stdpres(k500),50000.
       call mpi_abort(MPI_COMM_WORLD,1,ierr)
    endif

    if(size(grid%p700rv)>1) then
       ! Need a halo for winds in order to get vorticity and H point wind magnetudes:
       call HALO_EXCH(grid%V10,1,grid%U10,1,3,3)
    endif

    ! UPPER BOUND: MPRES

    ! Find mpres: the lowest pressure that we need to handle.  This is
    ! mostly for efficiency: we don't need to interpolate or relax
    ! anything higher in the atmosphere than the next pressure level
    ! above the domain-wide lowest surface pressure:
    call reduce(grid,pmin,grid%PINT,grid%LM,REDUCE_MIN,&
                ims,ime,jms,jme,kms,kme)

    ! FIXME: DON'T HANDLE ANYTHING ABOVE PMIN
    ! NOTE: MUST HANDLE TWO LEVELS ABOVE

    ! Step 1: calculate Tv, Q and Z on pressure levels using the same
    ! method as the NCEP Unified Post:
    call calculate_3D(grid,presTv,presZ,ground_mask,ground_level, &
         IDS,IDE,JDS,JDE,KDS,KDE, &
         IMS,IME,JMS,JME,KMS,KME, &
         IPS,IPE,JPS,JPE,KPS,KPE)

    ! Step 2: smooth Tv through an overrelaxation method:

    ! Modify the relax mask so that the outermost three rows and
    ! columns are always relaxed.  This is needed to overcome bad
    ! values fed in from the parent every timestep.  Setting the mask
    ! to true on the boundaries of the domain prevent them from being
    ! used as boundaries of the overrelaxation.  

    ! Some of the reasons for boundary issues: The parent and nest
    ! terrain will never match because the nest terrain is smoothed on
    ! the boundary, and the parent is not.  Also, the user may have
    ! set a different terrain data source for different domains, in
    ! which case you'll get an even worse mismatch.  Every time the
    ! nest moves, terrain changes on the leading and trailing edges of
    ! the nest.  That causes huge shocks when there are high mountains
    ! near the boundaries.  If you do a plot of 500mbar geopotential
    ! height, it looks like a piece of jello shaking every time the
    ! nest moves.  Also, there is some weirdness on the lateral
    ! boundaries of the outermost domain due to the mismatch between
    ! GFS terrain (which has its higher spectral components discarded)
    ! and the smoothed regional terrain.

    relaxmask=1

    ! Now loop over all vertical levels and relax them:
    do ipres=npres,1,-1
       ! In the inner regions (all but outermost row & col) set the
       ! relaxmask to the ground_mask:
       need_to_relax=0
       do j=max(jps,jds+1),min(jde-1,jpe)
          do i=max(ips,ids+1),min(ide-1,ipe)
             relaxmask(i,j)=ground_mask(i,j,ipres)
             if(relaxmask(i,j)/=0) need_to_relax=1
          enddo
       enddo

       ! If we do not need to relax any points, we are done.
       call max_integer(grid,need_to_relax)
       if(need_to_relax==0) then
! 38       format('end mslp relax loop at ',I0)
!          print 38,ipres
          exit
       endif

       ! Store Tv in relaxwork:
       do j=jps,jpe
          do i=ips,ipe
             relaxwork(i,j)=presTv(i,j,ipres)
          enddo
       enddo

       ! Overrelax:
       call relax4e(relaxwork,relaxmask,0.7,100, &
            IDS,IDE,JDS,JDE, &
            IMS,IME,JMS,JME, &
            IPS,IPE,JPS,JPE)

       ! Copy back the modified relaxation mask
       do j=jps,jpe
          do i=ips,ipe
             ground_mask(i,j,ipres)=relaxmask(i,j)
          enddo
       enddo

       ! Copy the relaxed values back to Tv:
       do j=jps,jpe
          do i=ips,ipe
             presTv(i,j,ipres)=relaxwork(i,j)
          enddo
       enddo
    end do

    ! Step 3: Solve for Z on interface levels (pressure space
    ! interface levels) using the hydrostatic equation.  Once Z=0 is
    ! reached, solve for Pmsl.
    call calculate_interP(presTv,presZ,grid%Z,Pmsl,grid%PINT, &
         grid%T(:,:,kde), grid%Q(:,:,kde), &
         ground_level,ground_mask,grid%fis, &
         IDS,IDE,JDS,JDE,KDS,KDE, &
         IMS,IME,JMS,JME,KMS,KME, &
         IPS,IPE,JPS,JPE,KPS,KPE)

    ! Copy the MSLP values back to the grid:
    do j=jps,jpe
       do i=ips,ipe
          grid%membrane_MSLP(i,j)=Pmsl(i,j)
       enddo
    enddo

    ! Smooth the membrane_mslp values:
    call smoothMSLP(grid,1,relaxmask,relaxwork, &
         IDS,IDE,JDS,JDE,KDS,KDE, &
         IMS,IME,JMS,JME,KMS,KME, &
         IPS,IPE,JPS,JPE,KPS,KPE)

    if(size(grid%p850z)>1) then
       ! Copy 700 and 850 mbar heights to their arrays:
       do j=jps,jpe
          do i=ips,ipe
             grid%p850z(i,j)=presZ(i,j,k850)
             grid%p700z(i,j)=presZ(i,j,k700)
          enddo
       enddo
    endif

  end subroutine quasipost_impl

  subroutine calculate_3D(grid,presTv,presZ,ground_mask,ground_level, &
       IDS,IDE,JDS,JDE,KDS,KDE, &
       IMS,IME,JMS,JME,KMS,KME, &
       IPS,IPE,JPS,JPE,KPS,KPE)
    use MODULE_SOLVER_INTERNAL_STATE, only : SOLVER_INTERNAL_STATE
    use mpi, only: MPI_COMM_WORLD,MPI_Abort
    use module_exchange, only: halo_exch
    implicit none

    type(SOLVER_INTERNAL_STATE), intent(inout) :: grid

    integer, intent(in) :: IDS,IDE,JDS,JDE,KDS,KDE
    integer, intent(in) :: IMS,IME,JMS,JME,KMS,KME
    integer, intent(in) :: IPS,IPE,JPS,JPE,KPS,KPE

    real, intent(inout) :: presTv(ims:ime,jms:jme,npres)
    real, intent(inout) :: presZ(ims:ime,jms:jme,npres)

    integer, intent(inout) :: ground_mask(ims:ime,jms:jme,npres)
    integer, intent(inout) :: ground_level(ims:ime,jms:jme)

    integer :: Tkdest(ims:ime,jms:jme), Zkdest(ims:ime,jms:jme), Zbottom(ims:ime,jms:jme)
    integer :: i,j,ks,kd,k
    real :: weight, TL,QL,PL, tempT, RHL, TVRL, TVRBLO, TBLO,QBLO

    integer,target, dimension(ims:ime,jms:jme) :: ks850,ks700,ks500
    real, target,dimension(ims:ime,jms:jme) :: dummy1,dummy2
    integer, pointer, dimension(:,:) :: ksX
    integer :: nanfound,ierr,iref,jref,b
    real, pointer, dimension(:,:) :: preswind,presrv,presu,presv

    real :: Pmass(ims:ime,jms:jme,kds:kde)
    real :: numsum,densum,modelP1,modelP2,pdiff,presQ,presT,ZL,QSAT, U1, V1, U2, V2, dudy1,dvdx1, dudy2,dvdx2
    character*255 :: message
    logical :: wantuv
    !call quasipost_message('calculate3D')
    call HALO_EXCH(grid%u10,1,grid%U,grid%LM,2,2)
    call HALO_EXCH(grid%v10,1,grid%V,grid%LM,2,2)

    ! ks: k in source (model level) array
    ! kd: k in destination (pressure level) array
    ground_level=0
    ground_mask=0
    Zkdest=1
    Tkdest=1

    ks850=0
    ks700=0
    ks500=0

    ! Interpolate geopotential height to post_stdpres pressure levels
    ! and create a temporary array with non-hydrostatic pressure
    ! (PINT) on model mass points:
    do ks=kds+1,kde
       do j=jps,jpe
          iZ: do i=ips,ipe
             Pmass(i,j,ks)=sqrt(grid%PINT(i,j,ks)*grid%PINT(i,j,ks+1))
          enddo iZ
       enddo
    enddo

    ! Interpolate temperature and specific humidity to post_stdpres
    ! pressure levels:
    do ks=kds+1,kde-1
       do j=jps,jpe
          iTQ: do i=ips,ipe
             kd=Tkdest(i,j)
             if(kd<=npres) then
                innerTQ: do while(kd<=npres)
                   if(.not.(post_stdpres(kd)<=Pmass(i,j,ks+1) &
                        .and. post_stdpres(kd)>=Pmass(i,j,ks))) then
                      cycle iTQ
                   endif
                   weight=log(post_stdpres(kd)/Pmass(i,j,ks))/log(Pmass(i,j,ks+1)/Pmass(i,j,ks))

                   presZ(i,j,kd)=weight*grid%Z(i,j,ks+1) + (1.-weight)*grid%Z(i,j,ks)

                   presT=weight*grid%T(i,j,ks+1) + (1.-weight)*grid%T(i,j,ks)
                   presQ=weight*grid%Q(i,j,ks+1) + (1.-weight)*grid%Q(i,j,ks)
                   presTv(i,j,kd)=presT*(1.+post_D608*presQ)

                   if(kd==k850) then
                      ks850(i,j)=ks
                   elseif(kd==k700) then
                      ks700(i,j)=ks
                   elseif(kd==k500) then
                      ks500(i,j)=ks
                   endif

103                format('interp ks=',I0,' kd=',I0,' presT(i=',I0,',j=',I0,',kd)=',F0.3, &
                        ' between T(i,j,ks-1)=',F0.3,' and T(i,j,ks)=', &
                        F0.3,' using weight=',F0.3)
                   !write(message,103) ks,kd,i,j,presT,grid%T(i,j,ks-1),grid%T(i,j,ks),weight
                   !call quasipost_message(message)
104                format(' Pmass(i,j,ks)=',F0.3,' Pmass(i,j,ks-1)=',F0.3,' post_stdpres(kd)=',F0.3)
                   !write(message,104) Pmass(i,j,ks),Pmass(i,j,ks-1),post_stdpres(kd)
                   !call quasipost_message(message)
                   if(weight<0 .or. weight>1) then
                      write(0,*) 'Bad weight: ',weight
                      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
                   endif
                   kd=kd+1
                   Tkdest(i,j)=kd
                   Zkdest(i,j)=kd
                end do innerTQ
             end if
          end do iTQ
       end do
    end do

   ! Interpolate to regions between the middle of the lowest mass
   ! level and the bottom of the atmosphere:
   do j=jps,jpe
      iTQ2: do i=ips,ipe
         kd=Zkdest(i,j)
         if(kd<=npres) then
            do while(kd<=npres)
               if(.not.(post_stdpres(kd)<=grid%PINT(i,j,kde+1) &
                    .and. post_stdpres(kd)>=Pmass(i,j,kde))) then
                  cycle iTQ2
               endif

               presT=grid%T(i,j,kde)
               presQ=grid%Q(i,j,kde)
               presTv(i,j,kd)=presT*(1.+post_D608*presQ)

               weight=log(post_stdpres(kd)/Pmass(i,j,kde))/log(grid%PINT(i,j,kde+1)/Pmass(i,j,kde))
               presZ(i,j,kd)=(1.-weight)*grid%Z(i,j,kde)+weight*grid%fis(i,j)/post_g

               kd=kd+1
               Tkdest(i,j)=kd
               Zkdest(i,j)=kd
            end do
         end if
      end do iTQ2
   end do

1234 format('grid size(',A,') = ',I0)
   !print 1234,   'grid%p700rv',size(grid%p700rv)
   !print 1234,   'grid%p700u',size(grid%p700u)

   ! do I need to calc. presu & presv? Yes.
   wantuv=.true. ! WAS: (grid%vortex_tracker == 7)

   ifwind: if(size(grid%p700rv)>1 .or. size(grid%p700u)>1) then
    ! Interpolate wind to H points on pressure levels, calculating
    ! horizontal wind vector magnitude and vertical component of
    ! vorticity.  Interpolate only to 700 and 850 mbar, except for U &
    ! V which are also interpolated to 500mbar.
    nullify(presu)
    nullify(presv)
    windloop: do k=0,2
       if(k==0) then
          ! Only need wind components at 500 mbar
          kd=k500
          ksX=>ks500
          preswind=>dummy1
          presrv=>dummy2
          if(wantuv) then
             presu=>grid%p500u
             presv=>grid%p500v
          endif
       elseif(k==1) then
          ksX=>ks700
          preswind=>grid%p700wind
          presrv=>grid%p700rv
          kd=k700
          if(wantuv) then
             presu=>grid%p700u
             presv=>grid%p700v
          endif
       elseif(k==2) then
          ksX=>ks850
          kd=k850
          preswind=>grid%p850wind
          presrv=>grid%p850rv 
          if(wantuv) then
             presu=>grid%p850u
             presv=>grid%p850v
          endif
       endif

      ! No wind on boundaries:
      if(jps<=jds) then
         do i=ips,ipe
            preswind(i,jds)=0
            presrv(i,jds)=0
         enddo
         if(wantuv) then
            do i=ips,ipe
               presu(i,jds)=0
               presv(i,jds)=0
            enddo
         endif
      endif
      if(jpe>=jde-1) then
         do i=ips,ipe
            preswind(i,jde-1)=0
            presrv(i,jde-1)=0
         enddo
         if(wantuv) then
            do i=ips,ipe
               presu(i,jde-1)=0
               presv(i,jde-1)=0
            enddo
         endif
      endif
      if(ips<=ids) then
         do j=jps,jpe
            preswind(ids,j)=0
            presrv(ids,j)=0
         enddo
         if(wantuv) then
            do j=jps,jpe
               presu(ids,j)=0
               presv(ids,j)=0
            enddo
         endif
      endif
      if(ipe>=ide-1) then
         do j=jps,jpe
            preswind(ide-1,j)=0
            presrv(ide-1,j)=0
         enddo
         if(wantuv) then
            do j=jps,jpe
               presu(ide-1,j)=0
               presv(ide-1,j)=0
            enddo
         endif
      endif

      ! Interpolate winds:
      do j=max(jps,jds+1),min(jde-2,jpe)
         do i=max(ips,ids+1),min(ide-2,ipe)
            ks=ksX(i,j)
            if(ks>1) then
               ! Interpolate between mass levels:
               weight=log(post_stdpres(kd)/Pmass(i,j,ks))/log(Pmass(i,j,ks+1)/Pmass(i,j,ks))

               U1=0.25*(grid%u(i,j-1,ks) + grid%u(i,j,ks) + grid%u(i-1,j,ks) + grid%u(i-1,j-1,ks))
               V1=0.25*(grid%v(i,j-1,ks) + grid%v(i,j,ks) + grid%v(i-1,j,ks) + grid%v(i-1,j-1,ks))
               U2=0.25*(grid%u(i,j-1,ks+1) + grid%u(i,j,ks+1) + grid%u(i-1,j,ks+1) + grid%u(i-1,j-1,ks+1))
               V2=0.25*(grid%v(i,j-1,ks+1) + grid%v(i,j,ks+1) + grid%v(i-1,j,ks+1) + grid%v(i-1,j-1,ks+1))
               
               dvdx1 = (grid%v(i+1,j,ks)-grid%v(i-1,j,ks))/(2.*grid%DXH(j))
               dudy1 = (grid%u(i,j+1,ks)-grid%u(i,j-1,ks))/(2.*grid%DYH)
               dvdx2 = (grid%v(i+1,j,ks+1)-grid%v(i-1,j,ks+1))/(2.*grid%DXH(j))
               dudy2 = (grid%u(i,j+1,ks+1)-grid%u(i,j-1,ks+1))/(2.*grid%DYH)

               if(wantuv) then
                  presu(i,j)=weight*u2+(1.-weight)*u1
                  presv(i,j)=weight*v2+(1.-weight)*v1
               endif
               preswind(i,j)=weight*sqrt(u2*u2+v2*v2) + (1.-weight)*sqrt(u1*u1+v1*v1)
               presrv(i,j)=(dvdx2-dudy2)*weight + (dvdx1-dudy1)*(1.-weight)
            elseif(post_stdpres(kd)>=Pmass(i,j,kds)) then
               ! At and below lowest mass level, use lowest model level winds
               ks=1
               U1=0.25*(grid%u(i,j-1,ks) + grid%u(i-1,j,ks) + grid%u(i,j,ks) + grid%u(i-1,j-1,ks))
               V1=0.25*(grid%v(i,j-1,ks) + grid%v(i-1,j,ks) + grid%v(i,j,ks) + grid%v(i-1,j-1,ks))
               
               dvdx1 = (grid%v(i+1,j,ks)-grid%v(i-1,j,ks))/(2.*grid%DXH(j))
               dudy1 = (grid%u(i,j+1,ks)-grid%u(i,j-1,ks))/(2.*grid%DYH)

               preswind(i,j)=sqrt(u1*u1 + v1*v1)
               presrv(i,j)=dvdx1-dudy1
               if(wantuv) then
                  presu(i,j)=u1
                  presv(i,j)=v1
               endif
            endif
         end do
      end do

      ! Copy wind-related fields to boundary points.
      do b=1,3
         if    (b==1) then ; i=ids   ; iref=ids+1
         elseif(b==2) then ; i=ide-1 ; iref=ide-2
         elseif(b==3) then ; i=ide   ; iref=ide-2
         endif
         if(i<=ipe .and. i>=ips) then
            do j=max(jps,jds+1),min(jde-2,jpe)
               if(wantuv) then
                  presu(i,j)=presu(iref,j)
                  presv(i,j)=presv(iref,j)
               endif
               presrv(i,j)=presrv(iref,j)
               preswind(i,j)=preswind(iref,j)
            enddo
         endif
      enddo

      do b=1,3
         if    (b==1) then ; j=jds   ; jref=jds+1
         elseif(b==2) then ; j=jde-1 ; jref=jde-2
         elseif(b==3) then ; j=jde   ; jref=jde-2
         endif
         if(j<=jpe .and. j>=jps) then
            do i=ips,ipe
               if(wantuv) then
                  presu(i,j)=presu(i,jref)
                  presv(i,j)=presv(i,jref)
               endif
               presrv(i,j)=presrv(i,jref)
               preswind(i,j)=preswind(i,jref)
            enddo
         endif
      enddo
   enddo windloop

   ! Calculate 10m wind magnitude and vorticity
   ! NOTE: u10 and v10 are already on H points
   nanfound=0
   do j=max(jps,jds+1),min(jpe,jde-1)
      do i=max(ips,ids+1),min(ipe,ide-1)
         grid%m10wind(i,j)=sqrt(grid%u10(i,j)*grid%u10(i,j) + grid%v10(i,j)*grid%v10(i,j))
         dvdx1 = (grid%v10(i+1,j)-grid%v10(i-1,j)) / (2*grid%DXH(j))
         dudy1 = (grid%u10(i,j+1)-grid%u10(i,j-1)) / (2*grid%DYH)
         grid%m10rv(i,j) = dvdx1 - dudy1
         if(grid%m10rv(i,j) == grid%m10rv(i,j)) then
            continue
         else
3088        format('NaN m10rv at i=',I0,' j=',I0,' dx=',F0.3,' dy=',F0.3)
            write(0,3088) i,j,grid%DXH(j),grid%DYH
            call MPI_Abort(MPI_COMM_WORLD,1,ierr)
3089        format('NaN m10rv at i=',I0,' j=',I0,': dvdx1=',F0.5,' dudy=',F0.5)
            write(0,3089) i,j,dvdx1,dudy1
            call MPI_Abort(MPI_COMM_WORLD,1,ierr)
            nanfound=1
         endif
      enddo
   enddo

   call max_integer(grid,nanfound)
   if(nanfound/=0) then
      write(0,*) 'ERROR: NaN m10rv seen; aborting.'
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
   endif
  endif ifwind

    do j=jps,jpe
       do i=ips,ipe
          ground_level(i,j)=min(Zkdest(i,j),Tkdest(i,j))
       enddo
    enddo

    do kd=1,npres
       do j=jps,jpe
          do i=ips,ipe
             if(kd>=ground_level(i,j)) then
                ground_mask(i,j,kd) = 1
             else
                ground_mask(i,j,kd) = 0
             endif
          enddo
       enddo
    enddo

    if(50>=ips .and. 50<=ipe .and. 50>=jps .and. 50<=jpe) then
       !print *,'Z(50,50) = ',grid%Z(50,50,:)
    endif

    ! Extrapolate below-ground temperature but not height.  Fill in
    ! badheight for height below ground.  
    jloop2: do j=jps,jpe
       iloop2: do i=ips,ipe
          if(ground_level(i,j)>npres) then
301          format('Extrap: i=',I0,' j=',I0,' NO EXTRAP: ground at ',I0)
             !write(message,301) i,j,ground_level(i,j)
             !call quasipost_message(message)
             cycle iloop2
          else
302          format('Extrap: i=',I0,' j=',I0,' extrap from ',F0.3,' ground at ',I0)
             !write(message,302) i,j,post_stdpres(ground_level(i,j)),ground_level(i,j)
             !call quasipost_message(message)
          endif
          kloop2: do kd=ground_level(i,j),npres
             ! Extrapolate first guess below-ground values using the
             ! exact same method used in the post.  Even the constants
             ! are copied from the post:
             PL=grid%PINT(I,J,kde-1)
             ZL=0.5*(grid%Z(I,J,kde-1)+grid%Z(I,J,kde))
             TL=0.5*(grid%T(I,J,kde-1)+grid%T(I,J,kde))
             QL=0.5*(grid%Q(I,J,kde-1)+grid%Q(I,J,kde))
             QSAT=post_PQ0/PL*EXP(post_A2*(TL-post_A3)/(TL-post_A4))
             !
             RHL=QL/QSAT
             !
             IF(RHL.GT.1.)THEN
                RHL=1.
                QL =RHL*QSAT
             ENDIF
             !
             IF(RHL.LT.post_RHmin)THEN
                RHL=post_RHmin
                QL =RHL*QSAT
             ENDIF
             !
             TVRL  =TL*(1.+post_D608*QL)
             TVRBLO=TVRL*(post_stdpres(kd)/PL)**post_RGAMOG
             TBLO  =TVRBLO/(1.+post_D608*QL)

             !QSAT=post_PQ0/post_stdpres(kd)*EXP(post_A2*(TBLO-post_A3)/(TBLO-post_A4))

             !QBLO =RHL*QSAT
             !presQ(i,j,kd)=AMAX1(post_smallQ,QBLO)

             presTv(i,j,kd)=TBLO

             ! Extrapolated virtual temperature:
             !presTv(i,j,kd)=TBLO*(1.+post_D608*QBLO)

             ! extrapolated temperature, with virtual part removed using extrapolated specific humidity:
             !presTv(i,j,kd)=TVRBLO/(1.+post_D608*QBLO)

             ! Below-ground Z is recalcluated after smoothing Tv.  We
             ! only fill in badval here:
             presZ(i,j,kd)=badheight

303          format('Extrap i=',I0,' j=',I0,' kd=',I0,' presTv=',F0.3,' presZ=',F0.3)
304          format('   TL=',F0.3,' QL=',F0.3,' ZL=',F0.3,' QSAT=',F0.3)
305          format('   TVRL=',F0.3,' TVRBLO=',F0.3,' TBLO=',F0.3,' RHL=',F0.3)
             !write(0,303) i,j,kd,presTv(i,j,kd),presZ(i,j,kd)
             !write(0,304) TL,QL,ZL,QSAT
             !write(0,305) TVRL,TVRBLO,TBLO,RHL
          enddo kloop2
       enddo iloop2
    enddo jloop2
  end subroutine calculate_3D

  subroutine calculate_interP( &
       presTv,presZ,modelZ,Pmsl,PINT,T1,Q1, &
       ground_level,ground_mask,fis, &
       IDS,IDE,JDS,JDE,KDS,KDE, &
       IMS,IME,JMS,JME,KMS,KME, &
       IPS,IPE,JPS,JPE,KPS,KPE)

    USE MPI, only: MPI_Abort,MPI_COMM_WORLD
    implicit none

    integer, intent(in) :: IDS,IDE,JDS,JDE,KDS,KDE
    integer, intent(in) :: IMS,IME,JMS,JME,KMS,KME
    integer, intent(in) :: IPS,IPE,JPS,JPE,KPS,KPE

    real, intent(in) :: PINT(ims:ime,jms:jme,kms:kme+1), modelZ(ims:ime,jms:jme,kms:kme)
    real, intent(in) :: T1(ims:ime,jms:jme,1)
    real, intent(in) :: Q1(ims:ime,jms:jme,1)

    real, intent(in) :: fis(ims:ime,jms:jme)
    real, intent(out) :: Pmsl(ims:ime,jms:jme)
    real, intent(inout) :: presTv(ims:ime,jms:jme,npres)
    real, intent(inout) :: presZ(ims:ime,jms:jme,npres)

    integer, intent(inout) :: ground_mask(ims:ime,jms:jme,npres)
    integer, intent(inout) :: ground_level(ims:ime,jms:jme)

    real :: Z,midTv,dZ,newZ,P,newP,TVRT,TLYR,DIS,oa,slope
    integer :: kp,ip,i,j,ierr

    ! What this code does:

    ! For every point where the surface is above Z=0, we start from
    ! the lowest above-ground pressure and integrate the hydrostatic
    ! equation downward to find P at Z=0.

    ! For points where the surface Z<=0 (surface is at or below sea
    ! level), we interpolate to get P at Z=0.


    ! STEP 1: extrapolate below-ground values
    do j=jps,jpe
       iloop: do i=ips,ipe
          !          nearground: if(modelZ(i,j,1)<50.0) then
          !             Pmsl(i,j)=pint1(i,j,1)
          !                method(i,j)=-30
          !          else
          if(ground_level(i,j)<npres+1) then
             kp=ground_level(i,j)-1
101          format('i=',I0,' j=',I0,' kp=',I0,' ground level =',I0)
             !write(0,101) i,j,kp,ground_level(i,j)
             if(kp<1) then
108             format('At ground_level(',I0,',',I0,')=',I0,', lowest model surface pressure is lower than second lowest standard pressure level.')
                write(0,108) i,j,ground_level(i,j)
                call MPI_Abort(MPI_COMM_WORLD,1,ierr)
             endif
             ! Ground is below lowest model level
             !newZ=presZ(i,j,kp)
             !newP=post_stdpres(kp)
             newZ=fis(i,j)/post_G
             newP=pint(i,j,kde+1)
             do ip=kp,npres-1
                P=newP
                Z=newZ
                !                midTv=0.5*(presTv(i,j,ip)+presTv(i,j,ip+1))
                midTv=presTv(i,j,ip+1)
                newP=post_stdpres(ip+1)
                dZ=post_Rd*midTv*alog(P/newP)/post_g
102             format('  make some Z at ip=',I0,': P=',F0.3,' newP=',F0.3)
1021            format('  Z=',F0.3,' midTv=',F0.3,' dZ=',F0.3)
                !write(0,102) ip,P,newP
                !write(0,1021) Z,midTv,dZ
                if(dZ>=0.) then
80881              format('ABORT: Bad dZ at ',I0,',',I0,': dZ=',F0.3,'>=0 Tv=',F0.3,' P=',F0.3,' newP=',F0.3,' kp=',I0)
                   write(0,80881) i,j,dZ,midTv,P,newP,kp
                   call MPI_Abort(MPI_COMM_WORLD,1,ierr)
                endif
                newZ=Z+dZ
                presZ(i,j,ip+1)=newZ
                if(newZ<=0) then
                   ! interpolate between Z and newZ
1022               format('  extrap using ',F0.3,'/exp(-',F0.3,'*',F0.3,'/(',F0.3,'*',F0.3,'))')
                   !write(0,1022) P,Z,post_G,post_RD,presTV(i,j,ip)


                   !Pmsl(i,j)=P/exp(-Z*post_G/(post_RD*presTv(i,j,ip)))
                   Pmsl(i,j)=(Z*newP-newZ*P)/(-dZ)
10221              format('  result: ',F0.3)
                   !write(0,10221) Pmsl(i,j)
!                   method(i,j)=ip
                   cycle iloop
                endif
             enddo
          endif
          ! If we get here, then Z=0 is below the lowest standard
          ! pressure level and we must extrapolate.

          !               if(pint1(i,j,1)>post_stdpres(npres) .and. modelZ(i,j,1)>0.)then
          !                  ! Model surface pressure is a higher pressure than the
          !                  ! highest standard pressure level.  Use the model
          !                  ! fields to extrapolate.
          !                  TVRT=T1(I,J,1)*(post_H1+post_D608*Q1(I,J,1))
          !                  !DIS=modelZ(I,J,kde-1)-modelZ(I,J,1)+0.5*modelZ(I,J,kde-1) ???
          !                  DIS=0.5*(modelZ(I,J,kde-1)+modelZ(I,J,1))
          !                  TLYR=TVRT-DIS*post_SLOPE*post_G*0.5
          !                  Pmsl(I,J)=PINT(I,J,1)*EXP((modelZ(I,J,1))*post_G/(post_RD*TLYR))
          ! ! 1023            format('  use model: TVRT=',F0.3,' DIS=',F0.3,' TLYR=',F0.3,' Pmsl=',F0.3)
          ! ! 1024            format('     result: ',F0.3,'*EXP(',F0.3,'/(',F0.3,'*',F0.3'))')
          ! !                 write(0,1023) TVRT,DIS,TLYR,Pmsl(i,j)
          ! !                 write(0,1024) PINT(I,J,1),modelZ(I,J,kde-1),post_RD,TLYR
          !                method(i,j)=-20
          !               ELSE
          ! Highest pressure level (post_stdpres(1)) has a
          ! higher pressure than the model surface pressure, so
          ! extrapolate using the pressure level values.
1025      format('  use npres: TLYR=',F0.3,' Pmsl=',F0.3)
1026      format('     result: ',F0.3,'/EXP(-',F0.3,'*',F0.3,'/(',F0.3,'*',F0.3,'))')
          TLYR=presTv(I,J,npres)-presZ(I,J,npres)*post_SLOPE*post_G*0.5
          Pmsl(I,J)=post_stdpres(npres)/EXP(-presZ(I,J,npres)*post_G/(post_RD*TLYR))
          !oa=0.5*post_SLOPE*post_g*presZ(i,j,npres)/TLYR
          !Pmsl(i,j)=post_stdpres(npres)*(1.-oa)**old_coef2
          !write(0,1025) TLYR,Pmsl(I,J)
          !write(0,1026) post_stdpres(npres),presZ(I,J,npres),post_G,post_RD,TLYR
!          method(i,j)=-10
          !             END IF
          !          endif nearground
       enddo iloop
    enddo
  end subroutine calculate_interP

  subroutine smoothMSLP(grid,iterations,relaxmask,relaxwork,  &
       IDS,IDE,JDS,JDE,KDS,KDE, &
       IMS,IME,JMS,JME,KMS,KME, &
       IPS,IPE,JPS,JPE,KPS,KPE)
    use module_relax4e
    use MODULE_SOLVER_INTERNAL_STATE, only : SOLVER_INTERNAL_STATE
    implicit none
    type(solver_internal_state), intent(inout) :: grid
    integer, intent(in) :: iterations
    integer :: relaxmask(ims:ime,jms:jme)
    real :: relaxwork(ims:ime,jms:jme)

    integer :: IDS,IDE,JDS,JDE,KDS,KDE
    integer :: IMS,IME,JMS,JME,KMS,KME
    integer :: IPS,IPE,JPS,JPE,KPS,KPE
    integer :: i,j

    do j=jps,jpe
       do i=ips,ipe
          relaxmask(i,j)=1
          relaxwork(i,j)=grid%membrane_mslp(i,j)
       enddo
    enddo

    call relax4e(relaxwork,relaxmask,0.5,iterations, &
         IDS,IDE,JDS,JDE, &
         IMS,IME,JMS,JME, &
         IPS,IPE,JPS,JPE)

    do j=jps,jpe
       do i=ips,ipe
          grid%membrane_mslp(i,j)=relaxwork(i,j)
       enddo
    enddo

  end subroutine smoothMSLP

end module module_quasipost
