module module_REDUCTION
  implicit none
  private
  integer, public, parameter :: FIND_MAX=1, FIND_MIN=2
  integer, public, parameter :: REDUCE_MAX=1, REDUCE_MIN=2, REDUCE_ADD=3, REDUCE_ADD_DOUBLE=4

  public :: reduce,find,max_integer,max_real,minloc_real,maxloc_real

  interface reduce
     module procedure reduce_3d_real2double
     module procedure reduce_2d_real2double
     module procedure reduce_3d_real
     module procedure reduce_2d_real
  end interface reduce

  interface find
     module procedure find_3d_real
     module procedure find_2d_real
     module procedure find_3d_integer
     module procedure find_2d_integer
  end interface find

  interface minloc_real
     module procedure minloc7_real ! value, full 3d location and grid index
     module procedure minloc4_real ! value and grid index
  end interface

  interface maxloc_real
     module procedure maxloc7_real ! value, full 3d location and grid index
     module procedure maxloc4_real ! value and grid index
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Array reductions without find
  subroutine reduce_3d_real2double(sis,dout,array,nk,method,&
                                   ims,ime,jms,jme,kms,kme)
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    integer, intent(in) :: ims,ime,jms,jme,kms,kme,nk
    real, intent(in) :: array(ims:ime,jms:jme,kms:kme)
    integer, intent(in) :: method
    real :: rout
    double precision,intent(out) :: dout
    call reduce_real_impl(array,rout,dout,&
         ims,ime,jms,jme,kms,kme, &
         sis%its,sis%ite,sis%jts,sis%jte,1,nk,&
         method,sis%MPI_COMM_COMP,sis%MYPE)
    if(method/=REDUCE_ADD_DOUBLE) &
         dout=rout
  end subroutine reduce_3d_real2double

  subroutine reduce_3d_real(sis,rout,array,nk,method,&
                            ims,ime,jms,jme,kms,kme)
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    integer, intent(in) :: ims,ime,jms,jme,kms,kme,nk
    real, intent(in) :: array(ims:ime,jms:jme,kms:kme)
    integer, intent(in) :: method
    real, intent(inout) :: rout
    double precision :: dout
    rout=0
    dout=0
    call reduce_real_impl(array,rout,dout,&
         ims,ime,jms,jme,kms,kme, &
         sis%its,sis%ite,sis%jts,sis%jte,1,nk,&
         method,sis%MPI_COMM_COMP,sis%MYPE)
    if(method==REDUCE_ADD_DOUBLE) &
         rout=dout
  end subroutine reduce_3d_real

  subroutine reduce_2d_real2double(sis,dout,array,method,&
                              ims,ime,jms,jme)
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    integer, intent(in) :: ims,ime,jms,jme
    real, intent(in) :: array(ims:ime,jms:jme)
    integer, intent(in) :: method
    real :: rout
    double precision,intent(out) :: dout
    call reduce_real_impl(array,rout,dout,&
         ims,ime,jms,jme,1,1, &
         sis%its,sis%ite,sis%jts,sis%jte,1,1,&
         method,sis%MPI_COMM_COMP,sis%MYPE)
    if(method/=REDUCE_ADD_DOUBLE) &
         dout=rout
  end subroutine reduce_2d_real2double

  subroutine reduce_2d_real(sis,rout,array,method,&
                            ims,ime,jms,jme)
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    integer, intent(in) :: ims,ime,jms,jme
    real, intent(in) :: array(ims:ime,jms:jme)
    integer, intent(in) :: method
    real, intent(out) :: rout
    double precision :: dout
    call reduce_real_impl(array,rout,dout,&
         ims,ime,jms,jme,1,1, &
         sis%its,sis%ite,sis%jts,sis%jte,1,1,&
         method,sis%MPI_COMM_COMP,sis%MYPE)
    if(method==REDUCE_ADD_DOUBLE) &
         rout=dout
  end subroutine reduce_2d_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Array reductions with find

  subroutine find_3d_real(sis,rfnd,iloc,jloc,kloc,rankloc,array,nk,method,&
                          ims,ime,jms,jme,kms,kme)
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    integer, intent(in) :: ims,ime,jms,jme,kms,kme,nk
    real, intent(in) :: array(ims:ime,jms:jme,kms:kme)
    integer, intent(in) :: method
    real, intent(out) :: rfnd
    integer, intent(out) :: iloc,jloc,kloc,rankloc
    call reduce_find_real_impl(array,rfnd,iloc,jloc,kloc,rankloc,&
         ims,ime,jms,jme,kms,kme, &
         sis%its,sis%ite,sis%jts,sis%jte,1,nk,&
         method,sis%MPI_COMM_COMP,sis%MYPE)
  end subroutine find_3d_real

  subroutine find_2d_real(sis,rfnd,iloc,jloc,rankloc,array,method,&
                          ims,ime,jms,jme)
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    integer, intent(in) :: ims,ime,jms,jme
    real, intent(in) :: array(ims:ime,jms:jme)
    integer, intent(in) :: method
    real, intent(out) :: rfnd
    integer, intent(out) :: iloc,jloc,rankloc
    integer :: kloc
    call reduce_find_real_impl(array,rfnd,iloc,jloc,kloc,rankloc,&
         ims,ime,jms,jme,1,1, &
         sis%its,sis%ite,sis%jts,sis%jte,1,1,&
         method,sis%MPI_COMM_COMP,sis%MYPE)
  end subroutine find_2d_real

  subroutine find_3d_integer(sis,ifnd,iloc,jloc,kloc,rankloc,array,nk,method,&
                          ims,ime,jms,jme,kms,kme)
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    integer, intent(in) :: ims,ime,jms,jme,kms,kme,nk
    integer, intent(in) :: array(ims:ime,jms:jme,kms:kme)
    integer, intent(in) :: method
    integer, intent(out) :: ifnd
    integer, intent(out) :: iloc,jloc,kloc,rankloc
    call reduce_find_integer_impl(array,ifnd,iloc,jloc,kloc,rankloc,&
         ims,ime,jms,jme,kms,kme, &
         sis%its,sis%ite,sis%jts,sis%jte,1,nk,&
         method,sis%MPI_COMM_COMP,sis%MYPE)
  end subroutine find_3d_integer

  subroutine find_2d_integer(sis,ifnd,iloc,jloc,rankloc,array,method,&
                          ims,ime,jms,jme)
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    integer, intent(in) :: ims,ime,jms,jme
    integer, intent(in) :: array(ims:ime,jms:jme)
    integer, intent(in) :: method
    integer, intent(out) :: ifnd
    integer, intent(out) :: iloc,jloc,rankloc
    integer :: kloc
    call reduce_find_integer_impl(array,ifnd,iloc,jloc,kloc,rankloc,&
         ims,ime,jms,jme,1,1, &
         sis%its,sis%ite,sis%jts,sis%jte,1,1,&
         method,sis%MPI_COMM_COMP,sis%MYPE)
  end subroutine find_2d_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Scalar reductions
  subroutine max_real(sis,val)
    use mpi
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    integer ierr
    real maxv,val
    call MPI_Allreduce(val,maxv,1,MPI_REAL,MPI_MAX,&
                       sis%MPI_COMM_COMP,ierr)
    val=maxv
  end subroutine max_real

  subroutine max_integer(sis,val)
    use mpi
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    integer ierr
    integer maxv,val
    call MPI_Allreduce(val,maxv,1,MPI_REAL,MPI_MAX,&
                       sis%MPI_COMM_COMP,ierr)
    val=maxv
  end subroutine max_integer

  subroutine minloc7_real(sis,val,lat,lon,z,idex,jdex)
    use mpi
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    REAL,intent(inout) :: val, lat, lon, z
    integer, intent(inout) :: idex,jdex
    INTEGER ierr, mrank
    REAL inreduce(2), outreduce(2), bcast(5)

    inreduce=(/ val, real(sis%MYPE) /)
    call MPI_Allreduce(inreduce,outreduce,1,MPI_2REAL,MPI_MINLOC,&
         sis%MPI_COMM_COMP,ierr)
    val=outreduce(1)
    mrank=outreduce(2)
    bcast=(/ lat,lon,z,real(idex),real(jdex) /)
    call MPI_Bcast(bcast,5,MPI_REAL,mrank,sis%MPI_COMM_COMP,ierr)
    lat=bcast(1)
    lon=bcast(2)
    z=bcast(3)
    idex=bcast(4)
    jdex=bcast(5)
  end subroutine minloc7_real

  subroutine maxloc7_real(sis,val,lat,lon,z,idex,jdex)
    use mpi
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    REAL,intent(inout) :: val, lat, lon, z
    integer, intent(inout) :: idex,jdex
    INTEGER ierr, mrank
    REAL inreduce(2), outreduce(2), bcast(5)

    inreduce=(/ val, real(sis%MYPE) /)
    call MPI_Allreduce(inreduce,outreduce,1,MPI_2REAL,MPI_MAXLOC,&
         sis%MPI_COMM_COMP,ierr)
    val=outreduce(1)
    mrank=outreduce(2)
    bcast=(/ lat,lon,z,real(idex),real(jdex) /)
    call MPI_Bcast(bcast,5,MPI_REAL,mrank,sis%MPI_COMM_COMP,ierr)
    lat=bcast(1)
    lon=bcast(2)
    z=bcast(3)
    idex=bcast(4)
    jdex=bcast(5)
  end subroutine maxloc7_real

  subroutine minloc4_real(sis,val,idex,jdex)
    use mpi
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    REAL,intent(inout) :: val
    integer, intent(inout) :: idex,jdex
    INTEGER ierr, mrank
    REAL inreduce(2), outreduce(2), bcast(2)

    inreduce=(/ val, real(sis%MYPE) /)
    call MPI_Allreduce(inreduce,outreduce,1,MPI_2REAL,MPI_MINLOC,&
         sis%MPI_COMM_COMP,ierr)
    val=outreduce(1)
    mrank=outreduce(2)
    bcast=(/ real(idex),real(jdex) /)
    call MPI_Bcast(bcast,2,MPI_REAL,mrank,sis%MPI_COMM_COMP,ierr)
    idex=bcast(1)
    jdex=bcast(2)
  end subroutine minloc4_real

  subroutine maxloc4_real(sis,val,idex,jdex)
    use mpi
    use module_solver_internal_state, only: solver_internal_state
    type(solver_internal_state), intent(in) :: sis
    REAL,intent(inout) :: val
    integer, intent(inout) :: idex,jdex
    INTEGER ierr, mrank
    REAL inreduce(2), outreduce(2), bcast(2)

    inreduce=(/ val, real(sis%MYPE) /)
    call MPI_Allreduce(inreduce,outreduce,1,MPI_2REAL,MPI_MAXLOC,&
         sis%MPI_COMM_COMP,ierr)
    val=outreduce(1)
    mrank=outreduce(2)
    bcast=(/ real(idex),real(jdex) /)
    call MPI_Bcast(bcast,2,MPI_REAL,mrank,sis%MPI_COMM_COMP,ierr)
    idex=bcast(1)
    jdex=bcast(2)
  end subroutine maxloc4_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Private implementation functions
  subroutine reduce_real_impl(dat,rfnd,dfnd, &
                         ims,ime,jms,jme,kms,kme,&
                         its,ite,jts,jte,kts,kte,&
                         method,comm,rank)
    use mpi
    integer, intent(in) :: method, comm, rank
    real, intent(in) :: dat(ims:ime,jms:jme,kms:kme)
    real, intent(inout) :: rfnd
    integer, intent(in) :: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte
    double precision, intent(out) :: dfnd

    real :: rfnd_mpi
    double precision :: dfnd_mpi
    integer :: i,j,k,ierr
    real :: rthread
    real :: inreduce(2),outreduce(2)
    double precision :: dthread
    logical :: found
    
    !$OMP PARALLEL PRIVATE(k,found)
    rthread=0
    dthread=0
    found=.false.
    if(method==REDUCE_ADD) then
       do k=kts,kte
          !$OMP DO PRIVATE(i,j) REDUCTION(+:rthread)
          do j=jts,jte
             do i=its,ite
                rthread=rthread+dat(i,j,k)
             enddo
          enddo
          !$OMP END DO
       enddo
       call MPI_Allreduce(rthread,rfnd_mpi,1,MPI_REAL,MPI_SUM,comm,ierr)
       dfnd=rfnd_mpi
       rfnd=rfnd_mpi
    elseif(method==REDUCE_ADD_DOUBLE) then
       do k=kts,kte
          !$OMP DO PRIVATE(i,j) REDUCTION(+:dthread)
          do j=jts,jte
             do i=its,ite
                dthread=dthread+dat(i,j,k)
             enddo
          enddo
          !$OMP END DO
       enddo
       call MPI_Allreduce(dthread,dfnd_mpi,1,MPI_DOUBLE_PRECISION,&
            MPI_SUM,comm,ierr)
       rfnd=real(dfnd_mpi)
       dfnd=dfnd_mpi
    elseif(method==REDUCE_MAX) then
       rthread=dat(its,jts,kts)
       do k=kts,kte
          !$OMP DO PRIVATE(i,j) REDUCTION(max:rthread)
          do j=jts,jte
             do i=its,ite
                rthread=max(rthread,dat(i,j,k))
             enddo
          enddo
          !$OMP END DO
       enddo
       call MPI_Allreduce(rthread,rfnd_mpi,1,MPI_REAL,MPI_MAX,comm,ierr)
       dfnd=rfnd_mpi
       rfnd=rfnd_mpi
    elseif(method==REDUCE_MIN) then
       rthread=dat(its,jts,kts)
       do k=kts,kte
          !$OMP DO PRIVATE(i,j) REDUCTION(min:rthread)
          do j=jts,jte
             do i=its,ite
                rthread=min(rthread,dat(i,j,k))
             enddo
          enddo
          !$OMP END DO
       enddo
       call MPI_Allreduce(rthread,rfnd_mpi,1,MPI_REAL,MPI_MIN,comm,ierr)
       dfnd=rfnd_mpi
       rfnd=rfnd_mpi
    endif
    !$OMP END PARALLEL
  end subroutine reduce_real_impl

  subroutine reduce_find_real_impl(dat,rfnd,iloc,jloc,kloc,rankloc,&
                         ims,ime,jms,jme,kms,kme,&
                         its,ite,jts,jte,kts,kte,&
                         method,comm,rank)
    use mpi
    integer, intent(in) :: method, comm, rank
    integer, intent(inout) :: iloc,jloc,kloc,rankloc
    real, intent(in) :: dat(ims:ime,jms:jme,kms:kme)
    real, intent(inout) :: rfnd
    integer, intent(in) :: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte

    integer :: i,j,k,iloc_local,jloc_local,kloc_local
    integer :: iloc_thread, jloc_thread,kloc_thread
    real :: rfnd_local, rfnd_thread
    real :: inreduce(2),outreduce(2)
    integer :: bcast(3), ierr

    if(method/=FIND_MAX .and. method/=FIND_MIN) then
       write(0,*) 'In module_REDUCTION reduce_find, method must be FIND_MAX or FIND_MIN.'
       call MPI_Abort(MPI_COMM_WORLD,2,ierr)
    endif

    iloc=-99 ; iloc_local=-99 ; iloc_thread=-99
    jloc=-99 ; jloc_local=-99 ; jloc_thread=-99
    kloc=-99 ; kloc_local=-99 ; kloc_thread=-99
    rfnd=0.  ; rfnd_local=0.  ; rfnd_thread=0.
    k=kts

    !$OMP PARALLEL FIRSTPRIVATE(k,iloc_thread,jloc_thread,kloc_thread,rfnd_thread)
    ifmax: if(method==FIND_MAX) then
       kmax: do k=kts,kte
          !$OMP DO PRIVATE(j,i)
          jmax: do j=jts,jte
             imax: do i=its,ite
                if(dat(i,j,k)>rfnd_thread .or. iloc_thread<0) then
                   rfnd_thread=dat(i,j,k)
                   iloc_thread=i
                   jloc_thread=j
                   kloc_thread=k
                endif
             enddo imax
          enddo jmax
          !$OMP END DO
       enddo kmax
       !$OMP CRITICAL
       if(iloc_local<0 .or. rfnd_thread>rfnd_local) then
          rfnd_local=rfnd_thread
          iloc_local=iloc_thread
          jloc_local=jloc_thread
          kloc_local=kloc_thread
       endif
       !$OMP END CRITICAL
    else ! FIND_MIN
       kmin: do k=kts,kte
          !$OMP DO PRIVATE(j,i)
          jmin: do j=jts,jte
             imin: do i=its,ite
                if(dat(i,j,k)<rfnd_thread .or. iloc_thread<0) then
                   rfnd_thread=dat(i,j,k)
                   iloc_thread=i
                   jloc_thread=j
                   kloc_thread=k
                endif
             enddo imin
          enddo jmin
          !$OMP END DO
       enddo kmin
       !$OMP CRITICAL
       if(iloc_local<0 .or. rfnd_thread<rfnd_local) then
          rfnd_local=rfnd_thread
          iloc_local=iloc_thread
          jloc_local=jloc_thread
          kloc_local=kloc_thread
       endif
       !$OMP END CRITICAL
    endif ifmax
    !$OMP END PARALLEL

    ! Find the maximum of maximums, and its location:
    inreduce=(/rfnd_local,real(rank)/)
    if(method==FIND_MAX) then
       call MPI_Allreduce(inreduce,outreduce,1,MPI_2REAL,MPI_MAXLOC,comm,ierr)
    else ! FIND_MIN
       call MPI_Allreduce(inreduce,outreduce,1,MPI_2REAL,MPI_MINLOC,comm,ierr)
    endif
    rfnd=outreduce(1)
    rankloc=outreduce(2)
    if(rank==rankloc) then
       bcast=(/ iloc_local, jloc_local, kloc_local /)
    endif
    call MPI_Bcast(bcast,3,MPI_INTEGER,rankloc,comm,ierr)
    iloc=bcast(1)
    jloc=bcast(2)
    kloc=bcast(3)
  end subroutine reduce_find_real_impl


  subroutine reduce_find_integer_impl(dat,ifnd,iloc,jloc,kloc,rankloc,&
                         ims,ime,jms,jme,kms,kme,&
                         its,ite,jts,jte,kts,kte,&
                         method,comm,rank)
    use mpi
    integer, intent(in) :: method, comm, rank
    integer, intent(inout) :: iloc,jloc,kloc,rankloc
    integer, intent(in) :: dat(ims:ime,jms:jme,kms:kme)
    integer, intent(inout) :: ifnd
    integer, intent(in) :: ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte

    integer :: i,j,k,iloc_local,jloc_local,kloc_local
    integer :: iloc_thread, jloc_thread,kloc_thread
    integer :: ifnd_local, ifnd_thread
    integer :: inreduce(2),outreduce(2)
    integer :: bcast(3), ierr
    if(method/=FIND_MAX .and. method/=FIND_MIN) then
       write(0,*) 'In module_REDUCTION reduce_find, method must be FIND_MAX or FIND_MIN.'
       call MPI_Abort(MPI_COMM_WORLD,2,ierr)
    endif

    iloc=-99 ; iloc_local=-99 ; iloc_thread=-99
    jloc=-99 ; jloc_local=-99 ; jloc_thread=-99
    kloc=-99 ; kloc_local=-99 ; kloc_thread=-99
    ifnd=0.  ; ifnd_local=0.  ; ifnd_thread=0.
    k=kts

    !$OMP PARALLEL FIRSTPRIVATE(k,iloc_thread,jloc_thread,kloc_thread,ifnd_thread)
    ifmax: if(method==FIND_MAX) then
       kmax: do k=kts,kte
          !$OMP DO PRIVATE(j,i)
          jmax: do j=jts,jte
             imax: do i=its,ite
                if(dat(i,j,k)>ifnd_thread .or. iloc_thread<0) then
                   ifnd_thread=dat(i,j,k)
                   iloc_thread=i
                   jloc_thread=j
                   kloc_thread=k
                endif
             enddo imax
          enddo jmax
          !$OMP END DO
       enddo kmax
       !$OMP CRITICAL
       if(iloc_local<0 .or. ifnd_thread>ifnd_local) then
          ifnd_local=ifnd_thread
          iloc_local=iloc_thread
          jloc_local=jloc_thread
          kloc_local=kloc_thread
       endif
       !$OMP END CRITICAL
    else ! FIND_MIN
       kmin: do k=kts,kte
          !$OMP DO PRIVATE(j,i)
          jmin: do j=jts,jte
             imin: do i=its,ite
                if(dat(i,j,k)<ifnd_thread .or. iloc_thread<0) then
                   ifnd_thread=dat(i,j,k)
                   iloc_thread=i
                   jloc_thread=j
                   kloc_thread=k
                endif
             enddo imin
          enddo jmin
          !$OMP END DO
       enddo kmin
       !$OMP CRITICAL
       if(iloc_local<0 .or. ifnd_thread<ifnd_local) then
          ifnd_local=ifnd_thread
          iloc_local=iloc_thread
          jloc_local=jloc_thread
          kloc_local=kloc_thread
       endif
       !$OMP END CRITICAL
    endif ifmax
    !$OMP END PARALLEL

    ! Find the maximum of maximums, and its location:
    inreduce=(/ifnd_local,rank/)
    if(method==FIND_MAX) then
       call MPI_Allreduce(inreduce,outreduce,1,MPI_2INTEGER,MPI_MAXLOC,comm,ierr)
    else
       call MPI_Allreduce(inreduce,outreduce,1,MPI_2INTEGER,MPI_MINLOC,comm,ierr)
    endif
    ifnd=outreduce(1)
    rankloc=outreduce(2)
    if(rank==rankloc) then
       bcast=(/ iloc_local, jloc_local, kloc_local /)
    endif
    call MPI_Bcast(bcast,3,MPI_INTEGER,rankloc,comm,ierr)
    iloc=bcast(1)
    jloc=bcast(2)
    kloc=bcast(3)
  end subroutine reduce_find_integer_impl
end module module_REDUCTION
