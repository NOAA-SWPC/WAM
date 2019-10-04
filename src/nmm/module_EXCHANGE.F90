!-----------------------------------------------------------------------
                        module module_exchange
!-----------------------------------------------------------------------
!
!***  MODULE_EXCHANGE contains the halo exchange routines.  There is a
!***  unique routine for every combination of 2-D and 3-D real array
!***  exchanges being done at one time.  Each subroutine name begins with
!***  "exch" which is then followed by a string of integers.  Each "2" 
!***  in the string indicates exchange being done for a 2-D array.  
!***  Similarly each "3" or "4" in the string indicates exchange being done
!***  for a 3-D or 4-D  array.
!***  Currently there are routines for these combinations:
!***
!***  2, 22, 222, 2222, 23, 223, 3, 33, 333, 3333, 4
!
!***  A generic interface exists so that all of the routines
!***  may be called with the name "halo_exch".  If new routines
!***  are added because new combinations are needed then also
!***  add the routine's name to the interface block.
!
!
!***  Buffer arrays are used during the exchange process.  Set the size
!***  below in the parameter ibufexch.  If an error occurs where the 
!***  MPI library indicates that the receive buffer is too small then
!***  increase the size of ibufexch.
!
!***  The 4-element IHANDLE array is used for the nonblocking requests
!***  for all the ISENDS/IRECVS and their MPI_WAITS.  Here is the key
!***  to their use:
!***
!***  IRECV/store from north --> IHANDLE(1)
!***  IRECV/store from south --> IHANDLE(2)
!***  ISEND to north --> IHANDLE(3)
!***  ISEND to south --> IHANDLE(4)
!***
!***  IRECV/store from west --> IHANDLE(1)
!***  IRECV/store from east --> IHANDLE(2)
!***  ISEND to east --> IHANDLE(3)
!***  ISEND to west --> IHANDLE(4)
!
!-----------------------------------------------------------------------
!
use mpi
!
use module_kinds
!
use module_my_domain_specs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      private
!
      public :: halo_exch
!
!-----------------------------------------------------------------------
!
      integer(kind=kint),parameter :: ibufexch=2500000
!
      real(kind=kfpt),dimension(ibufexch) :: buf0,buf1,buf2,buf3
!
!-----------------------------------------------------------------------
!
      interface halo_exch
        module procedure exch2
        module procedure exch22
        module procedure exch222
        module procedure exch2222
        module procedure exch23
        module procedure exch223
        module procedure exch3
        module procedure exch33
        module procedure exch333
        module procedure exch3333
        module procedure exch4
      end interface
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch2(arr1,ll1,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for a single 2-D array.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 2-D array (=1)
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
!
 arr1                   ! array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!***  The array of neighbors called my_neb is filled in subroutine
!***  decomp in module_dm_parallel.  Recall that the directional
!***  designations for my_neb are:
!
!***      north: 1
!***       east: 2
!***      south: 3
!***       west: 4
!***  northeast: 5
!***  southeast: 6
!***  southwest: 7
!***  northwest: 8
!
!***  If my_neb(n) holds the task ID of each neighbor.  If there is
!***  no neighbor due to the presence of a global boundary then the
!***  value of my_neb(n) in that direction is -1.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!
      if(my_neb(1)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j)
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j)
        enddo
        enddo
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1)=buf1(ic)
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1)=buf0(ic)
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j)
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j)
        enddo
        enddo
       call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf0(ic)
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf1(ic)
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch2
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch3(arr1,ll1,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for a single 3-D array.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
!
 arr1                   ! array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
       call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch3
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch23(arr1,ll1,arr2,ll2,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for 2-D and 3-D arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 2-D array (=1)
,ll2 &                  ! vertical dimension of 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1                   ! 2-D array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-D array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j)
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j)
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j)
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j)
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
       call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf0(ic)
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf1(ic)
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch23
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch222(arr1,ll1,arr2,ll2,arr3,ll3,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for three 2-D arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 2-D array (=1)
,ll2 &                  ! vertical dimension of 2nd 2-D array (=1)
,ll3 &                  ! vertical dimension of 3rd 2-D array (=1)
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1 &                 ! 2-D array whose haloes are exchanged
,arr2 &                 ! 2-D array whose haloes are exchanged
,arr3                   ! 2-D array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,jte-j)
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,jts+j)
        enddo
        enddo
!
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,j)
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,j)
        enddo
        enddo
!
       call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf0(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf0(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j)=buf0(ic)
        enddo
        enddo
!
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf1(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf1(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j)=buf1(ic)
        enddo
        enddo
!
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch222
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch2222(arr1,ll1,arr2,ll2,arr3,ll3,arr4,ll4           &
                         ,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange halos for four 2-D arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 2-D array (=1)
,ll2 &                  ! vertical dimension of 2nd 2-D array (=1)
,ll3 &                  ! vertical dimension of 3rd 2-D array (=1)
,ll4 &                  ! vertical dimension of 3rd 2-D array (=1)
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1 &                 ! 2-D array whose haloes are exchanged
,arr2 &                 ! 2-D array whose haloes are exchanged
,arr3 &                 ! 2-D array whose haloes are exchanged
,arr4                   ! 2-D array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,jte-j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr4(i,jte-j)
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,jts+j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr4(i,jts+j)
        enddo
        enddo
!
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr4(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr4(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr4(i,j)
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr4(i,j)
        enddo
        enddo
!
       call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf0(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf0(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j)=buf0(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr4(i,j)=buf0(ic)
        enddo
        enddo
!
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf1(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf1(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j)=buf1(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr4(i,j)=buf1(ic)
        enddo
        enddo
!
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch2222
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch22(arr1,ll1,arr2,ll2,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for two 2-D arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 2-D array (=1)
,ll2 &                  ! vertical dimension of 2nd 2-D array (=1)
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1 &                 ! 2-D array whose haloes are exchanged
,arr2                   ! 2-D array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j)
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j)
        enddo
        enddo
!
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j)
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j)
        enddo
        enddo
!
       call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf0(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf0(ic)
        enddo
        enddo
!
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf1(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf1(ic)
        enddo
        enddo
!
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch22
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch223(arr1,ll1,arr2,ll2,arr3,ll3,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for two 2-D arrays and one 3-D array.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 2-D array (=1)
,ll2 &                  ! vertical dimension of 2nd 2-D array (=1)
,ll3 &                  ! vertical dimension of 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1 &                 ! 2-D array whose haloes are exchanged
,arr2                   ! 2-D array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-D array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,jte-j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,jts+j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
       call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf0(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf0(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf1(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf1(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch223
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch22333(arr1,ll1,arr2,ll2,arr3,ll3,arr4,ll4,arr5,ll5 &
                          ,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for two 2-D arrays and three 3-D arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 2-D array (=1)
,ll2 &                  ! vertical dimension of 2nd 2-D array (=1)
,ll3 &                  ! vertical dimension of 1st 3-D array
,ll4 &                  ! vertical dimension of 2nd 3-D array
,ll5 &                  ! vertical dimension of 3rd 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1 &                 ! 2-D array whose haloes are exchanged
,arr2                   ! 2-D array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-D array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll4),intent(inout) :: &
 arr4                   ! 3-D array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll5),intent(inout) :: &
 arr5                   ! 3-D array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr4(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll5
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr5(i,jte-j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr4(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll5
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr5(i,jts+j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr4(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll5
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr5(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr4(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll5
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr5(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr4(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll5
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr5(i,j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr4(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll5
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr5(i,j,k)
        enddo
        enddo
        enddo
!
       call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf0(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf0(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr4(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll5
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr5(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf1(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf1(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr4(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll5
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr5(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch22333
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch33(arr1,ll1,arr2,ll2,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes two 3-D real arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 3-D array
,ll2 &                  ! vertical dimension of 2nd 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
 arr1                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-D array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
       call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch33
!
!--------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch333(arr1,ll1,arr2,ll2,arr3,ll3,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes three 3-D real arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 3-D array
,ll2 &                  ! vertical dimension of 2nd 3-D array
,ll3 &                  ! vertical dimension of 3rd 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
 arr1                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-D array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,jte-j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,jts+j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch333
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch3333(arr1,ll1,arr2,ll2,arr3,ll3,arr4,ll4,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes four 3-D real arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 3-D array
,ll2 &                  ! vertical dimension of 2nd 3-D array
,ll3 &                  ! vertical dimension of 3rd 3-D array
,ll4 &                  ! vertical dimension of 4th 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
 arr1                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll4),intent(inout) :: &
 arr4                   ! 3-D array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr4(i,jte-j,k)
        enddo
        enddo
        enddo
!
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr4(i,jts+j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr4(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr4(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr4(i,j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr4(i,j,k)
        enddo
        enddo
        enddo
        call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr4(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr4(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch3333
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch4(arr1,ll1,nl1,nstart,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for a single 4-D array.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of array
,nl1 &                  ! 4th dimension of array
,nstart &               ! index of the 4th dimension to start exchange
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1,nl1),intent(inout) :: &
!
 arr1                   ! array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  Local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend &
                           ,k,n
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  North/South
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
        do n=nstart,nl1
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j,k,n)
        enddo
        enddo
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(1),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
        do n=nstart,nl1
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j,k,n)
        enddo
        enddo
        enddo
        enddo
        call mpi_issend(buf3,ic,mpi_real,my_neb(3),mype &
                       ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
        do n=nstart,nl1
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1,k,n)=buf1(ic)
        enddo
        enddo
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do n=nstart,nl1
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1,k,n)=buf0(ic)
        enddo
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  East/West
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  Send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
        do n=nstart,nl1
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
        call mpi_issend(buf2,ic,mpi_real,my_neb(2),mype &
                       ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
        do n=nstart,nl1
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
       call mpi_issend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  Store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do n=nstart,nl1
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k,n)=buf0(ic)
        enddo
        enddo
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  Store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
        do n=nstart,nl1
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k,n)=buf1(ic)
        enddo
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch4
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      end module module_exchange
!
!-----------------------------------------------------------------------
