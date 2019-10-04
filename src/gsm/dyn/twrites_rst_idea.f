      subroutine twrites_rst_idea(fname,IOPROC,fhour,idate,
     &           si,ls_nodes,max_ls_nodes,trie_ls,trio_ls) 
!-------------------------------------------------------------------
!*** program log
!*** Dec, 2009 Jun Wang:  write spectral variables for restart
!*** Dec, 2010 Jun Wang:  change to nemsio library
!*** Feb, 2011 Henry Juang: use generic argument for spectral fit to mass_dp and ndsl
!    Jun 26 2014 S. Moorthi - added mpigathe8
!    May 11 2017 W. Yang    - For WAM restart run.
!-------------------------------------------------------------------
!
      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_coordinate_def	
      use namelist_dynamics_def
      use gfs_dyn_mpi_def
      use nemsio_module
!
      implicit none
!
      character(*),intent(in) :: fname
      integer,intent(in) :: ioproc
      real(kind=kind_evod),intent(in) :: fhour
      integer,intent(in) :: idate(4)
!
      integer,intent(in)             :: ls_nodes(ls_dim,nodes)
      integer,intent(in)             :: max_ls_nodes(nodes)
!
      real(kind=kind_evod),intent(in) :: trie_ls(len_trie_ls,2,lotls)
      real(kind=kind_evod),intent(in) :: trio_ls(len_trio_ls,2,lotls)
      real(kind=kind_evod) si(levp1)
!
!local variables:
!
      integer              ierr,j,k,lenrec
!
      integer i, jj, step
!
      real(kind=kind_mpi_r),allocatable :: trieo_ls_node (:,:,:)
      real(kind=kind_mpi_r),allocatable :: trieo_ls_nodes(:,:,:,:)
!
!
!---------------------------------------------------------------------
!
      call mpi_comm_size(MPI_COMM_ALL,i,ierr)
!
!-- allocate
      allocate ( trieo_ls_node  ( len_trie_ls_max+len_trio_ls_max,
     x                            2, lotls ) )
      trieo_ls_node = 0.0
!
      do k=1,lotls
        do j=1,len_trie_ls
          trieo_ls_node(j,1,k) = trie_ls(j,1,k)
          trieo_ls_node(j,2,k) = trie_ls(j,2,k)
        enddo
        do j=1,len_trio_ls
          jj = j+len_trie_ls_max
          trieo_ls_node(jj,1,k) = trio_ls(j,1,k)
          trieo_ls_node(jj,2,k) = trio_ls(j,2,k)
        enddo
      enddo
!
!-- collect data to ioproc
!-----------
      if ( me .eq. ioproc ) then
        write(0,*)'ALLOC PARMS TWRITE ',len_trie_ls_max+len_trio_ls_max,
     &      2, lotls, nodes,1
!
         allocate ( trieo_ls_nodes ( len_trie_ls_max+len_trio_ls_max,
     &                               2, lotls, nodes ),
     &              stat=ierr )
       else
         allocate (trieo_ls_nodes(1, 1, 1, 1), stat = ierr)
      endif
      if (ierr .ne. 0) then
        write (0,*) ' GWX trieo_ls_nodes allocate failed'
        call mpi_abort(mpi_comm_all,ierr,i)
       endif
!
      call mpi_barrier(MPI_COMM_ALL,ierr)
      if(nodes >1 )then
        lenrec = (len_trie_ls_max+len_trio_ls_max) * 2 * lotls
!
        call mpi_gather(trieo_ls_node , lenrec, MPI_R_MPI_R,
     &                  trieo_ls_nodes, lenrec, MPI_R_MPI_R,
     &                  ioproc, MPI_COMM_ALL, ierr)
        call mpi_barrier(MPI_COMM_ALL,ierr)
      else
       trieo_ls_nodes(:,:,:,1)=trieo_ls_node(:,:,:)
      endif
      deallocate ( trieo_ls_node  )
!
!-- write out data
      IF (me.eq.ioproc) THEN
       OPEN(1050, FILE=fname,FORM='unformatted')
       WRITE(1050) trieo_ls_nodes
       CLOSE(1050)
      endif   !me.eq.ioproc
      deallocate ( trieo_ls_nodes )
!!
      call mpi_barrier(MPI_COMM_ALL,ierr)

      return
      end
