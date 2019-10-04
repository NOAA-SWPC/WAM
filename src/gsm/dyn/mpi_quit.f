      SUBROUTINE mpi_quit(iret)
      use gfs_dyn_mpi_def
      implicit none
!
      integer iret
 
      write(0,*) 'CALL stop mpi_quit ',iret
!      CALL MPI_ABORT(MC_COMP,iret,info)
      CALL MPI_ABORT(MPI_COMM_WORLD,iret,info)
 
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
