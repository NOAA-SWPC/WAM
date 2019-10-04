      SUBROUTINE mpi_quit(iret)
      use mpi_def
      implicit none
!
      integer iret
 
      write(0,*) 'CALL stop mpi_quit ',iret
      CALL MPI_ABORT(MC_COMP,iret,info)
 
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
