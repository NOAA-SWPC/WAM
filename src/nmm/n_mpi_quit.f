      SUBROUTINE n_mpi_quit(iret)
      use n_mpi_def
      implicit none
!
      integer iret
 
      write(0,*) 'CALL stop mpi_quit ',iret
      CALL MPI_ABORT(MC_COMP,iret,info)
 
      RETURN
      END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
