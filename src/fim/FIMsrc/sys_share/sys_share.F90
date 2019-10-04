! Substitutes for non-standard system intrinsic subroutines.  
!JR Why is an empty module with a (also possibly empty) subroutine appended 
!JR to the end needed? Mac complained about empty file, so added a stub.

module module_sys_share
implicit none
end module module_sys_share

! Substitute for non-standard flush() intrinsic subroutine.  
!
! If your system does not support flush(), #define NO_FLUSH 
! and build this proxy.  

! if flush() is not supported...  
#ifdef NO_FLUSH
!SMS$IGNORE BEGIN
       subroutine flush(lun)
         implicit none
         integer,intent(in)::lun
!TODO:  this works on IBM, generalize if needed for other machines
         call flush_(lun)
       end subroutine flush
!SMS$IGNORE END
#endif

subroutine stub_to_satisfy_linkers_that_cant_handle_empty_doto_files
end subroutine stub_to_satisfy_linkers_that_cant_handle_empty_doto_files
