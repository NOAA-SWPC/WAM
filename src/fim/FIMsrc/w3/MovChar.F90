subroutine MovChar(KBUF,PFLD,LEN)
implicit none
integer    ,intent( IN) :: len
character*1,intent( IN) :: PFLD(LEN)
character*1,intent(OUT) :: kbuf(LEN)

  KBUF = PFLD

end subroutine MovChar
