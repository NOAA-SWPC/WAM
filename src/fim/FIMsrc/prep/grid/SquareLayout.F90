subroutine SquareLayout(nPEs,nip,perm,Rstart,Rend)
!This routine calculates an almost square almost equally sized layout.
!The number of processors must be divisible by 10.
!For an input processor count and number of points, this routine returns the 
!permutation array and the processor footprint start and end values.
!The algorithm is the same on each of the 10 rhombi.
!The algorithm divides a rhombus into columns.
!Each column has a possibly different number of rows.
!Each column is then divided equally into it's number of rows.
!The resulting layout has nearly square footprints of nearly equal size.
!For a square processor count this routine gives results identical to ij-block.
!Compared to ij-block for 150p G7, this algorithm reduced the difference in 
!  largest to smallest footprint size by a factor of three: 
!  from a 6% difference for ij-block to a 2% difference for this algorithm.
!
!Author: Jacques Middlecoff, November, 2009

implicit none
integer,intent(IN)  :: nPEs         ! Total number of processors
integer,intent(IN)  :: nip          ! Number of icosahedral points
integer,intent(OUT) :: perm(nip)    ! Mapping from ij order to desired order
INTEGER,intent(OUT) :: Rstart(nPEs) ! Global starting location for each PE
INTEGER,intent(OUT) :: Rend  (nPEs) ! Global ending location for each PE
integer,allocatable :: PEsPerCol(:)
integer             :: PEsPerRhombus
integer             :: PointsPerRhombus
integer             :: MinPEsPerCol,MaxPEsPerCol,PEsLeft
integer             :: col,ColsLeft,NumCols
integer             :: rhombus,start_gc
integer             :: PEstart,PEend
logical             :: FirstTime

 PEsPerRhombus    = nPEs/10
 if(10*PEsPerRhombus /= nPEs) then
   print*,'The number of processors must be divisible by 10, nPEs=',nPEs
   stop
 endif
 PointsPerRhombus = nip /10
 PEstart          = 1
 FirstTime        = .true.
 do rhombus = 1, 10
   start_gc        = (rhombus-1)*PointsPerRhombus+2
   MinPEsPerCol = sqrt(float(PEsPerRhombus))
   MaxPEsPerCol = sqrt(float(PEsPerRhombus)) + 1
   if(PEsPerRhombus-MinPEsPerCol**2<MaxPEsPerCol**2-PEsPerRhombus) then
      NumCols = MinPEsPerCol
   else
      NumCols = MaxPEsPerCol
   endif
   allocate(PEsPerCol(NumCols))
   PEsLeft=PEsPerRhombus
   do col = 1,NumCols
      ColsLeft = NumCols-col
      !See if the next column can have MaxPEsPerCol
      if(PEsLeft-MaxPEsPerCol < ColsLeft*MinPEsPerCol) then
        PEsPerCol(col) = MinPEsPerCol
      else
        PEsPerCol(col) = MaxPEsPerCol
      endif
      PEsLeft = PEsLeft - PEsPerCol(col)
   enddo
   if(PEsLeft /= 0) then
      print*,'Error in ProcessorLayout: PEsLeft = ',PEsLeft
      stop
   endif
   PEend = PEstart + PEsPerRhombus-1 !This assumes nPEs is divisible by 10
   call SquareDecomp(nip,nPEs,PEstart,NumCols,PEsPerCol,start_gc,  &
                     FirstTime,perm,Rstart,Rend)
   PEstart = PEend+1
   deallocate(PEsPerCol)
 enddo

end subroutine SquareLayout
