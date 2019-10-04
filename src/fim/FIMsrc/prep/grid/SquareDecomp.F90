subroutine SquareDecomp(nip,nPEs,PEstart,NumCols,PEsPerCol,start_gc, &
                        FirstTime,perm,Rstart,Rend)
!This routine decomposes a rhombus into NumRegions=sum(PEsPerCol) regions
!A region is a processor footprint.
!The regions are almost square and almost equally sized.
!The algorithm divides a rhombus into columns.
!Each column has a possibly different number of rows.
!Each column is then divided equally into it's number of rows.
!The resulting layout has nearly square regions of nearly equal size.
!For a square processor count this routine gives results identical to ij-block.
!This routine returns the permutation array and the region size for each region
!The default layout is ij on the rhombus

!Author: Jacques Middlecoff, November, 2009

implicit none
integer,intent(IN)    :: nip                 ! Number of icosahedral points
integer,intent(IN)    :: nPEs                ! Total number of processors
integer,intent(IN)    :: PEstart             ! Start index for Rstart and Rend
integer,intent(IN)    :: NumCols             ! Number of columns the rhombus is divided into
integer,intent(IN)    :: PEsPerCol(NumCols)  ! Number of regions per column
integer,intent(IN)    :: start_gc            ! Starting point for perm
logical,intent(INOUT) :: FirstTime           ! True means this is the first call
integer,intent(OUT)   :: perm(nip)           ! Mapping from ij order to desired order
INTEGER,intent(OUT)   :: Rstart(nPEs)        ! Global starting location for each region
INTEGER,intent(OUT)   :: Rend  (nPEs)        ! Global ending location for each region
integer               :: PointsPerRhombus    ! Number of points in the rhombus
integer               :: NumRegions          ! Number of regions to divide the rhombus into
integer               :: RhombusSide         ! Length of each side of the rhombus
integer               :: k                   ! index for the perm array
integer               :: mincol              ! index for the column with the minimum ratio
integer               :: PEcolStart,PEcolEnd ! Column start and end index for Rstart and Rend
integer               :: row,col,column      ! Row and column indexes
integer               :: ColStart,Cstart,Cend! Column starting and end indecies
integer               :: ColWidth(NumCols)
integer               :: PointsInThisCol
integer               :: ExtraColPoints
integer               :: ExtraPoint
real                  :: RegionSize          ! Number or points in each region
real                  :: RowWidth

 NumRegions       = sum(PEsPerCol)
 PointsPerRhombus = nip/10
 RhombusSide      = sqrt(Float(PointsPerRhombus))
 RegionSize       = float(PointsPerRhombus)/float(NumRegions)
 k                = start_gc
 ColStart         = start_gc-1

!Calculate column widths adding extra points to columns with minimum ratio
 do col = 1,NumCols
   RowWidth      = float(RhombusSide)/Float(PEsPerCol(col))
   ColWidth(col) = RegionSize/RowWidth
 enddo
 ExtraColPoints = RhombusSide - sum(ColWidth)
 if(ExtraColPoints < 0 .or. ExtraColPoints > NumCols) then
   Print*,'Error in SquareDecomp',ExtraColPoints,RhombusSide,sum(ColWidth),NumCols
   stop
 endif
 do ExtraPoint=1,ExtraColPoints
   mincol           = minloc(float(colWidth)/float(PEsPerCol),DIM=1)
   ColWidth(mincol) = ColWidth(mincol)+1
 enddo

!Fill perm column by column
 do column = 1,NumCols
   do row =1,RhombusSide
      do col=1,ColWidth(column)
        perm(k) = col + (row-1)*RhombusSide + ColStart
        k = k + 1
      enddo
   enddo
   ColStart = ColStart + ColWidth(column)
 enddo

!Calculate Rstart and Rend by equally dividing each column into nearly equal regions.
 PEcolStart = PEstart
 do col = 1,NumCols
   PointsInThisCol = RhombusSide*ColWidth(col)
   PEcolEnd        = PEcolStart + PEsPerCol(col)-1
   call GetRegions(PointsInThisCol,PEcolStart,PEcolEnd,nPEs,FirstTime,Rstart,Rend)
   FirstTime  = .false.
   PEcolStart = PEcolEnd+1
 enddo

 return
end subroutine SquareDecomp
