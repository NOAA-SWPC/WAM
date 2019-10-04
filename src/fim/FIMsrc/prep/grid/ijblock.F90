!======================================================
! This subroutine fills in sub-regions and cache blocks of a rhombus in IJ row major order.
! The regions follow each other in column major order because that is what SMS does.
! If no cache blocking, then RegionFactors=1 and BlockDim=RegionDim.
! Each processor has one region. The results are returned in perm.
! The default layout is ij on the rhombus and perm is a mapping from ij order on the rhombus to the desired order
! MAX_DECOMP_DIMS is hardcoded as 2
! Jacques Middlecoff   April, 2008
!======================================================
SUBROUTINE ijblock_curve(nip,nPEs,NB,PE,start_gc,Nregions,RhombusDim,Nblocks,RegionDim,BlockDim,Rstart,Rend,perm)
IMPLICIT NONE
INTEGER,intent(IN   ) :: nip                   !Number of Icosahedral points
INTEGER,intent(IN   ) :: nPEs                  !Number of processors
INTEGER,intent(IN   ) :: NB                    !Number of cache blocks per processor
INTEGER,intent(INOUT) :: PE                    !Processor number for saving start and end
INTEGER,intent(IN   ) :: start_gc              !Starting location in the globe
INTEGER,intent(IN   ) :: Nregions   (2)        !Number of regions in each dimension of the rhombus
INTEGER,intent(IN   ) :: RhombusDim (2)        !The dimension of the rhombus
INTEGER,intent(IN   ) :: Nblocks    (NB*nPEs,2)!Number of blocks in each dimension of each region
INTEGER,intent(IN   ) :: RegionDim  (NB*nPEs,2)!The dimensions of each region
INTEGER,intent(IN   ) :: BlockDim   (NB,nPEs,2)!The dimensions of each block
INTEGER,intent(  OUT) :: Rstart     (   nPEs  )!Global starting location for each PE
INTEGER,intent(  OUT) :: Rend       (   nPEs  )!Global ending location for each PE
INTEGER,intent(  OUT) :: perm       (nip)      !The calculated permutation

INTEGER :: k,rhf1,rhf2,rd1,rd2,ref1,ref2,bd1,bd2
INTEGER :: RowOffset,ColOffset,BrowOffset,BcolOffset,RegionStart

 k=start_gc
 ColOffset = 0
 DO rhf1 = 1,Nregions(1)
  RowOffset = 0
  DO rhf2 = 1,Nregions(2)
    PE = PE + 1
    Rstart(PE) = k
    RegionStart = start_gc+RowOffset+ColOffset
    BcolOffset = 0
    do ref1 = 1,Nblocks(rhf1,1)
      BrowOffset = 0
      do ref2 = 1, Nblocks(rhf2,2)
        do bd2 = 1, BlockDim(ref2,rhf2,2)
          do bd1 = 1, BlockDim(ref1,rhf1,1)
            perm(k)=RegionStart+bd1-1+(bd2-1)*RhombusDim(1)+BrowOffset+BcolOffset
            k=k+1
          enddo
        enddo
        BrowOffset = BrowOffset+RhombusDim(1)*BlockDim(ref2,rhf2,2)
      enddo
      BcolOffset = BcolOffset + BlockDim(ref1,rhf1,1)
    enddo
    RowOffset = RowOffset+RhombusDim(1)*RegionDim(rhf2,2)
    Rend(PE) = k-1
  END DO
  ColOffset = ColOffset + RegionDim(rhf1,1)
 END DO
END SUBROUTINE ijblock_curve
