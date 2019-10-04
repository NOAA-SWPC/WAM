subroutine GetRegions(BlockSize,PEstart,PEend,nPEs,FirstTIme,Rstart,Rend)
!This routine does a simple equally sized 1-D decomposition
!This routine divides a BlockSize sized block into nprocs equal sized regions.
!This routine returns the start and end of each region (processor footprint).

!Author: Jacques Middlecoff, September, 2009
!Added FirstTime so it will work section by section - Jacques Middlecoff, November 2009

implicit none
INTEGER,PARAMETER     :: MAX_DECOMP_DIMS = 2 !If this is changed ppp_factors,perm and ijblock must be changed
INTEGER,intent(IN )   :: BlockSize           !Number of icosahedral points in the block to be decomposed
integer,intent(IN)    :: PEstart,PEend       !Start and end index for Rstart and Rend
INTEGER,intent(IN )   :: nPEs                !Total number of processors
LOGICAL,intent(IN )   :: FirstTime           !true means this is the first call
INTEGER,intent(INOUT) :: Rstart(nPEs)        !Global starting location for each PE
INTEGER,intent(INOUT) :: Rend  (nPEs)        !Global ending location for each PE
INTEGER               :: RegDim(nPEs)        !Temporary for region dimensions
INTEGER               :: RD(MAX_DECOMP_DIMS) !Tempory for RegionDim
INTEGER               :: RF(MAX_DECOMP_DIMS) !Temporary for rhombus or region factors
INTEGER               :: nprocs              !Number of processors by which to divide BlockSIze
INTEGER               :: inc                 !Increment added to PEstart
INTEGER               :: PE                  !Processor number

 nprocs = PEend - PEstart+1
 RD(1)  = BlockSize
 RD(2)  = 1
 RF(1)  = nprocs
 RF(2)  = 1
 call ppp_regions(1,RD,RF,nprocs,RegDim(PEstart:PEend))
 if(FirstTime) then
   Rstart(PEstart) = 2
   Rend  (PEstart) = RegDim(PEstart)+1
   inc             = 1
 else
   inc             = 0
 endif
 do PE = PEstart+inc,PEend
   Rstart(PE) = Rend(PE-1) + 1
   Rend  (PE) = Rstart(PE) + RegDim(PE) - 1
 enddo
 return

end subroutine GetRegions
