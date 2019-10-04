subroutine IJblockLayout(nPEs,nip,NB,perm,Rstart,Rend)
!This routine creates a retangular block layout.
!The blocks are determined by the factors of the number of processors.
!The algorithm is the same on each of the 10 rhombi.
!The same blocking algorithm is reapplied to create NB cache sub-blocks
!If nPEs not divisable by 10 then decompose with PEsPerRhombus which is generally better than ij
!Created this routine using logic from perm.F90 - Jacques Middlecoff Nov. 2009

implicit none
integer,intent(IN)  :: nPEs                ! Total number of processors
integer,intent(IN)  :: nip                 ! Number of icosahedral points
INTEGER,intent(IN ) :: NB                  ! The number of cache blocks per processor.
integer,intent(OUT) :: perm(nip)           ! Mapping from ij order to desired order
INTEGER,intent(OUT) :: Rstart(nPEs)        ! Global starting location for each PE
INTEGER,intent(OUT) :: Rend  (nPEs)        ! Global ending location for each PE
INTEGER,PARAMETER   :: MAX_DECOMP_DIMS = 2 ! If this is changed ppp_factors GetRegions and ijblock must be changed
INTEGER             :: start_gc,i,rhombus,f1,f2
INTEGER             :: RhombusDim    (MAX_DECOMP_DIMS)         ! Number of points on rhombus side
INTEGER             :: RD            (MAX_DECOMP_DIMS)         ! Tempory for RegionDim
INTEGER             :: RF            (MAX_DECOMP_DIMS)         ! Temporary for rhombus or region factors
INTEGER             :: RhombusFactors(MAX_DECOMP_DIMS)         ! Factorization of PEsPerRhombus
INTEGER             :: RegDim        (nPEs)                    ! Temporary for region dimensions
INTEGER             :: BlkDim        (nPEs)                    ! Temporary for block dimensions
INTEGER             :: PointsPerRhombus
INTEGER             :: PEsPerRhombus
INTEGER             :: RegionFactors (NB*nPEs,MAX_DECOMP_DIMS) ! Factorization of PEsPerRhombus
INTEGER             :: RegionDim     (NB*nPEs,MAX_DECOMP_DIMS) ! Region dimensions
INTEGER             :: BlockDim      (NB,NPEs,MAX_DECOMP_DIMS) ! Cache block dimensions
INTEGER             :: PE                                      ! Processor number for start and end

 PEsPerRhombus    = nPEs/10
 PointsPerRhombus = nip /10
 RhombusDim(1)    = nint(sqrt(float(PointsPerRhombus)))
 RhombusDim(2)    = nint(sqrt(float(PointsPerRhombus)))
 call ppp_factors(PEsPerRhombus,RhombusDim,RF)
 call ppp_regions(2            ,RhombusDim,RF,PEsPerRhombus,RegDim)
 RegionDim(1:RF(1),1) = RegDim(      1:RF(1)      )
 RegionDim(1:RF(2),2) = RegDim(RF(1)+1:RF(1)+RF(2))
 RhombusFactors = RF
 do f1 = 1,RhombusFactors(1)
   RD(1) = RegionDim(f1,1)
   RD(2) = RegionDim( 1,2)
   call ppp_factors(NB,RD,RF)
   call ppp_regions(2,RD,RF,NB,BlockDim(:,f1,1))
   RegionFactors(f1,1) = RF(1)
 enddo
 DO f2 = 1,RhombusFactors(2)
   RD(1) = RegionDim( 1,1)
   RD(2) = RegionDim(f2,2)
   call ppp_factors(NB,RD,RF)
   call ppp_regions(2,RD,RF,NB,BlkDim)
   RegionFactors(f2,2) = RF(2)
   BlockDim(1:RF(2),f2,2) = BlkDim(RF(1)+1:RF(1)+RF(2))
 enddo
 PE = 0
 DO rhombus = 1, 10
   start_gc = (rhombus-1)*PointsPerRhombus+2
   call ijblock_curve(nip,nPEs,NB,PE,start_gc,   &
                      RhombusFactors,RhombusDim, &
                      RegionFactors, RegionDim,  &
                      BlockDim,Rstart,Rend,perm  )
 END DO 
 if(PE /= nPEs) then ! Must set Rstart and Rend corrrectly
   call GetRegions(nip-2,1,nPEs,nPEs,.true.,Rstart,Rend)
 endif
 return

end subroutine IJblockLayout
