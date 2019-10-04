!======================================================
!  This subroutine creates the permutation array for
!  the 2D filling curve. 
!  This subroutine also calculates the start and end location for each PE
!
!  Author: Ning Wang,   Oct., 2006
!  Added IJ per PE blocking and cache blocking - Jacques Middlecoff, May 2008
!  Moved calculation of start and end from GetDecomp to here J. Middlecoff, Sep 2009
!======================================================
SUBROUTINE mk_perm(nip, curve, nPEs, NB, Rstart, Rend) 
USE DataStru,only: perm,inv_perm
IMPLICIT NONE

INTEGER,PARAMETER   :: MAX_DECOMP_DIMS = 2 !If this is changed ppp_factors GetRegions and ijblock must be changed
INTEGER,intent(IN ) :: nip                 !Number of icosahedral points
INTEGER,intent(IN ) :: curve               !The type of space filling curve
INTEGER,intent(IN ) :: nPEs                !Number of processors
INTEGER,intent(IN ) :: NB                  !The number of cache blocks per processor.
INTEGER,intent(OUT) :: Rstart(nPEs)        !Global starting location for each PE
INTEGER,intent(OUT) :: Rend  (nPEs)        !Global ending location for each PE

INTEGER :: start_gc, i, rhombus, f1, f2
INTEGER :: RhombusDim    (MAX_DECOMP_DIMS)         ! Number of points on rhombus side
INTEGER :: RD            (MAX_DECOMP_DIMS)         ! Tempory for RegionDim
INTEGER :: RF            (MAX_DECOMP_DIMS)         ! Temporary for rhombus or region factors
INTEGER :: RhombusFactors(MAX_DECOMP_DIMS)         ! Factorization of PEsPerRhombus
INTEGER :: RegDim        (nPEs)                    ! Temporary for region dimensions
INTEGER :: BlkDim        (nPEs)                    ! Temporary for block dimensions
INTEGER :: PointsPerRhombus                        ! Number of points in each rhombus
INTEGER :: PEsPerRhombus                           ! Number of PEs allocated to each rhombus
INTEGER :: RegionFactors (NB*nPEs,MAX_DECOMP_DIMS) ! Factorization of PEsPerRhombus
INTEGER :: RegionDim     (NB*nPEs,MAX_DECOMP_DIMS) ! Region dimensions
INTEGER :: BlockDim      (NB,NPEs,MAX_DECOMP_DIMS) ! Cache block dimensions
INTEGER :: PE                                      ! Processor number for start and end
LOGICAL :: C3NotDivBy10                            ! For curve 3 is nPEs not divisible by 10?

 write (6,*) 'JR mk_perm allocating perm nip=', nip
 ALLOCATE(    perm(nip))
 ALLOCATE(inv_perm(nip))

 PointsPerRhombus = nip /10
 PEsPerRhombus    = nPEs/10
 C3NotDivBy10     = curve==3.and.10*PEsPerRhombus/=nPEs
 IF(curve==0.or.curve==1.or.PESperRhombus==0.or.C3NotDivBy10)THEN !IJ or Hilbert order
  IF (curve == 1) THEN                                            !Hilbert curve
    DO rhombus = 1, 10
      start_gc = (rhombus - 1) * PointsPerRhombus + 2
      CALL hilbert_curve(PointsPerRhombus, 1, start_gc)
    END DO
  ELSE                                                            !IJ order
    DO i = 1, nip
      perm(i) = i
    END DO
  ENDIF
  call GetRegions(nip-2,1,nPEs,nPEs,.true.,Rstart,Rend)
 ELSE IF (curve == 2) THEN                                        !IJ block order with cache blocking
  call IJblockLayout(nPEs,nip,NB,perm,Rstart,Rend)
 ELSE IF (curve == 3) THEN                                        !Square Layout - no cache blocking
  call SquareLayout(nPEs,nip,perm,Rstart,Rend)
 ELSE
  print*,'Error in perm.F90: Curve out of range, curve =',curve
  stop
 END IF
 !Add the 2 extra points at the poles to the first and last points
 perm(1)    = 1
 perm(nip)  = nip
 Rstart(1)  = Rstart(1)  - 1
 Rend(nPEs) = Rend(nPEs) + 1
 DO i = 1, nip
  inv_perm(perm(i)) = i
 END DO
END SUBROUTINE mk_perm


SUBROUTINE dealloc_perm_array
USE DataStru
IMPLICIT NONE

 write (6,*) 'JR dealloc_perm_array deallocating perm'
 DEALLOCATE(perm)
 DEALLOCATE(inv_perm)

END SUBROUTINE dealloc_perm_array

