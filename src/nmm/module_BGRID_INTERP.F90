!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      MODULE MODULE_BGRID_INTERP
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: V_TO_H_BGRID                                            &
               ,H_TO_V_BGRID
!
!-----------------------------------------------------------------------
!
      REAL :: BC_DUMMY=0.                                                  !<-- User-specified value for boundary mass points
                                                                           !    in regional mode
!
      interface V_TO_H_BGRID
        module procedure V_TO_H_BGRID2D
        module procedure V_TO_H_BGRID3D
      end interface V_TO_H_BGRID
!
      interface H_TO_V_BGRID
        module procedure H_TO_V_BGRID2D
        module procedure H_TO_V_BGRID3D
      end interface H_TO_V_BGRID
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE V_TO_H_BGRID2D(V_VALUE,IM,JM,GLOBAL,H_VALUE,V_SAVE)
!
!-----------------------------------------------------------------------
!***  PERFORM 4-POINT AVERAGING OF QUANTITIES ON B-GRID VELOCITY POINTS
!***  TO MASS POINTS.
!
!***  FOR THE REGIONAL MODE THE USER MUST SPECIFY THE DUMMY VALUE
!***  ON THE OUTPUT BOUNDARY SINCE INTERPOLATION IS NOT VALID THERE.
!***  THAT VALUE IS "BC_DUMMY" IN THE SPECIFICATION SECTION OF
!***  THIS MODULE.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!----------
!*** Input
!----------
!
      INTEGER                       ,INTENT(IN)  :: IM                  &  !<-- Full west-east array dimension
                                                   ,JM                     !<-- Full south-north array dimension
!
      REAL,DIMENSION(1:IM,1:JM),INTENT(IN)       :: V_VALUE                !<-- Input values on velocity points
!
      LOGICAL                       ,INTENT(IN)  :: GLOBAL                 !<-- Logical flag: True=>Global; False=>Regional
!
!-----------
!*** Output
!-----------
!
      REAL,DIMENSION(1:IM,1:JM),INTENT(OUT)      :: H_VALUE                !<-- Output values on mass (H) points
      REAL,DIMENSION(:,:),POINTER,INTENT(INOUT),OPTIONAL :: V_SAVE         !<-- Original values on boundary V points
!local vars 
      INTEGER  lm
      real,dimension(:,:,:),allocatable          :: tmp
!
      LM=1
      allocate(tmp(1:IM,1:JM,1))
      call  V_TO_H_BGRID3D(reshape(v_value,(/IM,JM,1/)),IM,JM,LM,GLOBAL, &
            tmp,v_save)
      H_VALUE(1:IM,1:JM)=TMP(1:IM,1:JM,1)
      deallocate(tmp)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE V_TO_H_BGRID2D
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE V_TO_H_BGRID3D(V_VALUE,IM,JM,LM,GLOBAL,H_VALUE,V_SAVE)
!
!-----------------------------------------------------------------------
!***  PERFORM 4-POINT AVERAGING OF QUANTITIES ON B-GRID VELOCITY POINTS
!***  TO MASS POINTS.
!
!***  FOR THE REGIONAL MODE THE USER MUST SPECIFY THE DUMMY VALUE
!***  ON THE OUTPUT BOUNDARY SINCE INTERPOLATION IS NOT VALID THERE.
!***  THAT VALUE IS "BC_DUMMY" IN THE SPECIFICATION SECTION OF
!***  THIS MODULE.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!----------
!*** Input
!----------
!
      INTEGER                       ,INTENT(IN)  :: IM                  &  !<-- Full west-east array dimension
                                                   ,JM                  &  !<-- Full south-north array dimension
                                                   ,LM                     !<-- Number of model layers
!
      REAL,DIMENSION(1:IM,1:JM,1:LM),INTENT(IN)  :: V_VALUE                !<-- Input values on velocity points
!
      LOGICAL                       ,INTENT(IN)  :: GLOBAL                 !<-- Logical flag: True=>Global; False=>Regional
!
!-----------
!*** Output
!-----------
!
      REAL,DIMENSION(1:IM,1:JM,1:LM),INTENT(OUT) :: H_VALUE                !<-- Output values on mass (H) points
      REAL,DIMENSION(:,:),POINTER,INTENT(INOUT),OPTIONAL :: V_SAVE         !<-- Original values on boundary V points
!
!----------
!*** Local
!----------
!
      INTEGER                                    :: I,J,L,N
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  LAYOUT IN I
!-----------------------------------------------------------------------
!
!                                   |
!                                   |<---- Integration western boundary
!                                   |
!               
!                          H(1)     H(2)     H(3)
!
!                              V(1)      V(2)     V(3)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!                     V(IM-3)  V(IM-2)   V(IM-1)  V(IM)                          
!                   
!                          H(IM-2)  H(IM-1)  H(IM)
!
!                                   |
! Integration eastern boundary ---->|
!                                   |
!
!-----------------------------------------------------------------------
!***  LAYOUT IN J
!
!***  MASS POINTS ACROSS THE POLE FROM EACH OTHER HAVE EQUAL VALUES.
!***  VELOCITY POINTS ACROSS THE POLE FROM EACH OTHER HAVE EQUAL
!***  VALUES OF OPPOSITE SIGN.  
!***  THE LOCATION ACROSS THE POLE FROM "I" IS "I_A=I+(IM-3)/2".
!***  IF I_A > IM THEN I_A=I_A-IM+3.
!***  IM IS ALWAYS ODD.
!
!***  THE SOUTH POLE COINCIDES WITH H POINTS AT J=2.
!***  H POINTS AT J=1 REFLECT "ACROSS THE POLE" VALUES FOR J=2.
!***  THE NORTH POLE COINCIDES WITH H POINTS AT J=JM-1.
!***  H POINTS AT J=JM REFLECT "ACROSS THE POLE" VALUES FOR J=JM-1.
!
!***  
!-----------------------------------------------------------------------
!
!           J=JM       V    V   V      <--- Phantom row; set equal to V at J=JM-1.
!
!        J=JM       H    H    H    H      <---  1 row "north" of pole reflects 1 row south
!                                                                  
!           J=JM-1     V    V    V     <--- Negative reflection of V across the pole at J=JM-2
!
!        J=JM-1     H    H    H    H      <---  NORTH POLE
!
!           J=JM-2     V    V    V
!
!        J=JM-2     H    H    H    H      <---  1 row south of pole
!
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
!
!          J=3      H    H    H    H      <---  1 row north of pole
!
!             J=2      V    V    V
!
!          J=2      H    H    H    H      <---  SOUTH POLE
!
!             J=1      V    V    V     <--- Negative reflection of V accross the pole at J=2
!
!          J=1      H    H    H    H      <---  1 row "south" of pole reflects 1 row north
!
!
!-----------------------------------------------------------------------
!***  PERFORM THE INTERPOLATION ON THE INTERNAL REGION OF THE
!***  OUTPUT ARRAY.  IT IS APPLICABLE TO BOTH REGIONAL AND
!***  GLOBAL MODES.
!-----------------------------------------------------------------------
!
      DO L=1,LM
      DO J=2,JM-1
      DO I=2,IM-1
!
        H_VALUE(I,J,L)=0.25*(V_VALUE(I-1,J-1,L)                         &
                            +V_VALUE(I-1,J  ,L)                         &
                            +V_VALUE(I  ,J-1,L)                         &
                            +V_VALUE(I  ,J  ,L) )
!
      ENDDO
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  FOR GLOBAL MODE THE OVERLAP VALUES ON VELOCITY POINTS WILL BE
!***  USED TO INTERPOLATE TO THE BOUNDARY MASS POINTS.
!-----------------------------------------------------------------------
!
      layers: DO L=1,LM
!
        mode: IF(GLOBAL)THEN
!
!--------------------------
!***  West/East Boundaries
!--------------------------
!
!***  General Case
!
          DO J=3,JM-2                                                     !<-- From 1 row north of S. Pole to 1 row south of N. Pole
!
            H_VALUE( 1,J,L)=0.25*(V_VALUE(IM-3,J-1,L)                   & !<-- West
                                 +V_VALUE(IM-2,J-1,L)                   &
                                 +V_VALUE(IM-3,J  ,L)                   &
                                 +V_VALUE(IM-2,J  ,L) )
!
            H_VALUE(IM,J,L)=0.25*(V_VALUE(2,J-1,L)                      & !<-- East
                                 +V_VALUE(3,J-1,L)                      &
                                 +V_VALUE(2,J  ,L)                      &
                                 +V_VALUE(3,J  ,L) )
!
          ENDDO
!
!***  One row south of South Pole reflects one row north of pole.
!
          H_VALUE( 1,1,L)=H_VALUE( 1,3,L)                                 !<-- West
          H_VALUE(IM,1,L)=H_VALUE(IM,3,L)                                 !<-- East
!
!***  H value at S. Pole averages V values to the north and
!***  the "south" which are the same.
!
          H_VALUE(1,2,L)=0.25*(V_VALUE(IM-3,2,L)                        & !<-- West
                              +V_VALUE(IM-2,2,L)                        &     
                              +V_VALUE(IM-3,2,L)                        &  
                              +V_VALUE(IM-2,2,L) )
!
          H_VALUE(IM,2,L)=0.25*(V_VALUE(2,2,L)                          & !<-- East
                               +V_VALUE(3,2,L)                          &     
                               +V_VALUE(2,2,L)                          &  
                               +V_VALUE(3,2,L) )
!
!***  H value at N. Pole averages V values to the south and
!***  the "north" which are the same.
!
          H_VALUE(1,JM-1,L)=0.25*(V_VALUE(IM-3,JM-2,L)                  & !<-- West
                                 +V_VALUE(IM-2,JM-2,L)                  & 
                                 +V_VALUE(IM-3,JM-2,L)                  & 
                                 +V_VALUE(IM-2,JM-2,L) )
!
          H_VALUE(IM,JM-1,L)=0.25*(V_VALUE(2,JM-2,L)                    & !<-- East
                                  +V_VALUE(3,JM-2,L)                    & 
                                  +V_VALUE(2,JM-2,L)                    & 
                                  +V_VALUE(3,JM-2,L) )
!
!***  One row north of North Pole reflects one row south of pole.
!
          H_VALUE( 1,JM,L)=H_VALUE( 1,JM-2,L)                             !<-- West
          H_VALUE(IM,JM,L)=H_VALUE(IM,JM-2,L)                             !<-- East
!
!----------------------------
!***  South/North Boundaries
!----------------------------
!
          DO I=1,IM
            H_VALUE(I, 1,L)=H_VALUE(I,   3,L)                             !<-- South
            H_VALUE(I,JM,L)=H_VALUE(I,JM-2,L)                             !<-- North
          ENDDO
!
!-----------------------------------------------------------------------
!***  FOR REGIONAL MODE JUST USE THE DUMMY VALUE FOR BOUNDARY POINTS.
!-----------------------------------------------------------------------
!
        ELSE mode
!
!----------------------------
!***  South/North Boundaries
!----------------------------
!
          DO I=1,IM
            H_VALUE(I, 1,L)=BC_DUMMY
            H_VALUE(I,JM,L)=BC_DUMMY
          ENDDO
!
!--------------------------
!***  West/East Boundaries
!--------------------------
!
          DO J=2,JM-1
            H_VALUE( 1,J,L)=BC_DUMMY
            H_VALUE(IM,J,L)=BC_DUMMY
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDIF mode
!
      ENDDO layers
!
!-----------------------------------------------------------------------
!***  IN THE REGIONAL MODE THE BOUNDARY VALUES ON VELOCITY POINTS
!***  WILL NOT BE RETRIEVABLE THROUGH REVERSE INTERPOLATION OF 
!***  THESE VALUES ON MASS POINTS SINCE THE BOUNDARY MASS POINTS
!***  DO NOT HAVE ACTUAL VALUES.  
!***  THEREFORE THE FOLLOWING SECTION SAVES THE BOUNDARY VELOCITY
!***  POINT VALUES THAT CAN THEN BE USED TO COMPLETE THE VELOCITY
!***  ARRAY AFTER A REVERSE INTERPOLATION IS DONE.
!
!***  ALL VELOCITY POINT VALUES AT I=IM and J=JM ARE MEANINGLESS
!***  BECAUSE THEY ARE OUTSIDE THE INTEGRATION DOMAIN.
!-----------------------------------------------------------------------
!
!jw      save: IF(.NOT.GLOBAL)THEN
      save: IF(.NOT.GLOBAL .and. present(v_save) )THEN
!
        DO L=1,LM
!
          N=0
!
!--------------------
!***  South Boundary
!--------------------
!
          DO I=1,IM-1
            N=N+1
            V_SAVE(N,L)=V_VALUE(I,1,L)
          ENDDO
!
!--------------------
!***  North Boundary
!--------------------
!
          DO I=1,IM-1
            N=N+1
            V_SAVE(N,L)=V_VALUE(I,JM-1,L)                                 !<-- Northernmost V points in integration domain
          ENDDO
!
!--------------------
!***  West Boundary
!--------------------
!
          DO J=2,JM-1
            N=N+1
            V_SAVE(N,L)=V_VALUE(1,J,L)
          ENDDO
!
!--------------------
!***  East Boundary
!--------------------
!
          DO J=2,JM-1
            N=N+1
            V_SAVE(N,L)=V_VALUE(IM-1,J,L)                                 !<-- Easternmost V points in integration domain
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDDO
!
      ENDIF save
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE V_TO_H_BGRID3D
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE H_TO_V_BGRID2D(H_VALUE,IM,JM,GLOBAL,V_SAVE,V_VALUE)
!
!-----------------------------------------------------------------------
!***  PERFORM 4-POINT AVERAGING OF QUANTITIES ON B-GRID MASS POINTS
!***  TO VELOCITY POINTS.  THIS ROUTINE IS INTENDED FOR USE IN
!***  INTERPOLATING VELOCITY VALUES BACK TO V POINTS AFTER HAVING
!***  BEEN INTERPOLATED ONTO H POINTS.  THE INPUT POINTER V_SAVE
!***  IS USED TO RECOVER THE VELOCITY VALUES ON THE BOUNDARY V POINTS
!***  WHEN IN REGIONAL MODE BECAUSE THEY CANNOT BE RECOVERED
!***  BY THE INTERPOLATION FROM H POINTS GIVEN THAT BOUNDARY H POINTS
!***  IN REGIONAL MODE COULD NOT BE OBTAINED BY INTERPOLATIONS FROM
!***  THE V POINTS.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!----------
!*** Input
!----------
!
      INTEGER                       ,INTENT(IN)  :: IM                  &  !<-- Full west-east array dimension
                                                   ,JM                     !<-- Full south-north array dimension
!
      REAL,DIMENSION(1:IM,1:JM),INTENT(IN)       :: H_VALUE                !<-- Input values on mass points
      REAL,DIMENSION(:,:),POINTER   ,INTENT(IN)  :: V_SAVE                 !<-- Saved values on boundary velocity points
!
      LOGICAL                       ,INTENT(IN)  :: GLOBAL                 !<-- Logical flag: True=>Global; False=>Regional
!
!-----------
!*** Output
!-----------
!
      REAL,DIMENSION(1:IM,1:JM),INTENT(OUT)      :: V_VALUE                !<-- Output values on mass (V) points
!
!----------
      INTEGER LM
      real,dimension(:,:,:),allocatable          :: tmp
!
!----------
      LM=1
      allocate(tmp(1:IM,1:JM,1))
      call H_TO_V_BGRID3D(reshape(H_VALUE,(/IM,JM,1/)),IM,JM,LM,GLOBAL,   &
        V_SAVE, tmp )
      V_VALUE(1:IM,1:JM)=tmp(1:IM,1:JM,1)
      deallocate(tmp)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE H_TO_V_BGRID2D
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE H_TO_V_BGRID3D(H_VALUE,IM,JM,LM,GLOBAL,V_SAVE,V_VALUE)
!
!-----------------------------------------------------------------------
!***  PERFORM 4-POINT AVERAGING OF QUANTITIES ON B-GRID MASS POINTS
!***  TO VELOCITY POINTS.  THIS ROUTINE IS INTENDED FOR USE IN
!***  INTERPOLATING VELOCITY VALUES BACK TO V POINTS AFTER HAVING
!***  BEEN INTERPOLATED ONTO H POINTS.  THE INPUT POINTER V_SAVE
!***  IS USED TO RECOVER THE VELOCITY VALUES ON THE BOUNDARY V POINTS
!***  WHEN IN REGIONAL MODE BECAUSE THEY CANNOT BE RECOVERED
!***  BY THE INTERPOLATION FROM H POINTS GIVEN THAT BOUNDARY H POINTS
!***  IN REGIONAL MODE COULD NOT BE OBTAINED BY INTERPOLATIONS FROM
!***  THE V POINTS.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!----------
!*** Input
!----------
!
      INTEGER                       ,INTENT(IN)  :: IM                  &  !<-- Full west-east array dimension
                                                   ,JM                  &  !<-- Full south-north array dimension
                                                   ,LM                     !<-- Number of model layers
!
      REAL,DIMENSION(1:IM,1:JM,1:LM),INTENT(IN)  :: H_VALUE                !<-- Input values on mass points
      REAL,DIMENSION(:,:),POINTER   ,INTENT(IN)  :: V_SAVE                 !<-- Saved values on boundary velocity points
!
      LOGICAL                       ,INTENT(IN)  :: GLOBAL                 !<-- Logical flag: True=>Global; False=>Regional
!
!-----------
!*** Output
!-----------
!
      REAL,DIMENSION(1:IM,1:JM,1:LM),INTENT(OUT) :: V_VALUE                !<-- Output values on mass (V) points
!
!----------
!*** Local
!----------
!
      INTEGER                                    :: I,J,L,N
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!***  NOTE THE NATURE OF THE BOUNDARY OVERLAP FOR THE GLOBAL MODE.
!
!***  H(1) and H(IM-2) COINCIDE.
!***  H(2) and H(IM-1) COINCIDE.
!***  H(3) and H(IM) COINCIDE.
!
!-----------------------------------------------------------------------
!
!                                               |
!                                               |<---- Integration western boundary for H
!                                               |
!
!                                                   |
!                                                   |<---- Integration western boundary for V
!                                                   |
!               
!                                         V(1)      V(2)      V(3)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!                              V(IM-3)    V(IM-2)   V(IM-1)   V(IM)
!
!                                         |
! Integration eastern boundary for V ---->|
!                                         |
!
!                                               |
!       Integration eastern boundary for H ---->|
!                                               |
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!           J=JM       V    V    V
!
!        J=JM       H    H    H    H      <---  1 row "north" of pole
!                                               is 1 row south.
!           J=JM-1     V    V    V
!
!        J=JM-1     H    H    H    H      <---  North Pole: Integration boundary for H
!
!           J=JM-2     V    V    V     <--- Integration boundary for V
!
!        J=JM-2     H    H    H    H      <---  1 row south of pole
!
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + 
!
!          J=3      H    H    H    H      <---  1 row north of pole
!
!             J=2      V    V    V     <--- Integration boundary for V
!
!          J=2      H    H    H    H      <---  South Pole:  Integration boundary for H
!
!             J=1      V    V    V     <--- Reflection of J=2
!
!          J=1      H    H    H    H      <---  1 row "south" of pole
!                                               is 1 row north.
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  PERFORM THE INTERPOLATION ON THE INTERNAL REGION OF THE
!***  OUTPUT ARRAY.  IT IS APPLICABLE TO BOTH REGIONAL AND
!***  GLOBAL MODES.
!***  WE CANNOT REACH IM-1 AND JM-1 ON V POINTS BECAUSE VALUES
!***  AT IM AND JM ON H POINTS CANNOT BE COMPUTED FROM V POINTS
!***  WHEN IN REGIONAL MODE.
!-----------------------------------------------------------------------
!
      DO L=1,LM
      DO J=2,JM-2
      DO I=2,IM-2
!
        V_VALUE(I,J,L)=0.25*(H_VALUE(I  ,J  ,L)                         &
                            +H_VALUE(I+1,J  ,L)                         &
                            +H_VALUE(I  ,J+1,L)                         &
                            +H_VALUE(I+1,J+1,L) )
!
      ENDDO
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  FOR GLOBAL MODE THE OVERLAP VALUES ON VELOCITY POINTS WILL BE
!***  USED TO INTERPOLATE TO THE BOUNDARY MASS POINTS.
!-----------------------------------------------------------------------
!
      layers: DO L=1,LM
!
        mode: IF(GLOBAL)THEN
!
!--------------------------
!***  West/East Boundaries
!--------------------------
!
          DO J=1,JM-1
!
            V_VALUE(   1,J,L)=0.25*(H_VALUE(1,J  ,L)                    & !<-- West
                                   +H_VALUE(2,J  ,L)                    &
                                   +H_VALUE(1,J+1,L)                    &
                                   +H_VALUE(2,J+1,L) )
!
            V_VALUE(IM-1,J,L)=0.25*(H_VALUE(2,J  ,L)                    & !<-- East - 1
                                   +H_VALUE(3,J  ,L)                    &
                                   +H_VALUE(2,J+1,L)                    &
                                   +H_VALUE(3,J+1,L) )
!
            V_VALUE(IM,J,L)=0.25*(H_VALUE(3,J  ,L)                      & !<-- East 
                                 +H_VALUE(4,J  ,L)                      &
                                 +H_VALUE(3,J+1,L)                      &
                                 +H_VALUE(4,J+1,L) )
!
          ENDDO
!
!----------------------------
!***  South/North Boundaries
!----------------------------
!
          DO I=1,IM-1
            V_VALUE(I,   1,L)=V_VALUE(I,   2,L)                            !<-- South
            V_VALUE(I,JM-1,L)=V_VALUE(I,JM-2,L)                            !<-- North-1 reflected around North Pole
            V_VALUE(I,JM  ,L)=V_VALUE(I,JM  ,L)                            !<-- North reflected around North Pole
          ENDDO
!
!-----------------------------------------------------------------------
!***  FOR REGIONAL MODE USE THE SAVED ORIGINAL BOUNDARY VALUES
!***  SINCE THOSE BOUNDARY VALUES CANNOT BE OBTAINED THROUGH
!***  INTEPOLATION FROM THE H POINTS.
!-----------------------------------------------------------------------
!
        ELSE mode
!
!--------------------
!***  South Boundary
!--------------------
!
          N=0
!
          DO I=1,IM-1
            V_VALUE(I, 1,L)=V_SAVE(N,L)
          ENDDO
!
!----------------------------
!***  North Boundary
!----------------------------
!
          DO I=1,IM-1
            V_VALUE(I,JM-1,L)=V_SAVE(N,L)                                 !<-- Northernmost V points in integration
          ENDDO
!
!--------------------------
!***  West Boundary
!--------------------------
!
          DO J=2,JM-1
            V_VALUE( 1,J,L)=V_SAVE(N,L)
          ENDDO
!
!--------------------------
!***  East Boundary
!--------------------------
!
          DO J=2,JM-1
            V_VALUE(IM-1,J,L)=V_SAVE(N,L)                                 !<-- Easternmost V points in integration
          ENDDO
!
!-----------------------------------------------------------------------
!***  ZERO OUT PHANTOM LOCATIONS.
!-----------------------------------------------------------------------
!
          DO J=1,JM
            V_VALUE(IM,J,L)=0.
          ENDDO
!
          DO I=1,IM
            V_VALUE(I,JM,L)=0.
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDIF mode
!
      ENDDO layers
!
!-----------------------------------------------------------------------
!***  ALL VELOCITY POINT VALUES AT I=IM and J=JM ARE MEANINGLESS
!***  BECAUSE THEY ARE OUTSIDE THE INTEGRATION DOMAIN.
!-----------------------------------------------------------------------
!
      END SUBROUTINE H_TO_V_BGRID3D
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_BGRID_INTERP
!
!-----------------------------------------------------------------------
