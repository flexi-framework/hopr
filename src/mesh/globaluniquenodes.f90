!============================================================================================================ xX =================
!        _____     _____    _______________    _______________   _______________            .xXXXXXXXx.       X
!       /    /)   /    /)  /    _____     /)  /    _____     /) /    _____     /)        .XXXXXXXXXXXXXXx  .XXXXx
!      /    //   /    //  /    /)___/    //  /    /)___/    // /    /)___/    //       .XXXXXXXXXXXXXXXXXXXXXXXXXx
!     /    //___/    //  /    //   /    //  /    //___/    // /    //___/    //      .XXXXXXXXXXXXXXXXXXXXXXX`
!    /    _____     //  /    //   /    //  /    __________// /    __      __//      .XX``XXXXXXXXXXXXXXXXX`
!   /    /)___/    //  /    //   /    //  /    /)_________) /    /)_|    |__)      XX`   `XXXXX`     .X`
!  /    //   /    //  /    //___/    //  /    //           /    //  |    |_       XX      XXX`      .`
! /____//   /____//  /______________//  /____//           /____//   |_____/)    ,X`      XXX`
! )____)    )____)   )______________)   )____)            )____)    )_____)   ,xX`     .XX`
!                                                                           xxX`      XXx
! Copyright (C) 2017 Claus-Dieter Munz <munz@iag.uni-stuttgart.de>
! This file is part of HOPR, a software for the generation of high-order meshes.
!
! HOPR is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! HOPR is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with HOPR. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "hopr.h"
MODULE MOD_GlobalUniqueNodes
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE GlobalUniqueNodes
  MODULE PROCEDURE GlobalUniqueNodes
END INTERFACE

PUBLIC::GlobalUniqueNodes
!===================================================================================================================================

CONTAINS

SUBROUTINE GlobalUniqueNodes(withOrientedOpt)
!===================================================================================================================================
! Eliminates multiple nodes, checks periodic boundary conditions and connects elements to their neighbours.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:tElem,tSide,tEdge,tNode,tNodePtr,FirstElem
USE MOD_Mesh_Vars, ONLY:N,deleteNode
USE MOD_Mesh_Vars, ONLY:SpaceQuandt
USE MOD_Mesh_Tools,ONLY:SetTempMarker
USE MOD_Mesh_Tolerances,ONLY:COMPAREPOINT
USE MOD_SpaceFillingCurve,ONLY:EVAL_MORTON,EVAL_MORTON_ARR
USE MOD_SortingTools,ONLY: Qsort1DoubleInt1Pint  
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,OPTIONAL,INTENT(IN) :: withOrientedOpt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER         :: Elem  ! ?
TYPE(tSide),POINTER         :: Side  ! ?
TYPE(tEdge),POINTER         :: Edge  ! ?
TYPE(tNode),POINTER         :: Node  ! ?
TYPE(tNodePtr),POINTER      :: Nodes(:)  ! ?
INTEGER                     :: i,iNode,jNode,NodeID  ! ?
INTEGER                     :: nTotalNodes,nDeletedNodes  ! ?
REAL                        :: tol   ! ?
REAL                        :: sLog2   ! ?
REAL                        :: box_min,box_max,box_sdx  ! ?
INTEGER                     :: box_nBits  ! ?
INTEGER                     :: maxStepsINVMAP  ! ?
INTEGER(KIND=8)             :: box_di,s_offset,maxIJK  ! ?
INTEGER(KIND=8)             :: s_minmax(8,2)  ! ?
INTEGER(KIND=8)             :: smin,smax  ! ?
INTEGER(KIND=8),ALLOCATABLE :: NodesIJK(:,:),SFCID(:)  ! ?
INTEGER,ALLOCATABLE         :: IDList(:)  ! ?
INTEGER                     :: nRanges  ! ?
INTEGER                     :: percent  ! ?
INTEGER                     :: lastNode,nextNode  ! ?
LOGICAL                     :: withOriented
LOGICAL,PARAMETER           :: T=.TRUE.
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'GLOBAL UNIQUE NODES ...'
sLog2=1./LOG(2.)
withOriented=.FALSE.
IF(PRESENT(withOrientedOpt)) withOriented=withOrientedOpt

! First step: set node marker=0 
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  CALL SetTempMarker(Elem,0,(/T,T,T,T,T,T,withOriented,T/))
  Elem=>Elem%nextElem
END DO

! Second step: set unique node marker and count 
NodeID=0
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  DO iNode=1,Elem%nNodes
    CALL SetCountNodeID(Elem%Node(iNode)%np%tmp,NodeID)
  END DO !iNodes
  IF(ASSOCIATED(Elem%CurvedNode))THEN
    DO iNode=1,Elem%nCurvedNodes
      CALL SetCountNodeID(Elem%curvedNode(iNode)%np%tmp,NodeID)
    END DO
  END IF
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    DO iNode=1,Side%nNodes
      CALL SetCountNodeID(Side%Node(iNode)%np%tmp,NodeID)
      IF(withOriented) CALL SetCountNodeID(Side%OrientedNode(iNode)%np%tmp,NodeID)
      IF(ASSOCIATED(Side%edge(iNode)%edp))THEN
        Edge=>Side%edge(iNode)%edp
        CALL SetCountNodeID(Edge%Node(1)%np%tmp,NodeID)
        CALL SetCountNodeID(Edge%Node(2)%np%tmp,NodeID)
        IF(ASSOCIATED(Edge%CurvedNode))THEN
          DO i=1,N+1
            CALL SetCountNodeID(Edge%curvedNode(i)%np%tmp,NodeID)
          END DO
        END IF
      END IF
    END DO
    DO iNode=1,Side%nCurvedNodes
      CALL SetCountNodeID(Side%curvedNode(iNode)%np%tmp,NodeID)
    END DO
    !periodic side only /= for connect!
    IF(Side%tmp2.GT.0)THEN
      !dummy side found
      DO iNode=1,Side%nNodes
        Node=>Side%connection%Node(iNode)%np
        CALL SetCountNodeID(Node%tmp,NodeID)
        Node%ind=SIGN(Node%ind,-1) ! ensure that ind of periodics is negative
      END DO
    END IF
    Side=>Side%nextElemSide
  END DO !associated(side)
  Elem=>Elem%nextElem
END DO !associated(elem)

nTotalNodes=NodeID

ALLOCATE(Nodes(nTotalNodes))
DO iNode=1,nTotalNodes
  NULLIFY(Nodes(iNode)%np)
END DO !iNode

!Associate nodelist
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  DO iNode=1,Elem%nNodes
    Nodes(Elem%Node(iNode)%np%tmp)%np=>Elem%Node(iNode)%np
  END DO !iNodes
  IF(ASSOCIATED(Elem%CurvedNode))THEN
    DO iNode=1,Elem%nCurvedNodes
      Nodes(Elem%curvedNode(iNode)%np%tmp)%np=>Elem%curvedNode(iNode)%np
    END DO
  END IF
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    DO iNode=1,Side%nNodes
      Nodes(Side%Node(iNode)%np%tmp)%np=>Side%Node(iNode)%np
      IF(withOriented) Nodes(Side%OrientedNode(iNode)%np%tmp)%np=>Side%OrientedNode(iNode)%np
      IF(ASSOCIATED(Side%edge(iNode)%edp))THEN
        Edge=>Side%edge(iNode)%edp
        Nodes(Edge%Node(1)%np%tmp)%np=>Edge%Node(1)%np
        Nodes(Edge%Node(2)%np%tmp)%np=>Edge%Node(2)%np
        IF(ASSOCIATED(Edge%CurvedNode))THEN
          DO i=1,N+1
            Nodes(Edge%curvedNode(i)%np%tmp)%np=>Edge%curvedNode(i)%np
          END DO
        END IF
      END IF
    END DO
    DO iNode=1,Side%nCurvedNodes
      Nodes(Side%curvedNode(iNode)%np%tmp)%np=>Side%curvedNode(iNode)%np
    END DO
    !periodic side only /= for connect!
    IF(Side%tmp2.GT.0)THEN
      !dummy side found
      DO iNode=1,Side%nNodes
        Nodes(Side%connection%Node(iNode)%np%tmp)%np=>Side%connection%Node(iNode)%np
      END DO
    END IF
    Side=>Side%nextElemSide
  END DO !associated(side)
  Elem=>Elem%nextElem
END DO !associated(elem)


! bounding cube (all sfc stuff is from routine spacefillingcurve.f90)
box_min=1.0E16
box_max=-1.0E16
DO iNode=1,nTotalNodes
  box_min=MIN(box_min,MINVAL(Nodes(iNode)%np%x(:)))
  box_max=MAX(box_max,MAXVAL(Nodes(iNode)%np%x(:)))
END DO !iNode
box_min=box_min-PP_MeshTolerance
box_max=box_max+PP_MeshTolerance
box_nbits = (bit_size(maxIJK)-1) / 3 
maxIJK = 2**box_nbits-1               ![0,2**box_nBits-1]
box_sdx = REAL(maxIJK)/(box_max-box_min)


ALLOCATE(NodesIJK(3,nTotalNodes),SFCID(nTotalNodes),IDList(nTotalNodes))
DO iNode=1,nTotalNodes
  NodesIJK(:,iNode)=NINT((Nodes(iNode)%np%x-box_min)*box_sdx)
  IDList(iNode)=iNode
END DO !iNode
! evaluate morton curve
CALL EVAL_MORTON_ARR(SFCID,NodesIJK,nTotalNodes,box_nBits)

CALL Qsort1DoubleInt1Pint(SFCID, IDList) !sort SFCID  and IDlist like SFCID!!

!resort by IDlist
NodesIJK(1,:)=NodesIJK(1,IDList(:))
NodesIJK(2,:)=NodesIJK(2,IDList(:))
NodesIJK(3,:)=NodesIJK(3,IDList(:))

DO iNode=1,nTotalNodes
  Nodes(IDList(iNode))%np%tmp=iNode
END DO !iNode
DEALLOCATE(IDList)

DO iNode=1,nTotalNodes
  NULLIFY(Nodes(iNode)%np)
END DO !iNode

!Associate nodelist again, now sorted by morton curve
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  DO iNode=1,Elem%nNodes
    Nodes(Elem%Node(iNode)%np%tmp)%np=>Elem%Node(iNode)%np
  END DO !iNodes
  IF(ASSOCIATED(Elem%CurvedNode))THEN
    DO iNode=1,Elem%nCurvedNodes
      Nodes(Elem%curvedNode(iNode)%np%tmp)%np=>Elem%curvedNode(iNode)%np
    END DO
  END IF
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    DO iNode=1,Side%nNodes
      Nodes(Side%Node(iNode)%np%tmp)%np=>Side%Node(iNode)%np
      IF(withOriented) Nodes(Side%OrientedNode(iNode)%np%tmp)%np=>Side%OrientedNode(iNode)%np
      IF(ASSOCIATED(Side%edge(iNode)%edp))THEN
        Edge=>Side%edge(iNode)%edp
        Nodes(Edge%Node(1)%np%tmp)%np=>Edge%Node(1)%np
        Nodes(Edge%Node(2)%np%tmp)%np=>Edge%Node(2)%np
        IF(ASSOCIATED(Edge%CurvedNode))THEN
          DO i=1,N+1
            Nodes(Edge%curvedNode(i)%np%tmp)%np=>Edge%curvedNode(i)%np
          END DO
        END IF
      END IF
    END DO
    DO iNode=1,Side%nCurvedNodes
      Nodes(Side%curvedNode(iNode)%np%tmp)%np=>Side%curvedNode(iNode)%np
    END DO
    !periodic side only /= for connect!
    IF(Side%tmp2.GT.0)THEN
      !dummy side found
      DO iNode=1,Side%nNodes
        Nodes(Side%connection%Node(iNode)%np%tmp)%np=>Side%connection%Node(iNode)%np
      END DO
    END IF
    Side=>Side%nextElemSide
  END DO !associated(side)
  Elem=>Elem%nextElem
END DO !associated(Elem)

DO iNode=1,nTotalNodes
  Nodes(iNode)%np%tmp=0
END DO !iNode
WRITE(*,*)'  All Nodes sorted...'
WRITE(*,*)'  Number of nodes to check: ',nTotalNodes
!================== Preparation done, now search ====================================
maxStepsINVMAP=INT(LOG(REAL(nTotalNodes))*sLog2)+1
nDeletedNodes=0

tol=SpaceQuandt*PP_MeshTolerance
 ! Size of tolerance gives a box_di for bisection (could be computed for each node seperately !!!)
!box_di  =MAX(0,box_nBits-FLOOR(LOG((box_max-box_min)/tol)*sLog2)) ! tol=2^x L=2^n => x=n-LOG(L/tol)/LOG(2) 
box_di  =MAX(0,CEILING(LOG(tol*box_sdx)*sLog2)) 
box_di=2**box_di
WRITE(*,*)'   size of tolerance box:',box_di
s_offset=box_di**3-1   !offset inside one box of size box_di, from the lower sfc index to to highest

percent=0
nextNode=1
DO iNode=1,nTotalNodes
  ! output of progress in %
  IF((nTotalNodes.GT.100000).AND.(MOD(iNode,(nTotalNodes/100)).EQ.0)) THEN
    percent=percent+1
    WRITE(0,'(I4,A23,A1)',ADVANCE='NO')percent, ' % of nodes evaluated...',ACHAR(13)
  END IF
  Node=>Nodes(iNode)%np
  IF(Node%tmp.GT.0) CYCLE ! node already checked
  Node%tmp=iNode !check this node

!WRITE(*,*)'============= iNODE=',iNode,'================'
!WRITE(*,'(A,4I)')'node (i,j,k,s)',NodesIJK(:,iNode),SFCID(iNode)
!WRITE(*,*)'==============='

  CALL FindBoxes(NodesIJK(:,iNode),box_di,box_nBits,maxIJK,s_minmax,smin,smax,nRanges)
  !next higher neighbor on SFC
  lastNode=nextNode
  DO nextNode=lastNode+1,nTotalNodes 
    IF(Nodes(nextNode)%np%tmp.EQ.0) EXIT ! not yet treated
  END DO
  nextNode=MIN(nextNode,nTotalNodes)
  IF(SFCID(nextNode).GT.smax) CYCLE 
  
  DO i=1,nRanges
     IF(s_minmax(i,1).EQ.-1)CYCLE
     IF(SFCID(iNode).GT.s_minmax(i,2)) CYCLE 
     IF(SFCID(nextNode).GT.s_minmax(i,2)) CYCLE 
     NodeID=INVMAP(s_minmax(i,1),nTotalNodes-(iNode-1),SFCID(iNode:nTotalNodes))
     IF(NodeID.EQ.-1) CYCLE !nothing found inside the box
     NodeID=MAX(nextNode,NodeID+iNode)
     DO jNode=NodeID,nTotalNodes 
       IF(SFCID(jNode).GT.s_minmax(i,2)) EXIT ! check if  > s_max
       IF(Nodes(jNode)%np%tmp.GT.0) CYCLE ! was already treated

       IF(COMPAREPOINT(Node%x,Nodes(jNode)%np%x,tol))THEN
         ! DOUBLE NODE FOUND
         nDeletedNodes=nDeletedNodes+1
         Nodes(jNode)%np%tmp=Node%tmp ! set pointer to unique node
         Nodes(jNode)%np%ind=MAX(Node%ind,Nodes(jNode)%np%ind) ! set unique node ID
         Node%ind=MAX(Node%ind,Nodes(jNode)%np%ind) ! set unique node ID
       END IF
       ! check next node
     END DO !NodeID +
  END DO

END DO !iNode
WRITE(*,*)' Number of deleted nodes',nDeletedNodes
WRITE(*,*)' Number of unique nodes',nTotalNodes-nDeletedNodes

!Associate nodelist back, now sorted by morton curve
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  DO iNode=1,Elem%nNodes
    Elem%Node(iNode)%np=>Nodes(Elem%Node(iNode)%np%tmp)%np
  END DO !iNodes
  IF(ASSOCIATED(Elem%CurvedNode))THEN
    DO iNode=1,Elem%nCurvedNodes
      Elem%curvedNode(iNode)%np=>Nodes(Elem%curvedNode(iNode)%np%tmp)%np
    END DO
  END IF
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    DO iNode=1,Side%nNodes
      Side%Node(iNode)%np=>Nodes(Side%Node(iNode)%np%tmp)%np
      IF(withOriented) Side%OrientedNode(iNode)%np=>Nodes(Side%OrientedNode(iNode)%np%tmp)%np
      IF(ASSOCIATED(Side%edge(iNode)%edp))THEN
        Edge=>Side%edge(iNode)%edp
        Edge%Node(1)%np=>Nodes(Edge%Node(1)%np%tmp)%np
        Edge%Node(2)%np=>Nodes(Edge%Node(2)%np%tmp)%np
        IF(ASSOCIATED(Edge%CurvedNode))THEN
          DO i=1,N+1
            Edge%curvedNode(i)%np=>Nodes(Edge%curvedNode(i)%np%tmp)%np
          END DO
        END IF
      END IF
    END DO
    DO iNode=1,Side%nCurvedNodes
      Side%curvedNode(iNode)%np=>Nodes(Side%curvedNode(iNode)%np%tmp)%np
    END DO
    !periodic side only /= for connect!
    IF(Side%tmp2.GT.0)THEN
      !dummy side found
      DO iNode=1,Side%nNodes
        Side%connection%Node(iNode)%np=>Nodes(Side%connection%Node(iNode)%np%tmp)%np
      END DO
    END IF
    Side=>Side%nextElemSide
  END DO !associated(side)
  Elem=>Elem%nextElem
END DO !associated(Elem)

DO iNode=1,nTotalNodes
  IF(Nodes(iNode)%np%tmp.NE.iNode) CALL deleteNode(Nodes(iNode)%np)
  NULLIFY(Nodes(iNode)%np)
END DO
DEALLOCATE(Nodes,NodesIJK,SFCID)
CALL Timer(.FALSE.)
END SUBROUTINE GlobalUniqueNodes


SUBROUTINE FindBoxes(IntCoord,box_di,box_nBits,maxIJK,s_minmax,smin,smax,nRanges)
!===================================================================================================================================
! finds the ranges of the spacefilling curve which correspond to a maximum of 8 boxes with a box size of 2*box_id,
! box_di (from node merging tolerance) defines the smallest box size. We now look at two levels below the octree  (4x4x4 box size)
! where node lies inside. We need a tolerance of one box  to find all possible nodes(makes 3x3x3 boxes), 
! but due to the octree, we always choose a 4x4x4 region, allowing a maximum of only 8 search boxes of size 2x2x2. 
! Due to the nature of the morton spacefilling curve, some boxes have  contiguous ranges and are merged. There are three cases:
! 1 box of 4x4x4 is used if the node is inside
!                    
! z       x0                  x1                 x2                  x3        
! ^  ___ ___ ___ ___    ___ ___ ___ ___     ___ ___ ___ ___    ___ ___ ___ ___ 
! | |   |   |   |   |  |   |   |   |   |   |   |   |   |   |  |   |   |   |   |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! | |   |   |   |   |  |   | . | . |   |   |   | . | . |   |  |   |   |   |   |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! | |   |   |   |   |  |   | . | . |   |   |   | . | . |   |  |   |   |   |   |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! | |   |   |   |   |  |   |   |   |   |   |   |   |   |   |  |   |   |   |   |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! *----------------------> y
! 2 boxes of 2x4x4,  if nodes are inside yz and +/-x 
!                    
! z       x0                  x1                 x2                  x3        
! ^  ___ ___ ___ ___    ___ ___ ___ ___     ___ ___ ___ ___    ___ ___ ___ ___ 
! | |   |   |   |   |  |   |   |   |   |   |   |   |   |   |  |   |   |   |   |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! | |   | . | . |   |  |   |   |   |   |   |   |   |   |   |  |   | . | . |   |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! | |   | . | . |   |  |   |   |   |   |   |   |   |   |   |  |   | . | . |   |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! | |   |   |   |   |  |   |   |   |   |   |   |   |   |   |  |   |   |   |   |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! *----------------------> y
! 4 boxes of 2x2x4, if nodes are inside z and +/-x +/- y
!                    
! z       x0                  x1                 x2                  x3        
! ^  ___ ___ ___ ___    ___ ___ ___ ___     ___ ___ ___ ___    ___ ___ ___ ___ 
! | |   |   |   |   |  |   |   |   |   |   |   |   |   |   |  |   |   |   |   |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! | | . |   |   | . |  |   |   |   |   |   |   |   |   |   |  | . |   |   | . |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! | | . |   |   | . |  |   |   |   |   |   |   |   |   |   |  | . |   |   | . |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! | |   |   |   |   |  |   |   |   |   |   |   |   |   |   |  |   |   |   |   |
! | |___|___|___|___|  |___|___|___|___|   |___|___|___|___|  |___|___|___|___|
! *----------------------> y
!
!or then 8 boxes of 2x2x2.
! This minimizes the number of overall spacefilling curve evaluations.
!===================================================================================================================================
! MODULES
USE MOD_SpaceFillingCurve,ONLY:EVAL_MORTON_ARR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)  :: IntCoord(3) ! ijk integer coordinates
INTEGER(KIND=8),INTENT(IN)  :: box_di      ! box size in integer counts
INTEGER,INTENT(IN)          :: box_nBits   ! number of bits of INTEGER(KIND=8)
INTEGER(KIND=8),INTENT(IN)  :: maxIJK      ! maximum domain size in integer counts (=2**box_nBits-1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(OUT) :: s_minmax(8,2)  ! sfc index ranges of the 27 boxes
INTEGER(KIND=8),INTENT(OUT) :: smin            ! minimum of s_minmax(:,1)
INTEGER(KIND=8),INTENT(OUT) :: smax            ! maximum of s_minmax(:,2)
INTEGER,INTENT(OUT)         :: nRanges          ! number of ranges 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER          :: i,j,k,l  ! ?
INTEGER          :: Ranges(3),start(3)  ! ?
LOGICAL          :: xi0(3),xi1(3)  ! ?
INTEGER(KIND=8)  :: IJK_min(3),IJK_tmp(3,8)  ! ?
INTEGER(KIND=8)  :: s_offset,box_di2,box_di4,dxi(3)
!===================================================================================================================================
box_di2=box_di+box_di
box_di4=box_di2+box_di2

IJK_min(:)=box_di4*(IntCoord(:)/(box_di4))
xi0(:)=(MOD(IntCoord(:)/box_di2,2).NE.0)
xi1(:)=(MOD(IntCoord(:)/box_di,2).NE.0)

Ranges=2
start=0
IF(.NOT.(xi0(1).NEQV.xi1(1))) start(1)=MERGE(1,-1,xi0(1))
IF(.NOT.(xi0(2).NEQV.xi1(2))) start(2)=MERGE(1,-1,xi0(2))
IF(.NOT.(xi0(3).NEQV.xi1(3))) start(3)=MERGE(1,-1,xi0(3))

IF(xi0(3).NEQV.xi1(3)) THEN
  Ranges(3)=1
  IF(xi0(2).NEQV.xi1(2)) THEN
    Ranges(2)=1
    IF(xi0(1).NEQV.xi1(1))THEN 
      Ranges(1)=1
    END IF
  END IF
END IF

dxi=box_di4/Ranges
nRanges=0
DO i=1,ranges(1); DO j=1,ranges(2); DO k=1,ranges(3)
  nRanges=nRanges+1
  !min
  IJK_tmp(1,nRanges)= IJK_min(1)+(start(1)+(i-1))*dxi(1)
  IJK_tmp(2,nRanges)= IJK_min(2)+(start(2)+(j-1))*dxi(2)
  IJK_tmp(3,nRanges)= IJK_min(3)+(start(3)+(k-1))*dxi(3)
END DO; END DO; END DO !i,j,k
s_offset = PRODUCT(dxi)-1


s_minmax(:,:)=-1
CALL EVAL_MORTON_ARR(s_minmax(1:nRanges,1),IJK_tmp(:,1:nRanges),nRanges,box_nBits)
s_minmax(1:nRanges,2)=s_minmax(1:nRanges,1)+s_offset
smin=MINVAL(s_minmax(1:nRanges,1))
smax=MAXVAL(s_minmax(1:nRanges,2))
!boundaries
DO l=1,nRanges
  IF((MINVAL(IJK_tmp(:,l)).LT.0).OR.(MAXVAL(IJK_tmp(:,l)).GT.maxIJK)) s_minmax(l,:)=-1
END DO

END SUBROUTINE FindBoxes

SUBROUTINE SetCountNodeID(NodeID_in,NodeID)
!===================================================================================================================================
! insert a new node id 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT) :: NodeID_in,NodeID  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF(NodeID_in.EQ.0)THEN
  NodeID=NodeID+1
  NodeID_in=NodeID
END IF
END SUBROUTINE SetCountNodeID

FUNCTION INVMAP(ID,nIDs,ArrID)
!===================================================================================================================================
! find the inverse Mapping of sfc index in a sorted list (a sorted array of unique NodeIDs), using bisection
! if Index is not in the range, -1 will be returned, gives back the first entry in the sorted list  which is <= ID 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8), INTENT(IN) :: ID            ! ID to search for
INTEGER, INTENT(IN)         :: nIDs          ! size of ArrID
INTEGER(KIND=8), INTENT(IN) :: ArrID(nIDs)   ! 1D array of IDs, SORTED!!
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                     :: INVMAP         ! position in arrID, where arrID(INVMAP)  <= ID 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: i,low,up,mid  ! ?
INTEGER                     :: maxSteps      ! =INT(LOG(REAL(nIDs))/LOG(2.))+1
!===================================================================================================================================
IF(ID.LE.ArrID(1))THEN
  INVMAP=1
  RETURN
ELSEIF(ID.EQ.ArrID(nIDs))THEN
  INVMAP=nIDs
  RETURN
ELSEIF(ID.GT.ArrID(nIDs))THEN
  INVMAP=-1 !nothing found!
  RETURN
END IF
low=1
up=nIDs
!bisection
maxSteps=INT(LOG(REAL(nIDs))*1.442695040888964)+1
DO i=1,maxSteps
  mid=(up-low)/2+low
  IF(ArrID(mid).GE.ID )THEN !seek in lower half 
    up=mid
  ELSE
    low=mid
  END IF
END DO
INVMAP=low
END FUNCTION INVMAP 

END MODULE MOD_GlobalUniqueNodes
