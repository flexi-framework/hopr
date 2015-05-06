#include "defines.f90"
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

SUBROUTINE GlobalUniqueNodes()
!===================================================================================================================================
! Eliminates multiple nodes, checks periodic boundary conditions and connects elements to their neighbours.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:tElem,tSide,tEdge,tNode,tNodePtr,FirstElem
USE MOD_Mesh_Vars, ONLY:N
USE MOD_Mesh_Vars,ONLY:SpaceQuandt
USE MOD_Mesh_Tolerances,ONLY:COMPAREPOINT
USE MOD_SpaceFillingCurve,ONLY:EVAL_MORTON,EVAL_MORTON_ARR
USE MOD_SortingTools,ONLY: Qsort1DoubleInt1Pint  
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
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
INTEGER(KIND=8)             :: s_minmax(27,2)  ! ?
INTEGER(KIND=8),ALLOCATABLE :: NodesIJK(:,:),SFCID(:)  ! ?
INTEGER,ALLOCATABLE         :: IDList(:)  ! ?
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'GLOBAL UNIQUE NODES ...'
sLog2=1./LOG(2.)
! First step: set node marker=0 
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  DO iNode=1,Elem%nNodes
    Elem%Node(iNode)%np%tmp = 0
  END DO !iNodes
  IF(ASSOCIATED(Elem%CurvedNode))THEN
    DO iNode=1,Elem%nCurvedNodes
      Elem%curvedNode(iNode)%np%tmp = 0
    END DO
  END IF
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    DO iNode=1,Side%nNodes
      Side%Node(iNode)%np%tmp=0
      IF(ASSOCIATED(Side%edge(iNode)%edp))THEN
        Edge=>Side%edge(iNode)%edp
        Edge%Node(1)%np%tmp=0
        Edge%Node(2)%np%tmp=0
        IF(ASSOCIATED(Edge%CurvedNode))THEN
          DO i=1,N+1
            Edge%curvedNode(i)%np%tmp=0
          END DO
        END IF
      END IF
    END DO
    DO iNode=1,Side%nCurvedNodes
      Side%curvedNode(iNode)%np%tmp=0
    END DO
    !periodic side, only /= for connect!
    IF(Side%tmp2.GT.0)THEN
      !dummy side found
      DO iNode=1,Side%nNodes
        Side%connection%Node(iNode)%np%tmp=0
      END DO
    END IF
    Side=>Side%nextElemSide
  END DO !associated(side)
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
        CALL SetCountNodeID(Side%connection%Node(iNode)%np%tmp,NodeID)
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
box_min=box_min-RealTolerance
box_max=box_max+RealTolerance
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

tol=SpaceQuandt*RealTolerance
 ! Size of tolerance gives a box_di for bisection (could be computed for each node seperately !!!)
!box_di  =MAX(0,box_nBits-FLOOR(LOG((box_max-box_min)/tol)*sLog2)) ! tol=2^x L=2^n => x=n-LOG(L/tol)/LOG(2) 
box_di  =MAX(0,CEILING(LOG(tol*box_sdx)*sLog2)) 
box_di=2**box_di
WRITE(*,*)'   size of tolerance box:',box_di
s_offset=box_di**3-1  !offset inside one box of size box_di, from the lower sfc index to to highest


DO iNode=1,nTotalNodes
  Node=>Nodes(iNode)%np
  IF(Node%tmp.GT.0) CYCLE ! node already checked
  Node%tmp=iNode !check this node

!WRITE(*,*)'============= iNODE=',iNode,'================'
!WRITE(*,'(A,4I)')'node (i,j,k,s)',NodesIJK(:,iNode),SFCID(iNode)
!WRITE(*,*)'==============='

  CALL FindBoxes(NodesIJK(:,iNode),box_di,s_offset,box_nBits,maxIJK,s_minmax)

  
  DO i=1,27
     IF(s_minmax(i,1).EQ.-1)CYCLE
     NodeID=INVMAP(s_minmax(i,1),nTotalNodes,SFCID,maxStepsINVMAP)
     IF(NodeID.EQ.-1) CYCLE !nothing found inside the box
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
  IF(Nodes(iNode)%np%tmp.NE.iNode) DEALLOCATE(Nodes(iNode)%np)
  NULLIFY(Nodes(iNode)%np)
END DO
DEALLOCATE(Nodes,NodesIJK,SFCID)
CALL Timer(.FALSE.)
END SUBROUTINE GlobalUniqueNodes

SUBROUTINE FindBoxes(IntCoord,box_di,s_offset,box_nBits,maxIJK,s_minmax)
!===================================================================================================================================
! finds the ranges of the spacefilling curve which correspond to the 27 boxes with a box size of box_id. 
! The list is sorted by the sfc index of the boxes, and contiguous boxes are also merged
!===================================================================================================================================
! MODULES
USE MOD_SortingTools,ONLY:Qsort1Int2double
USE MOD_SpaceFillingCurve,ONLY:EVAL_MORTON_ARR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)  :: IntCoord(3) ! ijk integer coordinates
INTEGER(KIND=8),INTENT(IN)  :: box_di      ! box size in integer counts
INTEGER(KIND=8),INTENT(IN)  :: s_offset    ! =box_di**3-1, range of sfc indizes in one box
INTEGER,INTENT(IN)          :: box_nBits   ! number of bits of INTEGER(KIND=8)
INTEGER(KIND=8),INTENT(IN)  :: maxIJK      ! maximum domain size in integer counts (=2**box_nBits-1)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(OUT) :: s_minmax(27,2)  ! sfc index ranges of the 27 boxes
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER          :: i,j,k,l  ! ?
INTEGER(KIND=8)  :: IJK_min(3),IJK_tmp(3,27)  ! ?
!===================================================================================================================================
!find the octant in which the node lies (by integer division)
IJK_min=box_di*(IntCoord(:)/box_di)
! find sfc index for all  27 octants
s_minmax=-1 !marker for nothing found
l=0
DO k=0,2; DO j=0,2; DO i=0,2
  l=l+1
  !min
  IJK_tmp(1,l)=IJK_min(1)+(i-1)*box_di
  IJK_tmp(2,l)=IJK_min(2)+(j-1)*box_di
  IJK_tmp(3,l)=IJK_min(3)+(k-1)*box_di
END DO; END DO; END DO !i,j,k
CALL EVAL_MORTON_ARR(s_minmax(:,1),IJK_tmp,27,box_nBits)
DO l=1,27
  IF((MINVAL(IJK_tmp(:,l)).LT.0).OR.(MAXVAL(IJK_tmp(:,l)).GT.maxIJK)) s_minmax(l,1)=-1
END DO

CALL Qsort1Int2double(s_minmax(1:27,1:1))
s_minmax(:,2)=s_minmax(:,1)+s_offset
! merge consecutive ranges
DO i=27,2,-1
  IF(s_minmax(i-1,1).EQ.-1) EXIT
  IF(s_minmax(i,1).EQ.s_minmax(i-1,2)+1)THEN
    s_minmax(i-1,2)=s_minmax(i,2)
    s_minmax(i,:)=-1
  END IF
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

FUNCTION INVMAP(ID,nIDs,ArrID,maxSteps)
!===================================================================================================================================
! find the inverse Mapping of sfc index in a sorted list (a sorted array of unique NodeIDs), using bisection
! if Index is not in the range, -1 will be returned, gives back the first entry in the sorted list  which is >= ID 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8), INTENT(IN) :: ID            ! ID to search for
INTEGER, INTENT(IN)         :: nIDs          ! size of ArrID
INTEGER(KIND=8), INTENT(IN) :: ArrID(nIDs)   ! 1D array of IDs, SORTED!!
INTEGER, INTENT(IN)         :: maxSteps      ! =INT(LOG(REAL(nIDs))/LOG(2.))+1, better to compute once and pass it 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                     :: INVMAP         ! position in arrID, where arrID(INVMAP)  >= ID 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: i,low,up,mid  ! ?
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
DO i=1,maxSteps
  mid=(up-low)/2+low
  IF(ID .EQ. ArrID(mid))THEN
    INVMAP=mid
    EXIT
  ELSEIF(ID .GT. ArrID(mid))THEN ! seek in upper half
    low=mid
  ELSE
    up=mid
  END IF
END DO
INVMAP=low
END FUNCTION INVMAP 

END MODULE MOD_GlobalUniqueNodes
