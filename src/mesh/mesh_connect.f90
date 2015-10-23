#include "defines.f90"
MODULE MOD_Mesh_Connect
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
INTERFACE Connect
  MODULE PROCEDURE Connect 
END INTERFACE

INTERFACE Connect2DMesh
  MODULE PROCEDURE Connect2DMesh
END INTERFACE

PUBLIC::Connect
PUBLIC::Connect2DMesh
!===================================================================================================================================

INTEGER :: next1(4,3:4)=RESHAPE((/2,3,1,0,2,3,4,1/),SHAPE(next1))
INTEGER :: prev1(4,3:4)=RESHAPE((/3,1,2,0,4,1,2,3/),SHAPE(prev1))
INTEGER :: next2(4,3:4)=RESHAPE((/3,1,2,0,3,4,1,2/),SHAPE(next2))

CONTAINS

SUBROUTINE Connect()
!===================================================================================================================================
! Eliminates multiple nodes, checks periodic boundary conditions and connects elements to their neighbours.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,        ONLY:tElem,tSide,FirstElem
USE MOD_Mesh_Vars,        ONLY:nVV,VV
USE MOD_Mesh_Vars,        ONLY:nNonConformingSides,nConformingSides,nInnerSides
USE MOD_Mesh_Vars,        ONLY:deleteNode,deleteBC
USE MOD_Mesh_Vars,        ONLY:getNewNode,getNewSide
USE MOD_Mesh_Basis,       ONLY:createSides,buildEdges
USE MOD_GlobalUniqueNodes,ONLY:GlobalUniqueNodes
USE MOD_Mesh_Tools,       ONLY:BCVisu
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER       :: Elem                                                       ! Local element pointers
TYPE(tSide),POINTER       :: Side                                                       ! Local side pointers
TYPE(tSide),POINTER       :: pSide  ! ?
INTEGER                   :: iNode,i  ! ?
INTEGER                   :: nInner(2),nPeriodic(2)  ! ?
INTEGER                   :: nBCSides,nTotalSides,nPeriodicSides,connectedSides  ! ?
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,*)'Mesh connect starts'

! Build element sides (side-node mapping) according to cgns standard
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  CALL createSides(Elem,.FALSE.)
  Elem=>Elem%nextElem
END DO

! Delete Null BCs and check periodic BCs
WRITE(UNIT_stdOut,'(A)')'Count sides, delete Null BCs and check periodic BCs...'
nBCSides=0
nTotalSides=0
nPeriodicSides=0
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    nTotalSides=nTotalSides+1
    IF (ASSOCIATED(Side%BC)) THEN
      nBCSides=nBCSides+1
      IF (Side%BC%BCType .EQ. 0) THEN
        CALL deleteBC(Side%BC)
        nBCSides=nBCSides-1
      ELSEIF (Side%BC%BCType .EQ. 100) THEN
        nBCSides=nBCSides-1
      ELSEIF (Side%BC%BCType .EQ. 1) THEN
        nPeriodicSides=nPeriodicSides+1
        IF (Side%BC%BCalphaInd.EQ.0) &
          CALL abort(__STAMP__, &
            'No displacement vector vv assigned to periodic BC')
        IF (abs(Side%BC%BCalphaInd).GT.nVV)   &
          CALL abort(__STAMP__, &
            'No defined displacement vector vv assigned to periodic BC')
      END IF
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

WRITE(UNIT_stdOut,'(A)')   '-----------------------------------'
WRITE(UNIT_stdOut,'(A,I12)')'number of sides          : ',nTotalSides
WRITE(UNIT_stdOut,'(A,I12)')'number of Inner sides    : ',nTotalSides-nBCSides
WRITE(UNIT_stdOut,'(A,I12)')'number of BC sides       : ',nBCSides
WRITE(UNIT_stdOut,'(A,I12)')'number of periodic sides : ',nPeriodicSides
WRITE(UNIT_stdOut,'(A)')   '-----------------------------------'

WRITE(UNIT_stdOut,'(A)')'Insert periodic sides...'
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    Side%tmp2=0
    IF (ASSOCIATED(Side%BC)) THEN
      IF (Side%BC%BCType .EQ. 1) THEN
        !build a dummy side as connection
        Side%tmp2=Side%BC%BCalphaInd  ! vv +/-1 as info and marker for dummy side
        IF(Side%tmp2.GT.0)THEN ! side with positive vv is the master and gets the dummy side
          CALL getNewSide(pSide,Side%nNodes)  
          DO iNode=1,Side%nNodes
            CALL getNewNode(pSide%Node(iNode)%np)
            pSide%Node(iNode)%np%x  =Side%Node(iNode)%np%x + VV(:,ABS(Side%tmp2)) 
            pSide%Node(iNode)%np%ind=Side%Node(iNode)%np%ind
          END DO
          pSide%elem=>side%elem
          side%connection=>pSide 
        END IF
      END IF
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

! Eliminate multiple nodes
WRITE(UNIT_stdOut,'(A)')'Eliminate multiple Nodes...'
CALL GlobalUniqueNodes()

!CALL buildEdges(oriented=.FALSE.)

WRITE(UNIT_stdOut,'(A)')'Connect Conforming inner and periodic sides...'
CALL ConnectMesh()
connectedSides=nConformingSides
WRITE(UNIT_stdOut,'(A)')'Connect Non-Conforming inner and periodic sides...'
DO i=1,10
  IF(connectedSides.NE.nInnerSides)THEN
    IF(i.GT.1) WRITE(UNIT_stdOut,'(A)')'Not all Non-Conforming could be connected in first run, try again...'
  ELSE
    EXIT
  END IF
  CALL NonconformConnectMesh()
  IF(nNonConformingSides.EQ.0) EXIT
  connectedSides=connectedSides+nNonConformingSides
  WRITE(UNIT_StdOut,*)'   --> ',connectedSides,' sides of ', nInnerSides,'  sides connected.'
END DO



! 4. Check connectivity
nInner=0     ! inner sides
nPeriodic=0  ! periodic sides
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    ! Reset values used for connect 3d
    Side%tmp=0
    Side%tmp2=0
    IF(Side%LocSide .LE. 0) &
      CALL abort(__STAMP__, &
        'Mesh Connect: Side%LocSide not set!')
    IF (ASSOCIATED(Side%BC)) THEN
      IF(Side%BC%BCType .EQ. 1) THEN
        nPeriodic(1)=nPeriodic(1)+1
        IF (.NOT. ASSOCIATED(Side%Connection))THEN
          IF(.NOT.ASSOCIATED(Side%MortarSide))THEN
            nPeriodic(2)=nPeriodic(2)+1
            Side%CurveIndex=SIGN(Side%BC%BCAlphaInd,-1) ! set visu marker
            DO iNode=1,Side%nNodes
              ERRWRITE(*,*)Side%Node(iNode)%np%x
            END DO
          END IF
        END IF
      END IF
      IF(Side%BC%BCType .EQ. 100) THEN
        STOP
        nInner(1)=nInner(1)+1
        IF(.NOT. ASSOCIATED(Side%Connection)) THEN
          IF(.NOT.ASSOCIATED(Side%MortarSide))THEN
            nInner(2)=nInner(2)+1
            DO iNode=1,Side%nNodes
              ERRWRITE(*,*)Side%Node(iNode)%np%x
            END DO
          END IF
        END IF
      END IF
    ELSE
      nInner(1)=nInner(1)+1
      IF(.NOT. ASSOCIATED(Side%Connection)) THEN
        IF(.NOT.ASSOCIATED(Side%MortarSide))THEN
          nInner(2)=nInner(2)+1
          DO iNode=1,Side%nNodes
            ERRWRITE(*,*)Side%Node(iNode)%np%x
          END DO
        END IF
      END IF
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO
IF(nInner(2)+nPeriodic(2) .GT. 0) THEN
  WRITE(*,*) 'Sides with missing connection found!'
  WRITE(*,*) 'Inner sides: ',nInner(2),' of ',nInner(1),' sides missing.'
  WRITE(*,*) 'Periodic sides ',nPeriodic(2),' of ',nPeriodic(1),' sides missing.'
  IF(nPeriodic(2).GT.0) CALL BCVisu()
  CALL abort(__STAMP__, &
    'Sides with missing connection found.')
END IF

WRITE(UNIT_stdOut,'(A,F0.3,A)')'Mesh Connect completed with success.  '
CALL Timer(.FALSE.)
END SUBROUTINE Connect
SUBROUTINE ConnectMesh()
!===================================================================================================================================
! Connect all sides which can be found by node association. Uses Quicksort 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:tElem,tSide,tSidePtr,FirstElem,nInnerSides,nConformingSides
USE MOD_Mesh_Vars, ONLY:deleteNode
USE MOD_SortingTools,ONLY:Qsort1Int,Qsort4Int,MSortNInt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER       :: Elem                                                       ! Local element pointers
TYPE(tSide),POINTER       :: Side,nSide                                               ! Local side pointers
TYPE(tSidePtr),POINTER    :: InnerSides(:)   ! ?
INTEGER,ALLOCATABLE       :: SideConnect(:,:)  ! ?
INTEGER                   :: iNode,iSide,fNode  ! ?
INTEGER                   :: counter   ! ?
INTEGER                   :: NodeID  ! ?
LOGICAL                   :: dominant  ! ?
LOGICAL                   :: isInner  ! ?
!===================================================================================================================================
!count inner Sides and set side IDs (side%tmp)
nInnerSides=0
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    isInner=.FALSE.
    IF(.NOT.ASSOCIATED(Side%BC))THEN
      isInner=.TRUE.
    ELSE
      IF(Side%BC%BCType.EQ.100) isInner=.TRUE.
    END IF
    IF(side%tmp2.LT.0) isInner=.TRUE.

    IF(isInner)THEN !no BC side or periodic slave side
      nInnerSides=nInnerSides+1
      Side%tmp=nInnerSides  !"SideID"
      DO iNode=1,Side%nNodes
        Side%Node(iNode)%np%tmp=0 
      END DO
    ELSE
      !adjust oriented nodes for BCs
      DO iNode=1,Side%nNodes
        Side%orientedNode(iNode)%np=>Side%Node(iNode)%np
      END DO
    END IF
    !periodic master side with dummy
    IF(Side%tmp2.GT.0)THEN
      nInnerSides=nInnerSides+1
      Side%tmp=nInnerSides  !"SideID"
      DO iNode=1,Side%nNodes
        Side%connection%Node(iNode)%np%tmp=0 !use dummy nodes! 
      END DO
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

ALLOCATE(InnerSides(nInnerSides))
DO iSide=1,nInnerSides
  NULLIFY(InnerSides(iSide)%sp)
END DO
ALLOCATE(SideConnect(5,nInnerSides)) 
SideConnect=0
! make a list of innersides and set unique node IDs (%np%tmp)
! and make a list of side node ids (SideConnect(iSide,:)1:4 --> nodeids, 5: SideID
NodeID=0
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    isInner=.FALSE.
    IF(.NOT.ASSOCIATED(Side%BC))THEN
      isInner=.TRUE.
    ELSE
      IF(Side%BC%BCType.EQ.100) isInner=.TRUE.
    END IF
    IF(side%tmp2.LT.0) isInner=.TRUE.

    IF(isInner)THEN !no BC side or periodic slave side
      iSide=Side%tmp
      !put into side list
      InnerSides(iSide)%sp=>Side
      !Set Unique node IDs
      DO iNode=1,Side%nNodes
        IF(Side%Node(iNode)%np%tmp.EQ.0)THEN
          NodeID=NodeID+1
          Side%Node(iNode)%np%tmp=NodeID
        END IF
      END DO
      ! save unique node IDs of iSide in SideConnect(iSide,1:4) and sort these
      DO iNode=1,Side%nNodes
        SideConnect(iNode,iSide)=Side%Node(iNode)%np%tmp
      END DO
      SideConnect(5,iSide)=iSide
      CALL Qsort1Int(SideConnect(1:4,iSide)) !Node IDs sorted, unique combination for one side (for trias and quads)
    END IF
    !periodic master side with dummy
    IF(Side%tmp2.GT.0)THEN
      iSide=Side%tmp
      !put into side list
      InnerSides(iSide)%sp=>Side
      !Set Unique node IDs of periodic side dummy !! 
      DO iNode=1,Side%nNodes
        IF(Side%connection%Node(iNode)%np%tmp.EQ.0)THEN
          NodeID=NodeID+1
          Side%connection%Node(iNode)%np%tmp=NodeID
        END IF
      END DO
      ! save unique node IDs of iSide in SideConnect(iSide,1:4) and sort these
      DO iNode=1,Side%nNodes
        SideConnect(iNode,iSide)=Side%Connection%Node(iNode)%np%tmp
      END DO
      SideConnect(5,iSide)=iSide
      CALL Qsort1Int(SideConnect(1:4,iSide)) !Node IDs sorted, unique combination for one side (for trias and quads)
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

CALL MSortNInt(SideConnect,4,1)

counter=0
DO iSide=1,nInnerSides-1
  Side=>InnerSides(SideConnect(5,iSide))%sp
  IF(ASSOCIATED(Side%Connection).AND.(Side%tmp2.EQ.0)) CYCLE ! already connected side
  ! now compare Side Node IDs
  IF(SUM(ABS(SideConnect(1:4,iSide+1)-SideConnect(1:4,iSide))).EQ.0)THEN
    ! Side connection found 
    counter=counter+2
    nSide=>InnerSides(SideConnect(5,iSide+1))%sp

    IF(Side%tmp2.EQ.0)THEN !normal inner side
      fNode=0
      DO iNode=1,nSide%nNodes
         IF(ASSOCIATED(Side%Node(1)%np,nSide%Node(iNode)%np)) fNode=iNode
      END DO
      IF(fNode.EQ.0)CALL abort(__STAMP__, & 
                    'Problem with OrientedNode !')
    ELSEIF(Side%tmp2.LT.0) THEN !only connect from periodic slave side to master side (->nSide has a dummy connection side)
      ! for adjustorientednodes by pointer association (no tolerance gedoens!) 
      fNode=0
      DO iNode=1,Side%nNodes
         IF(ASSOCIATED(Side%Node(1)%np,nSide%connection%Node(iNode)%np)) fNode=iNode
      END DO
      IF(fNode.EQ.0)CALL abort(__STAMP__, & 
                    'Problem with OrientedNode on periodic side!')
      !now delete dummy Side
      DO iNode=1,nSide%nNodes
        NULLIFY(nSide%connection%Node(iNode)%np)
      END DO
      NULLIFY(nSide%connection)
      Side%tmp2=0
      nSide%tmp2=0
    ELSE !Side is master periodic, nSide is slave
      fNode=0
      DO iNode=1,nSide%nNodes
         IF(ASSOCIATED(nSide%Node(1)%np,Side%connection%Node(iNode)%np)) fNode=iNode
      END DO
      IF(fNode.EQ.0)CALL abort(__STAMP__, & 
                    'Problem with OrientedNode on periodic side (master)!')
      !now delete dummy Side
      DO iNode=1,Side%nNodes
        NULLIFY(Side%connection%Node(iNode)%np)
      END DO
      NULLIFY(Side%connection)
      Side%tmp2=0
      nSide%tmp2=0
    END IF !check tmp2
    !set connection
    nSide%connection=>Side
    Side%connection=>nSide
    !adjustorientednodes 
    dominant=(Side%Elem%ind .LT. nSide%Elem%ind)
    IF(dominant) THEN
      DO iNode=1,nSide%nNodes
        Side%orientedNode(iNode)%np=>Side%Node(iNode)%np
        nSide%orientedNode(iNode)%np=>nSide%Node(fNode)%np
        fNode=fNode-1
        IF(fNode .LT.1) fNode=fNode+nSide%nNodes
      END DO
    ELSE
      DO iNode=1,nSide%nNodes
        Side%orientedNode(iNode)%np=>Side%Node(fNode)%np
        nSide%orientedNode(iNode)%np=>nSide%Node(iNode)%np
        fNode=fNode-1
        IF(fNode .LT.1) fNode=fNode+nSide%nNodes
      END DO
    END IF
  END IF
END DO !iSide

WRITE(UNIT_StdOut,*)'   --> ',counter,' conforming sides of ', nInnerSides,'  sides connected.'
nConformingSides=counter
END SUBROUTINE ConnectMesh



SUBROUTINE NonconformConnectMesh()
!===================================================================================================================================
! Connect all non-conforming sides. 
! This routine assumes that all conforming sides have already been connected and that all nodes in the mesh are unique.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:tElem,tSide,tSidePtr,tNode,FirstElem
USE MOD_Mesh_Vars, ONLY:nNonconformingSides,nInnerSides
USE MOD_Mesh_Vars, ONLY:deleteNode
USE MOD_SortingTools,ONLY:Qsort1Int,Qsort4Int,MSortNInt
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER       :: Elem                                                       ! Local element pointers
TYPE(tSide),POINTER       :: Side,dummySide
TYPE(tSide),POINTER       :: aSide,bSide,cSide,smallSide1,smallSide2,bigSide,tmpSide
TYPE(tSidePtr),POINTER    :: Sides(:)   ! ?
TYPE(tSidePtr)            :: quartett(4)   ! ?
TYPE(tNode),POINTER       :: node
INTEGER                   :: iNode,jNode,iSide,jSide,kSide,ind,nQuartett  ! ?
INTEGER                   :: aLocSide,bLocSide,counter  ! ?
INTEGER                   :: CGNSToCart(4)
INTEGER                   :: bigCorner(4),bigCornerNoSort(4)
INTEGER                   :: masterNode,slaveNode
INTEGER,ALLOCATABLE       :: edgeNodes(:,:),sideInds(:)
LOGICAL,ALLOCATABLE       :: SideDone(:),checkSide(:) ! ?
LOGICAL,ALLOCATABLE       :: foundEdge(:,:)
INTEGER                   :: iCheck,jCheck,nCheck
LOGICAL                   :: commonNode,check
LOGICAL                   :: aFoundEdge(4),bFoundEdge(4)
!===================================================================================================================================
!count inner Sides and set side IDs (side%tmp)

! here we assume that all conforming sides are already connected and all remaining sides are nonconforming
! first count and collect all unconnected sides
CGNSToCart=(/1,2,4,3/)

! count number of unconnected (potetially non-conforming) sides
nNonConformingSides=0
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%MortarType.NE.0)THEN
      Side=>Side%nextElemSide
      CYCLE
    END IF
    IF(.NOT.ASSOCIATED(Side%BC))THEN
      IF(.NOT.ASSOCIATED(Side%Connection))THEN
        nNonConformingSides=nNonConformingSides+1
        Side%tmp=nNonConformingSides
      END IF
    ELSE
      IF(.NOT.ASSOCIATED(Side%Connection).AND.((Side%BC%BCType.EQ.100).OR.(Side%BC%BCType.EQ.1)))THEN
        nNonConformingSides=nNonConformingSides+1
        Side%tmp=nNonConformingSides
      END IF
      IF(ASSOCIATED(Side%Connection).AND.Side%tmp2.NE.0)THEN
        nNonConformingSides=nNonConformingSides+1
        Side%Connection%tmp=nNonConformingSides
      END IF
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

! collect all non-conforming sides
ALLOCATE(Sides(nNonConformingSides))
ALLOCATE(SideDone(nNonConformingSides))
ALLOCATE(checkSide(nNonConformingSides))
ALLOCATE(edgeNodes(4,nNonConformingSides))
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%MortarType.NE.0)THEN
      Side=>Side%nextElemSide
      CYCLE
    END IF
    IF(.NOT.ASSOCIATED(Side%BC))THEN
      IF(.NOT.ASSOCIATED(Side%Connection))THEN
        Sides(Side%tmp)%sp=>Side
      END IF
    ELSE
      IF(.NOT.ASSOCIATED(Side%Connection).AND.((Side%BC%BCType.EQ.100).OR.(Side%BC%BCType.EQ.1)))THEN
        Sides(Side%tmp)%sp=>Side
      END IF
      IF(ASSOCIATED(Side%Connection).AND.Side%tmp2.NE.0)THEN
        Sides(Side%Connection%tmp)%sp=>Side%Connection
      END IF
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

DO iSide=1,nNonConformingSides
  aSide=>Sides(iSide)%sp
  DO iNode=1,aSide%nNodes
    edgeNodes(iNode,iSide)=aSide%Node(iNode)%np%ind
  END DO
END DO


! Nullify tmp and count number of adjacent sides for each nodes
DO iSide=1,nNonConformingSides
  aSide=>Sides(iSide)%sp
  DO iNode=1,aSide%nNodes
    aSide%Node(iNode)%np%tmp=0
  END DO
END DO
DO iSide=1,nNonConformingSides
  aSide=>Sides(iSide)%sp
  DO iNode=1,aSide%nNodes
    aSide%Node(iNode)%np%tmp=aSide%Node(iNode)%np%tmp+1
  END DO
END DO
DO iSide=1,nNonConformingSides
  aSide=>Sides(iSide)%sp
  NULLIFY(aSide%MortarSide)
  DO iNode=1,aSide%nNodes
    node=>aSide%Node(iNode)%np
    IF(.NOT.ASSOCIATED(node%firstNormal))THEN
      ALLOCATE(node%firstNormal)
      ALLOCATE(node%firstNormal%FaceID(node%tmp))
      node%tmp=0
    END IF
    node%tmp=node%tmp+1
    node%firstNormal%FaceID(node%tmp)=iSide
  END DO
END DO

! now connect 2->1: we need exactly 3 sides (two small one big)
! connected by 3 edges where no 2 edges may have common nodes
! |-----------------------| Edge 1
! | |-------------------| |
! | |                   | |
! | |                   | |
! | |-------------------| | Edge 2
! | |-------------------| |
! | |                   | |
! | |                   | |
! | |-------------------| | Edge 3
! |-----------------------|


SideDone=.FALSE.
counter=0
ALLOCATE(sideInds(nNonConformingSides))
ALLOCATE(foundEdge(4,nNonConformingSides))
sideInds=0
print*, nNonConformingSides
DO iSide=1,nNonConformingSides-2
  IF(SideDone(iSide)) CYCLE
  aSide=>Sides(iSide)%sp

  ! find all sides which share a common edge with aSide
  nCheck=0
  DO iNode=1,aSide%nNodes
    node=>aSide%node(iNode)%np
    DO kSide=1,node%tmp
      jSide=node%firstNormal%FaceID(kSide)
      IF(iSide.EQ.jSide)  CYCLE
      IF(ANY(sideInds(1:nCheck).EQ.jSide)) CYCLE
      IF(SideDone(jSide)) CYCLE
      CALL CommonEdge2(edgeNodes(:,iSide),edgeNodes(:,jSide),aSide%nNodes,Sides(jSide)%sp%nNodes,aFoundEdge)
      ! condition for 2->1, exactly one edge of two sides is identical
      IF(COUNT(aFoundEdge).EQ.1)THEN
        nCheck=nCheck+1
        foundEdge(:,nCheck)=aFoundEdge
        sideInds(nCheck)=jSide
      END IF
    END DO
  END DO
  IF(nCheck.LT.2) CYCLE

  DO iCheck=1,nCheck-1
    DO jCheck=iCheck+1,nCheck
      IF(SideDone(iSide)) CYCLE
      commonNode=.TRUE.
      ! check if the two edges found for side A are opposite edges (dont share a common node)
      DO iNode=1,aSide%nNodes
        IF(foundEdge(iNode,iCheck).AND.foundEdge(next2(iNode,aSide%nNodes),jCheck)) commonNode=.FALSE.
      END DO
      IF(.NOT.commonNode)THEN
        bSide=>Sides(sideInds(iCheck))%sp
        cSide=>Sides(sideInds(jCheck))%sp
        CALL CommonEdge2(edgeNodes(:,sideInds(iCheck)),edgeNodes(:,sideInds(jCheck)),bSide%nNodes,cSide%nNodes,bFoundEdge)
        IF(COUNT(bFoundEdge).NE.1) CYCLE ! if we have a 2->1 mortar interface a,b and c share one unique edge
        ! now we know that we have 3 common edges, we can thus identify small/big sides by checking if
        ! two elements of the connected sides share a common side (-> small sides)
        ! This is a 3D (!) check, from this point on there has to be a solution, if not the mesh is broken
        check=.FALSE.
        CALL CommonElementSide(aSide%Elem,bSide%Elem,aLocSide,bLocSide)
        IF(aLocSide.GT.0) THEN
          check=.TRUE.
          smallSide1=>aSide; smallSide2=>bSide; bigSide=>cSide
        END IF
        CALL CommonElementSide(aSide%Elem,cSide%Elem,aLocSide,bLocSide)
        IF(aLocSide.GT.0) THEN
          check=.TRUE.
          smallSide1=>aSide; smallSide2=>cSide; bigSide=>bSide
        END IF
        CALL CommonElementSide(bSide%Elem,cSide%Elem,aLocSide,bLocSide)
        IF(aLocSide.GT.0) THEN
          check=.TRUE.
          smallSide1=>bSide; smallSide2=>cSide; bigSide=>aSide
        END IF
        IF(.NOT.check) CYCLE  ! fuck, small elements are not yet connected, maybe they have another morter side, this is a problem
        
        ! check which edges of big and small side are identical to determine the mortar type (wheter mortar is in
        !  xi or eta direction), then set the connection
        CALL CommonEdge(bigSide,smallSide1,aFoundEdge)
        IF(aFoundEdge(1).OR.aFoundEdge(3))THEN
          ! mortar type 2
          bigSide%MortarType=2
          IF(aFoundEdge(3))THEN
            tmpSide=>smallSide1
            smallSide1=>smallSide2
            smallSide2=>tmpSide
          END IF
        ELSEIF(aFoundEdge(2).OR.aFoundEdge(4))THEN
          ! mortar type 3
          bigSide%MortarType=3
          IF(aFoundEdge(2))THEN
            tmpSide=>smallSide1
            smallSide1=>smallSide2
            smallSide2=>tmpSide
          END IF
        ELSE
          STOP 'ERROR: Mortar type could not be identified!'
        END IF

        bigSide%nMortars=2
      IF(ASSOCIATED(bigSide%MortarSide))THEN
        print*,iSide
        print*,aSide%tmp
        print*,bSide%tmp
        print*,cSide%tmp
        print*,bigSide%tmp
        print*,smallSide1%tmp
        print*,smallSide2%tmp
        print*,bigSide%MortarSide(1)%sp%tmp
        print*,bigSide%MortarSide(2)%sp%tmp
        print*,smallSide1%tmp
        print*,smallSide2%tmp
        STOP 'already assoc 2 1'
      END IF
      IF(ALLOCATED(bigSide%MortarSide)) STOP 'already allocated 2 1'
        ALLOCATE(bigSide%MortarSide(2))
        bigSide%MortarSide(1)%sp=>smallSide1
        bigSide%MortarSide(2)%sp=>smallSide2
        smallSide1%connection=>bigSide
        smallSide2%connection=>bigSide
        smallSide1%MortarType=-bigSide%MortarType
        smallSide2%MortarType=-bigSide%MortarType
        
        SideDone(bigSide%tmp)=.TRUE.
        SideDone(smallSide1%tmp)=.TRUE.
        SideDone(smallSide2%tmp)=.TRUE.
        counter=counter+3
      END IF
    END DO
  END DO
END DO
  print*,counter

! Find 4->1 interfaces: 
! We assume that only sides belonging to a 4-> mortar interface are left in the sidelist
! At 4->1 interfaces the center node has exactly 4 neighbour sides, while all other nodes have
! more neighbour sides. We want to find the nodes with exactly 4 neighbour sides
! This criterion is necessary but not sufficient

DO iSide=1,nNonConformingSides
  IF(SideDone(iSide)) CYCLE
  aSide=>Sides(iSide)%sp
  DO iNode=1,aSide%nNodes
    aSide%Node(iNode)%np%tmp=0
  END DO
  CALL Qsort1Int(edgeNodes(:,iSide))
END DO
DO iSide=1,nNonConformingSides
  IF(SideDone(iSide)) CYCLE
  aSide=>Sides(iSide)%sp
  DO iNode=1,aSide%nNodes
    aSide%Node(iNode)%np%tmp=aSide%Node(iNode)%np%tmp+1
  END DO
END DO

! start searching
DO iSide=1,nNonConformingSides
  aSide=>Sides(iSide)%sp
  IF(SideDone(iSide)) CYCLE

  ! find quartett
  DO iNode=1,aSide%nNodes
    node=>aSide%Node(iNode)%np
    nQuartett=0
    IF(node%tmp.EQ.4)THEN
      DO kSide=1,node%tmp
        jSide=node%firstNormal%FaceID(kSide)
        IF(SideDone(jSide)) EXIT
        nQuartett=nQuartett+1
        bSide=>Sides(jSide)%sp
        quartett(nQuartett)%sp=>bSide

        ! check if diagonal node is corner node of big, then node is midnode of big mortar interface
        DO jNode=1,bSide%nNodes
          IF(node%ind.EQ.bSide%node(jNode)%np%ind)THEN
            bigCorner(nQuartett)=bSide%Node(next2(jNode,bSide%nNodes))%np%ind
            EXIT
          END IF
        END DO
      END DO
    END IF
    IF(nQuartett.EQ.4) EXIT
  END DO
  IF(nQuartett.NE.4) CYCLE

  !! find 4 sides with same node ind as specified above and their corner nodes
  !! if bSide is a small side then the node at the opposite corner to ind is a big sides corner node
  !DO jNode=1,aSide%nNodes
  !  node=>aSide%node(jNode)%np
  !  IF(node%tmp.NE.4) CYCLE
  !  DO kSide=1,node%tmp
  !    jSide=node%firstNormal%FaceID(kSide)
  !    IF(iSide.EQ.jSide)  CYCLE
  !    IF(SideDone(jSide)) CYCLE
  !    bSide=>Sides(jSide)%sp
  !    !IF(bSide%tmp.EQ.iSide) CYCLE
  !    DO iNode=1,bSide%nNodes
  !      IF(bSide%Node(iNode)%np%ind.EQ.ind)THEN
  !        nQuartett=nQuartett+1
  !        quartett(nQuartett)%sp=>bSide
  !        bigCorner(nQuartett)=bSide%Node(next2(iNode,bSide%nNodes))%np%ind
  !        bSide%tmp=iSide ! marker
  !        EXIT
  !      END IF
  !    END DO
  !  END DO
  !END DO
  !IF(nQuartett.NE.4) THEN
  !  CALL Abort(__STAMP__,&
  !    'ERROR: Expected four adjacent sides to search node, but found',nQuartett)
  !END IF

  bigCornerNoSort=bigCorner
  CALL Qsort1Int(bigCorner)

  ! find big side belonging to quartett
  DO jSide=1,nNonConformingSides
    IF(SideDone(jSide)) CYCLE
    bSide=>Sides(jSide)%sp
    ! check potential big corner nodes of quartett against all nodes of bSide
    ! if nodes are identical then bSide is big side and quartett are small sides
    IF(ALL(bigCorner.EQ.edgeNodes(:,jSide)))THEN
      bSide%MortarType=1
      bSide%nMortars=4
      SideDone(jSide)=.TRUE.
      IF(ALLOCATED(bSide%MortarSide)) STOP 'already allocated 4 1'
      ALLOCATE(bSide%MortarSide(4))
      DO iNode=1,bSide%nNodes
        bSide%Node(iNode)%np%tmp=0
        DO jNode=1,bSide%nNodes
          quartett(iNode)%sp%Node(jNode)%np%tmp=0
        END DO
      END DO
      DO iNode=1,bSide%nNodes
        ! for type 1, small mortars are sorted on a cartesian grid (first xi, then eta)
        ! this means that e.g. the small side at node 3 of big side is stored in position 4 of mortar array
        jNode=-999
        DO jNode=1,bSide%nNodes
          IF(bigCornerNoSort(jNode).EQ.bSide%Node(iNode)%np%ind) EXIT
        END DO
        bSide%MortarSide(CGNSToCart(iNode))%sp=>quartett(jNode)%sp
        quartett(jNode)%sp%connection=>bSide
        quartett(jNode)%sp%MortarType=-1
        SideDone(quartett(jNode)%sp%tmp)=.TRUE.
      END DO
      counter=counter+5
    ELSE
      ! bSide is not a big side, reset tmp
      !DO iNode=1,aSide%nNodes
      !  IF(aSide%Node(iNode)%np%ind.EQ.ind) aSide%Node(iNode)%np%tmp=0
      !END DO
    END IF
  END DO
END DO

! now set the oriented nodes for mortar sides
! big sides are always master and small sides inherit order
! OrientedNode(iNode) of big side is identical to OrientedNode (iNode) of small
! side, if small side contains that node
DO iSide=1,nNonConformingSides
  aSide=>Sides(iSide)%sp
  IF(aSide%nMortars.LE.0) CYCLE  ! only check big mortar sides
  ! big side
  DO iNode=1,aSide%nNodes
    aSide%OrientedNode(iNode)%np=>aSide%Node(iNode)%np
  END DO
  ! small sides
  DO jSide=1,aSide%nMortars
    bSide=>aSide%MortarSide(jSide)%sp 
    commonNode=.FALSE.
    DO iNode=1,aSide%nNodes
      DO jNode=1,bSide%nNodes
        IF(ASSOCIATED(aSide%Node(iNode)%np,bSide%Node(jNode)%np))THEN
          masterNode=iNode
          slaveNode=jNode
          commonNode=.TRUE.
          EXIT
        END IF
      END DO
      IF(commonNode) EXIT
    END DO
    IF(.NOT.commonNode) STOP 'ERROR: no common node of big and small mortar sides found'
    DO iNode=1,aSide%nNodes
      bSide%orientedNode(masterNode)%np=>bSide%Node(slaveNode)%np
      masterNode=prev1(masterNode,aSide%nNodes)
      slaveNode=next1(slaveNode,aSide%nNodes)
    END DO
  END DO
END DO

! for periodic side connections, dummy sides have been introduced for periodic masters
! only these sides are connected by mortars
! remove temporary periodic dummy sides and update mortar connection to real (inner) sides
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  aSide=>Elem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    check=.TRUE.
    IF(aSide%tmp2.LE.0)                    check=.FALSE.
    IF(.NOT.ASSOCIATED(aSide%connection))THEN
                                           check=.FALSE.
    ELSE
      IF(aSide%connection%MortarType.EQ.0) check=.FALSE.
    END IF

    ! only check mortar sides which are periodic masters
    IF(.NOT.check)THEN
      aSide=>aSide%nextElemSide
      CYCLE ! normal side
    END IF

    dummySide=>aSide%connection
    IF(dummySide%nMortars.GT.0)THEN
      ! dummy side is mortar master, build mortar connection in aSide
      ALLOCATE(aSide%MortarSide(dummySide%nMortars))
      aSide%MortarType=dummySide%MortarType
      aSide%nMortars  =dummySide%nMortars
      DO iSide=1,aSide%nMortars
        aSide%MortarSide(iSide)%sp=>dummySide%MortarSide(iSide)%sp
        aSide%MortarSide(iSide)%sp%connection=>aSide
      END DO
      DO iNode=1,aSide%nNodes
        aSide%OrientedNode(iNode)%np=>aSide%Node(iNode)%np
      END DO
      NULLIFY(aSide%connection)
      DEALLOCATE(dummySide%node)
      DEALLOCATE(dummySide%orientedNode)
      DEALLOCATE(dummySide%MortarSide)
      DEALLOCATE(dummySide)
      aSide%tmp2=0
    ELSE
      ! dummy side is mortar slave, update connection of mortar master
      bigSide=>dummySide%connection
      IF(ASSOCIATED(dummySide,bigSide)) STOP 'ERROR: Periodic mortar slave has no connection to master!'
      DO iSide=1,bigSide%nMortars
        IF(ASSOCIATED(bigSide%MortarSide(iSide)%sp,dummySide))THEN
          bigSide%MortarSide(iSide)%sp=>aSide
          EXIT
        END IF
      END DO
      aSide%connection=>bigSide
      aSide%MortarType=-bigSide%MortarType
      DO iNode=1,aSide%nNodes
        DO jNode=1,aSide%nNodes
          IF(ASSOCIATED(dummySide%OrientedNode(iNode)%np,dummySide%Node(jNode)%np))THEN
            aSide%OrientedNode(iNode)%np=>aSide%Node(jNode)%np
            EXIT
          END IF
        END DO
      END DO
      DEALLOCATE(dummySide%node)
      DEALLOCATE(dummySide%orientedNode)
      DEALLOCATE(dummySide)
      aSide%tmp2=0
    END IF
    aSide=>aSide%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

! deallocate temporary side markers again
DO iSide=1,nNonConformingSides
  aSide=>Sides(iSide)%sp
  DO iNode=1,aSide%nNodes
    node=>aSide%Node(iNode)%np
    IF(ASSOCIATED(node%firstNormal)) DEALLOCATE(node%firstNormal)
  END DO
END DO

IF(counter.NE.nNonConformingSides) THEN
  WRITE(*,*) 'Warning: Number of remaining nonconforming sides:', nNonConformingSides
  WRITE(*,*) '         Number of found nonconforming sides:', counter
END IF

DEALLOCATE(SideDone)
DEALLOCATE(checkSide)
DEALLOCATE(sideInds)
DEALLOCATE(edgeNodes)
DEALLOCATE(Sides)
DEALLOCATE(foundEdge)
nNonConformingSides=counter

END SUBROUTINE NonconformConnectMesh


PURE SUBROUTINE CommonEdge2(aSide,bSide,nA,nB,aFoundEdge)
!===================================================================================================================================
! Check if two sides share a common edge just by node inds
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                   :: aSide(4),bSide(4)
INTEGER,INTENT(IN)                   :: nA,nB
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)                  :: aFoundEdge(4)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
LOGICAL                              :: aFoundNode(4)
INTEGER                              :: iNode,jNode  ! ?
!===================================================================================================================================
aFoundNode=.FALSE.
aFoundEdge=.FALSE.

DO iNode=1,nA
  DO jNode=1,nB
    IF(aSide(iNode).EQ.bSide(jNode))THEN
      aFoundNode(iNode)=.TRUE.
      EXIT
    END IF
  END DO
END DO

DO iNode=1,nA
  aFoundEdge(iNode)=(aFoundNode(iNode).AND.aFoundNode(next1(iNode,nA)))
END DO
END SUBROUTINE CommonEdge2
SUBROUTINE CommonEdge(aSide,bSide,aFoundEdge)
!===================================================================================================================================
! Check if two sides share a common edge just by node inds
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:tSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide),POINTER,INTENT(IN)       :: aSide,bSide
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)                  :: aFoundEdge(4)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
LOGICAL                              :: aFoundNode(4)
INTEGER                              :: iNode,jNode  ! ?
!===================================================================================================================================
aFoundNode=.FALSE.

DO iNode=1,aSide%nNodes
  DO jNode=1,bSide%nNodes
    IF(aSide%Node(iNode)%np%ind.EQ.bSide%Node(jNode)%np%ind)THEN
      aFoundNode(iNode)=.TRUE.
      EXIT
    END IF
  END DO
END DO

DO iNode=1,aSide%nNodes
  aFoundEdge(iNode)=(aFoundNode(iNode).AND.aFoundNode(next1(iNode,aSide%nNodes)))
END DO
END SUBROUTINE CommonEdge



SUBROUTINE CommonElementSide(aElem,bElem,aLocSide,bLocSide)
!===================================================================================================================================
! Check if two elements share one common side, ignore periodic sides
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:tSide,tElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(IN)       :: aElem,bElem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)                  :: aLocSide,bLocSide
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tSide),POINTER                  :: aSide,bSide
!===================================================================================================================================
aLocSide=0
aSide=>aElem%firstSide
DO WHILE(ASSOCIATED(aSide))
  aLocSide=aLocSide+1 
  bLocSide=0
  bSide=>bElem%firstSide
  DO WHILE(ASSOCIATED(bSide))
    ! TODO: Maybe dont rely on connection but use corners node inds instead
    bLocSide=bLocSide+1
    IF(ASSOCIATED(aSide%connection,bSide)) THEN
      IF     (ASSOCIATED(aSide%BC)) THEN
        IF(aSide%BC%BCType .NE. 1) RETURN ! no periodics
      ELSE
        RETURN
      END IF
    END IF
    IF(ASSOCIATED(bSide%connection,aSide)) THEN  ! double check is necessary in case one of the sides is a mortar side
      IF     (ASSOCIATED(bSide%BC)) THEN
        IF(bSide%BC%BCType .NE. 1) RETURN ! no periodics
      ELSE
        RETURN
      END IF
    END IF
    bSide=>bSide%nextElemSide
  END DO
  aSide=>aSide%nextElemSide
END DO
aLocSide=-1
bLocSide=-1

END SUBROUTINE CommonElementSide



SUBROUTINE Connect2DMesh(firstElem_in)
!===================================================================================================================================
! Connect all side edges (2D) which can be found by node association. Uses Quicksort 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:tElem,tSide,tSidePtr
USE MOD_SortingTools, ONLY:Qsort2Int
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(IN)       :: FirstElem_in   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                :: i,iSide,nEdges  ! ?
TYPE(tElem),POINTER    :: aElem   ! ?
TYPE(tSide),POINTER    :: aSide   ! ?
TYPE(tSidePtr),POINTER :: Edges(:)   ! ?
INTEGER,ALLOCATABLE    :: EdgeConnect(:,:)   ! ?
!===================================================================================================================================
nEdges=0
aElem=>firstElem_in
DO WHILE(ASSOCIATED(aElem))
  nEdges=nEdges+aElem%nNodes
  aElem=>aElem%nextElem
END DO  

!build up connect, new style!!
ALLOCATE(Edges(nEdges),EdgeConnect(nEdges,3))
DO i=1,nEdges
  NULLIFY(Edges(i)%sp)
END DO
!prepare connect
iSide=0
aElem=>firstElem_in
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    iSide=iSide+1
    Edges(iSide)%sp=>aSide
    EdgeConnect(iSide,1)=MIN(aSide%Node(1)%np%ind,aSide%Node(2)%np%ind)
    EdgeConnect(iSide,2)=MAX(aSide%Node(1)%np%ind,aSide%Node(2)%np%ind)
    EdgeConnect(iSide,3)=iSide
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO  
CALL Qsort2Int(EdgeConnect)
! only local 2D connect is needed
DO iSide=2,nEdges
  IF((EdgeConnect(iSide,1).EQ.EdgeConnect(iSide-1,1)).AND.(EdgeConnect(iSide,2).EQ.EdgeConnect(iSide-1,2)))THEN
    Edges(EdgeConnect(iSide,3))%sp%connection   => Edges(EdgeConnect(iSide-1,3))%sp
    Edges(EdgeConnect(iSide-1,3))%sp%connection => Edges(EdgeConnect(iSide,3))%sp
  END IF
END DO

DO i=1,nEdges
  NULLIFY(Edges(i)%sp)
END DO
DEALLOCATE(Edges,EdgeConnect)
END SUBROUTINE Connect2DMesh

END MODULE MOD_Mesh_Connect
