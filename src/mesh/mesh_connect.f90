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
USE MOD_Mesh_Basis,       ONLY:adjustOrientedNodes,createSides,buildEdges
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
INTEGER                   :: iNode  ! ?
INTEGER                   :: nInner(2),nPeriodic(2)  ! ?
INTEGER                   :: nBCSides,nTotalSides,nPeriodicSides  ! ?
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
            pSide%Node(iNode)%np%x=Side%Node(iNode)%np%x + VV(:,ABS(Side%tmp2)) 
          END DO
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
CALL NonconformConnectMesh()

WRITE(UNIT_StdOut,*)'   --> ',nConformingSides+nNonconformingSides,' sides of ', nInnerSides,'  sides connected.'


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
          nPeriodic(2)=nPeriodic(2)+1
          Side%CurveIndex=SIGN(Side%BC%BCAlphaInd,-1) ! set visu marker
          DO iNode=1,Side%nNodes
            ERRWRITE(*,*)Side%Node(iNode)%np%x
          END DO
        END IF
      END IF
      IF(Side%BC%BCType .EQ. 100) THEN
        nInner(1)=nInner(1)+1
        IF(.NOT. ASSOCIATED(Side%Connection)) THEN
          nInner(2)=nInner(2)+1
          DO iNode=1,Side%nNodes
            ERRWRITE(*,*)Side%Node(iNode)%np%x
          END DO
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
!CALL Qsort4Int(SideConnect)

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
! Connect all sides which can be found by node association. Uses Quicksort 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:tElem,tSide,tSidePtr,FirstElem
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
TYPE(tSide),POINTER       :: Side
TYPE(tSide),POINTER       :: aSide,bSide,cSide,smallSide1,smallSide2,bigSide
TYPE(tSidePtr),POINTER    :: InnerSides(:)   ! ?
INTEGER                   :: iNode,iSide,jSide,kSide  ! ?
INTEGER                   :: aLocSide,bLocSide,counter  ! ?
INTEGER                   :: next2(4)
LOGICAL,ALLOCATABLE       :: SideDone(:) ! ?
LOGICAL                   :: commonNode
LOGICAL                   :: aFoundEdge(4,2),bFoundEdge(4,2),cFoundEdge(4,2)
LOGICAL                   :: aFoundNode(4,2),bFoundNode(4,2),cFoundNode(4,2)
!===================================================================================================================================
!count inner Sides and set side IDs (side%tmp)

! here we assume that all conforming sides are already connected and all remaining sides are nonconforming
! first count and collect all unconnected sides
next2=(/3,4,1,2/)

nNonConformingSides=0
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(.NOT.ASSOCIATED(Side%Connection).AND..NOT.ASSOCIATED(Side%BC))THEN
      nNonConformingSides=nNonConformingSides+1
      Side%tmp=nNonConformingSides
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

ALLOCATE(InnerSides(nNonConformingSides))
ALLOCATE(SideDone(nNonConformingSides))

Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(.NOT.ASSOCIATED(Side%Connection).AND..NOT.ASSOCIATED(Side%BC)) &
      InnerSides(Side%tmp)%sp=>Side
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
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
DO iSide=1,nNonConformingSides
  IF(SideDone(iSide)) CYCLE
  aSide=>InnerSides(iSide)%sp
  DO jSide=iSide+1,nNonConformingSides
    IF(SideDone(jSide)) CYCLE
    bSide=>InnerSides(jSide)%sp
    CALL CommonNodeAndEdge(aSide,bSide,aFoundNode(:,1),bFoundNode(:,1),aFoundEdge(:,1),bFoundEdge(:,1))
    ! condition for 2->1, exactly one edge of two sides is identical
    IF(COUNT(aFoundEdge(:,1)).NE.1) CYCLE

    DO kSide=jSide+1,nNonConformingSides
      cSide=>InnerSides(kSide)%sp
      CALL CommonNodeAndEdge(aSide,cSide,aFoundNode(:,2),cFoundNode(:,1),aFoundEdge(:,2),cFoundEdge(:,1))
      IF(COUNT(aFoundEdge(:,2)).NE.1) CYCLE
      commonNode=.TRUE.
      DO iNode=1,aSide%nNodes
        IF(aFoundEdge(iNode,1).AND.aFoundEdge(next2(iNode),2)) commonNode=.FALSE.
      END DO
      IF(.NOT.commonNode)THEN
        CALL CommonNodeAndEdge(bSide,cSide,bFoundNode(:,2),cFoundNode(:,2),bFoundEdge(:,2),cFoundEdge(:,2))
        IF(COUNT(bFoundEdge(:,2)).NE.1) CYCLE
        ! now we know that we have 3 common edges, we can now identify small/big sides by checking if
        ! two elements of the connected sides share a common side (-> small sides)
        CALL CommonElementSide(aSide%Elem,bSide%Elem,aLocSide,bLocSide)
        IF(aLocSide.GT.0) THEN
          smallSide1=>aSide; smallSide2=>bSide; bigSide=>cSide
        END IF
        CALL CommonElementSide(aSide%Elem,cSide%Elem,aLocSide,bLocSide)
        IF(aLocSide.GT.0) THEN
          smallSide1=>aSide; smallSide2=>cSide; bigSide=>bSide
        END IF
        CALL CommonElementSide(bSide%Elem,cSide%Elem,aLocSide,bLocSide)
        IF(aLocSide.GT.0) THEN
          smallSide1=>bSide; smallSide2=>cSide; bigSide=>aSide
        END IF
        print*,aSide%tmp,bSide%tmp,cSide%tmp

        bigSide%nMortars=2
        ALLOCATE(bigSide%MortarSide(2))
        bigSide%MortarSide(1)%sp=>smallSide1
        bigSide%MortarSide(2)%sp=>smallSide2
        smallSide1%Connection=>bigSide
        smallSide2%Connection=>bigSide
        SideDone(aSide%tmp)=.TRUE.
        SideDone(bSide%tmp)=.TRUE.
        SideDone(cSide%tmp)=.TRUE.
        counter=counter+3
      END IF
    END DO
  END DO
END DO

IF(counter.NE.nNonConformingSides) THEN
  WRITE(*,*) 'Warning: Number of expected nonconforming sides:', nNonConformingSides
  WRITE(*,*) '         Number of found nonconforming sides:', counter
END IF

WRITE(UNIT_StdOut,*)'   --> ',counter,' nonconforming sides of ', nInnerSides,'  sides connected.'


DEALLOCATE(SideDone)
DEALLOCATE(InnerSides)

END SUBROUTINE NonconformConnectMesh


SUBROUTINE CommonNodeAndEdge(aSide,bSide,aFoundNode,bFoundNode,aFoundEdge,bFoundEdge)
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
LOGICAL,INTENT(OUT)                  :: aFoundEdge(4),bFoundEdge(4)
LOGICAL,INTENT(OUT)                  :: aFoundNode(4),bFoundNode(4)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                              :: iNode,jNode  ! ?
INTEGER                              :: nextNode(3:4,4)
!===================================================================================================================================
aFoundNode=.FALSE.
bFoundNode=.FALSE.
DO iNode=1,aSide%nNodes
  DO jNode=1,bSide%nNodes
    IF(aSide%Node(iNode)%np%ind.EQ.bSide%Node(jNode)%np%ind)THEN
      aFoundNode(iNode)=.TRUE.
      bFoundNode(jNode)=.TRUE.
      EXIT
    END IF
  END DO
END DO

nextNode(3,:)=(/2,3,1,0/)
nextNode(4,:)=(/2,3,4,1/)
aFoundEdge=.FALSE.
bFoundEdge=.FALSE.
DO iNode=1,aSide%nNodes
  IF(aFoundNode(iNode).AND.aFoundNode(nextNode(aSide%nNodes,iNode))) &
    aFoundEdge(iNode)=.TRUE.
END DO
DO iNode=1,bSide%nNodes
  IF(bFoundNode(iNode).AND.bFoundNode(nextNode(bSide%nNodes,iNode))) &
    bFoundEdge(iNode)=.TRUE.
END DO

END SUBROUTINE CommonNodeAndEdge



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
      IF (ASSOCIATED(aSide%BC)) THEN
        IF(aSide%BC%BCType .NE. 1) RETURN ! no periodics
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
