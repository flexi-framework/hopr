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
MODULE MOD_Readin_ANSA
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
INTERFACE ReadStar
  MODULE PROCEDURE ReadStar
END INTERFACE

INTERFACE readStar_split
  MODULE PROCEDURE readStar_split
END INTERFACE

PUBLIC::ReadStar,readStar_split
!===================================================================================================================================

CONTAINS
SUBROUTINE ReadStar()
!===================================================================================================================================
! Read mesh from ansa ascii file. Called by fillMesh.
! Read-in can be performed by just one processor
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tNodePtr
USE MOD_Mesh_Vars,ONLY:firstElem
USE MOD_Mesh_Vars,ONLY:MeshFileName
USE MOD_Mesh_Vars,ONLY:nUserDefinedBoundaries,BoundaryName,BoundaryType
USE MOD_Mesh_Basis,ONLY:createSides
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewNode
USE MOD_SortingTools,ONLY:Qsort1Int,Qsort4int
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! nMeshFiles                : Number of mesh files (INI-File)
! FileName(1)               : FileName of mesh file only 1
! nZones                    : Number of mesh zones (INI-File)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                :: os !openStatus
INTEGER                :: iNode,nNodes ! counter/number of Nodes
INTEGER                :: iElem,nElems,nElNodes ! counter/number of Elements
INTEGER                :: iZone,nZones,iBC,nBCs,i  ! ?
INTEGER                :: dummy1,dummy2,dummy3  ! ?
INTEGER                :: nodeInd   ! ?
LOGICAL,ALLOCATABLE    :: NodeBC(:,:)  ! ?
LOGICAL                :: found  ! ?
INTEGER,ALLOCATABLE    :: ZoneID(:),BCID(:),MapBC(:)  ! ?
REAL,ALLOCATABLE       :: nodeCoords(:,:)  ! ?
TYPE(tNodePtr),POINTER :: Nodes(:)  ! ?

TYPE(tElem),POINTER    :: aElem   ! ?
TYPE(tSide),POINTER    :: aSide   ! ?
INTEGER                :: conn(8),conn4(4) !connectivity list of eight nodes
INTEGER,ALLOCATABLE    :: BCSideBuffer(:,:)  ! ?
INTEGER,ALLOCATABLE    :: iDummyArray2(:,:),iDummyArray3(:)  ! ?
CHARACTER(LEN=100)     :: cdummy,cdummy2   ! ?
CHARACTER(LEN=200)     :: cDummy3  ! ?
INTEGER                :: NodeMap(5,8)   !for each Elemtype a special map
INTEGER                :: iBC_low,iBC_mid,iBC_up,nBCSides  ! ?
INTEGER                :: iBCSide,searchdir,maxSteps  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')'Reading Star(ANSA) mesh...'

!mappings found by Elem%nNodes-3
NodeMap=0
NodeMap(1,:)=(/1,2,3,5,0,0,0,0/)!tetra
NodeMap(2,:)=(/1,2,3,4,5,0,0,0/)!pyra, no mapping
NodeMap(3,:)=(/1,2,3,5,6,7,0,0/)!penta
NodeMap(5,:)=(/1,2,3,4,5,6,7,8/)!hexa

! start reading mesh
CALL openstarfile(TRIM(MeshFileName(1))//'.inp',106)
nZones=0
nBCs=0
DO
  READ(106,'(A100)',iostat=os) cdummy 
  IF(os.NE.0) EXIT !end of file
  IF(INDEX(TRIM(cdummy),'!').EQ.1)CYCLE !skip comments
  IF(INDEX(TRIM(cdummy),'ZONE_').NE.0) THEN
    nZones=nZones+1
    CYCLE
  END IF
  IF(INDEX(TRIM(cdummy),'BC_').NE.0) nBCs=nBCs+1
END DO
!  WRITE(*,*)'DEBUG number of Zones and BCs',nZones,nBCs
ALLOCATE(ZoneID(nZones),BCID(nBCs),MapBC(nBCs))
REWIND(106)
iZone=0
iBC=0
ZoneID=0
BCID=0
MapBC=0
DO
  READ(106,'(A100)',iostat=os) cdummy 
  IF(os.NE.0) EXIT !end of file
  IF(INDEX(TRIM(cdummy),'!').EQ.1) CYCLE
  IF(INDEX(TRIM(cdummy),'ZONE_').NE.0) THEN
    iZone=iZone+1
    READ(cdummy,*)cDummy2,ZoneID(iZone),cDummy3  !read ZoneID
    WRITE(*,*)'Zone found ',TRIM(cDummy3)
    CYCLE
  END IF
  IF(INDEX(TRIM(cdummy),'BC_').NE.0) THEN
    READ(cdummy,*)cDummy2,dummy1,cDummy3 !read BC ID
    found=.false.
    DO i=1,nUserDefinedBoundaries
      IF(INDEX(TRIM(cdummy3),TRIM(BoundaryName(i))).NE.0) THEN
        found=.true. 
        iBC=iBC+1
        BCID(iBC)=dummy1
        MapBC(iBC)=i
        WRITE(*,*)'BC found: ',TRIM(cDummy3)
        WRITE(*,*)'              -->  mapped to:',TRIM(BoundaryName(i))
        EXIT
      END IF
    END DO
    IF(.NOT.found) CALL abort(__STAMP__, &
                      'UserDefinedBoundary condition missing: '//TRIM(cdummy3),nUserDefinedBoundaries,999.)
  END IF
END DO
CLOSE(106) 
!read number of 3D Elements
CALL openstarfile(TRIM(MeshFileName(1))//'.cel',105)
os=0
nElems=0
DO                    ! ElemID, 8 Nodes, PID, ???
  READ(105,*,iostat=os) dummy1, conn(:),dummy2,dummy3
  IF(os.NE.0) EXIT !end of file
  IF(conn(5).EQ.0) CYCLE !2D Element
  nElems=nElems+1
  DO iZone=1,nZones
    IF(ZoneID(iZone).LT.0)CYCLE
    IF(ZoneID(iZone).EQ.dummy2)THEN
      ZoneID(iZone)=-ZoneID(iZone) ! mark used Zones with a minus sign
    END IF
  END DO
END DO
iZone=1
DO WHILE (iZone.LE.nZones)
  IF(ZoneID(iZone).GT.0)THEN
    ZoneID(iZone:nZones-1)=ZoneID(iZone+1:nZones)
    ZoneID(nZones)=0
    nZones=nZones-1
  END IF
  iZone=iZone+1
END DO
ZoneID=ABS(ZoneID) 
CLOSE(105)

! Count nodes
CALL openstarfile(TRIM(MeshFileName(1))//'.vrt',104)
os    =0 
nNodes=0
DO   
  READ(104,*,iostat=os)
  IF(os.NE.0) EXIT !end of file
  nNodes=nNodes+1
END DO
WRITE(UNIT_StdOut,*)' nNodes: ',nNodes
CLOSE(104)
! Now read node coordinates
ALLOCATE(nodeCoords(3,nNodes))
nodeCoords=0.
CALL openstarfile(TRIM(MeshFileName(1))//'.vrt',104)
DO iNode = 1,nNodes   ! READ nodes
  READ(104,*)dummy1,nodeCoords(:,iNode)
END DO
CLOSE(104)

! read Boundary Conditions and assign to nodes
ALLOCATE(NodeBC(nNodes,nBCs))
NodeBC=.false.
nBCSides=0
CALL openstarfile(TRIM(MeshFileName(1))//'.bnd',107)
DO                            ! FaceID, Nodes, BCID, ?, ??
  READ(107,*,iostat=os) dummy1,conn(1:4),dummy2,dummy3,cdummy 
  IF(os.NE.0) EXIT !end of file
  nBCSides=nBCSides+1
  found=.false.
  DO iBC=1,nBCs
    IF(BCID(iBC).EQ.dummy2) THEN
      found=.true.
      EXIT
    END IF
  END DO
  IF(.NOT.found) CALL abort(__STAMP__, &
     'a boundary condition was used without corresponding UserDefinedBoundary',dummy2,999.)
  DO iNode=1,4
    NodeBC(conn(iNode),iBC)=.true.
  END DO
END DO
CLOSE(107)
!!  !search for unused BCs
!!  DO iBC=1,nBCs
!!    found=.false.
!!    iNode=1
!!    DO WHILE((.NOT.found).AND. (iNode.LE.nNodes))
!!      found=found.OR.NodeBC(iNode,iBC)
!!      iNode=iNode+1
!!    END DO
!!    IF(.NOT.found)THEN
!!      NodeBC(:,iBC:nBCs-1)=NodeBC(:,iBC+1:nBCs) !unused BC
!!      NodeBC(:,nBCs)=.false.
!!      BCID(iBC:nBCs-1)=BCID(iBC+1:nBCs) !unused BC
!!      BCID(nBCs)=0
!!      MapBC(iBC:nBCs-1)=MapBC(iBC+1:nBCs) !unused BC
!!      MapBC(nBCs)=0
!!      nBCs=nBCs-1
!!    END IF
!!  END DO
!!  WRITE(*,*)'DEBUG,BCID',BCID
!!  WRITE(*,*)'DEBUG,MapBC',MapBC
!!  WRITE(*,*)'DEBUG,NodeBC',NodeBC
! node pointers!
ALLOCATE(Nodes(nNodes))
DO iNode=1,nNodes
  NULLIFY(Nodes(iNode)%np)
END DO
aElem=>firstElem ! READ elems
IF(ASSOCIATED(aElem))THEN
  DO WHILE(ASSOCIATED(aElem%nextElem))
    aElem%Ind=-1
    aElem=>aElem%nextElem
  END DO
  aElem%Ind=-1
END IF
ALLOCATE(iDummyArray2(nElems,9),iDummyArray3(nElems))
CALL openstarfile(TRIM(MeshFileName(1))//'.cel',105)
iElem=0
DO 
  READ(105,*,iostat=os) dummy1, conn(:),dummy2,dummy3
  IF(os.NE.0) EXIT !end of file
  IF(conn(5).EQ.0) CYCLE !2D Element
  iElem=iElem+1
  !find out number of element nodes
  nElNodes=3
  DO i=4,8
    IF(conn(i).NE.conn(i-1)) nElNodes=nElNodes+1
  END DO
  iDummyArray3(iElem)=nElNodes
  iDummyArray2(iElem,1:8)=conn
  iDummyArray2(iElem,9)=dummy2 !ZoneID
END DO
ALLOCATE (BCSideBuffer(1:nBCSides,1:5))
CALL openstarfile(TRIM(MeshFileName(1))//'.bnd',107)
iBC=0
DO                            ! FaceID, Nodes, BCID, ?, ??
  READ(107,*,iostat=os) dummy1,conn(1:4),dummy2,dummy3,cdummy 
  IF(os.NE.0) EXIT !end of file
  iBC=iBC+1
  IF(conn(3).EQ.conn(4)) conn(4)=0 ! tria
  CALL Qsort1Int(conn(1:4))
  BCSideBuffer(iBC,1:4)=conn(1:4)
  BCSideBuffer(iBC,5)  =dummy2
END DO
CLOSE(107)

CALL Qsort4Int(BCSideBuffer)

maxSteps=INT(LOG(REAL(nBCSides))/LOG(2.))+1  !max steps for bisection 

DO iElem=1,nElems
  ! build elem
  IF(.NOT. ASSOCIATED(firstElem)) THEN
    CALL getNewElem(aElem)
    firstElem=>aElem
  ELSE
    CALL getNewElem(aElem%nextElem)
    aElem%nextElem%prevElem=>aElem
    aElem=>aElem%nextElem
  END IF
  nElNodes=iDummyArray3(iElem)
  aElem%nNodes=nElNodes
  ALLOCATE(aElem%Node(nElNodes))
  DO iNode=1,nElNodes
    nodeInd=iDummyArray2(iElem,NodeMap(nElNodes-3,iNode))
    IF(nodeInd.GT.nNodes) CALL abort(__STAMP__, &
               'nodeInd>nNodes,node indizes not contiguous',nodeInd,REAL(nNodes))
    IF(.NOT.ASSOCIATED(Nodes(nodeInd)%np))THEN 
      CALL getNewNode(Nodes(nodeInd)%np)
    END IF
    Nodes(nodeInd)%np%ind=nodeInd
    Nodes(nodeInd)%np%x  =nodeCoords(:,nodeInd)
    aElem%Node(iNode)%np=>Nodes(nodeInd)%np
    Nodes(nodeInd)%np%refCount=Nodes(nodeInd)%np%refCount+1
  END DO
  CALL createSides(aElem,.TRUE.)
  aElem%ind=iElem
  aElem%zone=0
  DO iZone=1,nZones
    IF(ZoneID(iZone).EQ.iDummyArray2(iElem,9)) THEN 
      aElem%zone=iZone
      EXIT
    END IF
  END DO
  IF(aElem%zone.EQ.0) CALL abort(__STAMP__, &
               'Zone Number not found in readStar',iDummyArray2(iElem,9),999.)
  ! boundary conditions 
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    !all nodes of side must have a BCnode
    found=.TRUE.
    DO iNode=1,aSide%nNodes
      IF(.NOT.ANY(NodeBC(aSide%Node(iNode)%np%ind,:))) THEN !all are false
        found=.FALSE.
        EXIT
      END IF
    END DO !iNode
    IF (.NOT.found) THEN !side has one node which has no BC
      aSide=>aSide%nextElemSide
      CYCLE
    END IF
    ! now Side has BCnodes on all nodes (must not be a BC, if side is between two BC sides!!!)
    !look up side in BCSideBuffer
    conn4=0
    DO iNode=1,aSide%nNodes
       conn4(iNode)=aSide%Node(iNode)%np%ind
    END DO
    CALL Qsort1int(conn4)
    found=.FALSE.
    !find first node ID by bisection
    IF(conn4(1).EQ.BCSideBuffer(1,1))THEN
      iBCSide=1
      searchdir=1
      found=.TRUE.
    ELSEIF(conn4(1).EQ.BCSideBuffer(nBCSides,1))THEN
      iBCSide=nBCSides
      searchdir=-1
      found=.TRUE.
    ELSE
      !bisection for first node in conn4
      iBC_low=1
      iBC_up=nBCSides
      DO i=1,maxSteps
        iBC_mid=(iBC_up-iBC_low)/2+iBC_low
        IF(conn4(1) .EQ. BCSideBuffer(iBC_mid,1))THEN
          iBCSide=iBC_mid                     !index found!
          found=.TRUE.
          EXIT
        ELSEIF(conn4(1) .GT. BCSideBuffer(iBC_mid,1))THEN ! seek in upper half
          iBC_low=iBC_mid
        ELSE
          iBC_up=iBC_mid
        END IF
      END DO
    END IF ! find first nodeID
    IF(.NOT.found) THEN
      ! this is not a BC Side!!!
      aSide=>aSide%nextElemSide
      CYCLE
    END IF
    found=.FALSE.
    searchdir=1
    IF(conn4(2).LT.BCSideBuffer(iBCSide,2))THEN
      searchdir=-1
    ELSEIF(conn4(2).EQ.BCSideBuffer(iBCSide,2))THEN
      IF(conn4(3).LT.BCSideBuffer(iBCSide,3))THEN
        searchdir=-1
      ELSEIF(conn4(3).EQ.BCSideBuffer(iBCSide,3))THEN
        IF(conn4(4).LT.BCSideBuffer(iBCSide,4))THEN
          searchdir=-1
        ELSEIF(conn4(4).EQ.BCSideBuffer(iBCSide,4))THEN
          found=.TRUE.
        END IF
      END IF
    END IF
    ! search for the side in searchdir
    DO WHILE(.NOT.found)
      IF(SUM(ABS(conn4(1:4)-BCSideBuffer(iBCSide,1:4))).EQ.0)THEN
        found=.TRUE.
        EXIT
      ELSE
        iBCSide=iBCSide+searchdir
      END IF
      IF((iBCSide.GT.nBCSides).OR.(iBCSide.LT.1)) EXIT !bounds of array
      IF(conn4(1).NE.BCSideBuffer(iBCSide,1)) EXIT  !bounds of bi
    END DO  
    IF(.NOT.found) THEN
      ! this is not a BC Side!!!
      aSide=>aSide%nextElemSide
      CYCLE
    END IF
    !check BCID (and set iBC if found)
    found=.false.
    DO iBC=1,nBCs
      IF(BCID(iBC).EQ.BCSideBuffer(iBCSide,5)) THEN
        found=.true.
        EXIT
      END IF
    END DO !iBC 
    IF(.NOT.found) CALL abort(__STAMP__, &
       'a boundary condition was used without corresponding UserDefinedBoundary',iBCSide,999.)
    IF(found) THEN
      ALLOCATE(aSide%BC)
      i=MapBC(iBC) !iBC set from exit of upper do loop
      IF(i.GT.0)THEN
        aSide%BC%BCType     = BoundaryType(i,1)
        aSide%curveIndex    = BoundaryType(i,2)
        aSide%BC%BCState    = BoundaryType(i,3)
        aSide%BC%BCalphaInd = BoundaryType(i,4)
        aSide%BC%BCindex    = i
      ELSE
        CALL abort(__STAMP__, 'hier waere noch ein Problem bei den BCs...',999,999.)
      END IF
    END IF
    aSide=>aSide%nextElemSide
  END DO
END DO
DEALLOCATE (BCSideBuffer)
DEALLOCATE(nodeCoords)
DO iNode=1,nNodes
  NULLIFY(Nodes(iNode)%np)
END DO
DEALLOCATE(Nodes)
DEALLOCATE(iDummyArray2,iDummyArray3)

CALL Timer(.FALSE.)
WRITE(UNIT_stdOut,'(132("~"))')
END SUBROUTINE ReadStar


SUBROUTINE openstarfile(FileName,unit_in)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)       :: FileName  ! ?
INTEGER,INTENT(IN)                :: unit_in  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                :: os !openStatus
CHARACTER(LEN=255)     :: cdummy  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,*)'Opening Star file: ',FileName
OPEN(UNIT   = unit_in,                  &
     FILE   = FileName,             &
     STATUS = 'OLD',                &
     ACTION = 'READ',               &
     ACCESS = 'SEQUENTIAL',         &
     IOSTAT = os                    )
IF(os.NE.0) THEN  ! File Error
  CALL abort(__STAMP__, &
       'ERROR: cannot open star file: '//trim(FileName),999,999.)
END IF
!V4 or V3
READ(unit_in,*)cdummy
IF(INDEX(cdummy,'PROSTAR').NE.0)THEN !star V4
  CALL abort(__STAMP__, &
       'ERROR: star file of version V4, not implemented jet!',999,999.)
ELSE  !star V3
  REWIND(unit_in)
END IF
END SUBROUTINE openstarfile


SUBROUTINE readStar_split(firstElem_in,FileName)
!===================================================================================================================================
! Read splitted surface mesh from ansa ascii file. Called by fillMesh.
! we buld up a connected surface mesh of the splitted elements, which is then an element list starting
! at firstElem_in. 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tSidePtr,tNodePtr
USE MOD_Mesh_Vars,ONLY:ElemCount,SideCount
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewSide,getNewNode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(INOUT)           :: FirstElem_in  ! ?
CHARACTER(LEN=255),INTENT(IN) :: FileName  ! ?
! nMeshFiles                : Number of mesh fil  ! ?es (INI-File)
! MeshFileName(1)               : FileName of mesh file only 1  ! ?
! nZones                    : Number of mesh zones (INI-File)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!   
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                :: os !openStatus  ! ?
INTEGER                :: iNode,iSide,nNodes ! counter/number of Nodes
INTEGER                :: iElem,nElems,nElNodes,localnElems  ! ?
INTEGER                :: dummy1,dummy2,dummy3  ! ?
INTEGER                :: nodeInd  ! ?
INTEGER                :: conn(8) !connectivity list of eight nodes
INTEGER,ALLOCATABLE    :: ElemConnect(:,:),ElemnNodes(:)  ! ?
REAL                   :: x(3)  ! ?
REAL,ALLOCATABLE       :: nodeCoords(:,:)  ! ?
TYPE(tNodePtr),POINTER :: Nodes(:)  ! ?
TYPE(tElem),POINTER    :: aElem   ! ?
TYPE(tSide),POINTER    :: aSide   ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')'Reading Star(ANSA) mesh for splitted elements...'

!read number nodes
CALL openstarfile(TRIM(FileName)//'.vrt',104)
os=0
nNodes=0
DO 
  READ(104,*,iostat=os) dummy1,x
  IF(os.NE.0) EXIT !end of file
  nNodes=nNodes+1
END DO
CLOSE(104)

! read Elements
 !read number of 2D Elements
CALL openstarfile(TRIM(FileName)//'.cel',105)
os=0
nElems=0
DO                    ! ElemID, 8 Nodes, PID, ???
  READ(105,*,iostat=os) dummy1, conn(:),dummy2,dummy3
  IF(os.NE.0) EXIT !end of file
  IF(conn(5).NE.0) THEN
  STOP 'not surface elements!'
  END IF
  nElems=nElems+1
END DO
CLOSE(105)
WRITE(*,*)'total number of SplitElems',nElems 
ALLOCATE(nodeCoords(3,nNodes))
nodeCoords=0.
CALL openstarfile(TRIM(FileName)//'.vrt',104)
DO iNode = 1,nNodes   ! READ nodes
  READ(104,*)dummy1,nodeCoords(:,iNode)
END DO
CLOSE(104)

! node pointers!
ALLOCATE(Nodes(nNodes))
DO iNode=1,nNodes
  NULLIFY(Nodes(iNode)%np)
END DO
aElem=>firstElem_in ! READ elems
IF(ASSOCIATED(aElem))THEN
  DO WHILE(ASSOCIATED(aElem%nextElem))
    aElem%Ind=-1
    aElem=>aElem%nextElem
  END DO
  aElem%Ind=-1
END IF
ALLOCATE(ElemConnect(nElems,5),ElemnNodes(nElems))
CALL openstarfile(TRIM(FileName)//'.cel',105)
iElem=0
DO 
  READ(105,*,iostat=os) dummy1, conn(:),dummy2,dummy3
  IF(os.NE.0) EXIT !end of file
  IF(conn(5).NE.0) THEN
     CALL abort(__STAMP__, &
       'ERROR: split element list contains volumes!'//TRIM(FileName)//'.cel',999,999.)
  END IF      
  iElem=iElem+1
  !find out number of element nodes
  nElNodes=3
  IF(conn(4).NE.conn(3)) nElNodes=nElNodes+1
  ElemnNodes(iElem)=nElNodes
  ElemConnect(iElem,1:4)=conn(1:4)
  ElemConnect(iElem,5)=dummy2 !ZoneID
END DO
localnElems=0
DO iElem=1,nElems
  localnElems=localnElems+1
  IF(.NOT. ASSOCIATED(firstElem_in)) THEN
    CALL getNewElem(aElem)
    firstElem_in=>aElem
  ELSE
    CALL getNewElem(aElem%nextElem)
    Elemcount=elemcount-1 !undo global element count
    aElem%nextElem%prevElem=>aElem
    aElem=>aElem%nextElem
  END IF
  nElNodes=ElemnNodes(iElem)
  aElem%nNodes=nElNodes
  ALLOCATE(aElem%Node(nElNodes))
! WRITE(*,'(A,6I4)')'DEBUG,iElem,ElemConnect',iElem,ElemConnect(iElem,:)
  DO iNode=1,nElNodes
! WRITE(*,'(A,3F10.5)')'DEBUG,xnode',nodeCoords(:,ElemConnect(iElem,iNode))
  END DO
  DO iNode=1,nElNodes
    nodeInd=ElemConnect(iElem,iNode)
    IF(nodeInd .GT. nNodes) CALL abort(__STAMP__, &
       'ERROR: nodeInd>nNodes, node indizes are not contiguous!',nodeInd,REAL(nNodes)) 
    IF(.NOT.ASSOCIATED(Nodes(nodeInd)%np))THEN 
      CALL getNewNode(Nodes(nodeInd)%np)
      Nodes(nodeInd)%np%ind=nodeInd
      Nodes(nodeInd)%np%x  =nodeCoords(:,nodeInd)
    END IF
    aElem%Node(iNode)%np=>Nodes(nodeInd)%np
    Nodes(nodeInd)%np%refCount=Nodes(nodeInd)%np%refCount+1
  END DO
  aElem%ind=iElem
  ! buildSides
  DO iSide=1,nElNodes
    IF (iSide .EQ. 1) THEN
      CALL getNewSide(aElem%firstSide,2)
      sidecount=sidecount-1
      aSide=>aElem%firstSide
    ELSE
      CALL getNewSide(aSide%nextElemSide,2)  
      sidecount=sidecount-1
      aSide=>aSide%nextElemSide
    END IF
    aSide%LocSide=iSide
    aSide%Elem=>aElem
    DO iNode=1,2
      aSide%Node(iNode)%np=>aElem%Node(MOD(iSide+iNode-2,nElNodes)+1)%np
      aSide%Node(iNode)%np%refCount=aSide%Node(iNode)%np%refCount+1
    END DO
  END DO
END DO
WRITE(*,*)' SplitElems and Nodes: ',localnElems,nNodes
DEALLOCATE(nodeCoords)
DO iNode=1,nNodes
  NULLIFY(Nodes(iNode)%np)
END DO
DEALLOCATE(Nodes)
DEALLOCATE(ElemConnect,ElemnNodes)
CALL Timer(.FALSE.)
WRITE(UNIT_stdOut,'(132("~"))')
END SUBROUTINE readStar_split

END MODULE MOD_Readin_ANSA
