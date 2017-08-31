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
! Copyright (C) 2017  Florian Hindenlang <hindenlang@gmail.com>
! Copyright (C) 2015  Prof. Claus-Dieter Munz <munz@iag.uni-stuttgart.de>
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
MODULE MOD_Output_HDF5
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_IO_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE WriteMeshToHDF5
  MODULE PROCEDURE WriteMeshToHDF5
END INTERFACE

INTERFACE SpaceFillingCurve
  MODULE PROCEDURE SpaceFillingCurve
END INTERFACE

INTERFACE WriteAttribute
  MODULE PROCEDURE WriteAttribute
END INTERFACE

PUBLIC::WriteMeshToHDF5,SpaceFillingCurve,WriteAttribute
!===================================================================================================================================

CONTAINS
SUBROUTINE WriteMeshToHDF5(FileString)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:N
USE MOD_Output_Vars,ONLY:dosortIJK
USE MOD_Mesh_Vars,ONLY:nUserDefinedBoundaries,BoundaryName,BoundaryType
USE MOD_Mesh_Basis,ONLY:ISORIENTED
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: FileString  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER            :: Elem  ! ?
TYPE(tSide),POINTER            :: Side  ! ?
INTEGER                        :: ElemID,SideID,NodeID  ! ?
INTEGER                        :: locnSides
INTEGER                        :: iNode,i,iMortar
LOGICAL                        :: found
CHARACTER(LEN=26)              :: ElemTypeName(1:11)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')' WRITE DATA TO HDF5 FILE...'


!set all node and side indices =0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  DO i=1,Elem%nNodes
    Elem%Node(i)%np%ind=0
  END DO
  DO i=1,Elem%nCurvedNodes
    Elem%curvedNode(i)%np%ind=0
  END DO
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    Side%ind=0
    !CURVED
    DO i=1,side%nCurvedNodes
      side%curvedNode(i)%np%ind=0
    END DO
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO


! count Elements , unique sides and nodes are marked with ind=0
nNodeIDs=0 !number of unique nodeIDs
nSideIDs=0 !number of unique side IDs (side and side%connection have the same sideID)
nElems=0   !number of elements
nSides=0   !number of all sides
nNodes=0   !number of all nodes

Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  nElems=nElems+1
  ! Count nodes
  DO i=1,Elem%nNodes
    IF(Elem%Node(i)%np%ind.NE.0) CYCLE
    nNodeIDs=nNodeIDs+1
    Elem%Node(i)%np%ind=-88888  ! mark no MPI side
  END DO
  DO i=1,Elem%nCurvedNodes
    IF(Elem%CurvedNode(i)%np%ind.NE.0) CYCLE
    nNodeIDs=nNodeIDs+1
    Elem%CurvedNode(i)%np%ind=-88888
  END DO

  ! Count sides
  locnSides=0
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    locnSides = locnSides + 1
    IF(side%ind.EQ.0) THEN
      IF(Side%MortarType.EQ.0)THEN
        nSideIDs=nSideIDs+1
        Side%ind=-88888
        IF(ASSOCIATED(Side%connection))THEN      
          Side%connection%ind = -88888 ! count inner and periodic sides only once 
        END IF
      ELSEIF(Side%MortarType.GT.0)THEN
        locnSides = locnSides + Side%nMortars
        nSideIDs=nSideIDs+1
        Side%ind=-88888
        DO iMortar=1,Side%nMortars
          IF(Side%MortarSide(iMortar)%sp%ind.EQ.0)THEN
            nSideIDs=nSideIDs+1
            Side%MortarSide(iMortar)%sp%ind=-88888
          END IF 
        END DO !iMortar
      ELSE
        nSideIDs=nSideIDs+1
        Side%ind=-88888
      END IF
    END IF
    Side=>Side%nextElemSide
  END DO
  IF(N.EQ.1)THEN
    nNodes = nNodes+Elem%nNodes
  ELSE
    nNodes = nNodes+Elem%nCurvedNodes
  END IF
  nSides = nSides+locnSides
  Elem=>Elem%nextElem
END DO


!NOW CALLED IN FILLMESH!!
!! prepare sorting by space filling curve
!! NOTE: SpaceFillingcurve is not used, if existing hdf5 mesh is read in and the sorting should stay identical
!IF(useSpaceFillingCurve)THEN
!  CALL SpaceFillingCurve()
!END IF

IF(ALLOCATED(ElemBarycenters)) DEALLOCATE(ElemBarycenters)
ALLOCATE(ElemBarycenters(1:nElems,3))

!set unique nodes and Side Indices
ElemID=0
SideID=0
NodeID=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  ElemID=ElemID+1 
  Elem%ind=ElemID
  ElemBarycenters(ElemID,:)=0.
  DO i=1,Elem%nNodes
    ElemBarycenters(ElemID,:)=ElemBarycenters(ElemID,:)+Elem%Node(i)%np%x
    IF(Elem%Node(i)%np%ind.NE.-88888) CYCLE
    NodeID=NodeID+1
    Elem%Node(i)%np%ind=NodeID
  END DO
  ElemBarycenters(ElemID,:)=ElemBarycenters(ElemID,:)/REAL(Elem%nNodes)
  DO i=1,Elem%nCurvedNodes
    IF(Elem%CurvedNode(i)%np%ind.NE.-88888) CYCLE
    NodeID=NodeID+1
    Elem%CurvedNode(i)%np%ind=NodeID
  END DO

  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(side%ind.EQ.-88888) THEN  ! assign side ID 
      IF(Side%MortarType.EQ.0)THEN
        SideID=SideID+1
        Side%ind=SideID
        IF(ASSOCIATED(Side%connection))THEN      
          IF(Side%connection%ind.EQ.-88888) Side%connection%ind=SideID ! count inner and periodic sides only once 
        END IF
      ELSEIF(Side%MortarType.GT.0)THEN
        SideID=SideID+1
        Side%ind=SideID
        DO iMortar=1,Side%nMortars
          IF(Side%MortarSide(iMortar)%sp%ind.EQ.-88888)THEN
            SideID=SideID+1
            Side%MortarSide(iMortar)%sp%ind=SideID
          END IF 
        END DO !iMortar
      ELSE
        SideID=SideID+1
        Side%ind=SideID
      END IF
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO !Elem

IF(NodeID.NE.nNodeIDs) CALL abort(__STAMP__,&
                     'Sanity check: max(nodeID <> nNodeIDs!')
IF(SideID.NE.nSideIDs) CALL abort(__STAMP__,&
                     'Sanity check: max(sideID <> nSideIDs!')
IF(ElemID.NE.nElems) CALL abort(__STAMP__,&
                     'Sanity check: max(elemID <> nElems!')


!set Side Flip 
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    Side%flip=0
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%flip.GT.0)THEN
      Side=>Side%nextElemSide
      CYCLE
    END IF
    IF(.NOT.ISORIENTED(Side))THEN
      found=.FALSE.
      DO iNode=1,Side%nNodes
        IF(ASSOCIATED(Side%Node(iNode)%np,Side%OrientedNode(1)%np))THEN
          found=.TRUE.
          EXIT
        END IF
      END DO 
      IF(.NOT.found) STOP 'Flip not found'
      Side%flip=iNode
      IF(.NOT.ASSOCIATED(Side%connection)) CALL ABORT(__STAMP__, &
        'Side connection should be associated for non-oreinted side')
      IF (Side%connection%MortarType.LE.0) Side%connection%flip=iNode !flip is the same for the connection
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

CALL getMeshInfo() !allocates and fills ElemInfo,SideInfo,NodeInfo,NodeCoords

! Create the file
CALL OpenHDF5File(FileString,create=.TRUE.)  

!attributes 
CALL WriteAttribute(File_ID,'Version',1,RealScalar=1.0)
CALL WriteAttribute(File_ID,'Ngeo',1,IntScalar=N)
CALL WriteAttribute(File_ID,'nElems',1,IntScalar=nElems)
CALL WriteAttribute(File_ID,'nSides',1,IntScalar=nSides)
CALL WriteAttribute(File_ID,'nNodes',1,IntScalar=nNodes)
CALL WriteAttribute(File_ID,'nUniqueSides',1,IntScalar=nSideIDs)
CALL WriteAttribute(File_ID,'nUniqueNodes',1,IntScalar=nNodeIDs)

!WRITE ElemInfo,into (1,nElems)  
CALL WriteArrayToHDF5(File_ID,'ElemInfo',2,(/ElemInfoSize,nElems/),IntegerArray=ElemInfo)
DEALLOCATE(ElemInfo)

!WRITE SideInfo,into (1,nSides)   
CALL WriteArrayToHDF5(File_ID,'SideInfo',2,(/SideInfoSize,nSides/),IntegerArray=SideInfo)
DEALLOCATE(SideInfo)

! WRITE NodeCoords and NodeIDs
CALL WriteArrayToHDF5(File_ID,'NodeCoords',2,(/3,nNodes/),RealArray=NodeCoords)
CALL WriteArrayToHDF5(File_ID,'GlobalNodeIDs',1,(/nNodes/),IntegerArray=GlobalNodeIDs)
DEALLOCATE(GlobalNodeIDs)
DEALLOCATE(NodeCoords)

nBCs=nUserDefinedBoundaries
ALLOCATE(BCNames(nBCs))
ALLOCATE(BCType(4,nBCs))
DO i=1,nBCs
  BCNames(i)=BoundaryName(i) 
  BCType(:,i)=BoundaryType(i,:) 
END DO
! WRITE BC 
CALL WriteAttribute(File_ID,'nBCs',1,IntScalar=nBCs)

CALL WriteArrayToHDF5(File_ID,'BCNames',1,(/nBCs/),StrArray=BCNames)
CALL WriteArrayToHDF5(File_ID,'BCType',2,(/4,nBcs/),IntegerArray=BCType)

DEALLOCATE(BCNames)
DEALLOCATE(BCType)

!WRITE ElemWeight,into (1,nElems)  
CALL WriteArrayToHDF5(File_ID,'ElemBarycenters',2,(/3,nElems/),RealArray=TRANSPOSE(ElemBarycenters))
DEALLOCATE(ElemBarycenters)

CALL WriteArrayToHDF5(File_ID,'ElemWeight',1,(/nElems/),RealArray=ElemWeight)
DEALLOCATE(ElemWeight)

CALL WriteArrayToHDF5(File_ID,'ElemCounter',2,(/2,11/),IntegerArray=ElemCounter)
!WRITE(*,'(A)')'Mesh statistics:'
WRITE(*,'(A40)')   &
        '    ____________________________________'
WRITE(*,'(A40)') &
        '    #elements  .......... of type:      '
WRITE(*,'(A40)')   &
        '    ------------------------------------'
ElemTypeName(1:11)= (/' straight-edge Tetrahedra ', &
                      '        curved Tetrahedra ', &
                      '  planar-faced Prisms     ', &
                      ' straight-edge Prisms     ', &
                      '        curved Prisms     ', &
                      '  planar-faced Pyramids   ', &
                      ' straight-edge Pyramids   ', &
                      '        curved Pyramids   ', &
                      '  planar-faced Hexahedra  ', &
                      ' straight-edge Hexahedra  ', &
                      '        curved Hexahedra  '/)
DO i=1,11
IF(ElemCounter(2,i).GT.0) &
  WRITE(*,'(A4,I9,A26)')'    ',Elemcounter(2, i),ElemTypeName(i)
END DO !i=1,11
WRITE(*,'(A40)')   &
        '    ____________________________________'

IF(dosortIJK)THEN
  ! WRITE element ijk index (for postprocessing of structured/semistructured domains)
  CALL WriteArrayToHDF5(File_ID,'nElems_IJK',1,(/3/),IntegerArray=nElems_IJK)
  CALL WriteArrayToHDF5(File_ID,'Elem_IJK',2,(/3,nElems/),IntegerArray=TRANSPOSE(Elem_IJK))
  DEALLOCATE(Elem_IJK)
END IF

! Close the file.
CALL CloseHDF5File()
CALL Timer(.FALSE.)

END SUBROUTINE WriteMeshToHDF5


SUBROUTINE getMeshInfo()
!===================================================================================================================================
! Subroutine prepares ElemInfo,Sideinfo,Nodeinfo,NodeCoords arrays
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:N
USE MOD_Mesh_Basis,ONLY:ISORIENTED
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER            :: Elem  ! ?
TYPE(tSide),POINTER            :: Side  ! ?
INTEGER                        :: locnNodes,locnSides
INTEGER                        :: iNode,iSide,iElem,i,iMortar
TYPE(tSide),POINTER            :: aSide
!===================================================================================================================================
!fill ElementInfo. 
ALLOCATE(ElemInfo(ElemInfoSize,1:nElems))
ALLOCATE(ElemWeight(1:nElems))
ElemInfo=0
ElemWeight=1.  !equal weight for all elements
Elemcounter=0
Elemcounter(1,:)=(/104,204,105,115,205,106,116,206,108,118,208/)
iNode  = 0 
iSide  = 0
iElem  = 0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  iElem=iElem+1
  locnSides=nSidesElem(Elem%nNodes)

  aSide=>Elem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    locnSides = locnSides + aSide%nMortars
    aSide=>aSide%nextElemSide
  END DO 


  IF(N.EQ.1)THEN
    locnNodes=Elem%nNodes
  ELSE
    locnNodes=Elem%nCurvedNodes
  END IF
  SELECT CASE(elem%type)
    CASE(104) !linear tet
      elemcounter(2,1)=elemcounter(2,1)+1
    CASE(204) !spline tet
      elemcounter(2,2)=elemcounter(2,2)+1
    CASE(105) !linear pyr
      elemcounter(2,3)=elemcounter(2,3)+1
    CASE(115) !non-linear pyr
      elemcounter(2,4)=elemcounter(2,4)+1
    CASE(205) !spline pyr
      elemcounter(2,5)=elemcounter(2,5)+1
    CASE(106) !linear prism
      elemcounter(2,6)=elemcounter(2,6)+1
    CASE(116) !non-linear prism
      elemcounter(2,7)=elemcounter(2,7)+1
    CASE(206) !spline prism
      elemcounter(2,8)=elemcounter(2,8)+1
    CASE(108) !linear hex
      elemcounter(2,9)=elemcounter(2,9)+1
    CASE(118) !non-linear hex
      elemcounter(2,10)=elemcounter(2,10)+1
    CASE(208) !spline hex
      elemcounter(2,11)=elemcounter(2,11)+1
  END SELECT
  ElemInfo(ELEM_Type,iElem)         = Elem%Type        ! Element Type
  ElemInfo(ELEM_Zone,iElem)         = Elem%Zone        ! Zone Number
  ElemInfo(ELEM_FirstSideInd,iElem) = iSide            ! first index -1 in SideInfo
  ElemInfo(ELEM_LastSideInd,iElem)  = iSide+locnSides  ! last index in SideInfo
  ElemInfo(ELEM_FirstNodeInd,iElem) = iNode            ! first index -1 in NodeInfo
  ElemInfo(ELEM_LastNodeInd,iElem)  = iNode+locnNodes  ! last index in NodeInfo
  iNode = iNode + locnNodes
  iSide = iSide + locnSides
  Elem=>Elem%nextElem
END DO


!fill SideInfo
ALLOCATE(SideInfo(SideInfoSize,1:nSides)) 
SideInfo=0 
iSide=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    iSide=iSide+1
    !Side Tpye
    IF(Side%isCurved)THEN
      IF(N.GT.1)THEN
        SideInfo(SIDE_Type,iSide)=20+Side%nNodes         ! Side Type: NL tria/quad, 6/7
      ELSE
        IF(Side%nNodes.EQ.3)THEN
          SideInfo(SIDE_Type,iSide)=3
        ELSE  
          SideInfo(SIDE_Type,iSide)=10+Side%nNodes        ! Side Type: bilinear quad
        END IF
      END IF
    ELSE
      SideInfo(SIDE_Type,iSide)=Side%nNodes             ! Side Type: linear 3/4
    END IF
    !Side ID
    SideInfo(SIDE_ID,iSide)=Side%ind
    IF(.NOT.ISORIENTED(Side)) SideInfo(SIDE_ID,iSide)=-SideInfo(SIDE_ID,iSide)           
    !BCID
    IF (ASSOCIATED(Side%BC)) THEN
      SideInfo(SIDE_BCID,    iSide)= Side%BC%BCIndex
      IF(Side%BC%BCIndex.EQ.0) WRITE(*,*)'DEBUG, Warning, BC ind =0'
    ELSE 
      SideInfo(SIDE_BCID,    iSide)= 0
    END IF
      
    IF (Side%MortarType.GT.0) THEN ! Mortar master side (only implemented for Quad-sides!!!)
      IF(ASSOCIATED(Side%Connection)) CALL abort(__STAMP__,&
                                                 'Mortar master with connection is not allowed')
      IF(Side%flip.NE.0) STOP 'Problem with flip on mortar'
      SideInfo(SIDE_nbElemID,iSide)= -Side%MortarType
      SideInfo(SIDE_nbLocSide_flip,iSide)=0
      DO iMortar=1,Side%nMortars 
        iSide=iSide+1
        SideInfo(SIDE_Type,    iSide)= MERGE(104,204,N.EQ.1) 
        SideInfo(SIDE_ID,      iSide)= Side%MortarSide(iMortar)%sp%ind       ! small are always master
        SideInfo(SIDE_nbElemID,iSide)= Side%MortarSide(iMortar)%sp%Elem%ind  ! neighbour Element ID
        SideInfo(SIDE_nbLocSide_flip,iSide)=0
        IF (ASSOCIATED(Side%MortarSide(iMortar)%sp%BC)) THEN
          SideInfo(SIDE_BCID,    iSide)= Side%MortarSide(iMortar)%sp%BC%BCIndex
        ELSE 
          SideInfo(SIDE_BCID,    iSide)= 0
        END IF
      END DO
    ELSE IF (Side%MortarType.LT.0) THEN ! small Mortar side (only implemented for Quad-sides!!!)
       i=MERGE(104,204,N.EQ.1)          ! 104:bilinear,204:curved
       SideInfo(SIDE_Type,iSide)= -i    ! mark side with mortar neighbour
       SideInfo(SIDE_ID,  iSide)= -Side%ind ! neg. sideID for slaves
       IF(ASSOCIATED(Side%Connection))THEN
         SideInfo(SIDE_nbElemID,iSide)=Side%Connection%Elem%ind    ! neighbour Element ID
       END IF
       SideInfo(SIDE_nbLocSide_flip,iSide)=Side%flip
    ELSE
      !neighbor Element ID
      IF(ASSOCIATED(Side%Connection))THEN
        SideInfo(SIDE_nbElemID,iSide)=Side%Connection%Elem%ind                   ! Element ID of neighbor Element
        SideInfo(SIDE_nbLocSide_flip,iSide)=10*Side%connection%locSide+Side%connection%flip
      END IF
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

IF(iSide.NE.nSides) CALL abort(__STAMP__,&
                     'Sanity check: nSides not equal to total number of sides!')

!fill GlobalNodeID
ALLOCATE(NodeCoords(3,nNodes),GlobalNodeIDs(nNodes))
iNode=0
IF(N.EQ.1)THEN !linear
  Elem=>firstElem
  DO WHILE(ASSOCIATED(Elem))
    DO i=1,Elem%nNodes
      iNode=iNode+1
      NodeCoords(:,iNode)=Elem%Node(LinMap(i,Elem%nNodes))%np%x
      GlobalNodeIDs(iNode)=Elem%Node(LinMap(i,Elem%nNodes))%np%ind
    END DO
    Elem=>Elem%nextElem
  END DO
ELSE ! CurvedNodes
  Elem=>firstElem
  DO WHILE(ASSOCIATED(Elem))
    DO i=1,Elem%nCurvedNodes  !CurvedNodes already in IJK ordering
      iNode=iNode+1
      NodeCoords(:,iNode)=Elem%curvedNode(i)%np%x
      GlobalNodeIDs(iNode)=Elem%curvedNode(i)%np%ind
    END DO
    Elem=>Elem%nextElem
  END DO
END IF 
  
IF(iNode.NE.nNodes) CALL abort(__STAMP__,&
                     'Sanity check: nNodes not equal to total number of nodes!')

END SUBROUTINE getMeshinfo

SUBROUTINE spaceFillingCurve(nElems_in)
!===================================================================================================================================
! Subroutine prepares elementlist and barycentric coodrinates for element sorting by space filling curve and maps back to pointer
! srtructure
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElemPtr
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:MeshMode
USE MOD_Mesh_Vars,ONLY:AdaptedMesh
USE MOD_Output_Vars,ONLY:DebugVisu,dosortijk,sfc_boundbox
USE MOD_SpaceFillingCurve,ONLY:SortElemsBySpaceFillingCurve
USE MOD_sortIJK,ONLY:SortElemsByCoords
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nElems_in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: IDlist(1:nElems_in)  ! ?
TYPE(tElemPtr)                 :: Elems(1:nElems_in)  ! ?
REAL                           :: ElemBary(1:nElems_in,3)  ! ?
INTEGER                        :: ElemID  ! ?
INTEGER                        :: iNode 
!===================================================================================================================================

Elems(1)%ep=>firstElem
DO ElemID=2,nElems_in 
  Elems(ElemID)%ep=>Elems(ElemID-1)%ep%nextElem
END DO
DO ElemID=1,nElems_in
  IDList(ElemID)=ElemID 
  ElemBary(ElemID,:)=0.
  DO iNode=1,Elems(ElemID)%ep%nNodes
    ElemBary(ElemID,:)=ElemBary(ElemID,:)+Elems(ElemID)%ep%Node(iNode)%np%x
  END DO
  ElemBary(ElemID,:)=ElemBary(ElemID,:)/REAL(Elems(ElemID)%ep%nNodes)
END DO

IF((MeshMode.EQ.11).AND. (.NOT.AdaptedMesh))THEN 
  ! for Meshmode=11: if no splitting was done, this is a structured single block, elem_IJK already defined
  CALL SortElemsBySpaceFillingCurve(nElems_in,REAL(Elem_IJK),IDList,1) !use IJK for space filling curve
ELSE
  CALL SortElemsBySpaceFillingCurve(nElems_in,ElemBary,IDList,sfc_boundbox)
END IF

NULLIFY(Elems(IDlist(1))%ep%prevElem)
firstElem=>Elems(IDlist(1))%ep
DO ElemID=2,nElems_in 
  Elems(IDlist(ElemID-1))%ep%nextElem=>Elems(IDList(ElemID))%ep
  Elems(IDlist(ElemID))%ep%prevElem  =>Elems(IDList(ElemID-1))%ep
END DO
DO ElemID=1,nElems_in 
  Elems(IDlist(ElemID))%ep%ind=ElemID
END DO
IF(DebugVisu)THEN
  WRITE(*,*)'write space filling curve to sfc.dat'
  OPEN(UNIT=200,FILE='sfc.dat',STATUS='REPLACE')
  DO ElemID=1,nElems_in 
    WRITE(200,'(3E21.6)')ElemBary(IDlist(ElemID),:)
  END DO
  CLOSE(200)
END IF
IF(ALLOCATED(ElemBarycenters)) DEALLOCATE(ElemBarycenters)
ALLOCATE(ElemBarycenters(1:nElems_in,3))
DO ElemID=1,nElems_in 
  ElemBarycenters(ElemID,:)=ElemBary(IDlist(ElemID),:)
END DO
NULLIFY(Elems(IDlist(nElems_in))%ep%nextElem)
IF((MeshMode.EQ.11).AND.(.NOT.AdaptedMesh))THEN ! for Meshmode=11: structured single block, elem_IJK already defined
  !sort by spacefillingcurve
  Elem_IJK(:,1)=Elem_IJK(IDList(:),1)
  Elem_IJK(:,2)=Elem_IJK(IDList(:),2)
  Elem_IJK(:,3)=Elem_IJK(IDList(:),3)
ELSE
  IF(dosortijk) THEN
    ! do also directly the ijk coordinates of the elements
    IF(ALLOCATED(Elem_IJK)) DEALLOCATE(Elem_IJK)
    ALLOCATE(Elem_IJK(1:nElems_in,3))
    CALL SortElemsByCoords(nElems_in,ElemBarycenters(:,:),nElems_IJK,Elem_IJK)
  END IF
END IF
END SUBROUTINE spaceFillingCurve


! HFD5 STUFF
SUBROUTINE WriteArrayToHDF5(Loc_ID,ArrayName,Rank,nVal,RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to write Data to HDF5 format (ONLY FOR SINGLE USE)
!===================================================================================================================================
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T), INTENT(IN)     :: Loc_ID                                ! HDF5 file ID or dataset ID
CHARACTER(LEN=*), INTENT(IN)   :: ArrayName                             ! Name of the array
INTEGER,INTENT(IN)             :: Rank                                  ! number of dimensions
INTEGER,INTENT(IN)             :: nVal(Rank)                            ! dimensions of the array
REAL              ,DIMENSION(Rank),OPTIONAL,INTENT(IN) :: RealArray     ! Real array  
INTEGER           ,DIMENSION(Rank),OPTIONAL,INTENT(IN) :: IntegerArray  ! Integer array
CHARACTER(LEN=255),DIMENSION(Rank),OPTIONAL,INTENT(IN) :: StrArray      ! String array
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iError  ! ?
INTEGER(HID_T)                 :: PList_ID,DSet_ID,MemSpace,FileSpace,HDF5DataType  ! ?
INTEGER(HSIZE_T)               :: Dimsf(Rank)  ! ?
INTEGER(HSIZE_T)               :: offset(Rank)  ! offset of local data
!===================================================================================================================================
LOGWRITE(UNIT_stdOut,'(A,I1.1,A,A,A)')' WRITE ',Rank,'D ARRAY "',TRIM(ArrayName),'" TO HDF5 FILE...'

! Get global array size, always first dimension!!
Dimsf=nVal

! Create the data space for the  dataset.
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, FileSpace, iError)
! Create the dataset with default properties.
IF(PRESENT(RealArray))     HDF5DataType=H5T_NATIVE_DOUBLE
IF(PRESENT(IntegerArray))  HDF5DataType=H5T_NATIVE_INTEGER
IF(PRESENT(StrArray))THEN
  ! Create HDF5 datatype for the character array.
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, HDF5DataType, iError)
  SizeSet=255
  CALL H5TSET_SIZE_F(HDF5DataType, SizeSet, iError)
END IF

CALL H5DCREATE_F(Loc_ID, TRIM(ArrayName), HDF5DataType, filespace, DSet_ID, iError)
CALL H5SCLOSE_F(FileSpace, iError)

! Each process defines dataset in memory and writes it to the hyperslab in the file.
Dimsf=nVal  ! Now we need the local array size
offset=0  !ONLY SINGLE 
! Create the data space in the memory
IF(Dimsf(1) .NE. 0)THEN
  CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
ELSE
  CALL H5SCREATE_F(H5S_NULL_F,MemSpace,iError)
END IF
! Select hyperslab in the file.
CALL H5DGET_SPACE_F(DSet_id, FileSpace, iError)
IF(Dimsf(1) .NE. 0)THEN
  CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, Offset, Dimsf, iError)
ELSE
  CALL H5SSELECT_NONE_F(FileSpace,iError)
END IF

! Create property list for collective dataset write
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
!Write the dataset collectively.

IF(PRESENT(IntegerArray))THEN
  CALL H5DWRITE_F(DSet_ID,HDF5DataType,IntegerArray,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(RealArray))THEN
  CALL H5DWRITE_F(DSet_ID,HDF5DataType,RealArray   ,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF
IF(PRESENT(StrArray))THEN
  CALL H5DWRITE_F(DSet_ID,HDF5DataType,StrArray    ,Dimsf,iError,file_space_id=filespace,mem_space_id=memspace,xfer_prp=PList_ID)
END IF

! Close the property list.
CALL H5PCLOSE_F(PList_ID, iError)
! Close dataspaces.
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5SCLOSE_F(MemSpace, iError)
! Close the dataset.
CALL H5DCLOSE_F(DSet_ID, iError)

LOGWRITE(UNIT_stdOut,*)'...DONE!'

END SUBROUTINE WriteArrayToHDF5



!==================================================================================================================================
!> Subroutine to write Attributes to HDF5 format of a given Loc_ID, which can be the File_ID,datasetID,groupID. This must be opened
!> outside of the routine. If you directly want to write an attribute to a dataset, just provide the name of the dataset
!==================================================================================================================================
SUBROUTINE WriteAttribute(Loc_ID_in,AttribName,nVal,DataSetname,&
                          RealScalar,IntScalar,StrScalar,LogicalScalar, &
                          RealArray,IntArray,StrArray)
! MODULES
USE MOD_Globals
USE,INTRINSIC :: ISO_C_BINDING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER(HID_T)    ,INTENT(IN)           :: Loc_ID_in              !< Dataset ID (only if already open)
CHARACTER(LEN=*)  ,INTENT(IN)           :: AttribName             !< name of the attribute to be written
INTEGER           ,INTENT(IN)           :: nVal                   !< number of array entries if array is written
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL  :: DatasetName            !< name of the dataset created
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealScalar       !< real scalar
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntScalar        !< integer scalar
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL,TARGET :: StrScalar(1)     !< scalar string
LOGICAL           ,INTENT(IN),OPTIONAL        :: LogicalScalar    !< logical scalar
REAL              ,INTENT(IN),OPTIONAL,TARGET :: RealArray(nVal)  !< real array of length nVal
INTEGER           ,INTENT(IN),OPTIONAL,TARGET :: IntArray(nVal)   !< integer array of length nVal
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL,TARGET :: StrArray(nVal)   !< string array of length nVal
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: Rank
INTEGER(HID_T)                 :: DataSpace,Attr_ID,Loc_ID,Type_ID
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf
INTEGER(SIZE_T)                :: AttrLen
INTEGER,TARGET                 :: logtoint
TYPE(C_PTR)                    :: buf
!==================================================================================================================================
LOGWRITE(*,*)' WRITE ATTRIBUTE "',TRIM(AttribName),'" TO HDF5 FILE...'
IF(PRESENT(DataSetName))THEN
  ! Open dataset
  IF(TRIM(DataSetName).NE.'') CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
ELSE
  Loc_ID=Loc_ID_in
END IF
! Create scalar data space for the attribute.
Rank=1
Dimsf(:)=0 !???
Dimsf(1)=nVal
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, DataSpace, iError)
! Create the attribute for group Loc_ID.
IF(PRESENT(RealScalar)) Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(RealArray))  Type_ID=H5T_NATIVE_DOUBLE
IF(PRESENT(IntScalar))  Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(IntArray))   Type_ID=H5T_NATIVE_INTEGER
IF(PRESENT(LogicalScalar))THEN
  LogToInt=MERGE(1,0,LogicalScalar)
  Type_ID=H5T_NATIVE_INTEGER
END IF

! Create character string datatype for the attribute.
! For a attribute character, we have to build our own type with corresponding attribute length
IF(PRESENT(StrScalar))THEN
  AttrLen=LEN(StrScalar(1))
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, AttrLen, iError)
END IF
IF(PRESENT(StrArray))THEN
  AttrLen=LEN(StrArray(1))
  CALL H5TCOPY_F(H5T_NATIVE_CHARACTER, Type_ID, iError)
  CALL H5TSET_SIZE_F(Type_ID, AttrLen, iError)
ENDIF

CALL H5ACREATE_F(Loc_ID, TRIM(AttribName), Type_ID, DataSpace, Attr_ID, iError)
! Write the attribute data.
IF(PRESENT(RealArray))     buf=C_LOC(RealArray)
IF(PRESENT(RealScalar))    buf=C_LOC(RealScalar)
IF(PRESENT(IntArray))      buf=C_LOC(IntArray)
IF(PRESENT(IntScalar))     buf=C_LOC(IntScalar)
IF(PRESENT(LogicalScalar)) buf=C_LOC(LogToInt)
IF(PRESENT(StrScalar))     buf=C_LOC(StrScalar(1))
IF(PRESENT(StrArray))      buf=C_LOC(StrArray(1))
CALL H5AWRITE_F(Attr_ID, Type_ID, buf, iError)

! Close datatype
IF(PRESENT(StrScalar).OR.PRESENT(StrArray)) CALL H5TCLOSE_F(Type_ID, iError)
! Close dataspace
CALL H5SCLOSE_F(DataSpace, iError)
! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(Loc_ID.NE.Loc_ID_in)THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(*,*)'...DONE!'
END SUBROUTINE WriteAttribute


END MODULE MOD_Output_HDF5
