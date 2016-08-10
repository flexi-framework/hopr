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
MODULE MOD_Readin_GMSH
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tBC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadGMSH
  MODULE PROCEDURE ReadGMSH
END INTERFACE

PUBLIC::ReadGMSH
!===================================================================================================================================

TYPE tBCTemp
  INTEGER               :: nodeInds(4) !tri+quad (lazy)
  INTEGER               :: curveIndex
  TYPE(tBC),POINTER     :: BC
  TYPE(tBCTemp),POINTER :: nextBC
END TYPE

TYPE tBCTempPtr
  TYPE(tBCTemp),POINTER :: bp
END TYPE

INTEGER                 :: BCTagsFound(1337)

CONTAINS

FUNCTION GETNNODES(ElementType,bOrd)
!===================================================================================================================================
! Get nNodes from Element Type 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: ElementType  ! ?
INTEGER, INTENT(IN)                :: bOrd  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                            :: GETNNODES  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE(ElementType)
  CASE(3)
    GETNNODES=3
  CASE(6)
    GETNNODES=bOrd*(bOrd+1)/2
  CASE(4,5)
    GETNNODES=4
  CASE(7)
    GETNNODES=bOrd*bOrd
  CASE(104,204)
    GETNNODES=4
  CASE(105,115,205)
    GETNNODES=5
  CASE(106,116,206)
    GETNNODES=6
  CASE(108,118,208)
    GETNNODES=8
  END SELECT
END FUNCTION GETNNODES

SUBROUTINE readGMSH()
!===================================================================================================================================
! Read mesh from GMSH ascii or binary file. Called by fillMesh.
! Read-in can be performed by just one or all processors
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tElemPtr,tSide,tSidePtr,tNode,tNodePtr
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:nMeshFiles,MeshFileName
USE MOD_Mesh_Vars,ONLY:nUserDefinedBoundaries,BoundaryName
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewNode,getNewBC
USE MOD_Readin_GMSH_Vars
USE MOD_ReadinTools,ONLY:TRYREAD
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! nMeshFiles                : Number of mesh files (INI-File)
! MeshFileName(iFile)       : Filename of mesh file iFile (INI-File)
! nZones                    : Number of mesh zones (INI-File)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tNodePtr),POINTER :: Nodes(:)  ! ?
TYPE(tElemPtr),POINTER :: Elems(:)  ! ?
TYPE(tElem),POINTER    :: aElem  ! ?
TYPE(tSide),POINTER    :: aSide  ! ?
TYPE(tBCTemp),POINTER  :: aBCTemp  ! ?
TYPE(tBCTempPtr),POINTER :: BCList(:)  ! ?
INTEGER                :: os,iFile  ! ?
INTEGER                :: i,iElem,nElems,iNode,nNodes  ! ?
INTEGER                :: elemCount  ! ?
INTEGER                :: elemType,nTags  ! ?
INTEGER                :: tags(1337),nodeInds(1337)  ! ?
INTEGER                :: junk1  ! ?
INTEGER                :: minInd,tempNodeInds(4) ! three are enough
INTEGER                :: iBC,whichDim,BCInd
LOGICAL                :: isBCSide,BCFound(nUserDefinedBoundaries),found,s  ! ?
CHARACTER(LEN=255)     :: BCName
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')'Reading GMSH mesh...'
CALL buildTypes()
BCTagsFound=-1

! start reading mesh
DO iFile=1,nMeshFiles
  OPEN(UNIT   = 104,                  &
       FILE   = MeshFileName(iFile),  &
       STATUS = 'OLD',                &
       ACTION = 'READ',               &
       ACCESS = 'SEQUENTIAL',         &
       IOSTAT = os                    )
  WRITE(UNIT_stdOut,*)  'Reading mesh from ascii file: ',TRIM(MeshFileName(iFile))
  ! Read Version and other infos
  s=TRYREAD(104,'$MeshFormat')
  READ(104,*) ! format
  s=TRYREAD(104,'$EndMeshFormat')
  s=TRYREAD(104,'$PhysicalNames')
  READ(104,*) nBCs_GMSH
  ALLOCATE(MapBC(nBCs_GMSH))
  ALLOCATE(MapBCInd(nBCs_GMSH))
  MapBC=-1
  BCFound=.FALSE.
  DO iBC=1,nBCs_GMSH
    READ(104,*) whichDim, BCInd, BCName
    MapBCInd(iBC)=BCind
    IF(whichDim.EQ.2)THEN
      found=.FALSE.
      DO i=1,nUserDefinedBoundaries
        IF(INDEX(TRIM(BCName),TRIM(BoundaryName(i))).NE.0) THEN
          found=.TRUE. 
          BCFound(i)=.TRUE.
          MapBC(iBC)=i
          WRITE(*,*)'BC found: ',TRIM(BCName)
          WRITE(*,*)'              -->  mapped to:',TRIM(BoundaryName(i))
          EXIT
        END IF
      END DO
      IF(.NOT.found) CALL abort(__STAMP__, &
                    'UserDefinedBoundary condition missing: '//TRIM(BCName),nUserDefinedBoundaries)
    END IF
  END DO

  s=TRYREAD(104,'$EndPhysicalNames')
  s=TRYREAD(104,'$Nodes')
  READ(104,*) nNodes

  ! Read node coordinates
  ALLOCATE(Nodes(nNodes))
  DO iNode = 1, nNodes   ! READ nodes
    CALL GetNewNode(Nodes(iNode)%np,0)
    READ(104,*)Nodes(iNode)%np%ind,Nodes(iNode)%np%x !assume ordered nodes (1,2,3..)
    IF (iNode.NE.Nodes(iNode)%np%ind)&
      CALL abort(__STAMP__,&
                  'List of nodes in GMSH file has to be sorted (1,2,3,...)')
  END DO

  ALLOCATE(BCList(nNodes))
  DO i=1,nNodes
    NULLIFY(BCList(i)%bp)
  END DO

  s=TRYREAD(104,'$EndNodes')
  s=TRYREAD(104,'$Elements')
  READ(104,*) nElems !=vertices + lines + faces + elements (we want only elements)
  ALLOCATE(Elems(nElems))
  ! Read element connectivity
  elemCount=0
  DO iElem=1,nElems
    READ(104,*)junk1,elemType,nTags,tags(1:nTags),nodeInds(1:GMSH_TYPES(3,elemType))
    !IF (GMSH_TYPES(elemType,5) .NE. 1) & !not complete
    !  CALL abort(__STAMP__, &
    !  'The specified mesh contains element types, which are currently not supported',999,999.)
    SELECT CASE(GMSH_TYPES(6,elemType)) !2d=bc or 3d=element
    CASE(2)
      CALL addToBCs(BCList,elemType,nodeInds,nTags,tags)
    CASE(3)
      CALL buildElem(Elems(elemCount+1),elemCount,elemType,Nodes,nodeInds)
    END SELECT
  END DO ! nElems
  
  ! Assign Boundary Conditions
  DO iElem=1,elemCount
    aSide=>Elems(iElem)%ep%firstSide
    DO WHILE(ASSOCIATED(aSide))
      !get minimum index of side
      tempNodeInds=HUGE(1337)
      DO i=1,aSide%nNodes
        tempNodeInds(i)=aSide%node(i)%np%ind
      END DO
      minInd=MINVAL(tempNodeInds)
      ! get BC by minInd and check if primary nodes identical
      aBCTemp=>BCList(minInd)%bp
      isBCSide=.FALSE.
      DO WHILE(ASSOCIATED(aBCTemp))
        isBCSide=.TRUE.
        DO i=1,aSide%nNodes
          IF(COUNT(aBCTemp%nodeInds.EQ.tempNodeInds(i)).NE.1)THEN
            isBCSide=.FALSE.
            EXIT
          END IF
        END DO
        IF(isBCSide) EXIT
        aBCTemp=>aBCTemp%nextBC
      END DO
      IF(isBCSide)THEN
        aSide%curveIndex=aBCTemp%curveIndex
        aSide%BC=>aBCTemp%BC
      END IF
      aSide=>aSide%nextElemSide
    END DO
  END DO
  
  ! Build Pointer list from elements
  aElem=>Elems(1)%ep
  firstElem=>aElem
  DO iElem=2,elemCount
    aElem%nextElem=>Elems(iElem)%ep
    Elems(iElem)%ep%prevElem=>aElem
    aElem=>Elems(iElem)%ep
  END DO
  CLOSE(104)
! remember to delete all stuff (BCList etc.
  DO iNode=1,nNodes !throw away  
    IF (Nodes(iNode)%np%refCount.EQ.0) DEALLOCATE(Nodes(iNode)%np)
    DO WHILE(ASSOCIATED(BCList(iNode)%bp))
      aBCTemp=>BCList(iNode)%bp
      BCList(iNode)%bp=>aBCTemp%nextBC
      DEALLOCATE(aBCTemp)
    END DO
  END DO
  DEALLOCATE(Nodes,Elems,BCList)

  DO i=1,nUserDefinedBoundaries
    IF(.NOT.BCFound(i)) CALL abort(__STAMP__,&
            'One or more userdefined boundary conditions specified in ini file has not been found.',999,999.)
  END DO
END DO
WRITE(UNIT_stdOut,'(132("~"))')
END SUBROUTINE readGMSH

SUBROUTINE addToBCs(BCList,gmshElemType,nodeInds,nTags,tags)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:BoundaryType
USE MOD_Mesh_Vars,ONLY:getNewElem,tBC
USE MOD_Readin_GMSH_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)         :: gmshElemType  ! ?
INTEGER,INTENT(IN)         :: nodeInds(GMSH_TYPES(1,gmshElemType)) ! only take primary nodes of side
INTEGER,INTENT(IN)         :: nTags  ! ?
INTEGER,INTENT(IN)         :: tags(nTags)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tBCTempPtr),POINTER,INTENT(INOUT) :: BCList(:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                    :: i,iBC,minInd  ! ?
TYPE(tBCTemp),POINTER      :: aBCTemp  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------

IF(GMSH_TYPES(6,gmshElemType).NE.2) RETURN ! filter out volumes
IF(GMSH_TYPES(1,gmshElemType).LE.2) RETURN ! filter out lines 
! GMSH Tags: 1=physical group (aka BoundaryCondition), 2=geometric entitry(edge,face,vol) element belongs to, 3=mesh partition
! element belongs to, 4+=partition ids (negative partition = ghost cells)

iBC=-1
IF(tags(1).LE.MAXVAL(MapBCInd,1).AND.tags(1).GE.MINVAL(MapBCInd,1))THEN
  DO i=1,nBCs_GMSH
    IF(MapBCInd(i).EQ.tags(1))THEN
      iBC=MapBC(i)
      EXIT
    END IF
  END DO
END IF
IF(iBC.EQ.-1) CALL abort(__STAMP__,&
                         'Side with undefined BoundaryCondition found')

ALLOCATE(aBCTemp)
DO i=1,GMSH_TYPES(1,gmshElemType) !primary nodes
  aBCTemp%nodeInds(i)=nodeInds(i)
END DO
ALLOCATE(aBCTemp%BC)
! BCType
aBCTemp%BC%BCType    =BoundaryType(iBC,1)
! curveIndex
aBCTemp%curveIndex   =BoundaryType(iBC,2)  
! BCState
aBCTemp%BC%BCState   =BoundaryType(iBC,3)
! BCAlphaInd
aBCTemp%BC%BCAlphaInd=BoundaryType(iBC,4)
! BCInd
aBCTemp%BC%BCIndex   =iBC

minInd=MINVAL(nodeInds)
IF(ASSOCIATED(BCList(minInd)%bp))THEN
  aBCTemp%nextBC => BCList(minInd)%bp
  BCList(minInd)%bp => aBCTemp
ELSE
  BCList(minInd)%bp => aBCTemp
  NULLIFY(aBCTemp%nextBC)
END IF
END SUBROUTINE

SUBROUTINE buildElem(elem,elemCount,gmshElemType,Nodes,nodeInds)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Basis_Vars,ONLY:TetraMap,PyraMap,PrismMap,HexaMap
USE MOD_Mesh_Vars,ONLY:tElem,tElemPtr,tSide,tNode,tNodePtr
USE MOD_Mesh_Vars,ONLY:N
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewBC
USE MOD_Mesh_Vars,ONLY:useCurveds,rebuildCurveds
USE MOD_Mesh_Basis,ONLY:createSides
USE MOD_Readin_GMSH_Vars,ONLY:bOrd,getGMSHVolumeMapping,GMSH_TYPES
USE MOD_Readin_GMSH_Vars,ONLY:tetMapGMSH,pyrMapGMSH,priMapGMSH,hexMapGMSH
USE MOD_Readin_GMSH_Vars,ONLY:tetMapCGNSToGMSH,pyrMapCGNSToGMSH,priMapCGNSToGMSH,hexMapCGNSToGMSH
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)         :: gmshElemType  ! ?
INTEGER,INTENT(IN)         :: nodeInds(GMSH_TYPES(3,gmshElemType))  ! ?
TYPE(tNodePtr),POINTER,INTENT(IN)  :: Nodes(:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tElemPtr),INTENT(OUT) :: elem  ! ?
INTEGER, INTENT(INOUT)     :: elemCount  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                    :: i  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------

IF(GMSH_TYPES(6,gmshElemType).NE.3) RETURN ! no 3d element
IF (bOrd .EQ.0) THEN
  bOrd = GMSH_TYPES(4,gmshElemType)+1
  IF ((bOrd .NE. N+1).AND.useCurveds.AND..NOT.rebuildCurveds) & 
    CALL abort(__STAMP__,&
    'Mesh boundary order not equal to boundary order from ini file! Mesh order: ',N+1)
  CALL getGMSHVolumeMapping()
ELSE
  IF (bOrd.NE.GMSH_TYPES(4,gmshElemType)+1) &
    CALL abort(__STAMP__,&
    'All elements in the mesh are required to have the same boundary order!')
END IF

elemCount = elemCount + 1 ! element is valid, raise number of mesh elements by one
CALL getNewElem(elem%ep)
elem%ep%ind    = elemCount
elem%ep%Type   = GMSH_TYPES(2,gmshElemType)
elem%ep%nNodes = GETNNODES(GMSH_TYPES(2,gmshElemType),bOrd)
ALLOCATE(Elem%ep%Node(elem%ep%nNodes))

DO i=1,elem%ep%nNodes
  SELECT CASE(elem%ep%nNodes)
  CASE(4)
    elem%ep%node(i)%np => Nodes(nodeInds(TetMapCGNSToGMSH(i)))%np
  CASE(5)
    elem%ep%node(i)%np => Nodes(nodeInds(PyrMapCGNSToGMSH(i)))%np
  CASE(6)
    elem%ep%node(i)%np => Nodes(nodeInds(PriMapCGNSToGMSH(i)))%np
  CASE(8)
    elem%ep%node(i)%np => Nodes(nodeInds(HexMapCGNSToGMSH(i)))%np
  CASE DEFAULT
    STOP 'Unknown element type!'
  END SELECT
  elem%ep%node(i)%np%refCount = elem%ep%node(i)%np%refCount+1
END DO
CALL createSides(elem%ep,.TRUE.)

! assign curved elements if present, if curveds should be used and not rebuilt using our methods
IF(useCurveds .AND. (bOrd.GT.2) .AND.(.NOT.rebuildCurveds))THEN
  ALLOCATE(elem%ep%curvedNode(GMSH_TYPES(3,gmshElemType)))
  elem%ep%nCurvedNodes=GMSH_TYPES(3,gmshElemType)
  DO i=1,GMSH_TYPES(3,gmshElemType)
    SELECT CASE(elem%ep%nNodes)
    CASE(4)
      elem%ep%curvedNode(i)%np => Nodes(nodeInds(tetMapGMSH(tetraMap(i,1),tetraMap(i,2),tetraMap(i,3))))%np
    CASE(5)
      STOP 'High order pyramids not implemented yet for GMSH!'
      elem%ep%curvedNode(i)%np => Nodes(nodeInds(pyrMapGMSH(pyraMap(i,1),pyraMap(i,2),pyraMap(i,3))))%np
    CASE(6)
      STOP 'High order prisms not implemented yet for GMSH!'
      elem%ep%curvedNode(i)%np => Nodes(nodeInds(priMapGMSH(prismMap(i,1),prismMap(i,2),prismMap(i,3))))%np
    CASE(8)
      elem%ep%curvedNode(i)%np => Nodes(nodeInds(hexMapGMSH(hexaMap(i,1),hexaMap(i,2),hexaMap(i,3))))%np
    END SELECT
    !print*, tetraMap(i,1),tetraMap(i,2),tetraMap(i,3)
    !print*, elem%ep%curvedNode(i)%np%x
    !print*, tetMapGMSH(tetraMap(i,1),tetraMap(i,2),tetraMap(i,3))
    elem%ep%curvedNode(i)%np%refCount=elem%ep%curvedNode(i)%np%refCount+1
  END DO
END IF
END SUBROUTINE buildElem


END MODULE MOD_Readin_GMSH
