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
MODULE MOD_Readin_Gambit
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
INTERFACE ReadGambit
  MODULE PROCEDURE ReadGambit
END INTERFACE

PUBLIC::ReadGambit
!===================================================================================================================================

CONTAINS
SUBROUTINE ReadGambit()
!===================================================================================================================================
! Read mesh from gambit ascii or binary file. Called by fillMesh.
! Read-in can be performed by just one or all processors
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tElemPtr,tSide,tNode
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:MeshDim,n2dNodes,useBinary
USE MOD_Mesh_Vars,ONLY:nMeshFiles,MeshFileName
USE MOD_Mesh_Vars,ONLY:nUserDefinedBoundaries,BoundaryName,BoundaryType
USE MOD_Mesh_Basis,ONLY:createSides
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewNode,getNewBC
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! nMeshFiles                : Number of mesh files (INI-File)
! MeshFileName(iFile)           : Filename of mesh file iFile (INI-File)
! nZones                    : Number of mesh zones (INI-File)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tNode),POINTER    :: actualNode   ! ?
TYPE(tElem),POINTER    :: actualElem  ! ?
TYPE(tSide),POINTER    :: aSide   ! ?
TYPE(tElemPtr),POINTER :: elems(:)   ! ?
REAL,ALLOCATABLE       :: nodeCoords(:,:)  ! ?
INTEGER                :: i,nnodes,nElems,iFile,iGlZone,nElNodes  ! ?
INTEGER                :: iZone,nSet,dummy1,dummy2,iNode,iElem,dummy3,iSet,nZones  ! ?
INTEGER                :: BCtype,BCstate,curveIndex,BCalphaInd,nBCElems,iLine,nLine,nRest,os  ! ?
INTEGER                :: iDummyArray(10),strlen  ! ?
INTEGER                :: BCind  ! ?
INTEGER,ALLOCATABLE    :: iDummyArray2(:,:),iDummyArray3(:)  ! ?
CHARACTER(LEN=99)      :: formstr,BinaryFile  ! ?
CHARACTER(LEN=100)     :: cdummy,cdummy2   ! ?
CHARACTER(LEN=255)     :: strBC  ! ?
LOGICAL                :: foundBC  ! ?
INTEGER                :: HexNodeMap(8)  = (/3,2,7,6,4,1,8,5/)  ! ?
INTEGER                :: HexSideMap(6)  = (/1,2,6,4,3,5/)  ! ?
INTEGER                :: PyraNodeMap(5) = (/1,2,4,3,5/)  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')'Reading Gambit mesh...'

iGlZone=0
IF(meshDim .EQ. 2) n2dNodes=0  ! 2.5D mesh
! start reading mesh
DO iFile=1,nMeshFiles
  ! Open Gambit file
  IF(useBinary)THEN
    strlen=LEN_TRIM(MeshFileName(iFile))-3                                   ! Strip "neu"
    WRITE(formstr,'(a2,i2,a3)')'(a',strlen,',a)'
    WRITE(BinaryFile,formstr)MeshFileName(iFile),'halo'
    OPEN(UNIT   = 104,           &
         FILE   = BinaryFile,    &
         ACTION = 'READ',        &
         FORM   = 'UNFORMATTED', &
         IOSTAT = os             )
    WRITE(UNIT_stdOut,*)  'Reading mesh from binary file: ',TRIM(BinaryFile)
  ELSE
    OPEN(UNIT   = 104,                  &
         FILE   = MeshFileName(iFile), &
         STATUS = 'OLD',                &
         ACTION = 'READ',               &
         ACCESS = 'SEQUENTIAL',         &
         IOSTAT = os                    )
    WRITE(UNIT_stdOut,*)  'Reading mesh from ascii file: ',TRIM(MeshFileName(iFile))
  END IF ! useBinary
  ! Get number of nodes, elements, zone and BCs
  IF(useBinary)THEN
    READ(104)nNodes,nElems,nZones,nSet
  ELSE
    DO i = 1, 6                                                               ! Header weglesen
      READ(104,*)
    END DO
    READ(104,*)nNodes,nElems,nZones,nSet,dummy1,dummy2
  END IF ! useBinary
  IF(meshDim .EQ. 2) n2dNodes=n2dNodes+nNodes  ! 2.5D mesh
  ! Read node coordinates
  ALLOCATE(nodeCoords(3,nNodes))
  nodeCoords=0.
  IF(useBinary)THEN
    READ(104)nodeCoords(1:meshDim,:)
  ELSE
    DO i = 1, 2                                                               ! Header weglesen
      READ(104,*)
    END DO
    DO iNode = 1, nNodes   ! READ nodes
      READ(104,*)dummy1,nodeCoords(1:meshDim,iNode)
    END DO
    DO i = 1, 2                                                               ! Header weglesen
      READ(104,*)                                                             
    END DO
  END IF ! useBinary
  ! Read elements
  actualElem=>firstElem
  IF(ASSOCIATED(actualElem))THEN
    DO WHILE(ASSOCIATED(actualElem%nextElem))
      actualElem%Ind=-1
      actualElem=>actualElem%nextElem
    END DO
    actualElem%Ind=-1
  END IF
  ! Read element connectivity
  ALLOCATE(iDummyArray2(nElems,10),iDummyArray3(nElems))
  IF(useBinary)THEN 
    READ(104)iDummyArray2
    READ(104)iDummyArray3
  ELSE
    DO iElem=1,nElems
      READ(104,'(I8,I3,I3)', advance='no') dummy1,dummy2,nElNodes
      READ(104,*)iDummyArray2(iElem,1:nElNodes)
      iDummyArray3(iElem)=nElNodes
    END DO ! nElems
  END IF ! useBinary
  ! Build elements
  DO iElem=1,nElems
    IF(.NOT. ASSOCIATED(firstElem)) THEN
      CALL getNewElem(actualElem)
      firstElem=>actualElem
    ELSE
      CALL getNewElem(actualElem%nextElem)
      actualElem%nextElem%prevElem=>actualElem
      actualElem=>actualElem%nextElem
    END IF
    actualElem%nNodes=iDummyArray3(iElem)
    !CALL createSides(actualElem,iDummyArray3(iElem))
    ALLOCATE(actualElem%Node(actualElem%nnodes))
    DO iNode=1,actualElem%nnodes
      CALL getNewNode(actualNode,1)
      actualNode%x=nodeCoords(:,iDummyArray2(iElem,iNode))
      actualNode%ind=iDummyArray2(iElem,iNode)
      ! for some elements, gambit ordering is INCOMPATIBLE to standard (!!)
      SELECT CASE(actualElem%nnodes)
      CASE(5) ! Pyra
        actualElem%Node(PyraNodeMap(iNode))%np=>actualNode
      CASE(8) ! Hex
        actualElem%Node(HexNodeMap(iNode))%np=>actualNode
      CASE DEFAULT
        actualElem%Node(iNode)%np=>actualNode
      END SELECT
    END DO
    CALL createSides(actualElem,.TRUE.)
    actualElem%ind=iElem
  END DO
  DEALLOCATE(nodeCoords)
  DEALLOCATE(iDummyArray2,iDummyArray3)
  ALLOCATE(elems(nElems))
  DO iElem=1,nElems
    NULLIFY(elems(iElem)%ep)
  END DO
  actualElem=>firstElem
  DO WHILE(ASSOCIATED(actualElem))
    IF (actualElem%ind .GT.0) THEN 
      elems(actualElem%ind)%ep=>actualElem 
    END IF
    actualElem=>actualElem%nextElem
  END DO
  ! READ zones ("groups" in gambit)
  DO iZone = 1, nZones
    iGlZone=iGlZone+1
    IF (iGlZone .GT. nZones) &
      CALL abort(__STAMP__, &
       'inconsistent numer of Zones in readGambit',999,999.)
    IF(useBinary)THEN
      READ(104) dummy1
    ELSE
      DO i = 1, 2
        READ(104,*)                                                 ! Lese Header weg
      END DO
      READ(104,'(a28,i10,a40)') cdummy, dummy1, cdummy2
      DO i = 1, 2
        READ(104,*)                                                 ! Lese Header weg
      END DO
    END IF ! useBinary
    nLine  = dummy1/10
    nRest  = MOD(dummy1,10)
    ! Read zone elements for zone iZone
    DO iLine = 1, nLine
      READ(104,*) iDummyArray(1:10)
      DO i=1,10
        IF(ASSOCIATED(elems(idummyarray(i))%ep)) THEN
          elems(idummyarray(i))%ep%zone=iGlZone
        END IF
     END DO
    END DO
    ! If there is a rest, read last line
    IF(nRest.GT.0) THEN
      READ(104,*) iDummyArray(1:nRest)
      DO i=1,nRest
        IF(ASSOCIATED(elems(idummyarray(i))%ep)) THEN
          elems(idummyarray(i))%ep%zone=iGlZone
        END IF
      END DO
    END IF
  END DO !iZone
  ! READ Boundary Conditions ("sets" in gambit)
  DO iSet=1,nSet
    IF(useBinary)THEN
      READ(104)strBC
    ELSE
      DO i = 1, 2                                                                        ! Header weglesen
        READ(104,*)
      END DO
      READ(104,'(A255)') strBC                                                         ! Read line    
    END IF ! useBinary
    foundBC=.FALSE.
    DO i=1,nUserDefinedBoundaries
      dummy1=INDEX(TRIM(strBC),TRIM(BoundaryName(i)))
      IF(dummy1.NE.0) THEN
        foundBC=.TRUE. 
        BCType     = BoundaryType(i,1)
        curveIndex = BoundaryType(i,2)
        BCState    = BoundaryType(i,3)
        BCalphaInd = BoundaryType(i,4)
        BCind      = i
        WRITE(*,*)'BC found: ',TRIM(strBC)
        WRITE(*,*)'              -->  mapped to:',TRIM(BoundaryName(i))
        strBC=strBC(dummy1+LEN(TRIM(BoundaryName(i))):LEN(strBC))  ! First we need to cut off the boundary name...
        READ(strBC,'(I8,I8,I8,I8)')dummy1,nBCElems,dummy2,dummy3  ! ...before we can read the number of BC elements
        EXIT
      END IF
    END DO
    IF(.NOT.foundBC) CALL abort(__STAMP__, &
                      'UserDefinedBoundary condition missing: '//TRIM(strBC),nUserDefinedBoundaries,999.)
    ALLOCATE(iDummyArray2(nBCElems,3))
    IF(useBinary)THEN
      READ(104)iDummyArray2
    ELSE
      DO iElem = 1, nBCElems
        READ(104,*) iDummyArray2(iElem,1),iDummyArray2(iElem,2),iDummyArray2(iElem,3)  !ElementID, ??, Side ID (1..nsides)
      END DO
    END IF ! useBinary
    DO iElem = 1, nBCElems
      dummy1=iDummyArray2(iElem,1)
      dummy2=iDummyArray2(iElem,2)
      dummy3=iDummyArray2(iElem,3)
      IF(ASSOCIATED(elems(dummy1)%ep)) THEN
        actualElem=>elems(dummy1)%ep
        IF(meshDim .EQ. 2)THEN
          IF(actualElem%nNodes .EQ. 8)THEN  ! Quad -> Hexa
            SELECT CASE(dummy3)
            CASE(1)
              dummy3=2
            CASE(2)
              dummy3=3
            CASE(3)
              dummy3=4
            CASE(4)
              dummy3=5
            END SELECT
          END IF
        ELSE
          SELECT CASE(actualElem%nNodes)
          CASE(8)
            dummy3=HexSideMap(dummy3)
          END SELECT
        END IF
        aSide=>actualElem%firstSide
        DO i=1,dummy3-1
          aSide=>aSide%nextElemSide
        END DO
        CALL getNewBC(aSide%BC)
        aSide%BC%BCtype=BCtype !assign BC to side
        aSide%curveIndex=curveIndex
        aSide%BC%BCstate=BCstate
        aSide%BC%BCalphaInd=BCalphaInd
        aSide%BC%BCIndex=BCind
      END IF ! ASSOCIATED(elems(dummy1)%ep)
    END DO
    DEALLOCATE(iDummyArray2)
  END DO
  DEALLOCATE(elems)
  CLOSE(104)
END DO !iFile
CALL Timer(.FALSE.)
WRITE(UNIT_stdOut,'(132("~"))')
END SUBROUTINE readGambit

END MODULE MOD_Readin_Gambit
