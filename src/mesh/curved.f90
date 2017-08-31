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
MODULE MOD_Curved
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
INTERFACE reconstructNormals
  MODULE PROCEDURE reconstructNormals
END INTERFACE

INTERFACE deleteDuplicateNormals
  MODULE PROCEDURE deleteDuplicateNormals
END INTERFACE

INTERFACE getExactNormals
  MODULE PROCEDURE getExactNormals
END INTERFACE

INTERFACE readNormals
  MODULE PROCEDURE readNormals
END INTERFACE

INTERFACE checkNormals
  MODULE PROCEDURE checkNormals
END INTERFACE

INTERFACE getNewNormal
  MODULE PROCEDURE getNewNormal
END INTERFACE

INTERFACE ProjectToExactSurfaces
  MODULE PROCEDURE ProjectToExactSurfaces
END INTERFACE

INTERFACE create3DSplines
  MODULE PROCEDURE create3DSplines
END INTERFACE

INTERFACE buildCurvedElementsFromVolume
  MODULE PROCEDURE buildCurvedElementsFromVolume
END INTERFACE

INTERFACE buildCurvedElementsFromBoundarySides
  MODULE PROCEDURE buildCurvedElementsFromBoundarySides
END INTERFACE

INTERFACE curvedEdgesToSurf
  MODULE PROCEDURE curvedEdgesToSurf
END INTERFACE

INTERFACE curvedSurfacesToElem
  MODULE PROCEDURE curvedSurfacesToElem
END INTERFACE

INTERFACE splitToSpline
  MODULE PROCEDURE splitToSpline
END INTERFACE

INTERFACE RebuildMortarGeometry
  MODULE PROCEDURE RebuildMortarGeometry
END INTERFACE

PUBLIC::SplitToSpline,ReconstructNormals,getExactNormals,deleteDuplicateNormals,create3DSplines,curvEdedgesToSurf
PUBLIC::ProjectToExactSurfaces
PUBLIC::CurvedSurfacesToElem,readNormals,buildCurvedElementsFromVolume,buildCurvedElementsFromBoundarySides
PUBLIC::RebuildMortarGeometry
!===================================================================================================================================

CONTAINS


SUBROUTINE readNormals()
!===================================================================================================================================
! Read normal vectors from normal vector ascii file. Called by fillMesh.
! Read-in can be performed by just one processors
!-----------------------------------------------------------------------------------------------------------------------------------
! length of subroutine ~ xxx lines !
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tNode,tNormal,tNormalPtr
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:NormalVectFile
USE MOD_Mesh_Tolerances,ONLY:SAMEPOINT
USE MOD_Search
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
! normals on a Meshnode
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i,j,k,l,os,countAssignedNodes  ! ?
INTEGER                   :: Counter,nFaces,face(2)  ! ?
INTEGER                   :: nVects,NodeID,nTotalNodes,nPlanars,searchMeshNodes  ! ?
INTEGER,POINTER           :: tempFaceIDs(:)  ! ?
INTEGER,ALLOCATABLE       :: faceConnectivity(:,:),mergedFaces(:),DuplicatesInd(:)  ! ?
REAL                      :: actNodeCoords(3)  ! ?
REAL,POINTER              :: tempNormals(:,:)  ! ?
LOGICAL                   :: pointFound,useProjections,changeIndex,foundOnce  ! ?
LOGICAL,ALLOCATABLE       :: sameFaces(:,:)  ! ?
TYPE(tToObject),POINTER   :: ToObject !,nextToObject,searchToObject,NewToObject
TYPE(tSearchMesh),POINTER :: searchMesh  ! ?
TYPE(tNode),POINTER       :: aNode  ! ?
TYPE(tSide),POINTER       :: aSide  ! ?
TYPE(tElem),POINTER       :: aElem  ! ?
TYPE(tNormal),POINTER     :: aNormal,bNormal,firstNormal  ! ?
TYPE(tNormalPtr),POINTER  :: aNormals(:)  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'READING NORMAL VECTORS FROM FILE...',NormalVectFile
CALL Timer(.TRUE.)
OPEN(UNIT   = 150,                  &
   FILE   = NormalVectFile, &
   STATUS = 'OLD',                &
   ACTION = 'READ',               &
   ACCESS = 'SEQUENTIAL',         &
   IOSTAT = os                    )
IF(os.NE.0) THEN  ! File Error
  CALL abort(__STAMP__, &
         'ERROR: cannot open normal vector file: '//TRIM(NormalVectFile),999,999.)
END IF
LOGWRITE(UNIT_stdOut,*) 'Reading normals from ascii file: ',TRIM(NormalVectFile)
DO i = 1, 4                                                               ! Header weglesen
  READ(150,*)
END DO
READ(150,*)j
IF (j .EQ. 0) THEN
  useProjections=.FALSE.
ELSE
  useProjections=.TRUE.
END IF
! READ Face Connectivity
READ(150,*)
READ(150,*)nFaces
READ(150,*)
ALLOCATE(faceConnectivity(nFaces,nFaces),mergedFaces(nFaces),sameFaces(nFaces,nFaces))
faceConnectivity(:,:)=0
!pos numbers first in array (1,2,3..), negative numbers are last (..-3,-2,-1)
nPlanars=0 !planar faces have negative FaceIDs, curveds have positive
DO i=1, nFaces
  READ(150,'(I3)',ADVANCE='NO') j !number of entries in line
  READ(150,'(I6)',ADVANCE='NO') k !first FaceID
  IF (k .LT. 0) THEN
    nPlanars=nPlanars+1
    faceConnectivity(1,nFaces+1+k)=k   !actual Face ID  
    READ(150,*) faceConnectivity(2:j,nFaces+1+k) !Neighbor Face IDs of actual Face  
  ELSEIF (k .EQ. 0) THEN
    CALL abort(__STAMP__, &
         'ERROR: CAD-face connectivity invalid',999,999.)
  ELSE
    faceConnectivity(1,k)=k    
    READ(150,*) faceConnectivity(2:j,k)    
  END IF
END DO
DO i=1,nFaces !all entries have to be filled
  IF (faceConnectivity(1,i) .EQ. 0) THEN
    CALL abort(__STAMP__, &
         'ERROR: CAD-face connectivity invalid',999,999.)
  END IF
END DO

!use node%tmp as marker for unchecked elements
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO j=1, aElem%nNodes
    aElem%node(j)%np%tmp=0
  END DO
  aElem=>aElem%nextElem
END DO
!build local search mesh
searchmeshNodes=0
NULLIFY(searchMesh)
CALL getNewSearchMesh(searchMesh,.TRUE.)
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF(aSide%curveIndex.GT.0) THEN
      DO j=1,aSide%nNodes
        aNode=>aSide%node(j)%np
        IF(aNode%tmp .EQ. 0) THEN
          aNode%tmp=999 
          CALL insertNode(searchMesh,aNode)
          aNode%refCount=aNode%refCount-1
          searchMeshNodes=searchMeshNodes+1
        END IF
      END DO
    END IF
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO

READ(150,*)
READ(150,*)  !NNODES
READ(150,*)nTotalNodes
WRITE(UNIT_stdOut,*) 'Total number of curved surface nodes: ',nTotalNodes
READ(150,*)  !     NORMAL VECTORS
countAssignedNodes=0
!read Normalvectors in an Normalvector pointer list
DO i=1, nTotalNodes
  READ(150,*) NodeID,nVects,actNodeCoords
  !search if node is inside local search mesh
  ToObject=>getFirstToObject(searchMesh,.FALSE.,actNodeCoords)
  pointFound=.FALSE.
  DO WHILE (ASSOCIATED(ToObject))
    IF (SAMEPOINT(ToObject%Node%x,actNodeCoords,PP_MeshTolerance)) THEN
      pointFound=.TRUE.
      IF (useProjections) THEN
        READ(150,*) ToObject%Node%x !project the node to exact surface
      END IF
      ALLOCATE(tempNormals(3,nVects), tempFaceIDs(nVects))
      tempFaceIDs=0
      DO j = 1, nVects
        READ(150,*) tempFaceIDs(j),tempNormals(:,j)
      END DO
      !build normalvector pointer list
      NULLIFY(firstNormal)
      DO j=1, nVects
        IF(.NOT.ASSOCIATED(firstNormal))THEN
          CALL getNewNormal(aNormal, 1)
          firstNormal=>aNormal
        ELSE
          CALL getNewNormal(aNormal%nextNormal,1)
          aNormal%nextNormal%prevNormal=>aNormal
          aNormal=>aNormal%nextNormal
        END IF
        aNormal%normal(:)=tempNormals(:,j)
        aNormal%FaceID(:)=tempFaceIDs(j)
      END DO
      ToObject%Node%firstNormal=>firstNormal !write to node
      countAssignedNodes=countAssignedNodes+1
      ! delete ToObject from search mesh
      CALL deleteToObject(searchmesh%sm(searchMesh%actualInd(1),&
                                        searchMesh%actualInd(2),&
                                        searchMesh%actualInd(3))%ToObject,ToObject)
      DEALLOCATE(tempNormals, tempFaceIDs)
      EXIT
    END IF
    ToObject=>getNextToObject(searchMesh,.FALSE.)
  END DO
  IF (.NOT. pointFound) THEN
    IF (useProjections) READ(150,*) 
    DO j = 1, nVects
      READ(150,*) 
    END DO
  END IF
END DO !nTotalNodes
CLOSE(150)
IF (countAssignedNodes .NE. searchMeshNodes) THEN
  CALL abort(__STAMP__, &
           'ERROR: not all normals are assigned to local nodes in searchmesh',countAssignedNodes,REAL(searchMeshNodes))
END IF
CALL deleteDuplicateNormals
CALL checkNormals
!create list of merged faces
DO i=1,(nFaces-nPlanars)
  mergedFaces(i)=i
END DO
DO i=1,nPlanars
  mergedFaces(nFaces+1-i)=-i
END DO
sameFaces(:,:)=.FALSE.
DO i=1,nFaces
  sameFaces(i,i)=.TRUE.
END DO
ALLOCATE(aNormals(2))
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF(aSide%curveIndex.GT.0) THEN
      DO j=1,aSide%nNodes
        ! find an edge with only one normal at each node and which connects two different FaceIDs
        k=MOD(j,aSide%nNodes)+1
        aNormals(1)%np=>aSide%node(j)%np%firstNormal
        aNormals(2)%np=>aSide%node(k)%np%firstNormal
        IF (aNormals(1)%np%FaceID(1) .EQ. aNormals(2)%np%FaceID(1)) CYCLE !same faceID
        IF ((ASSOCIATED(aNormals(1)%np%nextNormal)) .OR. &
            (ASSOCIATED(aNormals(2)%np%nextNormal))) CYCLE !more than one normal at node
        IF ((SIZE(aNormals(1)%np%FaceID) .GT. 1) .OR. &
            (SIZE(aNormals(2)%np%FaceID) .GT. 1)) CYCLE !more than one face ID at node
        DO i=1,2
          IF (aNormals(i)%np%FaceID(1) .LT. 0) THEN
            face(i)=nFaces+1+aNormals(i)%np%FaceID(1)
          ELSE
            face(i)=aNormals(i)%np%FaceID(1)
          END IF
        END DO
        ! mark same Faces
        IF (sameFaces(face(1),face(2))) CYCLE            
        sameFaces(face(1),face(2))=.TRUE.
        sameFaces(face(2),face(1))=.TRUE.
        ! look if already another face is the same Face
        DO i=1,nFaces
          IF (sameFaces(face(1),i)) THEN
            sameFaces(i,face(2))=.TRUE.
            sameFaces(face(2),i)=.TRUE.
          END IF       
          IF (sameFaces(face(2),i)) THEN
            sameFaces(i,face(1))=.TRUE.
            sameFaces(face(1),i)=.TRUE.
          END IF
        END DO       
      END DO
    END IF
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO
!look for other sameFaces
DO i=1,nFaces
  DO j=1,nFaces-1
    IF (sameFaces(i,j)) THEN
      DO k=j+1,nFaces
        IF (sameFaces(i,k)) THEN
          sameFaces(j,k)=.TRUE.
          sameFaces(k,j)=.TRUE.
        END IF
      END DO
    END IF
  END DO
END DO
! determine for a group of sameFaces the mapping from the actual faceID to the biggest FaceID in the group
mergedFaces(:)=0
DO i=1,nFaces !
  DO j=(nFaces-nPlanars),1,-1 !curved Faces, backwards, so j is the biggest FaceID 
    IF (sameFaces(i,j)) THEN
      mergedFaces(i)=j
      EXIT
    END IF
  END DO
END DO
! planar faces are after the curved faces, but have a negative faceID
DO i=(nFaces-nPlanars),nFaces  ! planar Faces
  IF (mergedFaces(i) .NE. 0) CYCLE
  DO j=nFaces,(nFaces-nPlanars),-1
    IF (sameFaces(i,j)) THEN
      mergedFaces(i)=j-nFaces-1
      EXIT
    END IF
  END DO
END DO

!use node%tmp as marker for unchecked elements
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO j=1, aElem%nNodes
    aElem%node(j)%np%tmp=0
  END DO
  aElem=>aElem%nextElem
END DO
!apply new FaceIDs to nodes
i=0
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO j=1, aElem%nNodes
    IF (aElem%node(j)%np%tmp .NE. 0) CYCLE
    i=i+1
    aElem%node(j)%np%tmp=i       
    aNormal=>aElem%node(j)%np%firstNormal
    DO WHILE (ASSOCIATED(aNormal))
      DO k=1, SIZE(aNormal%FaceID)
        IF (aNormal%FaceID(k) .LT. 0) THEN !planar face
          aNormal%FaceID(k)=mergedFaces(nFaces+1+aNormal%FaceID(k)) !replace old FaceId with new one 
        ELSE
          aNormal%FaceID(k)=mergedFaces(aNormal%FaceID(k)) !replace old FaceId with new one
        END IF
      END DO
      aNormal=>aNormal%nextNormal
    END DO
  END DO
  aElem=>aElem%nextElem
END DO
!remove duplicate FaceIDs and merge normals
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO j=1, aElem%nNodes
    aElem%node(j)%np%tmp=0
  END DO
  aElem=>aElem%nextElem
END DO
i=0
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO j=1, aElem%nNodes
    IF (aElem%node(j)%np%tmp .NE. 0) CYCLE
    i=i+1
    aElem%node(j)%np%tmp=i      
    !remove duplicates inside each normal
    aNormal=>aElem%node(j)%np%firstNormal
    DO WHILE (ASSOCIATED(aNormal))
      Counter=0
      DO k=1,SIZE(aNormal%FaceID)-1
        IF(aNormal%FaceID(k).NE.0)THEN
          DO l=k+1, SIZE(aNormal%FaceID)
            IF (aNormal%FaceID(k) .EQ. aNormal%FaceID(l)) THEN
              aNormal%FaceID(l)=0
              Counter=Counter+1
            END IF
          END DO
        END IF
      END DO
      IF (Counter .GT. 0) THEN !duplicates found, reorder FaceID array
        l=0
        ALLOCATE(DuplicatesInd(SIZE(aNormal%FaceID)-Counter))
        DO k=1, (SIZE(aNormal%FaceID)-Counter)
          IF (aNormal%FaceID(k) .EQ. 0) THEN
            CYCLE
          ELSE
            l=l+1
            DuplicatesInd(l)=aNormal%FaceID(k)
          END IF
        END DO
        DEALLOCATE(aNormal%FaceID) 
        ALLOCATE(aNormal%FaceID(SIZE(duplicatesInd)))
        aNormal%FaceID=DuplicatesInd
        DEALLOCATE(duplicatesInd)
      END IF
      aNormal=>aNormal%nextNormal
    END DO
    !merge normals if duplicates are found in different normals
    aNormal=>aElem%node(j)%np%firstNormal
    DO WHILE (ASSOCIATED(aNormal))
      foundOnce=.FALSE.
      bNormal=>aNormal%nextNormal
      DO WHILE(ASSOCIATED(bNormal))
        Counter=0
        DO k=1,SIZE(aNormal%FaceID) !use aNormal as reference, set FaceIDs to 0 if found in bNormal
          DO l=1, SIZE(bNormal%FaceID)
            IF (aNormal%FaceID(k) .EQ. bNormal%FaceID(l)) THEN
              bNormal%FaceID(l)=0
              Counter=Counter+1
            END IF
          END DO
        END DO
        IF (Counter .GT. 0) THEN !duplicate FaceIDs found, bNormal is removed
          foundOnce=.TRUE.
          !copy FaceIDs of aNormal and bNormal into one array
          ALLOCATE(DuplicatesInd(SIZE(aNormal%FaceID)+SIZE(bNormal%FaceID)-Counter))
          DuplicatesInd(1:SIZE(aNormal%FaceID))=aNormal%FaceID(:)
          l=SIZE(aNormal%FaceID)
          DO k=1, SIZE(bNormal%FaceID) 
            IF (bNormal%FaceID(k) .EQ. 0) THEN
              CYCLE
            ELSE
              l=l+1
              DuplicatesInd(l)=bNormal%FaceID(k)
            END IF
          END DO
          DEALLOCATE(aNormal%FaceID) 
          ALLOCATE(aNormal%FaceID(SIZE(duplicatesInd)))
          aNormal%FaceID=DuplicatesInd
          DEALLOCATE(duplicatesInd)
          
          !average normals, remove bNormal and reorder pointers
          aNormal%normal=aNormal%normal+bNormal%normal
          IF (ASSOCIATED(bNormal%nextNormal)) THEN
            aNormal%nextNormal=>bNormal%nextNormal
            aNormal%nextNormal%prevNormal=>aNormal
            NULLIFY(bNormal%nextNormal)
          ELSE
            NULLIFY(aNormal%nextNormal)
          END IF
          NULLIFY(bNormal%prevNormal)
          DEALLOCATE(bNormal)
          bNormal=>aNormal
        END IF
        bNormal=>bNormal%nextNormal
      END DO
      IF (foundOnce) aNormal%normal(:)=aNormal%normal(:)/SQRT(SUM(aNormal%normal(:)*aNormal%normal(:)))
      aNormal=>aNormal%nextNormal
    END DO
  END DO
  aElem=>aElem%nextElem
END DO

! If each node of a side has same negative FaceID => CAD face is not curved => Side ist not curved! 
! => set aSide%curveIndex = 0 => no spline reconstruction => better performance (in work)
l=0
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    changeIndex=.FALSE.
    IF(aSide%curveIndex.GT.0) THEN
      aNormal=>aSide%node(1)%np%firstNormal
outer:DO WHILE(ASSOCIATED(aNormal)) !search neg FaceID in first node, check all normals. If no neg FaceIDs, side is curved!
        DO i=1,SIZE(aNormal%FaceID)
          IF (aNormal%FaceID(i) .GT. 0) CYCLE !loop till negative FaceID is found in aNormal. 
          DO j=2,aSide%nNodes !neg FaceID found, search in other nodes for that ID
            changeIndex=.FALSE.
            bNormal=>aSide%node(j)%np%firstNormal
            DO WHILE(ASSOCIATED(bNormal))
              DO k=1,SIZE(bNormal%FaceID)
                IF (bNormal%FaceID(k) .EQ. aNormal%FaceID(i)) changeIndex=.TRUE.
              END DO
              IF (changeIndex) EXIT
              bNormal=>bNormal%nextNormal
            END DO
            IF (.NOT. changeIndex) EXIT !one node does not have neg FaceID, don't search in other nodes
          END DO
          IF (changeIndex) EXIT outer !all nodes have same neg FaceID, side is not curved => don't continue with first Node
        END DO 
        aNormal=>aNormal%nextNormal
      END DO outer
      IF (changeIndex) THEN
        aSide%isCurved=.FALSE.
        aSide%curveIndex=0
        l=l+1
      END IF
    END IF !curveInd > 0
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO

CALL deleteSearchMesh(searchMesh)
DEALLOCATE(mergedFaces,sameFaces,faceConnectivity,aNormals)
LOGWRITE(UNIT_stdOut,*)'All normals have been assigned successfully!'
CALL Timer(.FALSE.)
END SUBROUTINE readNormals

SUBROUTINE checkNormals()
!===================================================================================================================================
! Count all normals
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tNormal
USE MOD_Mesh_Vars,ONLY:FirstElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER                :: nFaceIDs(5)      ! max Number of duplicate normals
INTEGER                :: globalnFaceIDs(5)      ! max Number of duplicate normals
INTEGER                :: nFaceID,j  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tNormal),POINTER  :: aNormal         ! New normal vectors object
TYPE(tElem),POINTER    :: aElem  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
!use nodeInd as marker for unchecked elements
nFaceIDs=0
globalnFaceIDs=0
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO j=1, aElem%nNodes
    aElem%node(j)%np%tmp=-1
  END DO
  aElem=>aElem%nextElem
END DO
!apply new FaceIDs to nodes
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO j=1, aElem%nNodes
    IF (aElem%node(j)%np%tmp .NE. -1) CYCLE
    aElem%node(j)%np%tmp=1       
    aNormal=>aElem%node(j)%np%firstNormal
    DO WHILE (ASSOCIATED(aNormal))
        nFaceID=min(5,SIZE(aNormal%FaceID))
        IF(nFaceID.LT.1) STOP 'nFaceID<1'
        nFaceIDs(nFaceID) = nFaceIDs(nFaceID)+1
      aNormal=>aNormal%nextNormal
    END DO
  END DO
  aElem=>aElem%nextElem
END DO
globalnFaceIDs=nFaceIDs
WRITE(*,*)'number of normals ...',globalnFaceIDs
END SUBROUTINE checkNormals

SUBROUTINE getNewNormal(aNormal,nFaceIDs)
!===================================================================================================================================
! Allocate and initialize new normalvector "aNormal"
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNormal
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                :: nFaceIDs      ! max Number of duplicate normals
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tNormal),POINTER,INTENT(OUT)  :: aNormal       ! New normal vectors object
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
ALLOCATE(aNormal)
ALLOCATE(aNormal%FaceID(nFaceIDs))
NULLIFY(aNormal%prevNormal)
NULLIFY(aNormal%nextNormal)
END SUBROUTINE getNewNormal


SUBROUTINE reconstructNormals()
!===================================================================================================================================
! create Normals on curved boundaries from the actual mesh, so just approximated normals! 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tNode,tNormal
USE MOD_Mesh_Vars,ONLY:FirstElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tSide),POINTER       :: aSide  ! ?
TYPE(tElem),POINTER       :: aElem  ! ?
TYPE(tNode),POINTER       :: Node  ! ?
TYPE(tNormal),POINTER     :: nv,nv2,nv3
INTEGER                   :: iNode,nn  ! ?
INTEGER                   :: prev1(4,3:4),next1(4,3:4)
REAL                      :: v1(3),v2(3)  ! ?
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'RECONSTRUCT NORMALS ... '

next1(:,3)=(/2,3,1,0/)
prev1(:,3)=(/3,1,2,0/)
next1(:,4)=(/2,3,4,1/)
prev1(:,4)=(/4,1,2,3/)

! compute the normals
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF (aSide%curveIndex.LE.0)THEN
      aSide=>aSide%nextElemSide
      CYCLE
    END IF
    nn=aSide%nNodes
    DO iNode=1,nn
      node=>aSide%Node(iNode)%np
      node%tmp=0
      v1=aSide%Node(next1(iNode,nn))%np%x-node%x
      v2=aSide%Node(prev1(iNode,nn))%np%x-node%x
      NULLIFY(nv)
      ALLOCATE(nv)
      nv%nextNormal=>node%firstNormal
      IF(ASSOCIATED(nv%nextNormal)) nv%nextNormal%prevNormal=>nv
      node%firstNormal=>nv
      ALLOCATE(nv%FaceID(1))
      nv%normal=cross(v1,v2)
      nv%FaceID(1)=aSide%curveIndex
    END DO
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO

! average normals by curveindex
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF (aSide%curveIndex.LE.0)THEN
      aSide=>aSide%nextElemSide
      CYCLE
    END IF
    nn=aSide%nNodes
    DO iNode=1,nn
      node=>aSide%Node(iNode)%np
      IF(node%tmp.NE.0)CYCLE
      node%tmp=999
      ! loop over all normals connected to a node and average them, deallocate averaged normals
      nv=>node%firstNormal
      DO WHILE(ASSOCIATED(nv))
        nv2=>nv%nextNormal
        DO WHILE(ASSOCIATED(nv2))
          IF(nv%FaceID(1).EQ.nv2%FaceID(1))THEN
            nv%normal=nv%normal+nv2%normal
            nv2%prevNormal%nextNormal=>nv2%nextNormal
            IF(ASSOCIATED(nv2%nextNormal)) nv2%nextNormal%prevNormal=>nv2%prevNormal
            nv3=>nv2
            nv2=>nv2%nextNormal
            DEALLOCATE(nv3)
            CYCLE
          END IF
          nv2=>nv2%nextNormal
        END DO
        nv%normal=nv%normal/NORM2(nv%normal)
        nv=>nv%nextNormal
      END DO
    END DO
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO

CALL Timer(.FALSE.)
END SUBROUTINE reconstructNormals


SUBROUTINE deleteDuplicateNormals()
!===================================================================================================================================
! ? 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tNode,tNormal
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:minNormalAngle
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
! normals on a Meshnode
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tNode),POINTER       :: aNode  ! ?
TYPE(tElem),POINTER       :: aElem  ! ?
TYPE(tSide),POINTER       :: aSide  ! ?
TYPE(tNormal),POINTER     :: aNormal,bNormal,cNormal  ! ?
INTEGER                   :: NVects,nFaceIDs  ! ?
REAL                      :: tempNormal(3),normalAngle  ! ?
INTEGER                   :: j  ! ?
INTEGER,ALLOCATABLE       :: TempFaceIDs(:)  ! ?
!===================================================================================================================================
!use nodeInd as marker for unchecked elements
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO j=1, aElem%nNodes
    aElem%node(j)%np%tmp=0
  END DO
  aElem=>aElem%nextElem
END DO
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF(aSide%curveIndex.GT.0) THEN
      DO j=1,aSide%nNodes
        aNode=>aSide%node(j)%np
        IF(aNode%tmp.EQ.-1)CYCLE
        aNode%tmp=-1
        aNormal=>aNode%firstNormal
        !count normals
        nVects=0
        DO WHILE(ASSOCIATED(aNormal))
          nVects=nVects+1
          aNormal=>aNormal%nextNormal
        END DO
        IF(nVects.GT.1)THEN
          !eliminate all double normals
          ALLOCATE(tempFaceIDs(nVects))
          aNormal=>aNode%firstNormal
          DO WHILE(ASSOCIATED(aNormal))
            tempFaceIDs=0
            tempFaceIDs(1)=aNormal%FaceID(1)
            tempNormal=aNormal%Normal
            nFaceIDs=1
            bNormal=>aNormal%nextNormal
            DO WHILE(ASSOCIATED(bNormal))
              normalAngle=SUM(aNormal%Normal*bNormal%Normal)
              IF (normalAngle .GT. minNormalAngle) THEN
                cNormal=>bNormal
                bNormal=>bNormal%prevNormal
                nFaceIDs=nFaceIDs+1
                tempFaceIDs(nFaceIDs)=cNormal%FaceID(1)
                tempNormal=tempNormal+cNormal%Normal
                !delete duplicate bNormal
                IF(ASSOCIATED(cNormal%nextNormal))THEN
                  cNormal%prevNormal%nextNormal=>cNormal%nextNormal
                  cNormal%nextNormal%prevNormal=>cNormal%prevNormal
                ELSE
                  NULLIFY(cNormal%prevNormal%nextNormal)
                END IF
                DEALLOCATE(cNormal)
              END IF
              bNormal=>bNormal%nextNormal
            END DO !bNormal
            IF(nFaceIDs.GT.1) THEN !duplicates!
              aNormal%Normal=tempNormal/SQRT(SUM(tempNormal*tempNormal))
              DEALLOCATE(aNormal%FaceID)
              ALLOCATE(aNormal%FaceID(nFaceIDs))
              aNormal%FaceID=tempFaceIDs(1:nFaceIDs)
            END IF
            aNormal=>aNormal%nextNormal
          END DO
          DEALLOCATE(tempFaceIDs)
        END IF !double Normals 
      END DO
    END IF !curveInd > 0
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO
END SUBROUTINE deleteDuplicateNormals

SUBROUTINE buildCurvedElementsFromVolume()
!===================================================================================================================================
! set surface curvednode pointers from existing curved volume 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,FirstElem,N
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER          :: Elem  ! ?
TYPE(tSide),POINTER          :: Side  ! ?
!===================================================================================================================================
IF(N .LE. 1) RETURN
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  ! check if element is really curved or only bilinear
  CALL curvedVolToSurf(Elem,onlyBoundarySides=.FALSE.)
  Side => Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    CALL curvedSurfToEdges(Side)
    Side => Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

END SUBROUTINE buildCurvedElementsFromVolume

SUBROUTINE buildCurvedElementsFromBoundarySides(nCurvedBoundaryLayers)
!===================================================================================================================================
! set surface curvednode pointers from existing curved volume 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,FirstElem,N,deleteNode
USE MOD_Basis_Vars,ONLY:MapSideToVol
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)           :: nCurvedBoundaryLayers
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER          :: Elem  ! ?
TYPE(tSide),POINTER          :: Side  ! ?
INTEGER                      :: iNode,iLayer  ! ?
LOGICAL                      :: onlyBL ! ?
!===================================================================================================================================
IF(N .LE. 1 .OR. nCurvedBoundaryLayers .LT. 0) RETURN

! mark layers, start with first layer at BC where either side or whole element can be curved
iLayer=MERGE(0,1,           nCurvedBoundaryLayers.EQ.0)
onlyBL=MERGE(.TRUE.,.FALSE.,nCurvedBoundaryLayers.EQ.0)
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Elem%tmp=-999
  DO iNode=1,Elem%nCurvedNodes
    Elem%CurvedNode(iNode)%np%tmp=-999
  END DO
  Elem=>Elem%nextElem
END DO !Elem


!mark nodes on curved surfac, iLayer=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side => Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(ASSOCIATED(Side%BC).AND.Side%curveIndex.GT.0)THEN
      DO iNode=1,MERGE((N+1)**2,(N+1)*(N+2)/2,(Side%nNodes.EQ.4))
        Elem%curvedNode(MapSideToVol(iNode,Side%locSide,Elem%nNodes))%np%tmp=0
      END DO !iNode
      Side%tmp=iLayer
      Elem%tmp=iLayer
      !!!Elem%zone=-iLayer !!! visualization hack
    ELSE
      Side%tmp=-999
    END IF
    Side => Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO !Elem
! now mark elements layer by layer
DO iLayer=1,nCurvedBoundaryLayers
  Elem=>firstElem
  DO WHILE(ASSOCIATED(Elem))
    IF(Elem%tmp.EQ.-999)THEN
      DO iNode=1,Elem%nCurvedNodes
        IF(Elem%CurvedNode(iNode)%np%tmp.EQ.iLayer-1) THEN
          Elem%tmp=iLayer
          !!!Elem%zone=-iLayer !!! visualization hack
          EXIT !do iNode loop
        END IF
      END DO !iNode
    END IF
    IF(Elem%tmp.EQ.iLayer)THEN
      !mark unmarked nodes
      DO iNode=1,Elem%nCurvedNodes
        IF(Elem%CurvedNode(iNode)%np%tmp.EQ.-999)THEN
          Elem%CurvedNode(iNode)%np%tmp=iLayer
        END IF
      END DO !iNode
    END IF
    Elem=>Elem%nextElem
  END DO !Elem
END DO !iLayer

!unmark nodes again
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  DO iNode=1,Elem%nCurvedNodes
    Elem%CurvedNode(iNode)%np%tmp=-999
  END DO
  DO iNode=1,Elem%nNodes
    Elem%Node(iNode)%np%tmp=999 ! never delete corner nodes
  END DO
  Elem=>Elem%nextElem
END DO !Elem

Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  IF(Elem%tmp.GE.0)   CALL curvedVolToSurf(Elem,onlyBoundarySides=onlyBL)
  Side => Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%tmp.GE.0 .OR. Elem%tmp.GT.0) CALL curvedSurfToEdges(Side)
    Side => Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

CALL curvedEdgesToSurf(keepExistingCurveds=.TRUE.)

! Mark all curved nodes which will not be deleted
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  IF(Elem%tmp.GT.0)THEN
    DO iNode=1,Elem%nCurvedNodes
      Elem%CurvedNode(iNode)%np%tmp=999
    END DO
  END IF
  Side => Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%tmp.GE.0)THEN
      DO iNode=1,Side%nNodes
        Side%Node(iNode)%np%tmp=999
      END DO
      DO iNode=1,Side%nCurvedNodes
        Side%CurvedNode(iNode)%np%tmp=999
      END DO
    END IF
    Side => Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  DO iNode=1,Elem%nCurvedNodes
    IF(Elem%CurvedNode(iNode)%np%tmp.NE.999)THEN
      CALL deleteNode(Elem%CurvedNode(iNode)%np)
    END IF
  END DO
  DEALLOCATE(Elem%CurvedNode)
  Elem%nCurvedNodes=0
  Elem=>Elem%nextElem
END DO
CALL curvedSurfacesToElem()

END SUBROUTINE buildCurvedElementsFromBoundarySides

SUBROUTINE curvedVolToSurf(Elem,onlyBoundarySides)
!===================================================================================================================================
! Curved elements are only created adjacent to boundary sides with a curveIndex > 0.
! Here only the surface information is taken for curving the elements.
! All other elements remain linear (except for elements adjacent to curved edges).
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tNodePtr
USE MOD_Mesh_Vars,ONLY:N
USE MOD_Basis_Vars,ONLY:MapSideToVol,TriaMapInv,QuadMapInv
!USE MOD_Basis_Vars,ONLY:TetraMapInv,PyraMapInv,PrismMapInv,HexaMapInv
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(IN)          :: Elem  ! ?
LOGICAL,INTENT(IN)                      :: onlyBoundarySides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tSide),POINTER          :: Side  ! ?
TYPE(tNodePtr),POINTER       :: refSide(:,:)  ! ?
INTEGER                      :: p,q  ! ?
!===================================================================================================================================
ALLOCATE(refSide(0:N,0:N))
DO q=0,N; DO p=0,N
  NULLIFY(refSide(p,q)%np)
END DO; END DO

Side=>Elem%firstSide
DO WHILE(ASSOCIATED(Side))
  IF(onlyBoundarySides)THEN
    IF(.NOT.ASSOCIATED(Side%BC).OR.Side%CurveIndex.LE.0)THEN
      Side=>Side%nextElemSide
      CYCLE
    END IF
  END IF
  ! Volume to reference side
  IF(Side%nNodes.EQ.4)THEN
    DO q=0,N
      DO p=0,N
        refSide(p,q)%np => Elem%curvedNode(MapSideToVol(QuadMapInv(p,q),Side%locSide,Elem%nNodes))%np
      END DO !p
    END DO !q
  ELSE
    DO q=0,N
      DO p=0,N-q
        refSide(p,q)%np => Elem%curvedNode(MapSideToVol(TriaMapInv(p,q),Side%locSide,Elem%nNodes))%np
      END DO !p
    END DO !q
  END IF !Side%nNodes == 4
  CALL referenceSideToFlipped(refSide,Side)
  Side=>Side%nextElemSide
END DO !iSide
DEALLOCATE(refSide)
END SUBROUTINE curvedVolToSurf


SUBROUTINE referenceSideToFlipped(refSide,Side)
!===================================================================================================================================
! set surface curvednode pointers from existing curved volume 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tSide,tNodePtr,N
USE MOD_Mesh_Basis,ONLY:getFlip
USE MOD_Basis_Vars,ONLY:TriaMapInv,QuadMapInv
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide),POINTER,INTENT(INOUT) :: Side  ! ?
TYPE(tNodePtr),POINTER,INTENT(IN) :: refSide(:,:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: flip,p,q  ! ?
!===================================================================================================================================
flip=getFlip(Side)

! Reference side to real side
SELECT CASE (Side%nNodes)
CASE(3)
  Side%nCurvedNodes=((N+1)*(N+2))/2
  ALLOCATE(Side%curvedNode(Side%nCurvedNodes))
  DO q=0,N 
    DO p=0,N-q
      SELECT CASE(flip)
      CASE(0)
        Side%curvedNode(TriaMapInv(p,q))%np => refSide(p,q)%np
      CASE(1)
        Side%curvedNode(TriaMapInv(p,q))%np => refSide(q,p)%np
      CASE(2)
        Side%curvedNode(TriaMapInv(p,q))%np => refSide(N-q-p,q)%np
      CASE(3)
        Side%curvedNode(TriaMapInv(p,q))%np => refSide(p,N-q-p)%np
      END SELECT  
      Side%curvedNode(TriaMapInv(p,q))%np%refCount=Side%curvedNode(TriaMapInv(p,q))%np%refCount+1
    END DO !p
  END DO !q
CASE(4)
  Side%nCurvedNodes=((N+1)**2)
  ALLOCATE(Side%curvedNode(Side%nCurvedNodes))
  DO q=0,N 
    DO p=0,N
      SELECT CASE(flip)
      CASE(0)
        Side%curvedNode(QuadMapInv(p,q))%np => refSide(p,q)%np
      CASE(1)
        Side%curvedNode(QuadMapInv(p,q))%np => refSide(q,p)%np
      CASE(2)
        Side%curvedNode(QuadMapInv(p,q))%np => refSide(N-p,q)%np
      CASE(3)
        Side%curvedNode(QuadMapInv(p,q))%np => refSide(N-q,N-p)%np
      CASE(4)
        Side%curvedNode(QuadMapInv(p,q))%np => refSide(p,N-q)%np
      END SELECT
      Side%curvedNode(QuadMapInv(p,q))%np%refCount=Side%curvedNode(QuadMapInv(p,q))%np%refCount+1
    END DO !p
  END DO !q
CASE DEFAULT
  CALL abort(__STAMP__,'Only triangular and quadrangular sides are supported.')
END SELECT
END SUBROUTINE referenceSideToFlipped


FUNCTION ElemIsCurved(Elem)
!===================================================================================================================================
! Check if element with curved nodes allocated is actually curved or only bi/trilinear.
! Get plane for all element sides and get distance to plane. If dist > tol elem is curved.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(IN) :: Elem  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                        :: elemIsCurved  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tSide),POINTER            :: Side  ! ?
REAL,DIMENSION(3)              :: v1,v2,n1  ! ?
REAL                           :: a,charLength  ! ?
INTEGER                        :: i  ! ?
!===================================================================================================================================
elemIsCurved=.FALSE.
Side=>Elem%firstSide
DO WHILE(ASSOCIATED(Side))
  v1=Side%node(2)%np%x-Side%node(1)%np%x
  v2=Side%node(3)%np%x-Side%node(1)%np%x
  n1=CROSS(v1,v2)
  n1=n1/SQRT(SUM(n1*n1)) !normal of linear side (by corner nodes)
  IF(SUM(v1*v1).GT.SUM(v2*v2))THEN !get characteristic length for cell
    charLength=SQRT(SUM(v1*v1))*0.000001
  ELSE
    charLength=SQRT(SUM(v2*v2))*0.000001
  END IF
  a=SUM(n1*Side%node(1)%np%x)
  DO i=1,Side%nCurvedNodes !dist curvednode/linear elem plane
    IF(ABS(SUM(n1*Side%curvedNode(i)%np%x)-a).GT.charLength)THEN
      elemIsCurved=.TRUE.
      RETURN
    END IF
  END DO
  Side=>Side%nextElemSide
END DO
END FUNCTION ElemIsCurved

SUBROUTINE curvedSurfToEdges(aSide)
!===================================================================================================================================
! set edge curvednode pointers from existing curved surface 
!===================================================================================================================================
! MODULES
USE MOD_Basis_Vars,ONLY:EdgeToTria,EdgeToQuad
USE MOD_Mesh_Vars,ONLY:tSide,tEdge
USE MOD_Mesh_Vars,ONLY:N
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide),POINTER,INTENT(IN)          :: aSide  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tEdge),POINTER          :: aEdge  ! ?
INTEGER                      :: j,iEdge,bOrd,nNodes  ! ?
INTEGER                      :: E2F(4,N+1)  ! ?
!===================================================================================================================================
bOrd=N+1
nNodes=aSide%nNodes
IF(nNodes.EQ.3)THEN
  E2F(1:3,:)=EdgeToTria
ELSE
  E2F=EdgeToQuad
END IF
DO iEdge=1,nNodes
  aEdge=>aSide%Edge(iEdge)%edp
  IF(.NOT.ASSOCIATED(aEdge%CurvedNode))THEN
    ALLOCATE(aEdge%CurvedNode(bOrd))
    aEdge%CurvedNode(1)%np=>aEdge%Node(1)%np
    aEdge%CurvedNode(bOrd)%np=>aEdge%Node(2)%np
    IF(aSide%EdgeOrientation(iEdge))THEN !oriented
      DO j=2,N
        aEdge%CurvedNode(j)%np=>aSide%CurvedNode(E2F(iEdge,j))%np
      END DO
    ELSE
      DO j=2,N
        aEdge%CurvedNode(bOrd+1-j)%np=>aSide%CurvedNode(E2F(iEdge,j))%np
      END DO
    END IF
  END IF
END DO !iEdge 
END SUBROUTINE curvedSurfToEdges


SUBROUTINE getExactNormals()
!===================================================================================================================================
! Evaluate Exact Normal Function and apply Normals to nodes
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tNode,tNormal
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:ExactNormals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: exactFunction,i  ! ?
TYPE(tElem),POINTER        :: aElem  ! ?
TYPE(tSide),POINTER        :: aSide  ! ?
TYPE(tNode),POINTER        :: aNode  ! ?
TYPE(tNormal),POINTER      :: aNormal  ! ?
LOGICAL                    :: bAlreadyThere  ! ?
!===================================================================================================================================
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    DO i=1,aSide%nNodes
      NULLIFY(aSide%node(i)%np%firstnormal)
    END DO
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO

aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF(aSide%curveIndex .GT. 0) THEN
      !find exact normal function assigned to curveindex
      exactFunction=0
      DO i=1, SIZE(exactNormals),2
        IF (aSide%CurveIndex .EQ. exactNormals(i)) THEN
          exactFunction = exactNormals(i+1)
          EXIT
        END IF
      END DO
      DO i=1, aSide%nNodes
        aNode=>aSide%node(i)%np
        aNormal=>aNode%firstNormal
        bAlreadyThere=.FALSE.
        DO WHILE(ASSOCIATED(aNormal)) !check if curveindex has already been evaluated
          IF (aNormal%FaceID(1) .EQ. aSide%CurveIndex) THEN
            bAlreadyThere= .TRUE.
            EXIT
          END IF
          aNormal=>aNormal%nextNormal
        END DO
        IF (.NOT. (bAlreadyThere)) THEN
          CALL getNewNormal(aNormal,1)
          CALL exactNormalFunction(aNode%x, exactFunction, aNormal%normal)
          aNormal%FaceID(1)=aSide%CurveIndex
          IF (ASSOCIATED(aNode%firstNormal)) THEN
            aNode%firstNormal%prevNormal=>aNormal
            aNormal%nextNormal=>aNode%firstNormal
          END IF
          aNode%firstNormal=>aNormal
        END IF
      END DO
    END IF
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO
END SUBROUTINE getExactNormals


SUBROUTINE exactNormalFunction(xpos,exactIndex,normal)
!===================================================================================================================================
! Calculate Normal from analytical geometry for CreateSplinePatches on a given point with xpos and CurveIndex
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)                       :: xpos(3)  ! ?
INTEGER,INTENT(IN)                    :: exactIndex  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)           :: normal(3)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(exactIndex)
CASE(0) !do nothing
CASE(1) !sphere, origin (0,0,0)
  normal(:)=xpos(:)
CASE(2) !cylinder, origin x=0.,y=0.
  normal(1)=xpos(1)
  normal(2)=xpos(2)
  normal(3)=0
CASE DEFAULT
  CALL abort(__STAMP__,  &
             'Exact normal function with specified index does not exist. Index:',exactIndex,999.)
END SELECT
normal(:)=normal(:)/SQRT(SUM(normal(:)*normal(:)))
END SUBROUTINE exactNormalFunction


SUBROUTINE ProjectToExactSurfaces()
!===================================================================================================================================
! Evaluate exact surface Function and change coordinates of surface nodes (+curved nodes)
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tNode,tNormal
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:nExactSurfFuncs,ExactSurfFunc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: exactFunction,i  ! ?
TYPE(tElem),POINTER        :: aElem  ! ?
TYPE(tSide),POINTER        :: aSide  ! ?
TYPE(tNode),POINTER        :: aNode  ! ?
REAL                       :: maxDist,xnew(3)  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'PROJECT SURFACE NODES TO EXACT SURF FUNCS... '
CALL Timer(.TRUE.)
!mark sides and nodes
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    aSide%tmp=0
    DO i=1,aSide%nNodes
      aSide%node(i)%np%tmp=0
    END DO
    DO i=1,aSide%nCurvedNodes
      aSide%curvedNode(i)%np%tmp=0
    END DO
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF(aSide%curveIndex .GT. 0) THEN
      exactFunction=0
      DO i=1, nExactSurfFuncs
        IF (aSide%CurveIndex .EQ. exactSurfFunc(2*(i-1)+1)) THEN
          exactFunction = exactSurfFunc(2*i)
          EXIT
        END IF
      END DO
      aSide%tmp=exactFunction 
      DO i=1,aSide%nNodes
        aSide%node(i)%np%tmp=1
      END DO
      DO i=1,aSide%nCurvedNodes
        aSide%curvedNode(i)%np%tmp=1
      END DO
    END IF
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO

maxDist=0.
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF(aSide%tmp .GT. 0) THEN
      DO i=1, aSide%nNodes
        aNode=>aSide%Node(i)%np
        IF (aNode%tmp .EQ.1) THEN
          CALL exactSurfaceFunction(aNode%x, aSide%tmp, xnew)
          maxDist=MAX(maxDist,SUM((xnew-aNode%x)**2))
          aNode%x=xnew
          aNode%tmp=0
        END IF
      END DO
      DO i=1, aSide%nCurvedNodes
        aNode=>aSide%CurvedNode(i)%np
        IF (aNode%tmp .EQ.1) THEN
          CALL exactSurfaceFunction(aNode%x, aSide%tmp, xnew)
          maxDist=MAX(maxDist,SUM((xnew-aNode%x)**2))
          aNode%x=xnew
          aNode%tmp=0
        END IF
      END DO
      aSide%tmp=0
    END IF
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO
maxDist=SQRT(maxDist)
WRITE(UNIT_stdOut,'(A,E11.4)')'maximum distance of projection points: ',maxDist
CALL Timer(.FALSE.)
END SUBROUTINE ProjectToExactSurfaces

SUBROUTINE exactSurfaceFunction(xold,exactFunction,xnew)
!===================================================================================================================================
! Calculate Normal from analytical geometry for CreateSplinePatches on a given point with xpos and CurveIndex
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)            :: xold(3)  ! ?
INTEGER,INTENT(IN)         :: exactFunction  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)           :: xnew(3)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: i,j  ! ?
REAL                       :: R,rloc  ! ?
REAL                       :: c(0:4),t,tt(8),xt,yt,dxt,dyt,ddyt,ft,dft,d  ! ?
!===================================================================================================================================
SELECT CASE(exactFunction)
CASE(0) !do nothing
CASE(1) !sphere, origin (0,0,0),radius 0.5
  R=0.5
  rloc=SQRT(SUM(xold(:)**2))
  xnew(:)=xold(:)*(R/rloc)
CASE(2) !cylinder, origin x=0.,y=0., radius=0.5
  R=0.5
  rloc=SQRT(SUM(xold(1:2)**2))
  xnew(1:2)=(xold(1:2))*R/rloc
  xnew(3)=xold(3)
CASE(11) !naca0012 profile,zero TE thickness,c=1,
         ! x=t^2,y=sign(yold)*0.12/0.2*(0.2969*t-0.1260*t**2-0.3516*t**4+0.2843*t**6-0.1036*t**8) 
  c(0:4)=(/0.2969,-0.1260,-0.3516,0.2843,-0.1036/)
  t=MIN(1.,SQRT(MAX(0.,xold(1)))) !start for newton iteration
  d=0.12/0.2
  !find t with minimum distance to xold (xt-x)^2+(yt-y)^2 -> min
  DO i=1,1000 
    tt(1)=t
    DO j=2,8
      tt(j)=tt(j-1)*t
    END DO
    xt  =tt(2)
    dxt =2.*t
    yt  =d*(c(0)*t +   c(1)*tt(2)+   c(2)*tt(4) +   c(3)*tt(6)+    c(4)*tt(8))
    dyt =d*(c(0)   +2.*c(1)*t    +4.*c(2)*tt(3) +6.*c(3)*tt(5)+ 8.*c(4)*tt(7))
    ddyt=d*(        2.*c(1)     +12.*c(2)*tt(2)+30.*c(3)*tt(4)+56.*c(4)*tt(6))
    ft=(yt-abs(xold(2)))*dyt+(xt-xold(1))*dxt
    IF(abs(ft).LT.1.0E-12)EXIT
    dft=(dyt**2+(yt-abs(xold(2)))*ddyt+dxt**2+(xt-xold(1))*2.)
    t=t-ft/dft
    t=MIN(1.,t)
  END DO
  xnew(1)=xt
  xnew(2)=yt*SIGN(1.,xold(2))
  xnew(3)=xold(3)
  IF(i.GE.1000) THEN 
    WRITE(UNIT_stdOut,'(A)')'  Warning: Newton iteration not converged in exactSurfaceFunction:'
    WRITE(UNIT_stdOut,'(A,3E21.11)')'    old position',xold
    WRITE(UNIT_stdOut,'(A,3E21.11)')'    new position',xnew
  END IF
CASE DEFAULT
  CALL abort(__STAMP__,  &
             'Exact surface function with specified index does not exist. Index:',exactFunction,999.)
END SELECT
END SUBROUTINE exactSurfaceFunction


SUBROUTINE create3DSplines()
!===================================================================================================================================
! Creates Spline Coefficients. 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tEdge
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:getNewNode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
     !  smoothInterface...can be used to smooth sharp corners and angles between two splines
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER       :: aElem  ! ?
TYPE(tSide),POINTER       :: aSide  ! ?
TYPE(tEdge),POINTER       :: aEdge  ! ?
INTEGER                   :: i  ! ?
REAL                      :: v(3,2)  ! ?
INTEGER                   :: normalCaseCount(2)  ! ?
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'CREATE CURVED EDGES FROM NORMALS ... '

normalCaseCount=0
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide =>aElem%FirstSide
  DO WHILE(ASSOCIATED(aSide))
    IF (aSide%curveIndex .EQ. 0) THEN
      aSide=>aSide%nextElemSide
      CYCLE
    END IF
    DO i=1,aSide%nNodes
      aEdge=>aSide%Edge(i)%edp
      IF(ASSOCIATED(aEdge%curvedNode))CYCLE
      ALLOCATE(aEdge%curvedNode(4))
      CALL getTangentialVectors(aSide,i,v,normalCaseCount) 
      aEdge%CurvedNode(1)%np=>aEdge%Node(1)%np
      aEdge%CurvedNode(4)%np=>aEdge%Node(2)%np
      CALL getNewNode(aEdge%CurvedNode(2)%np,1)
      CALL getNewNode(aEdge%CurvedNode(3)%np,1)
      IF(aSide%EdgeOrientation(i))THEN !edge oriented Node1->Node2
        aEdge%CurvedNode(2)%np%x=1./27.*(20.*aEdge%Node(1)%np%x+4.*v(:,1)+2.*v(:,2)+7.*aEdge%Node(2)%np%x)
        aEdge%CurvedNode(3)%np%x=1./27.*(7.*aEdge%Node(1)%np%x+2.*v(:,1)+4.*v(:,2)+20.*aEdge%Node(2)%np%x)
      ELSE
        aEdge%CurvedNode(2)%np%x=1./27.*(20.*aEdge%Node(1)%np%x+4.*v(:,2)+2.*v(:,1)+7.*aEdge%Node(2)%np%x)
        aEdge%CurvedNode(3)%np%x=1./27.*(7.*aEdge%Node(1)%np%x+2.*v(:,2)+4.*v(:,1)+20.*aEdge%Node(2)%np%x)
      END IF !Edge oriented
    END DO
    aSide%isCurved=.TRUE.
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO

IF (SUM(normalCaseCount) .GT. 0)  THEN
  ERRWRITE(UNIT_stdOut,*)'The normal assignment by FaceID has failed for some elements.'
  ERRWRITE(UNIT_stdOut,*)'Tangent calculation by projection:',normalCaseCount(1)
  ERRWRITE(UNIT_stdOut,*)'Tangent calculation by cross product:',normalCaseCount(2)
END IF

CALL Timer(.FALSE.)
END SUBROUTINE create3DSplines


SUBROUTINE getTangentialVectors(aSide,nodeInd,v,normalCaseCount)
!===================================================================================================================================
! Pick correct normal from normals in Node
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tSide,tNodePtr,tNormalPtr
USE MOD_Mesh_Basis,ONLY:isOriented
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide),POINTER,INTENT(IN)       :: aSide  ! ?
INTEGER,INTENT(IN)                   :: nodeInd  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                      :: v(3,2)  ! ?
INTEGER,INTENT(OUT)                   :: normalCaseCount(2)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i,j,c,n,iCounter(2),nNormals(2)  ! ?
INTEGER                   :: nextNodeInd,prevNodeInd  ! ?
REAL                      :: vTemp(3),approxNormal(3),sideVect(3,2)  ! ?
REAL                      :: AngleMin,AngleTmp,AngleLimit,length  ! ?
LOGICAL                   :: doSearch,tangentFound(2),bAlreadyFound(2),tryProjection  ! ?
TYPE(tNodePtr),POINTER    :: edgeNode(:)  ! ?
TYPE(tNormalPtr),POINTER  :: aNormal(:), foundNormals(:,:)  ! ?
!===================================================================================================================================
nextNodeInd=MOD(nodeInd,aSide%nNodes)+1
ALLOCATE(edgeNode(2),aNormal(2))
edgeNode(1)%np=>aSide%OrientedNode(nodeInd)%np
edgeNode(2)%np=>aSide%OrientedNode(nextNodeInd)%np
sideVect(:,1)=edgeNode(2)%np%x(:)-edgeNode(1)%np%x(:)
length=SQRT(SUM(sideVect(:,1)*sideVect(:,1)))
sideVect(:,2)=-1.*sideVect(:,1)
!special case: only one normal, calculate tangent directly without search
doSearch=.FALSE.
tangentFound(:)=.FALSE.
DO n=1,2
  aNormal(n)%np=>edgeNode(n)%np%firstNormal
  IF (.NOT. ASSOCIATED(aNormal(n)%np)) THEN
    ERRWRITE(*,'(132("~"))')
    ERRWRITE(*,*)'Warning: no normals exist for this node, using cell edge as tangent'
    ERRWRITE(*,*)'Coords:', edgeNode(n)%np%x(:)
    ERRWRITE(*,'(132("~"))')
    ERRWRITE(*,*)'No normals exist for node', edgeNode(n)%np%x(:)
    v(:,n) = sideVect(:,n)
    tangentFound(n)=.TRUE.
  ELSE
    IF (.NOT.(ASSOCIATED(aNormal(n)%np%nextNormal))) THEN
      v(:,n)=sideVect(:,n)-SUM(sideVect(:,n)*aNormal(n)%np%normal(:))*aNormal(n)%np%normal(:)
      tangentFound(n)=.TRUE.
    ELSE
      doSearch=.TRUE.
    END IF
  END IF
END DO
IF (.NOT. doSearch) THEN
  DEALLOCATE(edgeNode,aNormal)
  RETURN
END IF

!determine maxnormals
nNormals(:)=0
DO n=1,2
  aNormal(n)%np=>edgeNode(n)%np%firstNormal
  DO WHILE(ASSOCIATED(aNormal(n)%np))
    nNormals(n)=nNormals(n)+1
    aNormal(n)%np=>aNormal(n)%np%nextNormal
  END DO
END DO

c=MAX(nNormals(1), nNormals(2))
ALLOCATE(foundNormals(c,2))
! search normal pairs by FaceID
iCounter(:)=0
aNormal(1)%np=>edgeNode(1)%np%firstNormal
DO WHILE(ASSOCIATED(aNormal(1)%np))
  DO i=1, SIZE(aNormal(1)%np%FaceID,1)
    aNormal(2)%np=>edgeNode(2)%np%firstNormal
chk:DO WHILE(ASSOCIATED(aNormal(2)%np))
      DO j=1, SIZE(aNormal(2)%np%FaceID,1)
        IF (aNormal(1)%np%FaceID(i) .EQ. aNormal(2)%np%FaceID(j)) THEN
          bAlreadyFound(:)=.FALSE.
          DO n=1,2
            DO c=1, iCounter(n)
              IF (ASSOCIATED(foundNormals(c,n)%np, aNormal(n)%np)) THEN
                bAlreadyFound(n)=.TRUE.
              END IF
            END DO
            IF(.NOT.(bAlreadyFound(n))) THEN
              iCounter(n)=iCounter(n)+1
              foundNormals(iCounter(n),n)%np=>aNormal(n)%np
            END IF
          END DO
          EXIT chk
        END IF
      END DO
      aNormal(2)%np=>aNormal(2)%np%nextNormal
    END DO chk
  END DO
  aNormal(1)%np=>aNormal(1)%np%nextNormal
END DO

!calculate tangents with normals found in search
DO n=1,2
  IF (.NOT. tangentFound(n)) THEN
    SELECT CASE (iCounter(n))
      CASE (1) !one pair found: on face
        v(:,n)=sideVect(:,n)-SUM(sideVect(:,n)*foundNormals(1,n)%np%normal)*foundNormals(1,n)%np%normal
        tangentFound(n)=.TRUE.
      CASE (2) !two pairs found: on edge
        vTemp=cross(foundNormals(1,n)%np%normal,foundNormals(2,n)%np%normal)
        vTemp=vTemp/NORM2(vTemp)
        v(:,n)=vTemp(:)*SUM(vTemp*sideVect(:,n))

        tangentFound(n)=.TRUE.
    END SELECT
  END IF
END DO

!fallback: calculate normals if search did not succeed or more than two pairs have been found
DO n=1,2
  IF (.NOT. tangentFound(n)) THEN
    ERRWRITE(*,'(132("~"))')
    ERRWRITE(*,*)'Normal assignment by FaceID failed, using manual assignment'
    ERRWRITE(*,*)'Node1:nNormals', nNormals(1)
    ERRWRITE(*,*)'Node1:x', edgeNode(1)%np%x(:)
    aNormal(1)%np=>edgeNode(1)%np%firstNormal
    DO WHILE(ASSOCIATED(aNormal(1)%np))
      ERRWRITE(*,*) 'Node1:',  aNormal(1)%np%FaceID
      aNormal(1)%np=>aNormal(1)%np%nextNormal
    END DO
    ERRWRITE(*,*)'Node2:nNormals', nNormals(2)
    ERRWRITE(*,*)'Node2:x', edgeNode(2)%np%x(:)
    aNormal(2)%np=>edgeNode(2)%np%firstNormal
    DO WHILE(ASSOCIATED(aNormal(2)%np))
      ERRWRITE(*,*) 'Node2:', aNormal(2)%np%FaceID
      aNormal(2)%np=>aNormal(2)%np%nextNormal
    END DO
    ERRWRITE(*,'(132("~"))')
    
    prevNodeInd=nodeInd-1
    IF (prevNodeInd .EQ. 0) prevNodeInd=aSide%nNodes
    vTemp(:)=aSide%OrientedNode(prevNodeInd)%np%x(:)-aSide%OrientedNode(nodeInd)%np%x(:) !calculate normal of aSide
    approxNormal(:)=cross(sideVect(:,1),vTemp(:))
    approxNormal(:)=approxNormal(:)/SQRT(SUM(approxNormal(:)*approxNormal(:)))
    IF (.NOT. isOriented(aSide)) THEN
      approxNormal(:)=-1.*approxNormal(:)
    END IF
    IF (nNormals(MOD(n,2)+1) .GT. 1) THEN
      !if angle of crossp of normals is below limit then use as tangent
      IF (aSide%nNodes .EQ. 3) THEN
          AngleLimit = 0.99 !11.5 deg
      ELSE
          AngleLimit = 0.98 !30 deg
      END IF
      AngleMin=-1.
      NULLIFY(foundNormals(1,1)%np)
      aNormal(1)%np=>edgeNode(n)%np%firstNormal
      DO WHILE (ASSOCIATED(aNormal(1)%np%nextNormal))
       aNormal(2)%np=>aNormal(1)%np%nextNormal
       DO WHILE (ASSOCIATED(aNormal(2)%np))
         vTemp=cross(aNormal(1)%np%normal,aNormal(2)%np%normal)
         vTemp=vTemp/NORM2(vTemp)
         AngleTmp=ABS(SUM(vTemp(:)*approxNormal(:))) !angle to approxnormal not too big
         IF (AngleTmp .LE. 0.4) THEN 
           AngleTmp=ABS(SUM(vTemp(:)*sideVect(:,n)))/length ! result is cos of angle
           IF (AngleTmp .GT. AngleMin) THEN
             AngleMin=AngleTmp
             foundNormals(1,1)%np=>aNormal(1)%np
             foundNormals(1,2)%np=>aNormal(2)%np
           END IF
         END IF
         aNormal(2)%np=>aNormal(2)%np%nextNormal
       END DO
       aNormal(1)%np=>aNormal(1)%np%nextNormal
      END DO
      IF ((AngleMin .GT. AngleLimit) .AND. (ASSOCIATED(foundNormals(1,1)%np))) THEN !if small angle use cross product 
        normalCaseCount(2) = normalCaseCount(2)+1
        vTemp(:)=cross(foundNormals(1,1)%np%normal,foundNormals(1,2)%np%normal)
        vTemp(:)=vTemp(:)/SQRT(SUM(vTemp(:)*vTemp(:)))
        v(:,n)=vTemp(:)*SUM(vTemp(:)*sideVect(:,n))
        tryProjection=.FALSE.
      ELSE !if big angle use projection
        tryProjection=.TRUE. 
      END IF
    ELSE
      tryProjection=.TRUE.
    END IF
    
    IF (tryProjection) THEN
      normalCaseCount(1) = normalCaseCount(1)+1
      AngleMin=-1. !cos of angle 1=min, 0=max
      aNormal(n)%np=>edgeNode(n)%np%firstNormal
      DO WHILE(ASSOCIATED(aNormal(n)%np))
        AngleTmp=ABS(SUM(aNormal(n)%np%normal(:)*approxNormal(:)))
        IF (AngleTmp .GT. AngleMin) THEN
          AngleMin=AngleTmp
          foundNormals(1,n)%np=>aNormal(n)%np
        END IF
        aNormal(n)%np=>aNormal(n)%np%nextNormal
      END DO
      v(:,n)=sideVect(:,n)-SUM(sideVect(:,n)*foundNormals(1,n)%np%normal(:))*foundNormals(1,n)%np%normal(:)
    END IF
  END IF
END DO
DEALLOCATE(edgeNode,aNormal,foundNormals)
END SUBROUTINE getTangentialVectors


SUBROUTINE curvedEdgesToSurf(keepExistingCurveds)
!===================================================================================================================================
! Blend edges of a side to a curved side 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,FirstElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL                   :: keepExistingCurveds
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER       :: Elem  ! ?
TYPE(tSide),POINTER       :: Side  ! ?
INTEGER                   :: i,nAns  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,*)'CREATING CURVED SURFACES FROM CURVED EDGES ...'
CALL Timer(.TRUE.)

IF(keepExistingCurveds)THEN
  Elem=>firstElem
  DO WHILE(ASSOCIATED(Elem))
    Side=>Elem%firstSide
    DO WHILE(ASSOCIATED(Side))
      IF(.NOT.ASSOCIATED(Side%Connection)  .OR. &
         .NOT.ASSOCIATED(Side%curvedNode)  .OR. &
              ASSOCIATED(Side%BC)         ) THEN
        Side=>Side%nextElemSide
        CYCLE
      END IF

      IF(.NOT.ASSOCIATED(Side%Connection%CurvedNode))THEN
        nAns=Side%nCurvedNodes
        Side%Connection%nCurvedNodes=nAns
        ALLOCATE(Side%Connection%CurvedNode(nAns))
        DO i=1,nAns
          Side%Connection%CurvedNode(i)%np=>Side%CurvedNode(i)%np
        END DO
        Side%Connection%curveIndex=Side%curveIndex
      END IF
      Side=>Side%nextElemSide
    END DO
    Elem=>Elem%nextElem
  END DO
END IF

Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(ASSOCIATED(Side%curvedNode))THEN !already curved
      Side=>Side%nextElemSide
      CYCLE
    END IF
    ! compute the splines between the edges
    IF(Side%nNodes.EQ.3)THEN
      CALL curvedEdgesToTriaSurf(Side)
    ELSE 
      CALL curvedEdgesToQuadSurf(Side)
    END IF

    IF(Side%curveIndex.EQ.0)THEN
      DO i=1,Side%nNodes
        IF(ASSOCIATED(Side%Edge(i)%edp%curvedNode)) Side%curveIndex=-1
      END DO
    END IF

    !check connection, but for periodic BCs Splines are not copied!
    IF(ASSOCIATED(Side%Connection).AND.Side%MortarType.EQ.0)THEN
      IF(.NOT.ASSOCIATED(Side%BC)) THEN
        IF(ASSOCIATED(Side%Connection%CurvedNode).AND..NOT.keepExistingCurveds) THEN
          CALL ABORT(__STAMP__,'connection already curved...')
        END IF
        nAns=Side%nCurvedNodes
        Side%Connection%nCurvedNodes=nAns
        ALLOCATE(Side%Connection%CurvedNode(nAns))
        DO i=1,nAns
          Side%Connection%CurvedNode(i)%np=>Side%CurvedNode(i)%np
        END DO
        Side%Connection%curveIndex=Side%curveIndex
      END IF !Associated(Side%BC)
    END IF !connection
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

!check
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side =>Elem%FirstSide
  DO WHILE(ASSOCIATED(Side))
    IF(ASSOCIATED(Side%CurvedNode).AND.ASSOCIATED(Side%connection))THEN
      IF(.NOT.ASSOCIATED(Side%Connection%CurvedNode))&
          CALL ABORT(__STAMP__,'connection not curved')
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO
CALL Timer(.FALSE.)
END SUBROUTINE curvedEdgesToSurf


SUBROUTINE curvedEdgesToTriaSurf(Side)
!===================================================================================================================================
! apply curving from curved element edges (aElem%side) to element Side (Side). Boundary Splines on edges are just equidistant
! interpolations of the curve.
!===================================================================================================================================
! MODULES
USE MOD_Basis_Vars,ONLY:EdgeToTria,TriaMapInv
USE MOD_Mesh_Vars,ONLY:tSide,tEdge,tNodePtr,getNewNode,N
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide),POINTER,INTENT(IN)          :: Side  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tEdge),POINTER          :: Edge  ! ?
TYPE(tNodePtr),POINTER       :: xpos(:)  ! ?
INTEGER                      :: i,j,nAns,iEdge  ! ?
REAL,DIMENSION(3)            :: xa,xb,xc,xt,sN  ! ?
!===================================================================================================================================
sN=1./REAL(N)
! allocate Side 
nAns=(N+1)*(N+2)/2
Side%nCurvedNodes=nAns
ALLOCATE(Side%CurvedNode(nAns))
ALLOCATE(xpos(nAns))
DO i=1,nAns
  NULLIFY(xpos(i)%np)
END DO
! fill edge points 
DO iEdge=1,Side%nNodes
  Edge=>Side%Edge(iEdge)%edp
  IF(.NOT.ASSOCIATED(Edge%CurvedNode))THEN !create curved nodes for linear edge
    ALLOCATE(Edge%curvedNode(N+1))
    Edge%curvedNode(1  )%np=>Edge%Node(1)%np
    Edge%curvedNode(N+1)%np=>Edge%Node(2)%np
    DO i=2,N
      CALL getNewNode(Edge%curvedNode(i)%np)
      Edge%curvedNode(i)%np%x=REAL(i-1)*sN*(Edge%Node(2)%np%x-Edge%Node(1)%np%x)+Edge%Node(1)%np%x
    END DO
  END IF

  IF(Side%EdgeOrientation(iEdge))THEN !oriented 
    DO i=1,N
      xpos(EdgeToTria(iEdge,i))%np=>Edge%CurvedNode(i)%np
    END DO
  ELSE !edge not oriented
    DO i=1,N
      xpos(EdgeToTria(iEdge,i))%np=>Edge%CurvedNode(N+2-i)%np
    END DO
  END IF !oriented
END DO

! tria, blending of 3 coons
DO j=1,N-1
  DO i=1,N-j-1
    CALL getNewNode(xpos(TriaMapInv(i,j))%np)
    xa=(REAL(N-i-j)*xpos(TriaMapInv(i,  0))%np%x +REAL(j)*xpos(TriaMapInv(i,  N-i))%np%x)/REAL(N-i)
    xb=(REAL(N-i-j)*xpos(TriaMapInv(0,  j))%np%x +REAL(i)*xpos(TriaMapInv(N-j,j  ))%np%x)/REAL(N-j) 
    xc=(REAL(i)    *xpos(TriaMapInv(i+j,0))%np%x +REAL(j)*xpos(TriaMapInv(0,  i+j))%np%x)/REAL(i+j) 
    xt=(REAL(N-i-j)*xpos(TriaMapInv(0,  0))%np%x +REAL(i)*xpos(TriaMapInv(N,  0  ))%np%x+ &
         REAL(j)*xpos(TriaMapInv(0,N))%np%x)*sN
    xpos(TriaMapInv(i,j))%np%x=0.5*(xa+xb+xc-xt)
  END DO
END DO
! convert to curvedNode
DO i=1,nAns
  Side%CurvedNode(i)%np=>xpos(i)%np
  NULLIFY(xpos(i)%np)
END DO
DEALLOCATE(xpos) 
END SUBROUTINE curvedEdgesToTriaSurf


SUBROUTINE curvedEdgesToQuadSurf(aSide)
!===================================================================================================================================
! apply curving from curved element edges (aElem%side) to element Side (aside). Boundary Splines on edges are just equidistant
! interpolations of the curve.
!===================================================================================================================================
! MODULES
USE MOD_Basis_Vars,ONLY:EdgeToQuad,QuadMapInv
USE MOD_Mesh_Vars,ONLY:tSide,tEdge,tNodePtr,getNewNode,N
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide),POINTER,INTENT(IN)          :: aSide  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tEdge),POINTER          :: aEdge  ! ?
TYPE(tNodePtr),POINTER       :: xpos(:)  ! ?
INTEGER                      :: i,j,nAns,iEdge,nNodes  ! ?
REAL,DIMENSION(3)            :: xa,xb,xt,sN  ! ?
!===================================================================================================================================
sN=1./REAL(N)
nNodes=aSide%nNodes
! allocate aSide 
nAns=(N+1)**2
aSide%nCurvedNodes=nAns
ALLOCATE(aSide%CurvedNode(nAns))
ALLOCATE(xpos(nAns))
DO i=1,nAns
  NULLIFY(xpos(i)%np)
END DO

! fill edge points
DO iEdge=1,nNodes
  aEdge=>aSide%Edge(iEdge)%edp
  IF(.NOT.ASSOCIATED(aEdge%CurvedNode))THEN
    ! create curved nodes for linear edge
    ALLOCATE(aEdge%curvedNode(N+1))
    aEdge%curvedNode(1  )%np=>aEdge%Node(1)%np
    aEdge%curvedNode(N+1)%np=>aEdge%Node(2)%np
    DO i=2,N
      CALL getNewNode(aEdge%curvedNode(i)%np)
      aEdge%curvedNode(i)%np%x=REAL(i-1)*sN*(aEdge%Node(2)%np%x-aEdge%Node(1)%np%x)+aEdge%Node(1)%np%x
    END DO
  END IF
  
  IF(aSide%EdgeOrientation(iEdge))THEN !oriented 
    DO i=1,N
      xpos(EdgeToQuad(iEdge,i))%np=>aEdge%CurvedNode(i)%np
    END DO
  ELSE !edge not oriented
    DO i=1,N
      xpos(EdgeToQuad(iEdge,i))%np=>aEdge%CurvedNode(N+2-i)%np
    END DO
  END IF !oriented
END DO

!coons mapping
DO j=1,N-1
  DO i=1,N-1
    CALL getNewNode(xpos(QuadMapInv(i,j))%np)
    xa= ( REAL(N-j) *xpos(QuadMapInv(i,0))%np%x       +REAL(j)*xpos(QuadMapInv(i,N))%np%x)
    xb= ( REAL(N-i) *xpos(QuadMapInv(0,j))%np%x       +REAL(i)*xpos(QuadMapInv(N,j))%np%x)
    xt= ( REAL((N-i)*(N-j))*xpos(QuadMapInv(0,0))%np%x+ &
          REAL((N-i)*  (j))*xpos(QuadMapInv(0,N))%np%x+ &
          REAL(  (i)*(N-j))*xpos(QuadMapInv(N,0))%np%x+ &
          REAL(  (i)*  (j))*xpos(QuadMapInv(N,N))%np%x )*sN
    xpos(QuadMapInv(i,j))%np%x=(xa+xb-xt)*sN
  END DO
END DO
! convert to curvedNode
DO i=1,nAns
  aSide%CurvedNode(i)%np=>xpos(i)%np
  NULLIFY(xpos(i)%np)
END DO
DEALLOCATE(xpos) 
END SUBROUTINE curvedEdgesToQuadSurf

SUBROUTINE curvedSurfacesToElem()
!===================================================================================================================================
! Blend curved surfaces of an element to a curved volume (TODO: currently only hexas working) 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,FirstElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER       :: Elem  ! ?
!===================================================================================================================================
!RETURN !DEBUG
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,*)'CREATING CURVED ELEMENTS FROM CURVED SURFACES ...'
CALL Timer(.TRUE.)

! first build unique curved / bilinear sides
CALL curvedEdgesToSurf(keepExistingCurveds=.FALSE.)

Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  SELECT CASE(Elem%nNodes)
  CASE(8)
    CALL curvedSurfacesToHexa(Elem)
  END SELECT
  Elem=>Elem%nextElem
END DO
CALL Timer(.FALSE.)
END SUBROUTINE curvedSurfacesToElem

SUBROUTINE curvedSurfacesToHexa(Elem)
!===================================================================================================================================
! apply curving from curved element Surfaces (Elem%side) to element (Elem). Boundary Splines on edges are just equidistant
! interpolations of the curve.
!===================================================================================================================================
! MODULES
USE MOD_Basis_Vars,ONLY:QuadMapInv,HexaMapInv
USE MOD_Mesh_Vars, ONLY:tElem,tSide,N
USE MOD_Mesh_Vars, ONLY:getNewNode,tNodePtr
USE MOD_Mesh_Basis,ONLY:getFlip
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(IN)          :: Elem  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tSide),POINTER          :: Side  ! ?
TYPE(tNodePtr),POINTER       :: xpos(:,:,:),xSide(:,:)  ! ?
INTEGER                      :: i,j,k,nAns,p,q,flip  ! ?
REAL,DIMENSION(3)            :: xFaces,xEdges,xNodes,xi,rxi  ! ?
!===================================================================================================================================
! allocate 
nAns=(N+1)**3
Elem%nCurvedNodes=nAns
ALLOCATE(Elem%CurvedNode(nAns))
ALLOCATE(xpos(0:N,0:N,0:N),xSide(0:N,0:N))
DO i=0,N; DO j=0,N; DO k=0,N
  NULLIFY(xpos(i,j,k)%np)
END DO; END DO; END DO

! fill corners
xpos(0,0,0)%np=>Elem%Node(1)%np
xpos(N,0,0)%np=>Elem%Node(2)%np
xpos(N,N,0)%np=>Elem%Node(3)%np
xpos(0,N,0)%np=>Elem%Node(4)%np
xpos(0,0,N)%np=>Elem%Node(5)%np
xpos(N,0,N)%np=>Elem%Node(6)%np
xpos(N,N,N)%np=>Elem%Node(7)%np
xpos(0,N,N)%np=>Elem%Node(8)%np
! fill side points 
Side=>Elem%firstSide
DO WHILE(ASSOCIATED(Side)) 
  flip=getFlip(Side)

  DO q=0,N; DO p=0,N
    NULLIFY(xSide(p,q)%np)
    SELECT CASE(flip)
    CASE(0)
      Xside(p,q)%np=>Side%CurvedNode(QuadMapInv(p,q))%np
    CASE(1)
      Xside(p,q)%np=>Side%CurvedNode(QuadMapInv(q,p))%np
    CASE(2)
      Xside(p,q)%np=>Side%CurvedNode(QuadMapInv(N-p,q))%np
    CASE(3)
      Xside(p,q)%np=>Side%CurvedNode(QuadMapInv(N-q,N-p))%np
    CASE(4)
      Xside(p,q)%np=>Side%CurvedNode(QuadMapInv(p,N-q))%np
    END SELECT  
  END DO; END DO ! p,q

  DO q=0,N; DO p=0,N
    SELECT CASE(Side%LocSide)
    CASE(5) !XI_MINUS)
      IF(.NOT.ASSOCIATED(xpos(0,p,q)%np)) &
      xpos(0,p,q)%np => Xside(q,p)%np 
    CASE(3) !XI_PLUS)
      IF(.NOT.ASSOCIATED(xpos(N,p,q)%np)) &
      xpos(N,p,q)%np => Xside(p,q)%np
    CASE(2) !ETA_MINUS)
      IF(.NOT.ASSOCIATED(xpos(p,0,q)%np)) &
      xpos(p,0,q)%np => Xside(p,q)%np  
    CASE(4) !ETA_PLUS)
      IF(.NOT.ASSOCIATED(xpos(p,N,q)%np)) &
      xpos(p,N,q)%np => Xside(N-p,q)%np 
    CASE(1) !ZETA_MINUS)
      IF(.NOT.ASSOCIATED(xpos(p,q,0)%np)) &
      xpos(p,q,0)%np => Xside(q,p)%np
    CASE(6) !ZETA_PLUS)
      IF(.NOT.ASSOCIATED(xpos(p,q,N)%np)) &
      xpos(p,q,N)%np => Xside(p,q)%np
    END SELECT
  END DO; END DO ! p,q
  Side=>Side%nextElemSide
END DO

!coons mapping, inner Nodes
DO k=1,N-1; DO j=1,N-1; DO i=1,N-1
  CALL getNewNode(xpos(i,j,k)%np,0)
  xi     = REAL((/i,j,k/))/REAL(N)
  rxi    = 1.-xi 
  xFaces = rxi(1)*xpos( 0, j, k)%np%x + &
            xi(1)*xpos( N, j, k)%np%x + &
           rxi(2)*xpos( i, 0, k)%np%x + &
            xi(2)*xpos( i, N, k)%np%x + &
           rxi(3)*xpos( i, j, 0)%np%x + &
            xi(3)*xpos( i, j, N)%np%x 

  xEdges = rxi(1)*rxi(2)*xpos( 0, 0, k)%np%x + &
           rxi(1)* xi(2)*xpos( 0, N, k)%np%x + &
            xi(1)*rxi(2)*xpos( N, 0, k)%np%x + &
            xi(1)* xi(2)*xpos( N, N, k)%np%x + &
           rxi(1)*rxi(3)*xpos( 0, j, 0)%np%x + &
           rxi(1)* xi(3)*xpos( 0, j, N)%np%x + &
            xi(1)*rxi(3)*xpos( N, j, 0)%np%x + &
            xi(1)* xi(3)*xpos( N, j, N)%np%x + &
           rxi(2)*rxi(3)*xpos( i, 0, 0)%np%x + &
           rxi(2)* xi(3)*xpos( i, 0, N)%np%x + &
            xi(2)*rxi(3)*xpos( i, N, 0)%np%x + &
            xi(2)* xi(3)*xpos( i, N, N)%np%x

  xNodes = rxi(1)*rxi(2)*rxi(3)*xpos( 0, 0, 0)%np%x + &
           rxi(1)* xi(2)*rxi(3)*xpos( 0, N, 0)%np%x + &
            xi(1)*rxi(2)*rxi(3)*xpos( N, 0, 0)%np%x + &
            xi(1)* xi(2)*rxi(3)*xpos( N, N, 0)%np%x + &
           rxi(1)*rxi(2)* xi(3)*xpos( 0, 0, N)%np%x + &
           rxi(1)* xi(2)* xi(3)*xpos( 0, N, N)%np%x + &
            xi(1)*rxi(2)* xi(3)*xpos( N, 0, N)%np%x + &
            xi(1)* xi(2)* xi(3)*xpos( N, N, N)%np%x
  xpos(i,j,k)%np%x=(xFaces-xEdges+xNodes)
END DO; END DO; END DO

! convert to curvedNode
DO k=0,N; DO j=0,N; DO i=0,N
  Elem%CurvedNode(HexaMapInv(i,j,k))%np=>xpos(i,j,k)%np
  NULLIFY(xpos(i,j,k)%np)
END DO; END DO; END DO
DEALLOCATE(xpos,xside) 
END SUBROUTINE curvedSurfacesToHexa


SUBROUTINE SplitToSpline()
!===================================================================================================================================
! Finds splitted elements from ANSA and associate to Element. The splitElem is a surface element, aElem is the normal 3D element
! and aSide is the element side to be curved.
! 1) insert all split elements in searchmesh
! 2) for each element side having a curveind>0 => aSide
!    a) find a split element attached to the first node of the aSide. Until now it is only known that they share a node, we
!       have to find out if the other corners of the element side are corresponding too...
!    b) using the connections between the splitElements, we find all splitElements belonging to one side 
!       .p.e. for quads we have a (bOrd-1)*(bOrd-1) splitelement array (the 'left' side of each element is saved=>localSides)
!    c) if we can buld up the whole array (reachedcorner=T), we can check, if the other aSide Nodes share corners
!    d) then the points areinterpolated to the monomial spline basis of aSide by a inverted Vandermonde matrix (sVdMquad).
!       We suppose that the points in parameter space are equidistant, thus we can use only one VdM for quads and one for trias.
!
!===================================================================================================================================
! MODULES
USE MOD_Basis_Vars,ONLY:TriaMap,QuadMap
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tSidePtr,tNode,tNodePtr
USE MOD_Mesh_Vars,ONLY:FirstElem,FirstSplitElem
USE MOD_Mesh_Vars,ONLY:N
USE MOD_Mesh_Vars,ONLY:ConformConnect
USE MOD_Mesh_Vars,ONLY:SpaceQuandt
USE MOD_Mesh_Vars,ONLY:DeleteSide,DisconnectElem
USE MOD_Mesh_Tolerances,ONLY:SAMEPOINT
USE MOD_Search
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
! normals on a Meshnode
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tNode),POINTER       :: Nodekplus,Nodekminus,aNode  ! ?
TYPE(tSide),POINTER       :: aSide,aSplitSide  ! ?
TYPE(tElem),POINTER       :: aElem,aSplitElem  ! ?
TYPE(tNodePtr),POINTER    :: localNodes(:,:)  ! ?
Type(tSidePtr),POINTER    :: localSides(:,:)  ! ?
REAL                      :: tol  ! ?
REAL                      :: xmin(3),xmax(3),dxmax  ! ?
LOGICAL                   :: pointFound,sidefound,oriented,reachedcorner  ! ?
INTEGER                   :: i,j,k,l,bOrd,tria,nAns,nSplitElems,idx(3)  ! ?
TYPE(tToObject),POINTER   :: ToObject   ! ?
TYPE(tSearchMesh),POINTER :: searchMesh,redundantNodes  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'Bulding Splines from splitted ANSA elements...'
CALL Timer(.TRUE.)
!========================================== Init================================================================================
bOrd=N+1
ALLOCATE(localNodes(bOrd,bOrd))
ALLOCATE(localSides(bOrd-1,bOrd-1))
xmin=1.E14
xmax=-1.E14
dxmax=0.
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  sidefound=.false.
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF(aSide%curveIndex.GT.0) THEN
      sidefound=.TRUE.
      EXIT
    END IF
    aSide=>aSide%nextElemSide
  END DO
  IF(sidefound) THEN
    xmin=MIN(xmin,aElem%Node(1)%np%x)
    xmax=MAX(xmax,aElem%Node(1)%np%x)
    DO i=2,aElem%nNodes
      xmin=MIN(xmin,aElem%Node(i)%np%x)
      xmax=MAX(xmax,aElem%Node(i)%np%x)
      dxmax=MAX(dxmax,SUM((aElem%Node(i)%np%x-aElem%Node(1)%np%x)**2))
    END DO
  END IF !found
  aElem=>aElem%nextElem
END DO 
! enlarge box 
xmin=xmin-SQRT(dxmax)
xmax=xmax+SQRT(dxmax)
!build local search mesh of split sides
ConformConnect=.FALSE. ! damit nSearch nicht auf 1!!!
NULLIFY(searchMesh)
CALL getNewSearchMesh(searchMesh,.FALSE.,xmin,xmax)
NULLIFY(redundantNodes)
CALL getNewSearchMesh(redundantNodes,.FALSE.,xmin,xmax)
nSplitelems=0
aSplitElem=>firstSplitElem
DO WHILE(ASSOCIATED(aSplitElem))
 CALL getIdx(searchmesh,aSplitElem%Node(1)%np%x,idx)  ! Calculate search mesh indices of node
  IF(idxok(searchMesh,idx)) THEN
    nSplitElems=nSplitElems+1
    CALL insertNode(searchMesh,aSplitElem%Node(1)%np,aSplitElem)
  END IF
  aSplitElem=>aSplitElem%nextElem
END DO
WRITE(*,*)'number of splitElems in searchmesh',nSplitElems 
!========================================== Start search of split elements =====================================================

aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF(aSide%curveIndex.NE.0) THEN 
      aNode=>aSide%OrientedNode(1)%np
      !local element side tolerance
      tol=SQRT(SUM((aSide%OrientedNode(2)%np%x-aSide%OrientedNode(1)%np%x)**2))
      DO i=2,aSide%nNodes
        tol=min(tol,SQRT(SUM((aSide%OrientedNode(i)%np%x-aSide%OrientedNode(i-1)%np%x)**2)))
      END DO
      tol=0.1*tol/REAL(bOrd)/SpaceQuandt
!        WRITE(*,'(A,E21.10)')'DEBUG tol',tol
      ToObject=>getFirstToObject(searchMesh,.FALSE.,aNode%x)
      sidefound=.FALSE.
      DO WHILE (ASSOCIATED(ToObject))
        IF(ToObject%elem%ind.EQ.0)THEN
          ToObject=>getNextToObject(searchMesh,.FALSE.)
        ELSEIF(ToObject%elem%nNodes.NE.aSide%nNodes) THEN
          ToObject=>getNextToObject(searchMesh,.FALSE.)
        ELSE
          pointFound=.FALSE.
          aSplitElem=>ToObject%elem
          aSplitSide=>aSplitElem%firstSide
          DO k=1,aSplitElem%nNodes
            IF (SAMEPOINT(aSplitSide%Node(2)%np%x,aNode%x,tol)) THEN
              pointFound=.TRUE.
              EXIT 
            END IF
            aSplitSide=>aSplitSide%nextElemSide
          END DO
          IF(pointfound) THEN
            DO i=1,bOrd-1
              DO j=1,bOrd-1
                NULLIFY(localSides(i,j)%sp)
              END DO
            END DO
            !we suppose first that the element found is the corner element and check it later
            reachedcorner=.TRUE.
            localSides(1,1)%sp=>aSplitSide
            ! save left lying sides in local sides, bOrd=3 --> 4 sides
            !
            !  ^ j-dir 
            !  | 
            !  | 
            ! (4)---(+)---(3)
            !    |     |  
            ! (+)---(+)---(+)
            !    |     |  
            ! (1)---(+)---(2) --> i-dir
            !
            tria=0
            IF(aSplitElem%nNodes.EQ.3) tria=1
!              WRITE(*,'(A,6F10.5)')'DEBUG tria',tria 
            !now fill up the whole array
outer:      DO j=1,bOrd-1
              ! j direction, for j=1 do nothing
              IF(j.GT.1)THEN
                !
                !go upwards in j direction
                aSplitSide=>GoToSide(localSides(1,j-1)%sp,-1)
                DO l=1,1+tria  ! l=1,1 for quads, l=1,2 for trias
                  aSplitSide=>aSplitSide%Connection
                  IF(.NOT.ASSOCIATED(aSplitSide)) THEN
                    reachedcorner=.FALSE. 
                    EXIT !no corner found, because no connection (end of local domain)
                  END IF
                  aSplitSide=>GoToSide(aSplitSide,-1)
                END DO 
                localSides(1,j)%sp=>aSplitSide
              END IF
              DO i=2,bOrd-1-(j-1)*tria
                aSplitSide=>GoToSide(localSides(i-1,j)%sp,1)
                DO l=1,1+tria
                  aSplitSide=>GoToSide(aSplitSide,1)
                  !go to neighbor element
                  aSplitSide=>aSplitSide%connection
                  IF(.NOT.ASSOCIATED(aSplitSide)) THEN
                    reachedcorner=.FALSE. 
                    EXIT outer !no corner found, because no connection (end of local domain)
                  END IF
                END DO
                !go to neighbor element
                localSides(i,j)%sp=>aSplitSide
              END DO !i
            END DO outer !j
            IF(.NOT.reachedcorner)THEN
              WRITE(*,*)'DEBUG corner not reached, no mpi -> PROBLEM!'
            END IF ! corner not reached 
             
            IF(reachedcorner)THEN
              ! all local Sides are associated, check if localSides belong to aSide
              nodekminus=>localSides(1,bord-1)%sp%node(1)%np
              aSplitSide=>GoToSide(localsides(bord-1,1)%sp,1)
              nodekplus =>aSplitSide%node(2)%np
              IF((SAMEPOINT(aSide%OrientedNode(aSide%nNodes)%np%x,nodekminus%x,tol)).AND. &
                         (SAMEPOINT(aSide%OrientedNode(2)%np%x,nodekplus%x,tol))) THEN
                oriented=.TRUE.
                sidefound=.TRUE.
              ELSEIF((SAMEPOINT(aSide%OrientedNode(aSide%nNodes)%np%x,nodekplus%x,tol)).AND. &
                         (SAMEPOINT(aSide%OrientedNode(2)%np%x,nodekminus%x,tol))) THEN
                oriented=.FALSE.
                sidefound=.TRUE.
              END IF
              IF(sidefound)THEN
                !fill Node array!
                DO i=1,bOrd
                  DO j=1,bOrd
                    NULLIFY(localNodes(i,j)%np)
                  END DO
                END DO 
                DO j=1,bOrd-1
                  DO i=1,bOrd-1-(j-1)*tria
                    localNodes(i,j)%np=>localSides(i,j)%sp%Node(2)%np
                  END DO
                  aSplitSide=>GoToSide(localSides(bOrd-1-(j-1)*tria,j)%sp,1)
                  localNodes(bOrd-(j-1)*tria,j)%np=>aSplitSide%Node(2)%np
                END DO 
                IF(tria.EQ.0)THEN
                  DO i=1,bOrd-1
                    localNodes(i,bOrd)%np=>localSides(i,bOrd-1)%sp%Node(1)%np
                  END DO
                  aSplitSide=>GoToSide(localSides(bOrd-1,bOrd-1)%sp,2)
                  localNodes(bOrd,bOrd)%np=>aSplitSide%Node(2)%np
                ELSE
                  localNodes(1,bOrd)%np=>localSides(1,bOrd-1)%sp%Node(1)%np
                END IF

                !associate corner Nodes to element nodes exactly
                CALL insertNode(redundantNodes,localNodes(1,1)%np)
                CALL insertNode(redundantNodes,localNodes(bOrd,1)%np)
                CALL insertNode(redundantNodes,localNodes(1,bOrd)%np)
                localNodes(1,1)%np=>aSide%OrientedNode(1)%np
                IF(tria.EQ.0) THEN
                   CALL insertNode(redundantNodes,localNodes(bOrd,bOrd)%np)
                   localNodes(bord,bord)%np=>aSide%OrientedNode(3)%np
                END IF
                IF(oriented)THEN
                  localNodes(bord,1)%np=>aSide%OrientedNode(2)%np
                  localNodes(1,bord)%np=>aSide%OrientedNode(aSide%nNodes)%np
                ELSE
                  localNodes(1,bord)%np=>aSide%OrientedNode(2)%np
                  localNodes(bord,1)%np=>aSide%OrientedNode(aSide%nNodes)%np
                END IF
                EXIT ! searchloop
              END IF !side found
            END IF !reachedcorner
          END IF !point found
          ToObject=>getNextToObject(searchMesh,.FALSE.)
        END IF ! ToObject%elem%ind .NE.0
      END DO !ToObject
      IF(.NOT.sidefound) THEN
         ERRWRITE(*,*)'No associated side found for boundary order=',N+1
         DO i=1,aSide%nNodes
           ERRWRITE(*,'(A,I3,3E21.10)')'side corner node',i,aSide%OrientedNode(i)%np%x
         END DO
         CALL abort(__STAMP__, &
              'No associated side found for boundary order=',N+1)
      ELSE
        !convert to curvedNode pointers
        IF(aSide%nNodes.EQ.3) THEN
          nAns=bOrd*(bOrd+1)/2
          aSide%nCurvedNodes=nAns
          ALLOCATE(aSide%CurvedNode(nAns))
          aSide%isCurved=.TRUE.
          IF(oriented) THEN
            DO i=1,nAns
              aSide%CurvedNode(i)%np=>localNodes(TriaMap(i,1)+1,TriaMap(i,2)+1)%np
            END DO !i
          ELSE !not oriented
            DO i=1,nAns
              aSide%CurvedNode(i)%np=>localNodes(TriaMap(i,2)+1,TriaMap(i,1)+1)%np
            END DO !i
          END IF
        ELSE !quad
          nAns=bOrd*bOrd
          aSide%nCurvedNodes=nAns
          ALLOCATE(aSide%CurvedNode(nAns))
          aSide%isCurved=.TRUE.
          IF(oriented) THEN
            DO i=1,nAns
              aSide%CurvedNode(i)%np=>localNodes(QuadMap(i,1)+1,QuadMap(i,2)+1)%np
            END DO !i
          ELSE !not oriented
            DO i=1,nAns
              aSide%CurvedNode(i)%np=>localNodes(QuadMap(i,2)+1,QuadMap(i,1)+1)%np
            END DO !i
          END IF !oriented
        END IF !tri/quad
        CALL curvedSurfToEdges(aSide)
      END IF !not sidefound
    END IF !curveInd .NE.0
    aSide=>aSide%nextElemSide
  END DO
  aElem=>aElem%nextElem
END DO 
CALL deleteSearchMesh(SearchMesh,.FALSE.)
CALL deleteSearchMesh(redundantNodes,.TRUE.)
!delete all split elements and nodes 
aSplitElem=>firstSplitElem
DO WHILE (ASSOCIATED(aSplitElem)) 
  DO i=1,aSplitElem%nNodes
    aSplitElem%Node(i)%np%refcount=0
  END DO
  aSplitElem=>aSplitElem%nextElem
END DO
aSplitElem=>firstSplitElem
DO WHILE (ASSOCIATED(aSplitElem)) 
  DO i=1,aSplitElem%nNodes
    aSplitElem%Node(i)%np%refcount=aSplitElem%Node(i)%np%refcount+1
  END DO
  aSplitSide=>aSplitElem%firstSide
  DO WHILE(ASSOCIATED(aSplitSide))
    DO i=1,aSplitSide%nNodes
      aSplitSide%Node(i)%np%refcount=aSplitSide%Node(i)%np%refcount+1
    END DO
    aSplitSide=>aSplitSide%nextelemSide
  END DO
  aSplitElem=>aSplitElem%nextElem
END DO
i=0 !number of elements
l=0 !number of nodes
aSplitElem=>firstSplitElem
DO WHILE (ASSOCIATED(aSplitElem))
 ! Delete sides
  DO WHILE(ASSOCIATED(aSplitElem%firstSide))
    aSplitSide=>aSplitElem%firstSide
    CALL deleteSide(aSplitElem%firstSide,aSplitSide)
  END DO
  CALL DisconnectElem(firstSplitElem,aSplitElem)  ! Disconnect element from mesh
 ! Delete nodes
  IF(ASSOCIATED(aSplitElem%Node)) THEN
    DO j=1,aSplitElem%nNodes
      NULLIFY(aSplitElem%Node(j)%np)
      l=l+1
    END DO
    DEALLOCATE(aSplitElem%Node)
  END IF
  i=i+1
  DEALLOCATE(aSplitElem)
  aSplitElem=>firstSplitElem
END DO    
DEALLOCATE(localNodes,localSides)
WRITE(*,'(A,I8)')'nSplitElems deleted',i
WRITE(*,'(A,I8)')'nNodes of SplitElems deleted',l
CALL Timer(.FALSE.)
END SUBROUTINE SplitToSpline


FUNCTION GoToSide(Side,nJumps)
!===================================================================================================================================
! repointer side to the next nJumps 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)             :: nJumps  ! ?
TYPE(tSide),POINTER,INTENT(IN) :: Side  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tSide),POINTER            :: GoToSide  ! ?
! normals on a Meshnode
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i  ! ?
!===================================================================================================================================
GoToSide=>Side
IF(nJumps.GE.0)THEN
  DO i=1,nJumps
    GoToSide=>GoToSide%nextElemSide
    IF(.NOT.ASSOCIATED(GoToSide))GoToSide=>Side%Elem%firstSide
  END DO
ELSE
  DO i=1,Side%Elem%nNodes+nJumps
    GoToSide=>GoToSide%nextElemSide
    IF(.NOT.ASSOCIATED(GoToSide))GoToSide=>Side%Elem%firstSide
  END DO
END IF
END FUNCTION GoToSide 


SUBROUTINE RebuildMortarGeometry()
!===================================================================================================================================
! for curved mortarmeshes ensure that small mortar geometry is identical to big mortar geometry
!===================================================================================================================================
! MODULES
USE MOD_Globals  ,ONLY:abort
USE MOD_Mesh_Vars,ONLY:tElem,FirstElem,tSide,tEdge
USE MOD_Mesh_Vars,ONLY:M_0_2_T,M_0_1_T,N
USE MOD_Basis1D,  ONLY:GetMortarVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER       :: Elem  ! ?
TYPE(tSide),POINTER       :: Side  ! ?
TYPE(tEdge),POINTER       :: Edge  ! ?
INTEGER                   :: iEdge ! ?
LOGICAL                   :: periodicMortarFound ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,*)'REBUILDING CURVED MORTAR INTERFACES...'

ALLOCATE(M_0_1_T(0:N,0:N))
ALLOCATE(M_0_2_T(0:N,0:N))
CALL GetMortarVandermonde(N, M_0_1_T, M_0_2_T) 
M_0_1_T=TRANSPOSE(M_0_1_T)
M_0_2_T=TRANSPOSE(M_0_2_T)

periodicMortarFound=.FALSE.
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%MortarType.GT.0) THEN
      IF(ASSOCIATED(Side%BC))THEN
        IF(Side%BC%BCType.EQ.1) periodicMortarFound=.TRUE.
      END IF
      CALL MapBigSideToSmall(Side)
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%MortarType.GT.0)THEN
      DO iEdge=1,Side%nNodes
        Edge=>Side%Edge(iEdge)%edp
        IF(.NOT.ASSOCIATED(Edge%parentEdge).AND.ASSOCIATED(Edge%MortarEdge))THEN
          CALL MapBigEdgeToSmall(Edge)
        END IF
      END DO
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

! check if periodics are still watertight
IF(periodicMortarFound)THEN
  WRITE(UNIT_stdOut,'(132("!"))')
  WRITE(UNIT_stdOut,*) ' WARNING: Periodic mortar edges have been found in the mesh. Rebuilding these mortar edges is not yet implemented!'
  WRITE(UNIT_stdOut,*) '          The mesh may not be watertight along the edges of curved periodic mortars.'
  WRITE(UNIT_stdOut,'(132("!"))')
  !CALL abort(__STAMP__,&
  !  'Rebuilding curved periodic mortar edges is not yet implemented.')
END IF

WRITE(UNIT_stdOut,*)'...DONE.'
WRITE(UNIT_stdOut,'(132("~"))')
END SUBROUTINE RebuildMortarGeometry


SUBROUTINE MapBigSideToSmall(Side)
!===================================================================================================================================
! for curved mortarmeshes ensure that small mortar geometry is identical to big mortar geometry
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tSide,VV
USE MOD_Mesh_Vars,ONLY:M_0_1_T,M_0_2_T,N
USE MOD_Mesh_Basis,ONLY:PackGeo,UnpackGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide),POINTER,INTENT(IN)   :: Side  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: nMortars,p,q,l
REAL                             :: XGeo2DBig(   3,0:N,0:N)
REAL                             :: XGeo2DSmall4(3,0:N,0:N,4)
REAL                             :: XGeo2DSmall2(3,0:N,0:N,2)
REAL                             :: dir
!-----------------------------------------------------------------------------------------------------------------------------------
XGeo2DBig=0.
XGeo2DSmall4=0.
XGeo2DSmall2=0.
CALL PackGeo(N,Side,XGeo2DBig)

! In case of periodic BCs always add displacement vector. Pay attention to direction!
IF(ASSOCIATED(Side%BC))THEN
  IF(Side%BC%BCType.EQ.1.AND.Side%BC%BCalphaInd.NE.0)THEN
    dir=1.*SIGN(1,Side%BC%BCalphaInd)
    DO q=0,N; DO p=0,N
      XGeo2DBig(:,p,q)=XGeo2DBig(:,p,q) + dir*VV(:,ABS(Side%BC%BCalphaInd))
    END DO; END DO
  END IF
END IF

SELECT CASE(Side%MortarType)
CASE(1) !1->4
  nMortars=4
  !first in eta
  DO q=0,N; DO p=0,N
      XGeo2DSmall2(:,p,q,1)=                      M_0_1_T(0,q)*XGeo2DBig(:,p,0)
      XGeo2DSmall2(:,p,q,2)=                      M_0_2_T(0,q)*XGeo2DBig(:,p,0)
    DO l=1,N
      XGeo2DSmall2(:,p,q,1)=XGeo2DSmall2(:,p,q,1)+M_0_1_T(l,q)*XGeo2DBig(:,p,l)
      XGeo2DSmall2(:,p,q,2)=XGeo2DSmall2(:,p,q,2)+M_0_2_T(l,q)*XGeo2DBig(:,p,l)
    END DO
  END DO; END DO
  ! then in xi
  DO q=0,N; DO p=0,N
      XGeo2DSmall4(:,p,q,1)=                      M_0_1_T(0,p)*XGeo2DSmall2(:,0,q,1)
      XGeo2DSmall4(:,p,q,2)=                      M_0_2_T(0,p)*XGeo2DSmall2(:,0,q,1)
      XGeo2DSmall4(:,p,q,3)=                      M_0_1_T(0,p)*XGeo2DSmall2(:,0,q,2)
      XGeo2DSmall4(:,p,q,4)=                      M_0_2_T(0,p)*XGeo2DSmall2(:,0,q,2)
    DO l=1,N
      XGeo2DSmall4(:,p,q,1)=XGeo2DSmall4(:,p,q,1)+M_0_1_T(l,p)*XGeo2DSmall2(:,l,q,1)
      XGeo2DSmall4(:,p,q,2)=XGeo2DSmall4(:,p,q,2)+M_0_2_T(l,p)*XGeo2DSmall2(:,l,q,1)
      XGeo2DSmall4(:,p,q,3)=XGeo2DSmall4(:,p,q,3)+M_0_1_T(l,p)*XGeo2DSmall2(:,l,q,2)
      XGeo2DSmall4(:,p,q,4)=XGeo2DSmall4(:,p,q,4)+M_0_2_T(l,p)*XGeo2DSmall2(:,l,q,2)
    END DO !l=1,N
  END DO; END DO !p,q=0,N

CASE(2) !1->2 in eta
  nMortars=2
  DO q=0,N; DO p=0,N
      XGeo2DSmall4(:,p,q,1)=                      M_0_1_T(0,q)*XGeo2DBig(:,p,0)
      XGeo2DSmall4(:,p,q,2)=                      M_0_2_T(0,q)*XGeo2DBig(:,p,0)
    DO l=1,N
      XGeo2DSmall4(:,p,q,1)=XGeo2DSmall4(:,p,q,1)+M_0_1_T(l,q)*XGeo2DBig(:,p,l)
      XGeo2DSmall4(:,p,q,2)=XGeo2DSmall4(:,p,q,2)+M_0_2_T(l,q)*XGeo2DBig(:,p,l)
    END DO
  END DO; END DO

CASE(3) !1->2 in xi
  nMortars=2
  DO q=0,N; DO p=0,N
      XGeo2DSmall4(:,p,q,1)=                      M_0_1_T(0,p)*XGeo2DBig(:,0,q)
      XGeo2DSmall4(:,p,q,2)=                      M_0_2_T(0,p)*XGeo2DBig(:,0,q)
    DO l=1,N
      XGeo2DSmall4(:,p,q,1)=XGeo2DSmall4(:,p,q,1)+M_0_1_T(l,p)*XGeo2DBig(:,l,q)
      XGeo2DSmall4(:,p,q,2)=XGeo2DSmall4(:,p,q,2)+M_0_2_T(l,p)*XGeo2DBig(:,l,q)
    END DO
  END DO; END DO
END SELECT ! mortarType(SideID)

DO p=1,nMortars
  CALL UnpackGeo(N,XGeo2DSmall4(:,:,:,p),Side%MortarSide(p)%sp)
END DO

END SUBROUTINE MapBigSideToSmall


RECURSIVE SUBROUTINE MapBigEdgeToSmall(edge)
!===================================================================================================================================
! for curved mortarmeshes ensure that small mortar geometry is identical to big mortar geometry
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tEdge
USE MOD_Mesh_Vars,ONLY:M_0_1_T,M_0_2_T,N
USE MOD_Mesh_Basis,ONLY:PackGeo,UnpackGeo
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tEdge),POINTER,INTENT(IN)   :: Edge  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: p,l
REAL                             :: XGeo1DBig(  3,0:N)
REAL                             :: XGeo1DTmp(  3,0:N)
REAL                             :: XGeo1DSmall(3,0:N,2)
TYPE(tEdge),POINTER              :: smallEdge  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
XGeo1DBig=0.
XGeo1DSmall=0.
CALL PackGeo(N,Edge,XGeo1DBig)

! split edge into 2 in xi
DO p=0,N
    XGeo1DSmall(:,p,1)=                   M_0_1_T(0,p)*XGeo1DBig(:,0)
    XGeo1DSmall(:,p,2)=                   M_0_2_T(0,p)*XGeo1DBig(:,0)
  DO l=1,N
    XGeo1DSmall(:,p,1)=XGeo1DSmall(:,p,1)+M_0_1_T(l,p)*XGeo1DBig(:,l)
    XGeo1DSmall(:,p,2)=XGeo1DSmall(:,p,2)+M_0_2_T(l,p)*XGeo1DBig(:,l)
  END DO
END DO

DO p=1,2
  smallEdge=>Edge%MortarEdge(p)%edp
  IF(ASSOCIATED(Edge%Node(p)%np,smallEdge%Node(p)%np))THEN
  ELSEIF(ASSOCIATED(Edge%Node(p)%np,smallEdge%Node(3-p)%np))THEN
    XGeo1DTmp=XGeo1Dsmall(:,:,p)
    DO l=0,N
      XGeo1Dsmall(:,N-l,p)=XGeo1DTmp(:,l) 
    END DO
  ELSE 
    STOP "Error: Edges of mortar master and slave do not conform!"
  END IF
  CALL UnpackGeo(N,XGeo1DSmall(:,:,p),smallEdge)
  IF(ASSOCIATED(smallEdge%MortarEdge)) &
    CALL MapBigEdgeToSmall(smallEdge)
END DO

END SUBROUTINE MapBigEdgeToSmall


END MODULE MOD_Curved
