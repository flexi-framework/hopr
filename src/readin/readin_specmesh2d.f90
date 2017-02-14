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
MODULE MOD_Readin_SpecMesh2D
!==================================================================================================================================
! ?
!==================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part --------------------------------------------------------------------------------------------------------------------
! Public Part ---------------------------------------------------------------------------------------------------------------------
INTERFACE ReadSpecMesh2D
  MODULE PROCEDURE ReadSpecMesh2D
END INTERFACE

PUBLIC::ReadSpecMesh2D
!==================================================================================================================================

CONTAINS
SUBROUTINE ReadSpecMesh2D()
!==================================================================================================================================
! Read mesh from gambit ascii or binary file. Called by fillMesh.
! Read-in can be performed by just one or all processors
!==================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tElemPtr,tSide,tNode,tNodePtr
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:n2dNodes
USE MOD_Mesh_Vars,ONLY:nMeshFiles,MeshFileName
USE MOD_Mesh_Vars,ONLY:nUserDefinedBoundaries,BoundaryName,BoundaryType
USE MOD_Mesh_Vars,ONLY:N,BoundaryOrder,DZ,MeshDim
USE MOD_Mesh_Basis,ONLY:createSides
USE MOD_Mesh_Vars,ONLY:getNewNode
USE MOD_CurvedCartMesh,ONLY:getNewCurvedHexahedron
USE MOD_ChangeBasis,ONLY:ChangeBasis2D
USE MOD_Basis1D,ONLY:ChebyGaussLobNodesAndWeights
USE MOD_Basis1D,ONLY:BarycentricWeights
USE MOD_Basis1D,ONLY:InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! nMeshFiles                : Number of mesh files (INI-File)
! MeshFileName(iFile)           : Filename of mesh file iFile (INI-File)
! nZones                    : Number of mesh zones (INI-File)
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tNode),POINTER    :: Node   ! ?
TYPE(tElem),POINTER    :: Elem  ! ?
TYPE(tSide),POINTER    :: Side   ! ?
REAL,ALLOCATABLE       :: nodeCoords(:,:)  ! ?

!NEW
INTEGER                         :: os,iNode,iFile,iElem,i,j,k
INTEGER                         :: globalNodeID 
INTEGER                         :: fUnit = 108
INTEGER                         :: nodeIDs(4)
INTEGER                         :: curveFlags(4)
INTEGER                         :: nElems, nNodes, bCurveOrder
INTEGER                         :: nEdges
INTEGER                         :: EdgeMap(2,4)
REAL                            :: x(2),xStart(2), xEnd(2)
REAL,ALLOCATABLE                :: xCLedges(:,:,:)
REAL,ALLOCATABLE                :: xCLelem(:,:,:)
REAL,ALLOCATABLE                :: xElem(:,:,:)
REAL,ALLOCATABLE                :: xiCL(:),xiEq(:),wBaryCL(:),Vdm_CL_Eq(:,:)
CHARACTER(LEN=40)               :: BCnames(4)
TYPE(tNodePtr),POINTER          :: CurvedNode(:,:,:) 
!==================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')'Reading SpecMesh 2D mesh...'

EdgeMap(1:2,1)=(/1,2/)
EdgeMap(1:2,2)=(/2,3/)
EdgeMap(1:2,3)=(/3,4/)
EdgeMap(1:2,4)=(/4,1/)


!SET MESHDIM BACK TO 3, since we directly create the first layer of curved 3D elements (height DZ)
MeshDim=3

globalNodeID=0
n2dNodes=0  ! 2.5D mesh
! start reading mesh
DO iFile=1,nMeshFiles
  ! Open Gambit file
  OPEN(UNIT   = funit,                  &
       FILE   = MeshFileName(iFile), &
       STATUS = 'OLD',                &
       ACTION = 'READ',               &
       ACCESS = 'SEQUENTIAL',         &
       IOSTAT = os                    )
  IF ( os /= 0 )     THEN
     WRITE(UNIT_stdOut,*)  'Error opening file: ', TRIM(MeshFileName(iFile))
     STOP
  ELSE
    WRITE(UNIT_stdOut,*)  'Reading mesh from file: ',TRIM(MeshFileName(iFile))
  END IF
  !  -----------------------
  !  Read header information
  !  -----------------------
  ! Get number of nodes, elements and BoundaryOrder
  READ(fUnit,*) nNodes, nElems, bCurveOrder

  !ALLOCATE arrays
  ALLOCATE(nodeCoords(2,nNodes))

  IF(iFile.EQ.1) THEN
    N = bCurveOrder !init boundaryOrder for first mesh file
  ELSE
    IF(bCurveOrder.NE.N) STOP 'all meshfiles need to have the same bCurveOrder' !next file check
  END IF
  boundaryOrder=N+1

  IF(iFile.EQ.1) THEN
    ALLOCATE(xCLedges(2,0:N,4) )
    ALLOCATE(xCLelem(2,0:N,0:N) )
    ALLOCATE(xelem(2,0:N,0:N) )
    ALLOCATE(CurvedNode(0:N,0:N,0:N) )
    ALLOCATE( xiCL   (0:N) )
    ALLOCATE( xiEQ   (0:N) )
    ALLOCATE( wBaryCL(0:N) )
    ALLOCATE( Vdm_CL_EQ(0:N,0:N) )
    CALL ChebyGaussLobNodesAndWeights(N,xiCL)
    DO j = 0, N 
       xiEq(j) = 2.*REAL(j)/REAL(N)-1.
    END DO !j=0,N
    CALL BarycentricWeights(N,xiCL,wBaryCL)
    CALL InitializeVandermonde(N,N,wBaryCL,xiCL,xiEQ,Vdm_CL_Eq)
  END IF

  n2dNodes=n2dNodes+nNodes  
  ! Read node coordinates

  nodeCoords=0.
  DO iNode = 1, nNodes   ! READ nodes
    READ(funit,*)nodeCoords(1:2,iNode)
  END DO


  ! Read elements, go the last element
  Elem=>firstElem
  IF(ASSOCIATED(Elem))THEN
    DO WHILE(ASSOCIATED(Elem%nextElem))
      Elem%Ind=-1
      Elem=>Elem%nextElem
    END DO
    Elem%Ind=-1
  END IF

  ! -----------------------------------------
  ! Read elements: Elements have the format
  ! node1 node2 node3 node4
  ! b1 b2 b3 b4
  ! (=0 for straight side, 1 for curved)
  ! if curved boundaries, then for each:
  ! for j = 0 to boundaryOrder
  !    x_j  y_j
  ! next j
  ! bname1 bname2 bname3 bname4
  ! -----------------------------------------

  DO iElem = 1, nElems 
    READ( fUnit, * ) nodeIDs
    READ( fUnit, * ) curveFlags
    DO k = 1, 4 
      IF ( curveFlags(k) .EQ. 0 )     THEN !linear edge
        xStart = NodeCoords(1:2,nodeIDs(edgeMap(1,k)))
        xEnd   = NodeCoords(1:2,nodeIDs(edgeMap(2,k)))
        DO j=0,N
          xCLedges(:,j,k)=xstart+0.5*(xiCL(j)+1.)*(xEnd-xStart) !linear interpol. 
        END DO
      ELSE !curved edge
        DO j=0,N
          READ(fUnit,*) xCLedges(:,j,k)
        END DO
      END IF
    END DO
    xCLelem= TransfiniteEdgeToFace(N,xCLedges,xiCL)
    CALL ChangeBasis2D(2,N,N,Vdm_CL_EQ,XCLelem,xElem)
     
    DO k=0,N; DO j=0,N; DO i=0,N
      NULLIFY(CurvedNode(i,j,k)%np)
      CALL getNewNode(Node,1)
      globalNodeID=globalNodeID+1
      Node%ind=globalNodeID
      Node%x(1:2)=xElem(1:2,i,j)
      Node%x(3)  =REAL(k)/REAL(N)*DZ
      CurvedNode(i,j,k)%np=>Node     
    END DO; END DO; END DO !i,j,k 
    !puts a new curved element in the list
WRITE(*,*)'DEBUG,iELem',iElem
    CALL GetNewCurvedHexahedron(CurvedNode,N,iFile) !Zone=iFile


    READ( fUnit, * ) BCnames(1:4)

  END DO !iElem
     
  CLOSE(fUnit)
  DEALLOCATE(NodeCoords)
END DO !iFile
DEALLOCATE(xCLedges)
DEALLOCATE(xCLelem)
DEALLOCATE(xElem)
DEALLOCATE(CurvedNode)
DEALLOCATE( xiCL)
DEALLOCATE( xiEQ)
DEALLOCATE( wBaryCL)
DEALLOCATE( Vdm_CL_EQ)

CALL Timer(.FALSE.)
WRITE(UNIT_stdOut,'(132("~"))')
STOP
END SUBROUTINE readSpecMesh2D

FUNCTION TransfiniteEdgeToFace(N,xEdges,xi1D) RESULT(xFace)
!==================================================================================================================================
! transfinite map of 4 edges oriented counter-clockwise, to a quad
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N
REAL   ,INTENT(IN) :: xEdges(1:2,0:N,1:4)
REAL   ,INTENT(IN) :: xi1D(0:N) ![-1,1]
!----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL               :: xFace(1:2,0:N,0:N)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                   :: i,j,l
REAL,DIMENSION(2,0:N,0:N) :: xC,xF1,xF2
REAL                      :: xi,eta
!==================================================================================================================================
xC(:,0,0)=xEdges(:,0,1)
xC(:,N,0)=xEdges(:,0,2)
xC(:,N,N)=xEdges(:,0,3)
xC(:,0,N)=xEdges(:,0,4)
DO l=0,N
  xF1(:,  l,  0)=xEdges(:,l,1)
  xF1(:,N-l,  N)=xEdges(:,l,3)

  xF2(:,  0,N-l)=xEdges(:,l,4)
  xF2(:,  N,  l)=xEdges(:,l,2)
END DO
!linear corner interpolation
DO i=1,N-1
  xi=0.5*(xi1D(i)+1.)
  xC(:,i,0)=xi*(xC (:,N,0)-xC (:,0,0))+xC(:,0,0)
  xC(:,i,N)=xi*(xC (:,N,N)-xC (:,0,N))+xC(:,0,N)
END DO
DO j=1,N-1
  eta=0.5*(xi1D(j)+1.)
  DO i=0,N
    xC(:,i,j)=eta*(xC (:,i,N)-xC (:,i,0))+xC(:,i,0)
  END DO
END DO
DO j=1,N-1
  eta=0.5*(xi1D(j)+1.)
  DO i=0,N
    xF1(:,i,j)=eta*(xF1(:,i,N)-xF1(:,i,0))+xF1(:,i,0)
  END DO
END DO
DO i=1,N-1
  xi=0.5*(xi1D(i)+1.)
  DO j=0,N
    xF2(:,i,j)=xi*(xF2(:,N,j)-xF2(:,0,j))+xF2(:,0,j)
  END DO
END DO
xFace=xF1+xF2-xC

END FUNCTION TransfiniteEdgeToFace


END MODULE MOD_Readin_SpecMesh2D


