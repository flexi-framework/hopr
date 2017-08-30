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
USE MOD_Globals
USE MOD_Mesh_Vars     ,ONLY:tElem,tElemPtr,tSide,tNode,tNodePtr
USE MOD_Mesh_Vars     ,ONLY:FirstElem
USE MOD_Mesh_Vars     ,ONLY:n2dNodes
USE MOD_Mesh_Vars     ,ONLY:nMeshFiles,MeshFileName
USE MOD_Mesh_Vars     ,ONLY:nUserDefinedBoundaries,BoundaryName,BoundaryType
USE MOD_Mesh_Vars     ,ONLY:N,DZ,useCurveds
USE MOD_Mesh_Vars     ,ONLY:getNewNode,getNewBC
USE MOD_ChangeBasis   ,ONLY:ChangeBasis2D
USE MOD_CurvedCartMesh,ONLY:getNewCurvedHexahedron
USE MOD_CartMesh      ,ONLY:getNewHexahedron
USE MOD_Mesh_Basis    ,ONLY:createSides
USE MOD_Basis1D       ,ONLY:ChebyGaussLobNodesAndWeights
USE MOD_Basis1D       ,ONLY:BarycentricWeights
USE MOD_Basis1D       ,ONLY:InitializeVandermonde
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
INTEGER                         :: os,iNode,iBC,iFile,iElem,i,j,k
INTEGER                         :: globalNodeID 
INTEGER                         :: fUnit = 108
INTEGER                         :: nodeIDs(4)
INTEGER                         :: curveFlags(4)
INTEGER                         :: BC_ID(4)
INTEGER                         :: nElems, nNodes, Nspec
INTEGER                         :: EdgeMap(2,4)
INTEGER                         :: countBCs(nUserDefinedBoundaries)
INTEGER                         :: Nanalyze
REAL                            :: xStart(2), xEnd(2)
REAL                            :: maxDistErr
REAL                            :: v1(3),v2(3)
LOGICAL                         :: doAnalyze
REAL,ALLOCATABLE                :: SpecEdges(:,:,:),xCLedges(:,:,:)
REAL,ALLOCATABLE                :: xCLelem(:,:,:)
REAL,ALLOCATABLE                :: xElem(:,:,:)
REAL,ALLOCATABLE                :: xiCLNspec(:),wBaryCLNspec(:)
REAL,ALLOCATABLE                :: xiCLN(:),wBaryCLN(:),xiEqN(:),xiAnalyze(:)
REAL,ALLOCATABLE                :: Vdm_CLNspec_CLN(:,:),Vdm_CLN_EqN(:,:)
REAL,ALLOCATABLE                :: Vdm_CLNspec_analyze(:,:),Vdm_CLN_Analyze(:,:)
CHARACTER(LEN=40)               :: BCnames(4)
TYPE(tNodePtr),POINTER          :: CurvedNode(:,:,:) 
TYPE(tNodePtr),POINTER          :: CornerNode(:)
!==================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')'Reading SpecMesh 2D mesh...'

EdgeMap(1:2,1)=(/1,2/)
EdgeMap(1:2,2)=(/2,3/)
EdgeMap(1:2,3)=(/4,3/)
EdgeMap(1:2,4)=(/1,4/)
ALLOCATE(CornerNode(8))
NULLIFY(Elem)

maxDistErr=0.
doAnalyze=.FALSE.

globalNodeID=0
countBCs=0

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
  ! Get number of nodes, elements and degree of the boundary polynomials
  READ(fUnit,*) nNodes, nElems, Nspec

  !ALLOCATE arrays
  ALLOCATE(nodeCoords(2,nNodes))
  ALLOCATE( SpecEdges(2,0:Nspec,4) )

  IF(.NOT.useCurveds)THEN
      WRITE(UNIT_stdOut,*)'WARNING, useCurveds=F => we will only read linear elements!'
  ELSE
    IF(Nspec.NE.N) THEN
      WRITE(UNIT_stdOut,*)'WARNING, the SpecMesh geometry will be interpolated!'
      WRITE(UNIT_stdOut,*)'  ==> from degree ', Nspec, 'to given Ngeo= ',N
    END IF
    ALLOCATE( xCLEdges(2,0:N,4) )
    ALLOCATE( xCLelem(2,0:N,0:N) )
    ALLOCATE( xelem(2,0:N,0:N) )
    ALLOCATE( CurvedNode(0:N,0:N,0:N) )
    ALLOCATE( xiCLN  (0:N) )
    ALLOCATE( wBaryCLN(0:N) )
    ALLOCATE( Vdm_CLN_EQN(0:N,0:N) )
    
    CALL ChebyGaussLobNodesAndWeights(N,xiCLN)
    CALL BarycentricWeights(N,xiCLN,wBaryCLN)
    
    !change to equidistant points (used throughout hopr... :-( ) 
    ALLOCATE( xiEqN  (0:N) )
    DO j = 0, N 
       xiEqN(j) = 2.*REAL(j)/REAL(N)-1.
    END DO !j=0,N
    CALL InitializeVandermonde(N,N,wBaryCLN,xiCLN,xiEqN,Vdm_CLN_EqN)
    DEALLOCATE( xiEqN)
    
    !change of degree from CL Nspec => CL N
    IF(Nspec.NE.N)THEN
    
      Nanalyze=2*Nspec
      ALLOCATE( Vdm_CLNspec_CLN(0:N,0:Nspec) )
      ALLOCATE( Vdm_CLNspec_analyze(0:Nanalyze,0:Nspec))
      ALLOCATE( Vdm_CLN_analyze(0:Nanalyze,0:N))
      ALLOCATE( xiCLNspec  (0:Nspec) )
      ALLOCATE( wBaryCLNspec(0:Nspec) )
      ALLOCATE( xiAnalyze (0:Nanalyze) ) !for analyze
    
      CALL ChebyGaussLobNodesAndWeights(Nspec,xiCLNspec)
      CALL BarycentricWeights(Nspec,xiCLNspec,wBaryCLNspec)
    
      CALL InitializeVandermonde(Nspec,N,wBaryCLNspec,xiCLNspec,xiCLN,Vdm_CLNspec_CLN)
      !analyze interpolation error 
      IF(.NOT.doAnalyze) doAnalyze=(Nspec.GT.N)

      DO j = 0, Nanalyze 
         xiAnalyze(j) = 2.*REAL(j)/REAL(Nanalyze)-1.
      END DO !j=0,N
      CALL InitializeVandermonde(Nspec,Nanalyze,wBaryCLNspec,xiCLNspec,xiAnalyze,Vdm_CLNspec_Analyze)
      CALL InitializeVandermonde(    N,Nanalyze,wBaryCLN    ,xiCLN    ,xiAnalyze,Vdm_CLN_Analyze    )
    
      DEALLOCATE(xiAnalyze)
      DEALLOCATE(wBaryCLNspec)
    END IF !Nspec_old .NE. Nspec
  END IF !useCurved

  ! Read node coordinates

  nodeCoords=0.
  DO iNode = 1, nNodes   ! READ nodes
    READ(funit,*)nodeCoords(1:2,iNode)
  END DO



  DO iElem = 1, nElems 
    ! -----------------------------------------
    ! Read element: Elements have the format
    ! node1 node2 node3 node4
    ! b1 b2 b3 b4
    ! (=0 for straight side, 1 for curved)
    ! if curved boundaries, then for each:
    ! for j = 0 to boundaryOrder
    !    x_j  y_j
    ! next j
    ! bname1 bname2 bname3 bname4
    ! -----------------------------------------
    READ( fUnit, * ) nodeIDs
    READ( fUnit, * ) curveFlags

    !check first element orientation
    IF(iElem.EQ.1)THEN
      v1(1:2)=NodeCoords(1:2,nodeIDs(2))-NodeCoords(1:2,nodeIDs(1))
      v1(3)=0.
      v2(1:2)=NodeCoords(1:2,nodeIDs(4))-NodeCoords(1:2,nodeIDs(1))
      v2(3)=0.
      IF(SUM(CROSS(v1,v2)).LT.0.) THEN
        WRITE(UNIT_stdOut,*)'elements are oriented reversly, not implemented!',SUM(CROSS(v1,v2))
      END IF
    END IF !iElem=1
       
    SpecEdges=0.
    DO k=1,4 
      xStart = NodeCoords(1:2,nodeIDs(edgeMap(1,k)))
      xEnd   = NodeCoords(1:2,nodeIDs(edgeMap(2,k)))
      IF( curveFlags(k).NE.0)THEN
        DO j=0,Nspec
          READ(fUnit,*) SpecEdges(:,j,k)
        END DO !j=0,Nspec
        !WRITE(*,*)'DEBUG,check corner 0',k,SpecEdges(:,0,k),xStart,SUM((SpecEdges(:,0,k)-xStart)**2)
        !WRITE(*,*)'DEBUG,check corner N',k,SpecEdges(:,Nspec,k),xEnd,SUM((SpecEdges(:,Nspec,k)-xEnd)**2)
        !unique nodes at corners
        SpecEdges(:,0,k)=xStart
        SpecEdges(:,Nspec,k)=xEnd
      END IF
    END DO !k=1,4
    READ( fUnit, * ) BCnames(1:4)
    ! -----------------------------------------

    IF(useCurveds)THEN
      DO k=1,4 
        xCLedges(:,0,k) = NodeCoords(1:2,nodeIDs(edgeMap(1,k)))
        xCLedges(:,N,k) = NodeCoords(1:2,nodeIDs(edgeMap(2,k)))
        IF ( curveFlags(k) .EQ. 0 )     THEN !linear edge
          DO j=1,N-1
            xCLedges(:,j,k)=0.5*( (1.+xiCLN(j))*xCLedges(:,N,k)+ (1.-xiCLN(j))*xCLedges(:,0,k) ) !linear interpol. 
          END DO
        ELSE !curved edge, interpolate from Nspec => N
          IF(Nspec.NE.N)THEN
            xCLedges(1,1:N-1,k)=MATMUL(Vdm_CLNspec_CLN(1:N-1,:),SpecEdges(1,:,k))
            xCLedges(2,1:N-1,k)=MATMUL(Vdm_CLNspec_CLN(1:N-1,:),SpecEdges(2,:,k))
          ELSE
            xCLedges(:,1:N-1,k)=SpecEdges(:,1:Nspec-1,k)
          END IF
          !analyze interpolation error
          IF(Nspec.GT.N)THEN
            maxDistErr = MAX(maxDistErr, MAXVAL(                                                                           & 
                         SQRT( (MATMUL(Vdm_CLNspec_Analyze,SpecEdges(1,:,k))-MATMUL(Vdm_CLN_Analyze,xCLedges(1,:,k)))**2   &
                              +(MATMUL(Vdm_CLNspec_Analyze,SpecEdges(2,:,k))-MATMUL(Vdm_CLN_Analyze,xCLedges(2,:,k)))**2) )) 
          END IF !Nspec>N
        END IF !linear/curved edge
        !SAVE WAY
      END DO !k=1,4
      xCLelem= TransfiniteEdgeToFace(N,xCLedges,xiCLN)
      CALL ChangeBasis2D(2,N,N,Vdm_CLN_EqN,xCLelem,xElem)
     
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
      CALL GetNewCurvedHexahedron(CurvedNode,N,iFile) !Zone=iFile

    ELSE  !not useCurveds
      DO iNode=1,4
        NULLIFY(CornerNode(iNode)%np)
        CALL getNewNode(Node,1)
        globalNodeID=globalNodeID+1
        Node%ind=globalNodeID
        Node%x(1:2)=NodeCoords(1:2,nodeIDs(iNode))
        Node%x(3)  =0.
        CornerNode(iNode)%np=>Node
      END DO !inode=1,4
      DO iNode=5,8
        NULLIFY(CornerNode(iNode)%np)
        CALL getNewNode(Node,1)
        globalNodeID=globalNodeID+1
        Node%ind=globalNodeID
        Node%x(1:2)=CornerNode(iNode-4)%np%x(1:2)
        Node%x(3)  =DZ
        CornerNode(iNode)%np=>Node
      END DO !iNode=5,8
      CALL GetNewHexahedron(CornerNode)
      CALL CreateSides(FirstElem,.TRUE.)
    END IF !useCurveds
    !new element added always as the firstElem
    Elem=>FirstElem
    Elem%ind=-1 
  
    BC_ID=-1
    DO k=1,4
      IF(INDEX(TRIM(BCNames(k)),'---').NE.0)THEN
        ! side has no BC, do noting
        BC_ID(k)=0
      ELSE
        DO iBC=1,nUserDefinedBoundaries
          IF(INDEX(TRIM(ADJUSTL(BCnames(k))),TRIM(BoundaryName(iBC))) .GT. 0)THEN !  Found matching boundary name (can be only a part of the BCName)
            BC_ID(k)=iBC
            countBCs(iBC) = countBCs(iBC)+1
            IF(countBCs(iBC).EQ.1)THEN !information for first time found!
              WRITE(UNIT_stdout,*)'BC associated:',TRIM(BCNames(k)),'-->',TRIM(BoundaryName(iBC))
            END IF
            EXIT
          END IF
        END DO
        IF(BC_ID(k).EQ.-1) THEN
          WRITE(UNIT_stdOut,*)'ERROR - Could not find corresponding boundary definition of ',TRIM(BCnames(k))
          !STOP
        END IF
      END IF
    END DO !k=1,4
    !ASSOCIATE BC TO Side
    Side=>Elem%FirstSide !locSide 1
    Side=>Side%nextElemSide !go to locside 2
    DO k=1,4
      !associate locside 2...5 to BC_ID 1...4
      IF(BC_ID(k).GT.0)THEN
        iBC=BC_ID(k)
        CALL getNewBC(Side%BC)
        Side%BC%BCType    =BoundaryType(iBC,1)
        Side%CurveIndex   =BoundaryType(iBC,2)
        Side%BC%BCstate   =BoundaryType(iBC,3)
        Side%BC%BCalphaInd=BoundaryType(iBC,4)
        Side%BC%BCIndex   =iBC
      END IF !BC_ID(k)>0 (BC found)
      Side=>Side%nextElemSide
    END DO
  END DO !iElem
     
  CLOSE(fUnit)
  DEALLOCATE(NodeCoords)
  DEALLOCATE( SpecEdges )
  IF(useCurveds)THEN
    DEALLOCATE( xCLEdges )
    DEALLOCATE( xCLelem )
    DEALLOCATE( xelem )
    DO k=0,N; DO j=0,N; DO i=0,N
      NULLIFY(CurvedNode(i,j,k)%np)
    END DO; END DO; END DO !i,j,k
    DEALLOCATE( CurvedNode )
    DEALLOCATE( xiCLN )
    DEALLOCATE( wBaryCLN )
    DEALLOCATE( Vdm_CLN_EQN )
    IF(Nspec.NE.N)THEN
      DEALLOCATE( Vdm_CLNspec_CLN )
      DEALLOCATE( Vdm_CLNspec_analyze )
      DEALLOCATE( Vdm_CLN_analyze )
      DEALLOCATE( xiCLNspec )
    END IF !Nspec /= N
  END IF!useCurveds
END DO !iFile


n2DNodes=globalNodeID !total number of points

DO iNode=1,8
  NULLIFY(CornerNode(iNode)%np)
END DO!iNode=1,8
DEALLOCATE( CornerNode )

IF(doAnalyze)THEN
  WRITE(UNIT_stdOut,*) 'maximum error (distance) from interpolation Nspec => N: ',maxDistErr
END IF
CALL Timer(.FALSE.)
WRITE(UNIT_stdOut,'(132("~"))')
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
INTEGER                   :: i,j
REAL,DIMENSION(2,0:N,0:N) :: xC,xF1,xF2
REAL                      :: xi,eta
!==================================================================================================================================
xC(:,0,0)=xEdges(:,0,1)
xC(:,N,0)=xEdges(:,N,1)
xC(:,0,N)=xEdges(:,0,3)
xC(:,N,N)=xEdges(:,N,3)

xF1(:,  :,  0)=xEdges(:,:,1)
xF1(:,  :,  N)=xEdges(:,:,3)

xF2(:,  0,  :)=xEdges(:,:,4)
xF2(:,  N,  :)=xEdges(:,:,2)

!linear corner interpolation
DO i=1,N-1
  xi=0.5*(xi1D(i)+1.)
  xC(:,i,0)=xi*xC (:,N,0)+(1.-xi)*xC (:,0,0)
  xC(:,i,N)=xi*xC (:,N,N)+(1.-xi)*xC (:,0,N)
END DO
DO j=1,N-1
  eta=0.5*(xi1D(j)+1.)
  DO i=0,N
    xC(:,i,j)=eta*xC (:,i,N)+(1.-eta)*xC (:,i,0)
  END DO !i=0,N
END DO! j=1,N-1

DO j=1,N-1
  eta=0.5*(xi1D(j)+1.)
  DO i=0,N
    xF1(:,i,j)=eta*xF1(:,i,N)+(1.-eta)*xF1(:,i,0)
  END DO
END DO

DO i=1,N-1
  xi=0.5*(xi1D(i)+1.)
  DO j=0,N
    xF2(:,i,j)=xi*xF2(:,N,j)+(1.-xi)*xF2(:,0,j)
  END DO
END DO
xFace=xF1+xF2-xC

END FUNCTION TransfiniteEdgeToFace


END MODULE MOD_Readin_SpecMesh2D


