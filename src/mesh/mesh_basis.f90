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
MODULE MOD_Mesh_Basis
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tEdge
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ElemGeometry
  MODULE PROCEDURE ElemGeometry
END INTERFACE

INTERFACE FindElemTypes
  MODULE PROCEDURE FindElemTypes
END INTERFACE

INTERFACE getNewHexa
  MODULE PROCEDURE getNewHexa
END INTERFACE

INTERFACE CreateSides
 MODULE PROCEDURE CreateSides
END INTERFACE

!INTERFACE AdjustOrientedNodes
  !MODULE PROCEDURE AdjustOrientedNodes
!END INTERFACE

INTERFACE GetBoundaryIndex
  MODULE PROCEDURE GetBoundaryIndex
END INTERFACE

INTERFACE BuildEdges
  MODULE PROCEDURE BuildEdges
END INTERFACE

INTERFACE FlushMesh
  MODULE PROCEDURE FlushMesh
END INTERFACE

INTERFACE assignBC
  MODULE PROCEDURE assignBC
END INTERFACE

INTERFACE isOriented
   MODULE PROCEDURE isOriented
END INTERFACE

INTERFACE getFlip
   MODULE PROCEDURE getFlip
END INTERFACE

INTERFACE PACKGEO
   MODULE PROCEDURE PACK1D
   MODULE PROCEDURE PACK2D
   MODULE PROCEDURE PACK3D
END INTERFACE

INTERFACE UNPACKGEO
   MODULE PROCEDURE UNPACK1D
   MODULE PROCEDURE UNPACK2D
   MODULE PROCEDURE UNPACK3D
END INTERFACE


PUBLIC::ElemGeometry
PUBLIC::getNewHexa
PUBLIC::CreateSides
!PUBLIC::AdjustOrientedNodes
PUBLIC::GetBoundaryIndex
PUBLIC::BuildEdges
PUBLIC::FlushMesh
PUBLIC::assignBC
PUBLIC::isOriented
PUBLIC::getFlip
PUBLIC::FindElemTypes
PUBLIC::PackGeo
PUBLIC::UnpackGeo

!===================================================================================================================================
  !---------------------------------------------------------------------------!
CONTAINS
SUBROUTINE ElemGeometry(Elem,TrafoOpt,TrafoInvOpt)
!===================================================================================================================================
! Compute the (linear) transformation for each element to its reference element.
! This linear transformation starts always from the first element node x_1!!
! Tranfsormation matrix (jacobi matrix)+Inverse of this mapping and Jacobi determinand of this mapping are computed.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tSidePtr,tBC,tNode,N,jacobianTolerance
USE MOD_Basis_Vars,ONLY:TetraMapInv,PrismMapInv,PyraMapInv,HexaMapInv
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(INOUT) :: Elem
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT),OPTIONAL      :: TrafoOpt(3,3)    ! see below
REAL,INTENT(OUT),OPTIONAL      :: TrafoInvOpt(3,3)   ! ? 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tSide),POINTER :: Side  ! ?
TYPE(tNode),POINTER :: p  ! ?
TYPE(tBC),POINTER   :: bc   ! ?
TYPE(tSidePtr)      :: sides(6)  ! ?
INTEGER             :: i,j,k   ! ?
INTEGER             :: CurvInd   ! ?
REAL,DIMENSION(3)   :: v1,v2,v3  ! ?
INTEGER,POINTER     :: map(:,:,:)  ! ?
REAL                :: Trafo(3,3)    ! element tranformation used in curved and nodal
REAL                :: TrafoInv(3,3) ! inverse element tranformation used in curved and nodal
!===================================================================================================================================

IF(Elem%nNodes .LT. 4)&
  CALL abort(__STAMP__,'Less than 4 Nodes for an Elem')
DO i=1,Elem%nNodes
  IF (.NOT. ASSOCIATED(Elem%Node(i)%np))&
    CALL abort(__STAMP__,'Elem node not existing in Elemgeometry')
END DO
IF(.NOT. ASSOCIATED(Elem%firstSide))&
  CALL abort(__STAMP__,'Elem has no sides in Elemgeometry')

!Determine linear mapping where 3 or 4 points are mapped to the unit triangle or tetrahedra
SELECT CASE(Elem%nNodes)
CASE(4) ! Tetraedron
  v1=Elem%Node(2)%np%x-Elem%Node(1)%np%x
  v2=Elem%Node(3)%np%x-Elem%Node(1)%np%x
  v3=Elem%Node(4)%np%x-Elem%Node(1)%np%x
CASE(6) ! Prism
  v1=Elem%Node(2)%np%x-Elem%Node(1)%np%x
  v2=Elem%Node(3)%np%x-Elem%Node(1)%np%x
  v3=Elem%Node(4)%np%x-Elem%Node(1)%np%x
CASE(5) ! Pyramid
  v1=Elem%Node(2)%np%x-Elem%Node(1)%np%x
  v2=Elem%Node(4)%np%x-Elem%Node(1)%np%x
  v3=Elem%Node(5)%np%x-Elem%Node(1)%np%x
CASE(8) ! Hexaedron
  v1=Elem%Node(2)%np%x-Elem%Node(1)%np%x
  v2=Elem%Node(4)%np%x-Elem%Node(1)%np%x
  v3=Elem%Node(5)%np%x-Elem%Node(1)%np%x
CASE DEFAULT
  v1=Elem%Node(2)%np%x-Elem%Node(1)%np%x
  v2=Elem%Node(Elem%firstSide%nNodes)%np%x-Elem%Node(1)%np%x
  v3=Elem%Node(Elem%nNodes)%np%x-Elem%Node(1)%np%x
END SELECT
Trafo(:,1)=v1
Trafo(:,2)=v2
Trafo(:,3)=v3
! compute also the inverse transformation
CALL INV33(Trafo,TrafoInv,Elem%detT)
IF(Elem%detT .LE. jacobianTolerance)THEN
  Trafo(:,1)=v2
  Trafo(:,2)=v1
  CALL INV33(Trafo,TrafoInv,Elem%detT)
  IF(Elem%detT .LE. jacobianTolerance)&
    CALL abort(__STAMP__,'Element with null-negative detT found! Elem%ind: ',Elem%ind)

  ! Element has wrong orientation rotate it into right-handed system
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    sides(Side%LocSide)%sp=>Side
    Side=>Side%nextElemSide
  END DO

  SELECT CASE(Elem%nNodes)
  CASE(4) ! Tetraedron
    p=>Elem%Node(2)%np; Elem%Node(2)%np=>Elem%Node(3)%np; Elem%Node(3)%np=>p
    bc=>sides(2)%sp%BC; sides(2)%sp%BC=>sides(4)%sp%BC  ; sides(4)%sp%BC=>bc
    CurvInd=sides(2)%sp%CurveIndex; sides(2)%sp%CurveIndex=sides(4)%sp%CurveIndex  ; sides(4)%sp%CurveIndex=CurvInd
    map=>TetraMapInv
  CASE(6) ! Prism
    p=>Elem%Node(2)%np; Elem%Node(2)%np=>Elem%Node(3)%np; Elem%Node(3)%np=>p
    p=>Elem%Node(5)%np; Elem%Node(5)%np=>Elem%Node(6)%np; Elem%Node(6)%np=>p
    bc=>sides(1)%sp%BC; sides(1)%sp%BC=>sides(3)%sp%BC  ; sides(3)%sp%BC=>bc
    CurvInd=sides(1)%sp%CurveIndex; sides(1)%sp%CurveIndex=sides(3)%sp%CurveIndex ; sides(3)%sp%CurveIndex=CurvInd
    map=>PrismMapInv
  CASE(5) ! Pyramid
    p=>Elem%Node(2)%np; Elem%Node(2)%np=>Elem%Node(4)%np; Elem%Node(4)%np=>p
    bc=>sides(2)%sp%BC; sides(2)%sp%BC=>sides(5)%sp%BC  ; sides(5)%sp%BC=>bc
    CurvInd=sides(2)%sp%CurveIndex; sides(2)%sp%CurveIndex=sides(5)%sp%CurveIndex  ; sides(5)%sp%CurveIndex=CurvInd
    bc=>sides(3)%sp%BC; sides(3)%sp%BC=>sides(4)%sp%BC  ; sides(4)%sp%BC=>bc
    CurvInd=sides(3)%sp%CurveIndex; sides(3)%sp%CurveIndex=sides(4)%sp%CurveIndex  ; sides(4)%sp%CurveIndex=CurvInd
    map=>PyraMapInv
  CASE(8) ! Hexaedron
    p=>Elem%Node(2)%np; Elem%Node(2)%np=>Elem%Node(4)%np; Elem%Node(4)%np=>p
    p=>Elem%Node(6)%np; Elem%Node(6)%np=>Elem%Node(8)%np; Elem%Node(8)%np=>p
    bc=>sides(2)%sp%BC; sides(2)%sp%BC=>sides(5)%sp%BC  ; sides(5)%sp%BC=>bc
    CurvInd=sides(2)%sp%CurveIndex; sides(2)%sp%CurveIndex=sides(5)%sp%CurveIndex  ; sides(5)%sp%CurveIndex=CurvInd
    bc=>sides(3)%sp%BC; sides(3)%sp%BC=>sides(4)%sp%BC  ; sides(4)%sp%BC=>bc
    CurvInd=sides(3)%sp%CurveIndex; sides(3)%sp%CurveIndex=sides(4)%sp%CurveIndex  ; sides(4)%sp%CurveIndex=CurvInd
    map=>HexaMapInv
  CASE DEFAULT
    CALL Abort(__STAMP__,'Unknown elem type. nNodes:',Elem%nNodes)
  END SELECT
  IF(ASSOCIATED(Elem%CurvedNode))THEN
    DO k=0,N; DO j=0,N; DO i=j+1,N
      IF(map(i,j,k).EQ.0) CYCLE
      p=>Elem%CurvedNode(map(i,j,k))%np
      Elem%CurvedNode(map(i,j,k))%np=>Elem%CurvedNode(map(j,i,k))%np
      Elem%CurvedNode(map(j,i,k))%np=>p
    END DO; END DO; END DO
  END IF
  CALL CreateSides(Elem,.FALSE.)
END IF
IF(PRESENT(TrafoOpt)) TrafoOpt=Trafo
IF(PRESENT(TrafoInvOpt)) TrafoInvOpt=TrafoInv

END SUBROUTINE elemGeometry

SUBROUTINE INV33(M,MInv,detM)
!===================================================================================================================================
! Computes the inverse of a 3x3 matrix
!===================================================================================================================================
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)     :: M(3,3)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: MInv(3,3),detM  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
detM =   M(1,1)*M(2,2)*M(3,3)  &
       - M(1,1)*M(2,3)*M(3,2)  &
       - M(1,2)*M(2,1)*M(3,3)  &
       + M(1,2)*M(2,3)*M(3,1)  &
       + M(1,3)*M(2,1)*M(3,2)  &
       - M(1,3)*M(2,2)*M(3,1)

IF(ABS(detM).LE.1.E-12)THEN
   MInv = 0.
   RETURN
END IF

MInv(1,1) =  (M(2,2)*M(3,3)-M(2,3)*M(3,2))
MInv(2,1) = -(M(2,1)*M(3,3)-M(2,3)*M(3,1))
MInv(3,1) =  (M(2,1)*M(3,2)-M(2,2)*M(3,1))
MInv(1,2) = -(M(1,2)*M(3,3)-M(1,3)*M(3,2))
MInv(2,2) =  (M(1,1)*M(3,3)-M(1,3)*M(3,1))
MInv(3,2) = -(M(1,1)*M(3,2)-M(1,2)*M(3,1))
MInv(1,3) =  (M(1,2)*M(2,3)-M(1,3)*M(2,2))
MInv(2,3) = -(M(1,1)*M(2,3)-M(1,3)*M(2,1))
MInv(3,3) =  (M(1,1)*M(2,2)-M(1,2)*M(2,1))
MInv=MInv/detM

END SUBROUTINE INV33

SUBROUTINE FindElemTypes()
!===================================================================================================================================
! determine the element type, for memory efficiency (involved tolerance).
! Requires existing linear mapping M1 for the Element "aElem" and mapped boundary splines
! 3/4             = triangle/quadrangle with linear sides
! 5               = bilinear quadrangle
! 6/7             = triangle/quadrangle with curved sides
! 104/105/106/108 = tetra/pyramid/prism/hexaeder  with linear sides
!     115/116/118 = pyramid/prism/hexaeder        with bilinear sides   
! 204/205/206/208 = tetra/pyramid/prism/hexaeder with curved sides
! >1000           = polygon or polyeder 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,FirstElem
USE MOD_Mesh_Tolerances,ONLY:SAMEPOINT
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER :: Elem     ! %nNodes,%Trafo,%node,  is used. %type is set
TYPE(tSide),POINTER :: Side     ! %curvedNode are checked, %iscurved is set =true if curved
LOGICAL             :: elemCurved,sideCurved  ! ?
REAL                :: Trafo(3,3)  ! transformation to reference element
!===================================================================================================================================
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  SELECT CASE(Elem%nNodes)
  CASE(4)
    Elem%Type=104
  CASE(5)
    Elem%Type=105
  CASE(6)
    Elem%Type=106
  CASE(8)
    Elem%Type=108
  CASE DEFAULT
    Elem%Type=1000+Elem%nNodes
  END SELECT
  
  elemCurved=.FALSE.
  IF(ASSOCIATED(Elem%CurvedNode)) elemCurved=.TRUE.
  
  sideCurved=.FALSE.
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(ASSOCIATED(Side%CurvedNode)) THEN
      Side%isCurved=.TRUE.
      sideCurved=.TRUE.
    END IF
    Side=>Side%nextElemSide
  END DO
  
  !IF(sideCurved.AND.(.NOT.elemCurved)) CALL Abort(__STAMP__,&
  !  'Sanity check: Element with uncurved volume but curved sides found. ElemInd:',Elem%Ind)
  !
  !IF(.NOT.sideCurved.AND.elemCurved)&
  !  WRITE(*,*) 'WARNING: Element is curved without any side beeing curved. ElemInd:',Elem%ind
  
  IF ((Elem%Type .EQ. 104) .AND. (sideCurved)) Elem%Type=204
  IF ((Elem%Type .EQ. 105) .AND. (sideCurved)) Elem%Type=205
  IF ((Elem%Type .EQ. 106) .AND. (sideCurved)) Elem%Type=206
  IF ((Elem%Type .EQ. 108) .AND. (sideCurved)) Elem%Type=208
  CALL elemGeometry(Elem,Trafo)
  SELECT CASE(Elem%Type)
  CASE(105)
    IF(.NOT.SAMEPOINT(Elem%Node(1)%np%x+MATMUL(Trafo,(/1.,1.,0./)),Elem%Node(3)%np%x)) Elem%Type=115
  CASE(106)
    IF(.NOT.SAMEPOINT(Elem%Node(1)%np%x+MATMUL(Trafo,(/1.,0.,1./)),Elem%Node(5)%np%x)) Elem%Type=116
    IF(.NOT.SAMEPOINT(Elem%Node(1)%np%x+MATMUL(Trafo,(/0.,1.,1./)),Elem%Node(6)%np%x)) Elem%Type=116
  CASE(108)
    IF(.NOT.SAMEPOINT(Elem%Node(1)%np%x+MATMUL(Trafo,(/0.,1.,1./)),Elem%Node(8)%np%x)) Elem%Type=118
    IF(.NOT.SAMEPOINT(Elem%Node(1)%np%x+MATMUL(Trafo,(/1.,1.,0./)),Elem%Node(3)%np%x)) Elem%Type=118
    IF(.NOT.SAMEPOINT(Elem%Node(1)%np%x+MATMUL(Trafo,(/1.,0.,1./)),Elem%Node(6)%np%x)) Elem%Type=118
    IF(.NOT.SAMEPOINT(Elem%Node(1)%np%x+MATMUL(Trafo,(/1.,1.,1./)),Elem%Node(7)%np%x)) Elem%Type=118
  END SELECT

  Elem=>Elem%nextElem
END DO

END SUBROUTINE findElemTypes

SUBROUTINE getNewHexa(Elem,zone,node1,node2,node3,node4,node5,node6,node7,node8)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:getNewElem
USE MOD_Mesh_Vars,ONLY:tElem,tNode
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(INOUT)            :: Elem  ! ?
INTEGER, INTENT(IN)            :: zone  ! ?
TYPE(tNode),POINTER,INTENT(IN)            :: node1   ! ?
TYPE(tNode),POINTER,INTENT(IN)            :: node2   ! ?
TYPE(tNode),POINTER,INTENT(IN)            :: node3   ! ?
TYPE(tNode),POINTER,INTENT(IN)            :: node4   ! ?
TYPE(tNode),POINTER,INTENT(IN)            :: node5   ! ?
TYPE(tNode),POINTER,INTENT(IN)            :: node6   ! ?
TYPE(tNode),POINTER,INTENT(IN)            :: node7   ! ?
TYPE(tNode),POINTER,INTENT(IN)            :: node8   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
  CALL GetNewElem(elem)
  Elem%zone=zone
  Elem%nNodes=8
  ALLOCATE(Elem%Node(8))
  Elem%Node(1)%np=>Node1
  Elem%Node(2)%np=>Node2
  Elem%Node(3)%np=>Node3
  Elem%Node(4)%np=>Node4
  Elem%Node(5)%np=>Node5
  Elem%Node(6)%np=>Node6
  Elem%Node(7)%np=>Node7
  Elem%Node(8)%np=>Node8
  CALL CreateSides(Elem,.TRUE.)
END SUBROUTINE getNewHexa


SUBROUTINE CreateSides(Elem,buildSides)
!===================================================================================================================================
! Creates the pointer Structure of the Element Sides, with respect to number of Side points in CGNS standard
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tNodePtr
USE MOD_Mesh_Vars,ONLY:MeshDim
USE MOD_Mesh_Vars,ONLY:DZ
USE MOD_Mesh_Vars,ONLY:nNodesElemSideMapping,ElemSideMapping
USE MOD_Mesh_Vars,ONLY:getNewSide,getNewNode
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(INOUT) :: Elem        ! pointer to Element
LOGICAL,INTENT(IN)                :: buildSides  ! determines if Sides should also be build 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tSide),POINTER :: Side  ! ?
TYPE(tNodePtr)      :: TempNodeArray(Elem%nNodes)  ! ?
INTEGER             :: iSide,iNode  ! ?
!===================================================================================================================================
IF((MeshDim .EQ. 2) .AND. buildSides)THEN
  ! Convert to 3D element
  DO iNode=1,Elem%nNodes
    TempNodeArray(iNode)%np => Elem%Node(iNode)%np
  END DO
  DEALLOCATE(Elem%Node)
  ALLOCATE(Elem%Node(2*Elem%nNodes))
  DO iNode=1,Elem%nNodes
    Elem%Node(iNode)%np     => TempNodeArray(iNode)%np
  END DO
  DO iNode=Elem%nNodes+1,2*Elem%nNodes
    CALL getNewNode(Elem%Node(iNode)%np,1)
    Elem%Node(iNode)%np%x   = Elem%Node(iNode-Elem%nNodes)%np%x + (/0.,0.,DZ/)
    ! Same node ind: this is needed to identify boundary sides. Unique node inds are created in fill25DMesh
    Elem%Node(iNode)%NP%Ind = Elem%Node(iNode-Elem%nNodes)%NP%Ind
  END DO
  Elem%nNodes=2*Elem%nNodes
END IF

IF((Elem%nNodes.LT.4) .OR. (Elem%nNodes.GT.8)) THEN
  WRITE(*,*)'Elem%nNodes: ',Elem%nNodes
  DO iNode=1, Elem%nNodes
    WRITE(*,*)'Node= ',iNode,' x= ',Elem%Node(iNode)%np%x
  END DO
  CALL abort(__STAMP__, 'Unknown 3D Elem%nNodes in CreateSides')
ENDIF

iSide=1
DO
  IF(iSide .EQ. 1) THEN
    IF(buildSides) CALL getNewSide(Elem%firstSide,nNodesElemSideMapping(Elem%nNodes,iSide))   !sidenNodesmapping(iSide))
    Side=>Elem%firstSide
  ELSE
    IF(buildSides)THEN
      CALL getNewSide(Side%nextElemSide,nNodesElemSideMapping(Elem%nNodes,iSide))   !sidenNodesmapping(iSide))
    END IF
    Side=>Side%nextElemSide
  END IF
  Side%LocSide=iSide
  Side%Elem=>Elem
  DO iNode=1,Side%nNodes
    IF(ElemSideMapping(Elem%nNodes,iSide,iNode).EQ.0) then
      WRITE(*,*)Elem%nNodes,iSide,iNode,Side%nNodes
      CALL abort(__STAMP__, &
        'Error in ElemSideMapping in CreateSides')
    END IF
    Side%Node(iNode)%np=>Elem%Node(ElemSideMapping(Elem%nNodes,iSide,iNode))%np
    IF(buildSides) Side%Node(iNode)%np%refCount=Side%Node(iNode)%np%refCount+1
  END DO
  iSide=iSide+1

  IF(iSide.GT.6) EXIT
  IF(nNodesElemSideMapping(Elem%nNodes,iSide).EQ.0) EXIT
END DO

END SUBROUTINE CreateSides


!SUBROUTINE AdjustOrientedNodes(Side,countRef)
!!==================================================================================================================================
!! Nodes of a Side and the nodes of its neighbor side have to be oriented in the
!!   same manner. Therefor we have the structure nodes (which is sorted according
!!   to CGNS standards in local element system) and the Oriented nodes which are
!!   equally oriented for Side and Side%connection
!!==================================================================================================================================
!! MODULES
!USE MOD_Mesh_Vars,ONLY:tSide
!USE MOD_Mesh_Vars,ONLY:VV
!USE MOD_Mesh_Tolerances,ONLY:SAMEPOINT
!! IMPLICIT VARIABLE HANDLING
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT VARIABLES
!! VV...periodic displacement vector for Cartmesh generator
!! CALC%spaceQuandt...used to secure geometric operations
!! RealTolerance is defined in hopr.h and used to secure geometric operations
!TYPE(tSide),POINTER,INTENT(IN) :: Side     ! pointer to the actual considered Side
!LOGICAL,INTENT(IN)             :: countRef ! determines if the Node%countref counter has to be updated
!!----------------------------------------------------------------------------------------------------------------------------------
!! OUTPUT VARIABLES
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES 
!TYPE(tSide),POINTER :: nSide   ! ?
!INTEGER             :: iNode,fNode,deriv(2)  ! ?
!REAL                :: VV_loc(3)  ! ?
!LOGICAL             :: isPeriodic,dominant   ! ?
!!==================================================================================================================================
!deriv=0
!IF (countRef) THEN
  !DO iNode=1,Side%nNodes
    !Side%Node(iNode)%np%refCount=Side%Node(iNode)%np%refCount+1
  !END DO
!END IF
!nSide=>Side%Connection
!IF(.NOT. ASSOCIATED(nSide)) THEN
  !DO iNode=1,Side%nNodes
    !Side%orientedNode(iNode)%np=>Side%Node(iNode)%np
  !END DO
!ELSE
  !isPeriodic=.FALSE.
  !IF(ASSOCIATED(Side%BC)) THEN
    !IF(.NOT. ASSOCIATED(nSide%BC)) &
      !CALL abort(__STAMP__, &
      !'periodic Side error')
    !IF(Side%BC%BCType.EQ.1) THEN
      !IF(nSide%BC%BCType.NE.1) &
        !CALL abort(__STAMP__, &
        !'periodic Side error')
      !isPeriodic=.TRUE.
      !VV_loc=VV(:,abs(Side%BC%BCalphaInd))*sign(1,Side%BC%BCalphaInd)
    !END IF
  !END IF
  !IF (countRef) THEN
    !DO iNode=1,nSide%nNodes
      !nSide%Node(iNode)%np%refCount=nSide%Node(iNode)%np%refCount+1
    !END DO
  !END IF
  !fNode=0
  !DO iNode=1,Side%nNodes
    !IF(isPeriodic) THEN
      !IF(SAMEPOINT(Side%Node(1)%np%x+VV_loc,nSide%Node(iNode)%np%x)) THEN
        !fNode=iNode
      !END IF
    !ELSE
      !IF(ASSOCIATED(Side%Node(1)%np,nSide%Node(iNode)%np)) THEN
        !fNode=iNode
      !END IF
    !END IF
  !END DO
  !dominant=.FALSE.
  !IF(Side%Elem%ind .LT. nSide%Elem%ind) dominant=.TRUE.
  !IF(dominant) THEN
    !DO iNode=1,nSide%nNodes
      !Side%orientedNode(iNode)%np=>Side%Node(iNode)%np
      !nSide%orientedNode(iNode)%np=>nSide%Node(fNode)%np
      !fNode=fNode-1
      !IF(fNode .LT.1) fNode=fNode+nSide%nNodes
    !END DO
  !ELSE
    !DO iNode=1,nSide%nNodes
      !Side%orientedNode(iNode)%np=>Side%Node(fNode)%np
      !nSide%orientedNode(iNode)%np=>nSide%Node(iNode)%np
      !fNode=fNode-1
      !IF(fNode .LT.1) fNode=fNode+nSide%nNodes
    !END DO
  !END IF
!END IF
!END SUBROUTINE AdjustOrientedNodes



FUNCTION GetBoundaryIndex(BCString)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:nUserDefinedBoundaries,BoundaryName
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)     :: BCString  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                         :: GetBoundaryIndex  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i  ! ?
CHARACTER(LEN=255)              :: BCName  ! ?
!===================================================================================================================================
GetBoundaryIndex=-1                                                                          ! Used as error flag
BCName=ADJUSTL(BCString)                                                                     ! Remove leading blanks
DO i=1,nUserDefinedBoundaries
  IF(INDEX(TRIM(BCName),TRIM(BoundaryName(i))) .GT. 0)THEN !  Found matching boundary name (can be only a part of the BCName)
    GetBoundaryIndex=i
    WRITE(UNIT_stdout,*)'BC associated:',TRIM(BCName),'-->',TRIM(BoundaryName(i))
    EXIT
  END IF
END DO
END FUNCTION GetBoundaryIndex

SUBROUTINE buildEdges()
!===================================================================================================================================
! Create Edge datastructure, each edge is unique, and has a pointer from each side and from the node with the lower index.
! on the node, a list beginning with node%firstEdge is build up. On the Element sides, a edge pointer array Edge(1:nNodes) is 
! filled, together with their orientation inside the side. Very important: OrientedNodes are used!!!! 
! If the edge is oriented, it goes from orientedNode(i)-> orientedNode(i+1), and 
! If the edge is not  oriented, it goes from orientedNode(i+1)-> orientedNode(i)
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tEdge,tNode,tEdgePtr
USE MOD_Mesh_Vars,ONLY:firstElem
USE MOD_Mesh_Vars,ONLY:GetNewEdge
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER          :: aElem  ! ?
TYPE(tSide),POINTER          :: aSide,bSide   ! ?
TYPE(tEdge),POINTER          :: aEdge,bEdge  ! ?
TYPE(tEdgePtr)               :: smallEdges(4)  ! ?
TYPE(tNode),POINTER          :: aNode,bNode  ! ?
INTEGER                      :: iSide,jSide,iEdge,jEdge,kEdge,iNode,iPlus,nSides,EdgeInd,nNodes  ! ?
INTEGER                      :: indA(2),indB(2,4),indTmp(2)
INTEGER                      :: edgeCount  ! ?
LOGICAL                      :: edgeFound  ! ?
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'BUILD EDGES ...'

! count unique corner nodes
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO iNode=1,aElem%nNodes
    aElem%Node(iNode)%np%tmp=0
  END DO
  aElem=>aElem%nextElem
END DO !! ELEMS!!
nNodes=0
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO iNode=1,aElem%nNodes
    IF(aElem%Node(iNode)%np%tmp.EQ.0) nNodes=nNodes+1
    aElem%Node(iNode)%np%tmp=nNodes
  END DO
  aElem=>aElem%nextElem
END DO !! ELEMS!!

EdgeInd=0
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  SELECT CASE(aElem%nNodes)
  CASE(8)
    nSides=6
  CASE(6)
    nSides=5
  CASE(5)
    nSides=5
  CASE(4)
    nSides=4
  END SELECT
  aSide=>aElem%firstSide
  DO iSide=1,nSides    !!SIDES!!***********
    DO iEdge=1,aSide%nNodes     !!EDGES!! nNodes=nEdges**************
      iPlus=iEdge+1
      IF(iEdge.EQ.aSide%nNodes) iPlus=1
      IF((aSide%OrientedNode(iEdge)%np%ind).GT.(aSide%OrientedNode(iPlus)%np%ind)) THEN
        aNode=>aSide%OrientedNode(iPlus)%np
        bNode=>aSide%OrientedNode(iEdge)%np
        aSide%edgeOrientation(iEdge)=.FALSE.
      ELSEIF((aSide%OrientedNode(iEdge)%np%ind).LT.(aSide%OrientedNode(iPlus)%np%ind))THEN
        aNode=>aSide%OrientedNode(iEdge)%np
        bNode=>aSide%OrientedNode(iPlus)%np
        aSide%edgeOrientation(iEdge)=.TRUE.
      ELSE 
        WRITE(*,*) 'Problem with node%ind in buildEdges'
        WRITE(*,*) 'node IDs',aSide%OrientedNode(iEdge)%np%ind,aSide%OrientedNode(iPlus)%np%ind
        WRITE(*,*) 'node1%x',aSide%OrientedNode(iEdge)%np%x
        WRITE(*,*) 'node2%x',aSide%OrientedNode(iPlus)%np%x
        
        STOP
      END IF

      edgeFound=.FALSE.
      aEdge=>aNode%firstEdge
      DO WHILE (ASSOCIATED(aEdge))
        IF (aEdge%Node(2)%np%ind .EQ. bNode%ind) THEN
          edgeFound=.TRUE.
          EXIT
        END IF
        aEdge=>aEdge%nextEdge
      END DO
      IF (.NOT.edgeFound) THEN
        CALL getNewEdge(aEdge,aNode,bNode)
        EdgeInd=EdgeInd+1
        IF (ASSOCIATED(aNode%firstEdge)) THEN
          aEdge%nextEdge=>aNode%firstEdge 
        END IF
        aNode%firstEdge=>aEdge
      END IF 
      !WRITE(*,*)'DEBUG a',aNode%ind,'b',bNode%ind,edgeFound,aEdge%ind
      aSide%Edge(iEdge)%edp=>aEdge
    END DO !!EDGES!!***************
    aSide=>aSide%nextElemSide
  END DO !!SIDES!!**************
  aElem=>aElem%nextElem
END DO !! ELEMS!!

! in case of nonconforming meshes, build nonconforming edge connectivity
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  SELECT CASE(aElem%nNodes)
  CASE(8)
    nSides=6
  CASE(6)
    nSides=5
  CASE(5)
    nSides=5
  CASE(4)
    nSides=4
  END SELECT
  aSide=>aElem%firstSide
  DO iSide=1,nSides
    IF(aSide%nMortars.LE.0)THEN  ! only check big mortar sides
      aSide=>aSide%nextElemSide
      CYCLE
    END IF
    IF(ASSOCIATED(aSide%BC))THEN ! dont build periodic curveds
      IF(aSide%BC%BCType.EQ.1)THEN
        print*,'Warning: Periodic mortar edges are treated like conforming edges!'
        aSide=>aSide%nextElemSide
        CYCLE
      END IF
    END IF

    ! check all edges of a side
    DO iEdge=1,aSide%nNodes
      aEdge=>aSide%Edge(iEdge)%edp
      IF(ASSOCIATED(aEdge%MortarEdge)) CYCLE ! already built

      ! search for all small edges that share at least one node with big edge
      ! 0-2 edges found: this should NEVER happen
      ! 3   edges found: the conforming case, occurring for the non-mortar edges of 2 -> 1 side mortars: ignore
      ! 4   edges found: the one we search for and the only case which should occur, if we're using 4 -> mortars.
      !                Two of these small edges belong to the mortar.
      ! >4  edges found: this should NEVER happen

      indA(1)=aEdge%Node(1)%np%ind
      indA(2)=aEdge%Node(2)%np%ind
      edgeCount=0
      DO jSide=1,aSide%nMortars
        bSide=>aSide%MortarSide(jSide)%sp
        DO jEdge=1,bSide%nNodes
          bEdge=>bSide%Edge(jEdge)%edp
          indTmp(1)=bEdge%Node(1)%np%ind
          indTmp(2)=bEdge%Node(2)%np%ind
          IF(ANY(indA(1).EQ.indTmp).OR.ANY(indA(2).EQ.indTmp))THEN
            edgeCount=edgeCount+1
            indB(:,edgeCount)=indTmp
            smallEdges(edgeCount)%edp=>bEdge
          END IF
        END DO
      END DO

      IF(edgeCount.EQ.3) CYCLE
      IF(edgeCount.NE.4) THEN
        STOP 'Mismatch of neighbour edge count of non-conforming edges.'
      END IF

      ! now select the 2 edges of the found 4 edges, sharing a common node
      DO jEdge=1,3
        DO kEdge=jEdge+1,4
          IF(ANY(indB(1,jEdge).EQ.indB(:,kEdge)).OR.ANY(indB(2,jEdge).EQ.indB(:,kEdge)))THEN
            IF(ANY(indA(1).EQ.indB(:,jEdge)).AND.ANY(indA(2).EQ.indB(:,kEdge))) THEN
              ALLOCATE(aEdge%MortarEdge(2))
              aEdge%MortarEdge(1)%edp=>smallEdges(jEdge)%edp
              aEdge%MortarEdge(2)%edp=>smallEdges(kEdge)%edp
              smallEdges(jEdge)%edp%parentEdge=>aEdge
              smallEdges(kEdge)%edp%parentEdge=>aEdge
            ELSEIF(ANY(indA(1).EQ.indB(:,kEdge)).AND.ANY(indA(2).EQ.indB(:,jEdge))) THEN
              ALLOCATE(aEdge%MortarEdge(2))
              aEdge%MortarEdge(2)%edp=>smallEdges(jEdge)%edp
              aEdge%MortarEdge(1)%edp=>smallEdges(kEdge)%edp
              smallEdges(jEdge)%edp%parentEdge=>aEdge
              smallEdges(kEdge)%edp%parentEdge=>aEdge
            END IF
          END IF
        END DO
      END DO
    END DO
    aSide=>aSide%nextElemSide
  END DO !!SIDES!!**************
  aElem=>aElem%nextElem
END DO !! ELEMS!!

CALL timer(.FALSE.)
END SUBROUTINE buildEdges


SUBROUTINE FlushMesh()
!===================================================================================================================================
! Delete mesh geometry
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem
USE MOD_Mesh_Vars,ONLY:firstElem
USE MOD_Mesh_Vars,ONLY:DeleteElem
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER :: Elem 
!===================================================================================================================================
DO WHILE(ASSOCIATED(FirstElem))
   Elem=>FirstElem
   CALL DeleteElem(FirstElem,Elem)  ! Delete elements
END DO
END SUBROUTINE flushMesh


SUBROUTINE assignBC(BCcopy,BCorig)
!===================================================================================================================================
! Makes a copy of boundary condition "BCorig"
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tBC
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tBC),POINTER,INTENT(IN) :: BCorig ! Original boundary condition
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tBC),POINTER,INTENT(OUT) :: BCcopy ! Copy of boundary condition "BCorig"
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF (ASSOCIATED(BCorig)) THEN
  IF (.NOT. ASSOCIATED(BCcopy)) THEN
    ALLOCATE(BCcopy)
  END IF
  BCcopy%BCtype     = BCorig%BCtype
  BCcopy%BCstate    = BCorig%BCstate
  BCcopy%BCalphaInd = BCorig%BCalphaInd
ELSE
  NULLIFY(BCcopy)  ! No boundary condition
END IF
END SUBROUTINE assignBC


FUNCTION isOriented(Side)
!===================================================================================================================================
! Check orientation of "Side" -> sign of normal vector
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tSide
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide),POINTER,INTENT(IN) :: Side         ! ? 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL             :: isOriented   ! ? 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
isOriented=.FALSE.
IF ((ASSOCIATED(Side%Node(1)%np, Side%orientedNode(1)%np)) .AND. &
    (ASSOCIATED(Side%Node(2)%np, Side%orientedNode(2)%np))) THEN
      isOriented=.TRUE.
END IF
END FUNCTION isOriented


FUNCTION getFlip(Side)
!===================================================================================================================================
! Get flip of side
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tSide
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide),POINTER,INTENT(IN) :: Side         ! ? 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                        :: getFlip   ! ? 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                        :: i
!===================================================================================================================================
getFlip=-1
IF(isOriented(Side))THEN !oriented side
  getFlip=0
ELSE !not oriented
  DO i=1,Side%nNodes
    IF(ASSOCIATED(Side%Node(i)%np,Side%OrientedNode(1)%np))THEN
      getFlip=i 
      EXIT
    END IF
  END DO
END IF
END FUNCTION getFlip


SUBROUTINE Pack1D(Ngeo,edge,data_out) 
!===================================================================================================================================
! description
!===================================================================================================================================
! MODULES  
USE MOD_Mesh_Vars,ONLY:tEdge
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)             :: Ngeo
TYPE(tEdge),POINTER,INTENT(IN) :: edge
REAL,INTENT(OUT)               :: data_out(3,0:Ngeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNgeo 
!===================================================================================================================================
IF(NGeo.GT.1)THEN
  DO iNgeo=0,Ngeo
    data_out(:,iNgeo) = edge%curvedNode(iNgeo+1)%np%x
  END DO 
ELSE 
  data_out(:,0) = edge%Node(1)%np%x
  data_out(:,1) = edge%Node(2)%np%x
END IF
END SUBROUTINE Pack1D

SUBROUTINE Pack2D(Ngeo,side,data_out) 
!===================================================================================================================================
! description
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Mesh_Vars,ONLY:tSide
USE MOD_Basis_Vars ,ONLY:QuadMapInv
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: Ngeo
TYPE(tSide),POINTER,INTENT(IN) :: side
REAL,INTENT(OUT)    :: data_out(3,0:Ngeo,0:Ngeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNgeo,jNgeo,i1D 
!===================================================================================================================================
IF(Ngeo.GT.1)THEN
  DO jNgeo=0,Ngeo
    DO iNgeo=0,Ngeo
      i1D = QuadMapInv(iNgeo,jNgeo)
      data_out(:,iNgeo,jNgeo) = side%curvedNode(i1D)%NP%x
    END DO 
  END DO 
ELSE
  data_out(:,0,0) = side%OrientedNode(1)%NP%x
  data_out(:,1,0) = side%OrientedNode(2)%NP%x
  data_out(:,1,1) = side%OrientedNode(3)%NP%x
  data_out(:,0,1) = side%OrientedNode(4)%NP%x
END IF
END SUBROUTINE Pack2D

SUBROUTINE Pack3D(Ngeo,elem,data_out) 
!===================================================================================================================================
! description
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Mesh_Vars,ONLY:tElem
USE MOD_Basis_Vars ,ONLY:HexaMapInv
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: Ngeo
TYPE(tElem),POINTER,INTENT(IN) :: elem
REAL,INTENT(OUT)    :: data_out(3,0:Ngeo,0:Ngeo,0:Ngeo)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNgeo,jNgeo,kNgeo,i1D 
!===================================================================================================================================
IF(Ngeo.GT.1)THEN
  DO kNgeo=0,Ngeo
    DO jNgeo=0,Ngeo
      DO iNgeo=0,Ngeo
        i1D = HexaMapInv(iNgeo,jNgeo,kNgeo)
        data_out(:,iNgeo,jNgeo,kNgeo) = elem%curvedNode(i1D)%NP%x
      END DO 
    END DO 
  END DO 
ELSE
  data_out(:,0,0,0) = elem%Node(1)%NP%x
  data_out(:,1,0,0) = elem%Node(2)%NP%x
  data_out(:,1,1,0) = elem%Node(3)%NP%x
  data_out(:,0,1,0) = elem%Node(4)%NP%x
  data_out(:,0,0,1) = elem%Node(5)%NP%x
  data_out(:,1,0,1) = elem%Node(6)%NP%x
  data_out(:,1,1,1) = elem%Node(7)%NP%x
  data_out(:,0,1,1) = elem%Node(8)%NP%x
END IF
END SUBROUTINE Pack3D

SUBROUTINE Unpack1D(Ngeo,data_in,edge) 
!===================================================================================================================================
! description
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Mesh_Vars,ONLY:tEdge
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: Ngeo
REAL,INTENT(IN)     :: data_in(3,0:Ngeo)
TYPE(tEdge),POINTER,INTENT(IN) :: edge
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNgeo 
!===================================================================================================================================
IF(NGeo.GT.1)THEN
  DO iNgeo=0,Ngeo
    edge%curvedNode(iNgeo+1)%NP%x = data_in(:,iNgeo)
  END DO 
ELSE
  edge%Node(1)%NP%x = data_in(:,0)
  edge%Node(2)%NP%x = data_in(:,1)
END IF
END SUBROUTINE Unpack1D

SUBROUTINE Unpack2D(Ngeo,data_in,side) 
!===================================================================================================================================
! description
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Mesh_Vars,ONLY:tSide
USE MOD_Basis_Vars ,ONLY:QuadMapInv
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: Ngeo
REAL,INTENT(IN)     :: data_in(3,0:Ngeo,0:Ngeo)
TYPE(tSide),POINTER,INTENT(IN) :: side
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNgeo,jNgeo,i1D
!===================================================================================================================================
IF(NGeo.GT.1)THEN
  DO jNgeo=0,Ngeo
    DO iNgeo=0,Ngeo
      i1D = QuadMapInv(iNgeo,jNgeo)
      side%curvedNode(i1D)%NP%x = data_in(:,iNgeo,jNgeo)
    END DO 
  END DO 
ELSE
  side%OrientedNode(1)%NP%x = data_in(:,0,0)
  side%OrientedNode(2)%NP%x = data_in(:,1,0)
  side%OrientedNode(3)%NP%x = data_in(:,1,1)
  side%OrientedNode(4)%NP%x = data_in(:,0,1)
END IF
END SUBROUTINE Unpack2D

SUBROUTINE Unpack3D(Ngeo,data_in,elem) 
!===================================================================================================================================
! description
!===================================================================================================================================
! MODULES                                                                                                                          !
USE MOD_Mesh_Vars,ONLY:tElem
USE MOD_Basis_Vars ,ONLY:HexaMapInv
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------!
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN)  :: Ngeo
REAL,INTENT(IN)     :: data_in(3,0:Ngeo,0:Ngeo,0:Ngeo)
TYPE(tElem),POINTER,INTENT(IN) :: elem
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iNgeo,jNgeo,kNgeo,i1D
!===================================================================================================================================
IF(NGeo.GT.1)THEN
  DO kNgeo=0,Ngeo
    DO jNgeo=0,Ngeo
      DO iNgeo=0,Ngeo
        i1D = HexaMapInv(iNgeo,jNgeo,kNgeo)
        elem%curvedNode(i1D)%NP%x = data_in(:,iNgeo,jNgeo,kNgeo)
      END DO 
    END DO 
  END DO 
ELSE
  elem%Node(1)%NP%x = data_in(:,0,0,0)
  elem%Node(2)%NP%x = data_in(:,1,0,0)
  elem%Node(3)%NP%x = data_in(:,1,1,0)
  elem%Node(4)%NP%x = data_in(:,0,1,0)
  elem%Node(5)%NP%x = data_in(:,0,0,1)
  elem%Node(6)%NP%x = data_in(:,1,0,1)
  elem%Node(7)%NP%x = data_in(:,1,1,1)
  elem%Node(8)%NP%x = data_in(:,0,1,1)
END IF
END SUBROUTINE Unpack3D

END MODULE MOD_Mesh_Basis
