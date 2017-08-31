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
MODULE MOD_CartMesh
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
TYPE tCartesianMesh ! provides data structure for Cartesian mesh
  REAL                                ::       Corner(3,8)            ! Corner nodes of zone
  REAL                                ::       l0(3)                  ! first length (+/-) = direction
  REAL                                ::       factor(3)              ! strech factor (+/-) = direction
  INTEGER                             ::       BCIndex(6)             ! index in UserdefinedBoundaries
  INTEGER                             ::       elemType               ! Element type in zone
  INTEGER                             ::       nElems(3)              ! number of elements in zone
  INTEGER                             ::       meshTemplate           ! defines which template is used to build the cartmesh
 ! especially for tetras
END TYPE tCartesianMesh

TYPE tCartesianMeshPtr
  TYPE(tCartesianMesh),POINTER        ::       CM                     ! cartesian mesh pointer
END TYPE tCartesianMeshPtr

! Public Part ----------------------------------------------------------------------------------------------------------------------
TYPE(tCartesianMeshPtr),POINTER       ::       CartMeshes(:)          ! cartMeshes(nZones) pointer to Cartesian mesh information

PUBLIC::tCartesianMesh,CartMeshes,tCartesianMeshPtr

INTERFACE CartesianMesh
  MODULE PROCEDURE CartesianMesh
END INTERFACE

INTERFACE GetNewHexahedron
  MODULE PROCEDURE GetNewHexahedron
END INTERFACE

PUBLIC::CartesianMesh, GetNewHexahedron
!===================================================================================================================================

CONTAINS
SUBROUTINE GetNewTetrahedron(CornerNode,cartMesh,l,m,n,ind)
!===================================================================================================================================
! Build new tetrahedron for cartesian mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewNode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNodePtr),INTENT(IN)       :: CornerNode(8)  ! corner nodes of box
TYPE(tCartesianMesh),INTENT(IN) :: cartMesh       ! cartmesh in use
INTEGER,INTENT(IN)              :: l,m,n          ! position in mesh
INTEGER,INTENT(INOUT)           :: ind            ! max node ind
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: aElem   ! ?
INTEGER                         :: tetMap(4,24)
INTEGER                         :: i,nTets
TYPE(tNodePtr)                  :: node(14)
!===================================================================================================================================
! Three different tet strategies are implemented
! 1. strategy: 6 tets per box, all tets have same volume and angle, not periodic but isotropic
! 2. strategy: 6 tets per box, split hex into two prisms and each prism into 3 tets, periodic but strongly anisotropic
! 3. strategy: 5 tets per box (minimum number of tets), 4 small tets at the edges one big in the middle, not periodic not isotropic
! 4. strategy: 24 tets per box: one tet per edge (12), 4 in the center and 8 at the edges, fully symmetric and periodic mesh
! 5. strategy: 12 tets per box: one pyramid per face, each of them splitted into to tets
! Corner node ordering:
! 1. (-,-,-)
! 2. (+,-,-)
! 3. (+,+,-)
! 4. (-,+,-)
! 5. (-,-,+)
! 6. (+,-,+)
! 7. (+,+,+)
! 8. (-,+,+)

DO i=1,8
  node(i)%np=>CornerNode(i)%np
END DO

SELECT CASE(cartmesh%meshTemplate)
CASE(1)
  nTets=6
  tetMap(:,1)=(/1,3,4,5/)
  tetMap(:,2)=(/1,2,3,5/)
  tetMap(:,3)=(/3,5,7,8/)
  tetMap(:,4)=(/3,5,6,7/)
  tetMap(:,5)=(/2,3,5,6/)
  tetMap(:,6)=(/3,4,5,8/)
CASE(2)
  nTets=6
  tetMap(:,1)=(/1,2,4,5/)
  tetMap(:,2)=(/2,5,6,8/)
  tetMap(:,3)=(/2,4,5,8/)
  tetMap(:,4)=(/2,3,6,8/)
  tetMap(:,5)=(/2,3,4,8/)
  tetMap(:,6)=(/3,6,7,8/)
CASE(3)
  nTets=24
  ! center nodes of the box faces are required
  ! 9.  (-,o,o) XI_MINUS
  ! 10. (+,o,o) XI_PLUS
  ! 11. (o,-,o) ETA_MINUS
  ! 12. (o,+,o) ETA_PLUS
  ! 13. (o,o,-) ZETA_MINUS
  ! 14. (o,o,+) ZETA_PLUS

  CALL getNewNode(node(9 )%np,ind=ind)
  CALL getNewNode(node(10)%np,ind=ind)
  CALL getNewNode(node(11)%np,ind=ind)
  CALL getNewNode(node(12)%np,ind=ind)
  CALL getNewNode(node(13)%np,ind=ind)
  CALL getNewNode(node(14)%np,ind=ind)
  node(9 )%np%x=(node(1)%np%x+node(4)%np%x+node(5)%np%x+node(8)%np%x)/4.
  node(10)%np%x=(node(2)%np%x+node(3)%np%x+node(6)%np%x+node(7)%np%x)/4.
  node(11)%np%x=(node(1)%np%x+node(2)%np%x+node(5)%np%x+node(6)%np%x)/4.
  node(12)%np%x=(node(3)%np%x+node(4)%np%x+node(7)%np%x+node(8)%np%x)/4.
  node(13)%np%x=(node(1)%np%x+node(2)%np%x+node(3)%np%x+node(4)%np%x)/4.
  node(14)%np%x=(node(5)%np%x+node(6)%np%x+node(7)%np%x+node(8)%np%x)/4.
  IF(l.EQ.1)                  node(9 )%np%tmp=50000  !xi minus
  IF(l.EQ.cartmesh%nElems(1)) node(10)%np%tmp=300    !xi plus
  IF(m.EQ.1)                  node(11)%np%tmp=20     !eta minus
  IF(m.EQ.cartmesh%nElems(2)) node(12)%np%tmp=4000   !eta plus
  IF(n.EQ.1)                  node(13)%np%tmp=1      !zeta minus
  IF(n.EQ.cartmesh%nElems(3)) node(14)%np%tmp=600000 !zeta plus

  ! first the edges
  tetMap(:,1 )=(/1,2,13,11/)
  tetMap(:,2 )=(/2,3,13,10/)
  tetMap(:,3 )=(/3,4,13,12/)
  tetMap(:,4 )=(/4,1,13,9/)
  tetMap(:,5 )=(/5,6,14,11/)
  tetMap(:,6 )=(/6,7,14,10/)
  tetMap(:,7 )=(/7,8,14,12/)
  tetMap(:,8 )=(/8,5,14,9/)
  tetMap(:,9 )=(/1,5,11,9/)
  tetMap(:,10)=(/2,6,11,10/)
  tetMap(:,11)=(/3,7,10,12/)
  tetMap(:,12)=(/4,8,12,9/)
  ! now split the 2 pyras in the middle into 4 tets
  tetMap(:,13)=(/9,10,11,13/)
  tetMap(:,14)=(/9,10,11,14/)
  tetMap(:,15)=(/9,10,12,13/)
  tetMap(:,16)=(/9,10,12,14/)
  ! now the connect edges to 8 pyra faces (corner node+adjacent center nodes)
  tetMap(:,17)=(/1,9, 11,13/)
  tetMap(:,18)=(/2,10,11,13/)
  tetMap(:,19)=(/3,10,12,13/)
  tetMap(:,20)=(/4,9 ,12,13/)
  tetMap(:,21)=(/5,9 ,11,14/)
  tetMap(:,22)=(/6,10,11,14/)
  tetMap(:,23)=(/7,10,12,14/)
  tetMap(:,24)=(/8,9 ,12,14/)
CASE(4)   ! template is not periodic only use for single element
  nTets=5
  tetMap(:,1)=(/1,2,4,5/)
  tetMap(:,2)=(/2,3,4,7/)
  tetMap(:,3)=(/5,6,7,2/)
  tetMap(:,4)=(/7,8,5,4/)
  tetMap(:,5)=(/2,4,5,7/)
  IF(ANY(cartmesh%nElems.GT.1)) STOP 'The selected mesh template is not periodic and can only be used for a single element.'
CASE(5)
  nTets=12
  !Center node of box is required
  CALL getNewNode(node(9)%np,ind=ind)
  node(9 )%np%x=(node(1)%np%x+node(2)%np%x+node(3)%np%x+node(4)%np%x+node(5)%np%x+node(6)%np%x+node(7)%np%x+node(8)%np%x)/8.

  !x_minus
  tetMap(:,1) =(/1,2,6,9/)
  tetMap(:,2) =(/1,5,6,9/)
  !y_plus
  tetMap(:,3) =(/2,3,7,9/)
  tetMap(:,4) =(/2,6,7,9/)
  !x_plus
  tetMap(:,5) =(/3,4,7,9/)
  tetMap(:,6) =(/4,7,8,9/)
  !y_minus
  tetMap(:,7) =(/1,5,8,9/)
  tetMap(:,8) =(/1,4,8,9/)
  !z_plus
  tetMap(:,9) =(/5,6,8,9/)
  tetMap(:,10)=(/6,7,8,9/)
  !z_minus
  tetMap(:,11)=(/1,2,4,9/)
  tetMap(:,12)=(/2,3,4,9/)


CASE DEFAULT
  STOP 'The selected mesh template does not exist for tetrahedra.'
END SELECT

DO i=1,nTets
  IF(i.EQ.1)THEN
    CALL GetNewElem(aElem)
  ELSE
    CALL GetNewElem(aElem%nextElem)
    aElem%nextElem%prevElem => aElem
    aElem                   => aElem%nextElem
  END IF
  aElem%nNodes=4
  ALLOCATE(aElem%Node(aElem%nNodes))
  aElem%Node(1)%np=>node(tetMap(1,i))%np
  aElem%Node(2)%np=>node(tetMap(2,i))%np
  aElem%Node(3)%np=>node(tetMap(3,i))%np
  aElem%Node(4)%np=>node(tetMap(4,i))%np
END DO

! Add elements to list
IF(.NOT.ASSOCIATED(FirstElem))THEN
  FirstElem=>aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
ELSE
  aElem%nextElem          => FirstElem
  aElem%nextElem%prevElem => aElem
  FirstElem          => aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
END IF
NULLIFY(aElem)

END SUBROUTINE GetNewTetrahedron


SUBROUTINE GetNewPyramid(CornerNode)
!===================================================================================================================================
! Build new pyramid for cartesian mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem,tNode
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewNode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNodePtr),INTENT(IN)                  :: CornerNode(8)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: aElem  ! ?
TYPE(tNode),POINTER             :: CenterNode  ! ?
REAL                            :: xm(3)   ! ?
INTEGER                         :: i  ! ?
!===================================================================================================================================
! Get center node
xm=0.
DO i=1,8
  xm=xm+CornerNode(i)%np%x
END DO
xm=0.125*xm
CALL GetNewNode(CenterNode)
CenterNode%x=xm

! 1. Pyramid
CALL GetNewElem(aElem)
aElem%nNodes=5
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(1)%NP
aElem%Node(2)%NP=>CornerNode(2)%NP
aElem%Node(3)%NP=>CornerNode(3)%NP
aElem%Node(4)%NP=>CornerNode(4)%NP
aElem%Node(5)%NP=>CenterNode
! 2. Pyramid
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=5
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(1)%NP
aElem%Node(2)%NP=>CornerNode(5)%NP
aElem%Node(3)%NP=>CornerNode(6)%NP
aElem%Node(4)%NP=>CornerNode(2)%NP
aElem%Node(5)%NP=>CenterNode
! 3. Pyramid
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=5
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(2)%NP
aElem%Node(2)%NP=>CornerNode(6)%NP
aElem%Node(3)%NP=>CornerNode(7)%NP
aElem%Node(4)%NP=>CornerNode(3)%NP
aElem%Node(5)%NP=>CenterNode
! 4. Pyramid
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=5
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(1)%NP
aElem%Node(2)%NP=>CornerNode(4)%NP
aElem%Node(3)%NP=>CornerNode(8)%NP
aElem%Node(4)%NP=>CornerNode(5)%NP
aElem%Node(5)%NP=>CenterNode
! 5. Pyramid
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=5
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(5)%NP
aElem%Node(2)%NP=>CornerNode(8)%NP
aElem%Node(3)%NP=>CornerNode(7)%NP
aElem%Node(4)%NP=>CornerNode(6)%NP
aElem%Node(5)%NP=>CenterNode
! 6. Pyramid
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=5
ALLOCATE(aElem%node(aElem%nnodes))
aElem%Node(1)%NP=>CornerNode(7)%NP
aElem%Node(2)%NP=>CornerNode(8)%NP
aElem%Node(3)%NP=>CornerNode(4)%NP
aElem%Node(4)%NP=>CornerNode(3)%NP
aElem%Node(5)%NP=>CenterNode

! Add elements to list
IF(.NOT.ASSOCIATED(FirstElem))THEN
  FirstElem=>aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
ELSE
  aElem%nextElem          => FirstElem
  aElem%nextElem%prevElem => aElem
  FirstElem          => aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
END IF
NULLIFY(aElem)

END SUBROUTINE GetNewPyramid


SUBROUTINE GetNewPrism(CornerNode)
!===================================================================================================================================
! Build new prism for cartesian mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:getNewElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNodePtr),INTENT(IN)                  :: CornerNode(8)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: aElem   ! ?
!===================================================================================================================================
CALL getNewElem(aElem)
aElem%nNodes=6
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(1)%NP
aElem%Node(2)%NP=>CornerNode(4)%NP
aElem%Node(3)%NP=>CornerNode(5)%NP
aElem%Node(4)%NP=>CornerNode(2)%NP
aElem%Node(5)%NP=>CornerNode(3)%NP
aElem%Node(6)%NP=>CornerNode(6)%NP
CALL getNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=6
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(4)%NP
aElem%Node(2)%NP=>CornerNode(8)%NP
aElem%Node(3)%NP=>CornerNode(5)%NP
aElem%Node(4)%NP=>CornerNode(3)%NP
aElem%Node(5)%NP=>CornerNode(7)%NP
aElem%Node(6)%NP=>CornerNode(6)%NP

! Add elements to list
IF(.NOT.ASSOCIATED(FirstElem))THEN
  FirstElem=>aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
ELSE
  aElem%nextElem          => FirstElem
  aElem%nextElem%prevElem => aElem
  FirstElem          => aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
END IF
NULLIFY(aElem)

END SUBROUTINE GetNewPrism


SUBROUTINE GetNewHexahedron(CornerNode)
!===================================================================================================================================
! Build new hexahedron for cartesian mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:getNewElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNodePtr),INTENT(IN)                  :: CornerNode(8)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: aElem   ! ?
INTEGER                         :: i  ! ?
!===================================================================================================================================
CALL getNewElem(aElem)
aElem%nNodes=8
ALLOCATE(aElem%Node(aElem%nNodes))
DO i=1,8
  aElem%Node(i)%NP=>CornerNode(i)%NP
END DO

! Add elements to list
IF(.NOT.ASSOCIATED(FirstElem))THEN
  FirstElem=>aElem
ELSE
  aElem%nextElem          => FirstElem
  aElem%nextElem%prevElem => aElem
  FirstElem          => aElem
END IF
NULLIFY(aElem)
END SUBROUTINE GetNewHexahedron


SUBROUTINE CartesianMesh()
!===================================================================================================================================
! Builds cartesian mesh. Called by fillMesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:nZones,BoundaryType
USE MOD_Mesh_Vars,ONLY:getNewSide,getNewNode,getNewBC
USE MOD_Mesh_Vars,ONLY:deleteSide,deleteNode
USE MOD_Mesh_Vars,ONLY:InnerElemStretch,BoundaryOrder
USE MOD_Mesh_Basis,ONLY:CreateSides
USE MOD_CurvedCartMesh,ONLY:GetNewCurvedHexahedron
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! Some parameters set in INI-File:
! nZones                    : Number of mesh zones
! CartMeshes(iCm)%cm%nElems : Number of elements in zone icm
! CartMesh%Corner                : Corner points of zone
! CartMesh%elemType              : Element type in zone
! CartMesh%BCtype                : Type of boundary condition for zone boundaries
! CartMesh%BCstate               : Boundary state corresponding to boundary type
! CartMesh%BCalphaInd            : Boundary vector index for periodic zone boundaries
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tCartesianMesh),POINTER :: CartMesh  ! ?
TYPE(tNodePtr)               :: CornerNode(8)   ! ?
TYPE(tNodePtr),POINTER       :: Mnodes(:,:,:)   ! ?
TYPE(tElem),POINTER          :: aElem   ! ?
TYPE(tSide),POINTER          :: aSide  ! ?
TYPE(tNodePtr),POINTER       :: CurvedNode(:,:,:)  ! ?
REAL                         :: x1(3),x2(3),e1(3)  ! ?
REAL                         :: Co(3,8) ,dx(3),fac(3)   ! ?
REAL                         :: F, dF  ! ?
INTEGER                      :: iZone,i_Dim,iSide,iter  ! ?
INTEGER                      :: i,j,k,l,m,n  ! ?
INTEGER                      :: nElems(3),NodeInd  ! ?
INTEGER                      :: ne(3)  ! ?
INTEGER                      :: Ngeo
LOGICAL                      :: onBnd  ! ?
CHARACTER(LEN=32)            :: formatstr
!===================================================================================================================================
WRITE(UNIT_stdOut,*)'Building cartesian mesh'
NodeInd=0
DO iZone=1,nZones
  WRITE(UNIT_stdOut,*)'  build zone no.:' ,iZone
  CartMesh=>CartMeshes(iZone)%CM
  nElems=CartMesh%nElems
  ! Stretching elements with cartmesh%factor and/or carmesh%l0
  DO i_Dim=1,3
    SELECT CASE(i_Dim)
    CASE(1)
      e1(i_Dim)=CartMesh%Corner(i_Dim,2)-CartMesh%Corner(i_Dim,1) ! e1 container for cartmesh-dimensions
    CASE(2)
      e1(i_Dim)=CartMesh%Corner(i_Dim,4)-CartMesh%Corner(i_Dim,1)
    CASE(3)
      e1(i_Dim)=CartMesh%Corner(i_Dim,5)-CartMesh%Corner(i_Dim,1)
    END SELECT
    IF(ABS(CartMesh%l0(i_Dim)) .LT. PP_RealTolerance ) THEN !l0 = 0 = inactive
      IF(ABS(CartMesh%factor(i_Dim)) .LT. PP_RealTolerance ) THEN !fac= 0 = equidistant
        fac(i_Dim)=1.       
      ELSE ! stretched, (nElem,f) given
        fac(i_Dim)=(ABS(CartMesh%factor(i_Dim)))**(SIGN(1.,CartMesh%factor(i_Dim))) ! sign for direction 
      END IF
    ELSE !l0 active
      dx(i_Dim)=e1(i_Dim)/ABS(cartMesh%l0(i_Dim)) ! l/l0
      IF(dx(i_Dim) .LT. (1.-PP_RealTolerance)) THEN ! l0 > l
        !CALL abort(__STAMP__,  & 
         stop  'stretching error, length l0 longer than grid region, in direction '  !,i_Dim,999.)
      END IF
      IF(ABS(CartMesh%factor(i_Dim)) .LT. PP_RealTolerance ) THEN ! fac=0 , (nElem,l0) given, fac calculated
        ! 
      ELSE ! (factor, l0) given, change nElems
        fac(i_Dim)=(ABS(CartMesh%factor(i_Dim)))**(SIGN(1.,CartMesh%factor(i_Dim)*CartMesh%l0(i_Dim))) ! sign for direction 
        IF(fac(i_Dim) .NE. 1.) THEN
          CartMesh%nElems(i_Dim)=NINT(LOG(1.-dx(i_Dim)*(1.-fac(i_Dim))) / LOG(fac(i_Dim))) !nearest Integer 
        ELSE
          CartMesh%nElems(i_Dim)=NINT(dx(i_Dim))
        END IF
        IF (CartMesh%nElems(i_Dim) .LT. 1) CartMesh%nElems(i_Dim)=1
        WRITE(UNIT_stdOut,*)'   -Element number in dir',i_Dim,'changed from', &
                              nElems(i_Dim),'to',CartMesh%nElems(i_Dim) 
        nElems(i_Dim)=CartMesh%nElems(i_Dim) 
      END IF
      IF(nElems(i_Dim).EQ. 1) THEN 
        fac(i_Dim)=1.0          
      ELSEIF(nElems(i_Dim).EQ. 2) THEN 
        fac(i_Dim)=dx(i_Dim)-1.
      ELSE !nElems > 2
        fac(i_Dim)=dx(i_Dim)/nElems(i_Dim) !start value for Newton iteration
        IF (abs(fac(i_Dim)-1.) .GT. PP_RealTolerance) THEN ! NEWTON iteration, only if not equidistant case
          F=1.
          dF=1.
          iter=0
          DO WHILE ((ABS(F) .GT. PP_RealTolerance) .AND. (ABS(F/dF) .GT. PP_RealTolerance) .AND. (iter .LT. 1000))
            F = fac(i_Dim)**nElems(i_Dim) + dx(i_Dim)*(1.-fac(i_Dim)) -1. ! non-linear function
            dF= nElems(i_Dim)*fac(i_Dim)**(nElems(i_Dim)-1) -dx(i_Dim)  !dF/dfac
            fac(i_Dim)= fac(i_Dim) - F/dF
            iter=iter+1
          END DO
          IF(iter.GT.1000) STOP 'Newton iteration for computing the stretching function has failed.'
          fac(i_Dim)=fac(i_Dim)**SIGN(1.,CartMesh%l0(i_Dim)) ! sign for direction
        END IF
      END IF
      WRITE(UNIT_stdOut,*)'   -stretching factor in dir',i_Dim,'is now', fac(i_Dim)  
    END IF

    IF( ABS((nElems(i_Dim)-1.)*LOG(fac(i_Dim))/LOG(10.)) .GE. 4. )                      &
      !CALL abort(__STAMP__, &
       stop  'stretching error, length ratio > 1.0E4 in direction '  !,i_Dim,999.)
    IF(ABS(fac(i_Dim)-1.) .GT. PP_RealTolerance ) THEN
      dx(i_Dim)=(1.-fac(i_Dim))/(1.-fac(i_Dim)**nElems(i_Dim)) ! first length to start
    ELSE !equidistant case
      dx(i_Dim)=1./ nElems(i_Dim)
    END IF
  END DO  !i_Dim=1,3
  IF(ABS(SUM(fac(:))-3.).GT. PP_RealTolerance) THEN
    WRITE(UNIT_stdOut,*) '   -STRETCHING INFORMATION'
    WRITE(formatstr,'(A5,I1,A8)')'(A19,',3,'(I9,5X))'
    WRITE(UNIT_stdOut,formatstr) '      -nElems(:) : ', CartMesh%nElems(:)
    WRITE(formatstr,'(A5,I1,A11)')'(A19,',3,'(G13.7,1X))'
    WRITE(UNIT_stdOut,formatstr) '      -l0(:)     : ', e1(:)*dx(:)
    WRITE(UNIT_stdOut,formatstr) '      -factor(:) : ', fac(:)
  END IF

  IF(InnerElemStretch.AND.CartMesh%ElemType.EQ.108)THEN
    Ngeo=BoundaryOrder-1
    ALLOCATE(CurvedNode(0:Ngeo,0:Ngeo,0:Ngeo))
  ELSE
    Ngeo=1
  END IF

  ne=CartMesh%nElems*Ngeo
  ALLOCATE(Mnodes(0:ne(1),0:ne(2),0:ne(3)))
  Co=CartMesh%Corner
  !recompute factor for high order mesh
  IF(Ngeo.NE.1) fac=fac**(1./REAL(Ngeo))
  DO i_dim=1,3
    IF(fac(i_dim).NE.1.)THEN
      dx(i_dim)=dx(i_dim)*(1.-fac(i_dim))/(1.-fac(i_dim)**Ngeo)
    ELSE
      dx(i_dim)=dx(i_dim)/REAL(Ngeo)
    END IF
  END DO
  ! Build nodes
  e1(:)=dx(:)/fac(:) !container for dx
  !start position x2 in unit square grid
  x2(:)=0.
  ! calculate node postions, x1(:) coordinates in unit square grid
  x1(1)=x2(1)
  dx(1)=e1(1)
  DO l=0,ne(1)
!    x1(1)=REAL(l)/REAL(nSplit(1))
    x1(2)=x2(2)
    dx(2)=e1(2) 
    DO m=0,ne(2)
!      x1(2)=REAL(m)/REAL(nSplit(2))
      x1(3)=x2(3)
      dx(3)=e1(3)
      DO n=0,ne(3)
!        x1(3)=REAL(n)/REAL(nSplit(3))
        CALL GetNewNode(Mnodes(l,m,n)%np,ind=NodeInd)
        Mnodes(l,m,n)%np%x=Co(:,1)+x1(1)*(Co(:,2)-Co(:,1))+x1(2)*(Co(:,4)-Co(:,1))+x1(3)*(Co(:,5)-Co(:,1))+ &
                         (Co(:,3)-Co(:,4)-Co(:,2)+Co(:,1))*x1(1)*x1(2)+(Co(:,6)-Co(:,5)-Co(:,2)+Co(:,1))*x1(1)*x1(3)+&
                            (Co(:,8)-Co(:,4)-Co(:,5)+Co(:,1))*x1(2)*x1(3)+(Co(:,7)-Co(:,1)+Co(:,4)+Co(:,5)+Co(:,2)-  &
                                  Co(:,8)-Co(:,3)-Co(:,6))*x1(1)*x1(2)*x1(3)
        dx(3)=dx(3)*fac(3)
        x1(3)=x1(3)+dx(3)
      END DO
      dx(2)=dx(2)*fac(2)
      x1(2)=x1(2)+dx(2)
    END DO
    dx(1)=dx(1)*fac(1)
    x1(1)=x1(1)+dx(1)
  END DO

  DO l=0,ne(1)-1,Ngeo
    DO m=0,ne(2)-1,Ngeo
      DO n=0,ne(3)-1,Ngeo
        CornerNode(1)%np=>Mnodes(l     ,m     ,n     )%np
        CornerNode(2)%np=>Mnodes(l+Ngeo,m     ,n     )%np
        CornerNode(3)%np=>Mnodes(l+Ngeo,m+Ngeo,n     )%np
        CornerNode(4)%np=>Mnodes(l     ,m+Ngeo,n     )%np
        CornerNode(5)%np=>Mnodes(l     ,m     ,n+Ngeo)%np
        CornerNode(6)%np=>Mnodes(l+Ngeo,m     ,n+Ngeo)%np
        CornerNode(7)%np=>Mnodes(l+Ngeo,m+Ngeo,n+Ngeo)%np
        CornerNode(8)%np=>Mnodes(l     ,m+Ngeo,n+Ngeo)%np
        SELECT CASE(CartMesh%ElemType)
        CASE(104) ! Tetrahedron
          CALL GetNewTetrahedron(CornerNode,CartMesh,l+1,m+1,n+1,ind=NodeInd)
        CASE(105) ! Pyramid
          CALL GetNewPyramid(CornerNode)
        CASE(106) ! Dreiecks-prismon
          CALL GetNewPrism(CornerNode)
        CASE(108) ! Hexaeder
          IF(InnerElemStretch)THEN
            DO i=0,Ngeo; DO j=0,Ngeo; DO k=0,Ngeo
              CurvedNode(i,j,k)%np=>Mnodes(l+i,m+j,n+k)%np
            END DO; END DO; END DO; 
            CALL GetNewCurvedHexahedron(CurvedNode,Ngeo,iZone)
          ELSE 
            CALL GetNewHexahedron(CornerNode)
          END IF
        CASE DEFAULT
          CALL abort(__STAMP__,&
            'The specified element type is not known. Valid types: 104,105,106,108',CartMesh%ElemType)
        END SELECT
      END DO !n
    END DO !m
  END DO !l
  !for Boundary Conditions
  DO l=0,ne(1),Ngeo
    DO m=0,ne(2),Ngeo
      DO n=0,ne(3),Ngeo
        Mnodes(l,m,n)%np%tmp=0 
        IF(n.EQ.0) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+1 !zeta minus
        IF(m.EQ.0) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+20 !eta minus
        IF(l.EQ.ne(1) ) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+300 !xi plus
        IF(m.EQ.ne(2) ) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+4000  !eta plus
        IF(l.EQ.0) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+50000 !xi minus
        IF(n.EQ.ne(3) ) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+600000 !zeta plus
      END DO !n
    END DO !m
  END DO !l

  ! Set element zone
  aElem=>FirstElem
  DO WHILE(ASSOCIATED(aElem))
    IF(aElem%Zone .EQ. 0)THEN
      aElem%Zone=iZone
      CALL CreateSides(aElem,.TRUE.)
    END IF
    aElem=>aElem%nextElem
  END DO

  aElem=>FirstElem
  DO WHILE(ASSOCIATED(aElem))
    IF(aElem%Zone .EQ. iZone)THEN
      aSide=>aElem%firstSide
      DO WHILE(ASSOCIATED(aSide))
        DO iSide=1,6
          onBnd=.TRUE.
          DO i=1, aSide%nNodes
            IF((MOD(aSide%Node(i)%np%tmp,10**iSide)/(10**(iSide-1))).NE.iSide)THEN
              onBnd=.FALSE.
              EXIT ! Loop
            END IF
          END DO
          IF(CartMesh%BCIndex(iSide).EQ.0) onBnd=.FALSE.
          IF(onBnd)THEN
            CALL getNewBC(aSide%BC)
            aSide%BC%BCIndex    = CartMesh%BCIndex(iSide)
            aSide%BC%BCType     = BoundaryType(aSide%BC%BCIndex,1)
            aSide%curveIndex    = BoundaryType(aSide%BC%BCIndex,2)
            aSide%BC%BCstate    = BoundaryType(aSide%BC%BCIndex,3)
            aSide%BC%BCalphaInd = BoundaryType(aSide%BC%BCIndex,4)
          END IF
        END DO  ! iSide=1,6
        aSide=>aSide%nextElemSide
      END DO ! WHILE(ASSOCIATED(aSide))
    END IF  ! aElem%Zone .EQ. iZone
    aElem=>aElem%nextElem
  END DO ! WHILE(ASSOCIATED(aElem))
  DEALLOCATE(Mnodes)
END DO ! iZone
END SUBROUTINE CartesianMesh

END MODULE MOD_CartMesh
