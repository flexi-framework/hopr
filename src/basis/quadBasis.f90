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
MODULE MOD_QuadBasis
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE getQuadBasis
  MODULE PROCEDURE getQuadBasis
END INTERFACE

INTERFACE getBasisMappingQuad
  MODULE PROCEDURE getBasisMappingQuad
END INTERFACE

PUBLIC::getQuadBasis,getBasisMappingQuad
!===================================================================================================================================

CONTAINS
SUBROUTINE getQuadBasis(Deg,nNodes1D,Vdm_visu,D_visu)
!===================================================================================================================================
! given the degree of the orthogonal basis and the number of 1D nodes, equidistant nodes in the tetrahedron are generated and 
! we compute the Vadndermonde matrix and the Vandermondematrix of the gradient of the basis function
!===================================================================================================================================
! MODULES
USE MOD_Globals, ONLY: GetInverse
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: Deg  ! ?
INTEGER, INTENT(IN) :: nNodes1D  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,ALLOCATABLE,INTENT(OUT)        :: Vdm_visu(:,:)  ! ?
REAL,ALLOCATABLE,INTENT(OUT)        :: D_visu(:,:,:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,ALLOCATABLE        :: VdmInv(:,:)  ! ?
REAL,ALLOCATABLE        :: r1D(:)            ! equidistant 1D Lobatto nodes in [-1,1]
INTEGER,ALLOCATABLE     :: bMap(:,:)         ! basis mapping iAns=>i,j,k
INTEGER,ALLOCATABLE     :: nodeMap(:,:)      ! mapping for equidistant nodes iNode=>i,j,k
INTEGER                 :: nNodes  ! ?
INTEGER                 :: nAns  ! ?
INTEGER                 :: i    ! ?
!===================================================================================================================================
!nodal Vandermonde
ALLOCATE(r1D(0:Deg))
!equidistant nodes
DO i=0,Deg
  r1D(i)=-1.+2.*REAL(i)/REAL(Deg)
END DO
CALL getBasisMappingQuad(Deg,nAns,bMap)
nNodes=nAns
ALLOCATE(VdmInv(nNodes,nAns))
CALL VandermondeQuad(Deg+1,Deg,bMap,bMap,r1D,VdmInv)
VdmInv=GetInverse(nAns,VdmInv)
DEALLOCATE(r1D)

!visualization
ALLOCATE(r1D(0:nNodes1D-1))
!equidistant nodes
DO i=0,nNodes1D-1
  r1D(i)=-1.+2.*REAL(i)/REAL(nNodes1D-1)
END DO
CALL getBasisMappingQuad(nNodes1D-1,nNodes,nodeMap)
ALLOCATE(Vdm_visu(nNodes,nAns),D_visu(nNodes,nAns,2))

CALL VandermondeQuad(nNodes1D,Deg,nodeMap,bMap,r1D,VdM_visu)
CALL GradVandermondeQuad(nNodes1D,Deg,nodeMap,bMap,r1D,D_visu)

Vdm_visu=MATMUL(Vdm_visu,VdmInv)
D_visu(:,:,1)=MATMUL(D_visu(:,:,1),VdmInv)
D_visu(:,:,2)=MATMUL(D_visu(:,:,2),VdmInv)
DEALLOCATE(bMap,nodeMap) 
END SUBROUTINE getQuadBasis


SUBROUTINE getBasisMappingQuad(Deg,nAns,bMap,bMapInv)
!===================================================================================================================================
! mapping from iAns -> i,j  in [0,Deg], can be used for nodeMap too: CALL getBasisMapping(nNodes1D-1,nodeMap) 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)         :: Deg  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)        :: nAns  ! ?
INTEGER,ALLOCATABLE,INTENT(OUT)         :: bMap(:,:)  ! ?
INTEGER,ALLOCATABLE,OPTIONAL,INTENT(OUT):: bMapInv(:,:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                     :: iAns,i,j,d  ! ?
!===================================================================================================================================
nAns=(Deg+1)**2
ALLOCATE(bMap(nAns,2))
!iAns=0
!DO j=0,Deg
!  DO i=0,Deg
!    iAns=iAns+1
!    bMap(iAns,:)=(/i,j/)
!  END DO
!END DO 
iAns=0
DO d=1,2*Deg+1
  DO i=0,Deg
    DO j=0,Deg
      IF(i+j+1 .EQ. d) THEN
        iAns=iAns+1
        bMap(iAns,:)=(/i,j/)
      END IF
    END DO
  END DO
END DO
IF(PRESENT(bMapInv))THEN
  ALLOCATE(bMapInv(0:Deg,0:Deg))
  bMapInv=0
  DO iAns=1,nAns
    bMapInv(bMap(iAns,1),bMap(iAns,2))=iAns
  END DO
END IF
END SUBROUTINE getBasisMappingQuad


SUBROUTINE VandermondeQuad(nNodes1D,Deg,nodeMap,bMap,r1D,VdM)
!===================================================================================================================================
! given the degree of the polynomial and the 1D nodes, nodes in the quadrangle are generated and 
! we compute the 2D Vandermonde matrix. The polynomial basis is tensor-product orthonormal on the reference square [-1,1]^2
!===================================================================================================================================
! MODULES
USE MOD_basis1D, ONLY : Vandermonde1D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: nNodes1D                      ! number of nodes in 1D 
INTEGER, INTENT(IN) :: Deg                           ! polynomial degree in 1D
INTEGER, INTENT(IN) :: nodeMap(nNodes1D**2,2)        ! node mapping iNode=> i,j in [0,nNodes1D-1]
INTEGER, INTENT(IN) :: bMap((Deg+1)**2,2)            ! basis mapping iAns=> i,j in [0,Deg]
REAL,    INTENT(IN) :: r1D(0:nNodes1D-1)             ! one dimensional node coordinates in reference space [1-,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,   INTENT(OUT) :: Vdm(nNodes1D**2,(Deg+1)**2)   ! 3D Vandermondematrix
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                 :: nNodes  ! ?
INTEGER                 :: nAns  ! ?
INTEGER                 :: iNode,iAns  ! ?
REAL                    :: Vdm1D(0:nNodes1D-1,0:Deg)  ! ?
!===================================================================================================================================
CALL Vandermonde1D(nNodes1D,Deg,r1D,VdM1D)
nAns=(Deg+1)**2
nNodes=nNodes1D**2
DO iAns=1,nAns
  DO iNode=1,nNodes
    VdM(iNode,iAns)=Vdm1D(nodeMap(iNode,1),bMap(iAns,1))*Vdm1D(nodeMap(iNode,2),bMap(iAns,2))
  END DO
END DO
END SUBROUTINE VandermondeQuad


SUBROUTINE GradVandermondeQuad(nNodes1D,Deg,nodeMap,bMap,r1D,GradVdM)
!===================================================================================================================================
! given the degree of the polynomial and the 1D nodes, nodes in the quadrangle are generated and 
! we compute the 2D Gradient Vandermonde matrix. The polynomial basis is tensor-product orthonormal on the reference square [-1,1]^2
!===================================================================================================================================
! MODULES
USE MOD_basis1D, ONLY : Vandermonde1D,GradVandermonde1D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: nNodes1D                      ! number of nodes in 1D 
INTEGER, INTENT(IN) :: Deg                           ! polynomial degree in 1D
INTEGER, INTENT(IN) :: nodeMap(nNodes1D**2,2)        ! node mapping iNode=> i,j in [0,nNodes1D-1] 
INTEGER, INTENT(IN) :: bMap((Deg+1)**2,2)            ! basis mapping iAns=> i,j in [0,Deg] 
REAL,    INTENT(IN) :: r1D(0:nNodes1D-1)             ! one dimensional node coordinates in reference space [1-,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,   INTENT(OUT) :: GradVdm(nNodes1D**2,(Deg+1)**2,2)   ! 3D Vandermondematrix
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                 :: nNodes  ! ?
INTEGER                 :: nAns  ! ?
INTEGER                 :: iNode,iAns  ! ?
REAL                    :: Vdm1D(0:nNodes1D-1,0:Deg),GradVdm1D(0:nNodes1D-1,0:Deg)  ! ?
!===================================================================================================================================
CALL Vandermonde1D(nNodes1D,Deg,r1D,VdM1D)
CALL GradVandermonde1D(nNodes1D,Deg,r1D,GradVdM1D)
nAns=(Deg+1)**2
nNodes=nNodes1D**2
DO iAns=1,nAns
  DO iNode=1,nNodes
    GradVdM(iNode,iAns,1)=GradVdm1D(nodeMap(iNode,1),bMap(iAns,1))*Vdm1D(nodeMap(iNode,2),bMap(iAns,2))
    GradVdM(iNode,iAns,2)=Vdm1D(nodeMap(iNode,1),bMap(iAns,1))*GradVdm1D(nodeMap(iNode,2),bMap(iAns,2))
  END DO
END DO
END SUBROUTINE GradVandermondeQuad

END MODULE MOD_QuadBasis
