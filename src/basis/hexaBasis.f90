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
MODULE MOD_HexaBasis
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
INTERFACE getHexaBasis
  MODULE PROCEDURE getHexaBasis
END INTERFACE

INTERFACE getBasisMappingHexa
  MODULE PROCEDURE getBasisMappingHexa
END INTERFACE

PUBLIC::getHexaBasis,getBasisMappingHexa
!===================================================================================================================================

CONTAINS
SUBROUTINE getHexaBasis(Deg,nNodes1D,Vdm_visu,D_visu)
!===================================================================================================================================
! given the degree of the orthogonal basis and the number of 1D nodes, equidistant nodes in the tetrahedron are generated and 
! we compute the Vadndermonde matrix and the Vandermondematrix of the gradient of the basis function
!===================================================================================================================================
! MODULES
USE MOD_Basis1D, ONLY:BarycentricWeights,InitializeVandermonde,PolynomialDerivativeMatrix
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: Deg  ! ?
INTEGER, INTENT(IN) :: nNodes1D  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,ALLOCATABLE,INTENT(OUT)        :: Vdm_visu(:,:) ! ?
REAL,ALLOCATABLE,INTENT(OUT)        :: D_visu(:,:,:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,ALLOCATABLE        :: Vdm1D(:,:),D0(:,:),D_1D(:,:)  ! ?
REAL,ALLOCATABLE        :: r1D(:),wBary(:)    ! equidistant 1D Lobatto nodes in [-1,1]
REAL,ALLOCATABLE        :: r1Dvisu(:)    ! equidistant 1D Lobatto nodes in [-1,1]
INTEGER,ALLOCATABLE     :: bMap(:,:)         ! basis mapping iAns=>i,j,k
INTEGER,ALLOCATABLE     :: nodeMap(:,:)      ! mapping for equidistant nodes iNode=>i,j,k
INTEGER                 :: nNodes  ! ?
INTEGER                 :: nAns  ! ?
INTEGER                 :: i  ! ?
!===================================================================================================================================
ALLOCATE(r1D(0:Deg),wBary(0:Deg))
!equidistant nodes of the mapping
DO i=0,Deg
  r1D(i)=-1.+2.*REAL(i)/REAL(Deg)
END DO
CALL BarycentricWeights(Deg,r1D,wBary)

!visualization
ALLOCATE(r1Dvisu(0:nNodes1D-1))
!equidistant nodes for output
DO i=0,nNodes1D-1
  r1Dvisu(i)=-1.+2.*REAL(i)/REAL(nNodes1D-1)
END DO

CALL getBasisMappingHexa(Deg,nAns,bMap)

CALL getBasisMappingHexa(nNodes1D-1,nNodes,nodeMap)
ALLOCATE(Vdm_visu(nNodes,nAns),D_visu(nNodes,nAns,3))

ALLOCATE(Vdm1D(0:nNodes1D-1,0:Deg),D0(0:Deg,0:Deg))
CALL InitializeVandermonde(Deg,nNodes1D-1,wBary,r1D,r1Dvisu,Vdm1D)
CALL PolynomialDerivativeMatrix(Deg,r1D,D0)
ALLOCATE(D_1D(0:nNodes1D-1,0:Deg))
D_1D=MATMUL(Vdm1D,D0)


!nodal 3D Vandermonde
CALL VandermondeHexa(nNodes1D,Deg,nodeMap,bMap,Vdm1D,VdM_visu)

!nodal 3D Dmatrix
CALL GradVandermondeHexa(nNodes1D,Deg,nodeMap,bMap,Vdm1D,D_1D,D_visu)

DEALLOCATE(bMap,nodeMap) 
END SUBROUTINE getHexaBasis


SUBROUTINE getBasisMappingHexa(Deg,nAns,bMap,bMapInv)
!===================================================================================================================================
! mapping from iAns -> i,j  in [0,Deg], can be used for nodeMap too: CALL getBasisMapping(nNodes1D-1,nodeMap) 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: Deg  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)         :: nAns  ! ?
INTEGER,ALLOCATABLE,INTENT(OUT)          :: bMap(:,:)  ! ?
INTEGER,ALLOCATABLE,OPTIONAL,INTENT(OUT) :: bMapInv(:,:,:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                     :: iAns,i,j,k  ! ?
!===================================================================================================================================
nAns=(Deg+1)**3
ALLOCATE(bMap(nAns,3))
iAns=0
DO k=0,Deg
  DO j=0,Deg
    DO i=0,Deg
      iAns=iAns+1
      bMap(iAns,:)=(/i,j,k/)
    END DO
  END DO 
END DO 
IF(PRESENT(bMapInv))THEN
  ALLOCATE(bMapInv(0:Deg,0:Deg,0:Deg))
  bMapInv=0
  iAns=0
  DO k=0,Deg
    DO j=0,Deg
      DO i=0,Deg
        iAns=iAns+1
        bMapInv(i,j,k)=iAns
      END DO
    END DO 
  END DO 
END IF
END SUBROUTINE getBasisMappingHexa


SUBROUTINE VandermondeHexa(nNodes1D,Deg,nodeMap,bMap,Vdm1D,VdM)
!===================================================================================================================================
! given the degree of the polynomial and the 1D nodes, nodes in the hexahedron are generated and 
! we compute the 3D Vandermonde matrix. The polynomial basis is tensor-product orthonormal on the reference hexa [-1,1]^3
!===================================================================================================================================
! MODULES
USE MOD_Basis1D, ONLY : Vandermonde1D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: nNodes1D                      ! number of nodes in 1D 
INTEGER, INTENT(IN) :: Deg                           ! polynomial degree in 1D
INTEGER, INTENT(IN) :: nodeMap(nNodes1D**3,3)        ! node mapping iNode=> i,j,k in [0,nNodes1D-1]
INTEGER, INTENT(IN) :: bMap((Deg+1)**3,3)            ! basis mapping iAns=> i,j,k in [0,Deg]
REAL,    INTENT(IN) :: Vdm1D(0:nNodes1D-1,0:Deg)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,   INTENT(OUT) :: Vdm(nNodes1D**3,(Deg+1)**3)   ! 3D Vandermondematrix
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                 :: nNodes  ! ?
INTEGER                 :: nAns  ! ?
INTEGER                 :: iNode,iAns  ! ?
!===================================================================================================================================
nAns=(Deg+1)**3
nNodes=nNodes1D**3
DO iAns=1,nAns
  DO iNode=1,nNodes
    VdM(iNode,iAns)=Vdm1D(nodeMap(iNode,1),bMap(iAns,1))* &
                    Vdm1D(nodeMap(iNode,2),bMap(iAns,2))* &
                    Vdm1D(nodeMap(iNode,3),bMap(iAns,3))
  END DO
END DO
END SUBROUTINE VandermondeHexa



SUBROUTINE GradVandermondeHexa(nNodes1D,Deg,nodeMap,bMap,Vdm1D,D_1D,GradVdM)
!===================================================================================================================================
! given the degree of the polynomial and the 1D nodes, nodes in the hexahedron are generated and 
! we compute the 3D Gradient Vandermonde matrix. The polynomial basis is tensor-product orthonormal on the reference hexa [-1,1]^3
!===================================================================================================================================
! MODULES
USE MOD_Basis1D, ONLY : Vandermonde1D,GradVandermonde1D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: nNodes1D                      ! number of nodes in 1D 
INTEGER, INTENT(IN) :: Deg                           ! polynomial degree in 1D
INTEGER, INTENT(IN) :: nodeMap(nNodes1D**3,3)        ! node mapping iNode=> i,j,k in [0,nNodes1D-1]
INTEGER, INTENT(IN) :: bMap((Deg+1)**3,3)            ! basis mapping iAns=> i,j,k in [0,Deg]
REAL,    INTENT(IN) :: Vdm1D(0:nNodes1D-1,0:Deg)  ! ?
REAL,    INTENT(IN) :: D_1D(0:nNodes1D-1,0:Deg)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,   INTENT(OUT) :: GradVdm(nNodes1D**3,(Deg+1)**3,3)   ! 3D Vandermondematrix
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                 :: nNodes  ! ?
INTEGER                 :: nAns  ! ?
INTEGER                 :: iNode,iAns  ! ?
!===================================================================================================================================
nAns=(Deg+1)**3
nNodes=nNodes1D**3
DO iAns=1,nAns
  DO iNode=1,nNodes
    GradVdM(iNode,iAns,1)=  D_1D(nodeMap(iNode,1),bMap(iAns,1))* &
                           Vdm1D(nodeMap(iNode,2),bMap(iAns,2))* &
                           Vdm1D(nodeMap(iNode,3),bMap(iAns,3))
    GradVdM(iNode,iAns,2)=  Vdm1D(nodeMap(iNode,1),bMap(iAns,1))* &
                             D_1D(nodeMap(iNode,2),bMap(iAns,2))* &
                            Vdm1D(nodeMap(iNode,3),bMap(iAns,3))
    GradVdM(iNode,iAns,3)=  Vdm1D(nodeMap(iNode,1),bMap(iAns,1))* &
                            Vdm1D(nodeMap(iNode,2),bMap(iAns,2))* &
                             D_1D(nodeMap(iNode,3),bMap(iAns,3))
  END DO
END DO
END SUBROUTINE GradVandermondeHexa

END MODULE MOD_hexaBasis
