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
MODULE MOD_triaBasis
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
INTERFACE getTriaBasis
  MODULE PROCEDURE getTriaBasis
END INTERFACE

INTERFACE getBasisMappingTria
  MODULE PROCEDURE getBasisMappingTria
END INTERFACE

PUBLIC::getTriaBasis,getBasisMappingTria
!===================================================================================================================================

CONTAINS
SUBROUTINE getTriaBasis(Deg,nNodes1D,Vdm_visu,D_visu)
!===================================================================================================================================
! given the degree of the orthogonal basis and the number of 1D nodes, equidistant nodes in the tetrahedron are generated and 
! we compute the Vadndermonde matrix and the Vandermondematrix of the gradient of the basis function
!===================================================================================================================================
! MODULES
USE MOD_Globals  ,ONLY: GetInverse
USE MOD_quadBasis,ONLY: getBasisMappingQuad
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
REAL,ALLOCATABLE        :: VdmInv(:,:)  ! ?
REAL,ALLOCATABLE        :: r(:),s(:)         ! coordinates in tetrahedra reference space  r,s,t in [-1,1]
REAL,ALLOCATABLE        :: r1D(:)            ! equidistant 1D Lobatto nodes in [-1,1]
INTEGER,ALLOCATABLE     :: bMap(:,:)         ! basis mapping iAns=>i,j,k
INTEGER,ALLOCATABLE     :: nodeMap(:,:)      ! mapping for equidistant nodes iNode=>i,j,k
INTEGER                 :: nNodes  ! ?
INTEGER                 :: nAns  ! ?
INTEGER                 :: i  ! ?
!===================================================================================================================================
!nodal Vandermonde
ALLOCATE(r1D(0:Deg))
!equidistant nodes
DO i=0,Deg
  r1D(i)=-1.+2.*REAL(i)/REAL(Deg)
END DO
CALL getBasisMappingTria(Deg,nAns,bMap)
nNodes=nAns
ALLOCATE(VdmInv(nNodes,nAns))
ALLOCATE(r(nNodes),s(nNodes))
r(:)=r1D(bMap(:,1))
s(:)=r1D(bMap(:,2))
CALL VandermondeTria(nNodes,nAns,bMap,r,s,VdMinv)
VdmInv=GetInverse(nAns,VdmInv)
DEALLOCATE(r1D,r,s)

!visualization
ALLOCATE(r1D(0:nNodes1D-1))
!equidistant nodes
DO i=0,nNodes1D-1
  r1D(i)=-1.+2.*REAL(i)/REAL(nNodes1D-1)
END DO
CALL getBasisMappingQuad(nNodes1D-1,nNodes,nodeMap) !used for visualization, tensorproduct
ALLOCATE(Vdm_visu(nNodes,nAns),D_visu(nNodes,nAns,2))
ALLOCATE(r(nNodes),s(nNodes))
r(:)=r1D(nodeMap(:,1))
s(:)=r1D(nodeMap(:,2))
r=0.5*(1+r)*(1-s)-1.  ! trafo tensor-product points to triangle

CALL VandermondeTria(nNodes,nAns,bMap,r,s,VdM_visu)
CALL GradVandermondeTria(nNodes,nAns,bMap,r,s,D_visu)
Vdm_visu=MATMUL(Vdm_visu,VdmInv)
D_visu(:,:,1)=MATMUL(D_visu(:,:,1),VdmInv)
D_visu(:,:,2)=MATMUL(D_visu(:,:,2),VdmInv)
DEALLOCATE(bMap,nodeMap) 
END SUBROUTINE getTriabasis


SUBROUTINE getBasisMappingTria(Deg,nAns,bMap,bMapInv)
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
INTEGER,ALLOCATABLE, INTENT(OUT)         :: bMap(:,:)  ! ?
INTEGER,ALLOCATABLE,OPTIONAL, INTENT(OUT):: bMapInv(:,:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                     :: iAns,i,j,d  ! ?
!===================================================================================================================================
nAns=(Deg+1)*(Deg+2)/2
ALLOCATE(bMap(nAns,2))
iAns=0
DO d=1,Deg+1
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
  iAns=0
  DO d=1,Deg+1
    DO i=0,Deg
      DO j=0,Deg
        IF(i+j+1 .EQ. d) THEN
          iAns=iAns+1
          bMapInv(i,j)=iAns
        END IF
      END DO
    END DO 
  END DO 
END IF 
END SUBROUTINE getBasisMappingTria


SUBROUTINE rs2abTria(nNodes,r,s,a,b)
!===================================================================================================================================
! Transforms reference coordinates of the triangle (r,s) to tensorial coordinates (a,b)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: nNodes  ! ?
REAL,INTENT(IN)              :: r(nNodes),s(nNodes)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: a(nNodes),b(nNodes)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                     :: iNode  ! ?
!===================================================================================================================================
DO iNode=1,nNodes
  IF(ABS(1.-s(iNode)).LT.PP_RealTolerance)THEN
    a(iNode)=-1.
  ELSE
    a(iNode)=2.*(1.+r(iNode))/(1.-s(iNode))-1.
  END IF
END DO !iNode
b=s 
END SUBROUTINE rs2abTria


SUBROUTINE VandermondeTria(nNodes,nAns,bMap,r,s,VdM)
!===================================================================================================================================
! For a given vector of nNodes in reference coordinates (r,s) of the triangle, and nAns=(Deg+1)*(Deg+2)/2 basis functions, we
! compute the 2D Vandermonde matrix. The polynomial basis is orthonormal with respect to the reference triangle r,s>=-1, r+s<=1 
!===================================================================================================================================
! MODULES
USE MOD_Basis1D, ONLY : JacobiP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: nNodes,nAns  ! ?
INTEGER, INTENT(IN)          :: bMap(nAns,2)   ! basis mapping iAns=>i,j  i,j in [0,Deg]
REAL,INTENT(IN)              :: r(nNodes),s(nNodes)   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: VdM(nNodes,nAns)   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(nNodes)      :: a,b,Pa,Pb  ! ?
INTEGER                     :: iAns,i,j  ! ?
REAL                        :: norm  ! ?
!===================================================================================================================================
norm=SQRT(2.)
CALL rs2abTria(nNodes,r,s,a,b) 
DO iAns=1,nAns
  i=bMap(iAns,1)
  j=bMap(iAns,2)
  ! compute Vandermonde
  CALL JacobiP(nNodes,a,         0, 0,i,Pa)
  CALL JacobiP(nNodes,b,     2*i+1, 0,j,Pb)
  VdM(:,iAns)=norm*(Pa*Pb) !orthonormal
  IF(i.GT.0) VdM(:,iAns)=VdM(:,iAns)*((1.-b)**i)
END DO 
END SUBROUTINE VandermondeTria


SUBROUTINE GradVandermondeTria(nNodes,nAns,bMap,r,s,gradVdM)
!===================================================================================================================================
! For a given vector of nNodes in reference coordinates (r,s) of the triangle, and nAns=(Deg+1)*(Deg+2)/2 basis functions, we
! compute the 2D Gradient Vandermonde matrix. The polynomial basis is orthonormal with respect to reference triangle r,s>=-1, r+s<=1
!===================================================================================================================================
! MODULES
USE MOD_Basis1D, ONLY : JacobiP,GradJacobiP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: nNodes,nAns  ! ?
INTEGER, INTENT(IN)          :: bMap(nAns,2)   ! basis mapping iAns=>i,j  i,j in [0,Deg]
REAL,INTENT(IN)              :: r(nNodes),s(nNodes)   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: gradVdM(nNodes,nAns,2)   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(nNodes)      :: a,b,Pa,Pb,dPa,dPb,dPhi_da,dPhi_db  ! ?
INTEGER                     :: iAns,i,j  ! ?
REAL                        :: norm  ! ?
!===================================================================================================================================
norm=SQRT(2.)
CALL rs2abTria(nNodes,r,s,a,b) 
DO iAns=1,nAns
  i=bMap(iAns,1)
  j=bMap(iAns,2)
  
  CALL JacobiP(nNodes,a,         0, 0,i,Pa)
  CALL JacobiP(nNodes,b,     2*i+1, 0,j,Pb)
  CALL GradJacobiP(nNodes,a,         0, 0,i,dPa)
  CALL GradJacobiP(nNodes,b,     2*i+1, 0,j,dPb)
  !dphi/da*(1-b)^-1
  dPhi_da=norm*(dPa*Pb) 
  IF(i.GT.1) dPhi_da=dPhi_da*((1.-b)**(i-1))
  !dphi/db
  dPhi_db=dPb
  IF(i.GT.0) dPhi_db=dPhi_db*(1.-b)**i
  IF(i.EQ.1) dPhi_db=dPhi_db-Pb
  IF(i.GT.1) dPhi_db=dPhi_db-i*(1.-b)**(i-1)*Pb
  dPhi_db=norm*Pa*dPhi_db 
  ! r-derivative
  gradVdM(:,iAns,1)=dPhi_da
  ! s-derivative
  gradVdM(:,iAns,2)=(1.+a)*dPhi_da+dPhi_db
END DO 
END SUBROUTINE GradVandermondeTria

END MODULE MOD_triaBasis
