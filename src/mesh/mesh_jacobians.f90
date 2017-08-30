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
MODULE MOD_Mesh_Jacobians
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
INTERFACE CheckJacobians
  MODULE PROCEDURE CheckJacobians
END INTERFACE

PUBLIC:: CheckJacobians
!===================================================================================================================================

CONTAINS
SUBROUTINE CheckJacobians()
!===================================================================================================================================
! Solve the Jacobian determinant and sharp boundaries for J_min and J_max
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars ,ONLY: FirstElem,nMeshElems
USE MOD_Mesh_Vars ,ONLY: tElem
USE MOD_Mesh_Vars ,ONLY: N,negativeJacobians
USE MOD_Mesh_Vars ,ONLY: JacobianTolerance
USE MOD_Mesh_Basis,ONLY: PackGeo
USE MOD_Basis_Vars,ONLY:nAnalyze
USE MOD_Basis1D   ,ONLY: LegGaussLobNodesAndWeights
USE MOD_Basis1D   ,ONLY: PolynomialDerivativeMatrix
USE MOD_Basis1D   ,ONLY: BarycentricWeights
USE MOD_Basis1D   ,ONLY: InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iElem,i,iNode  ! ?
REAL                 :: maxJac,minJac  ! ?
REAL                 :: xBary(3)
REAL                 :: scaledJac(nMeshElems)  ! ?
INTEGER              :: scaledJacStat(0:10)  ! ?
TYPE(tElem),POINTER  :: aElem  ! ?
REAL                 :: xEq(0:N),wBaryEq(0:N)
REAL                 :: xGL(0:N),wBaryGL(0:N),DGL(0:N,0:N),VdmEqToGL(0:N,0:N),D_EqToGL(0:N,0:N)
REAL                 :: xAP(0:nAnalyze)
REAL                 :: VdmGLToAP(0:nAnalyze,0:N)
REAL                 :: xGeo(3,0:N,0:N,0:N)
REAL                 :: Jac(0:nAnalyze,0:nAnalyze,0:nAnalyze)
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)') 'Check Elements Jacobian...'

!==== 1D basis
DO i=0,N
  xEq(i)=2.*REAL(i)/REAL(N) -1.
END DO !i
! take derivatives on Gauss-Lobatto points
! interpolate Eq (N)-> GL (N), derivative on GL (N) 

CALL LegGaussLobNodesAndWeights(N,xGL)
CALL PolynomialDerivativeMatrix(N,xGL,DGL)
CALL BarycentricWeights(N,xEq,wBaryEq)
CALL InitializeVandermonde(N,N,wBaryEq,xEq,xGL,VdmEqtoGL)

D_EqToGL=MATMUL(DGL,VdmEqToGL)

!interpolate derivatives on GL (N) to nAnalyze points
DO i=0,nAnalyze
  xAP(i)=2.*REAL(i)/REAL(nAnalyze) -1.
END DO !i
CALL BarycentricWeights(N,xGL,wBaryGL)
CALL InitializeVandermonde(N,nAnalyze,wBaryGL,xGL,xAP,VdmGLtoAP)

!=== 1D basis finished

scaledJac=1.
iElem=0
! Allocate Curved Node Coordinates, Array for Jacobians_Determinants J, Bezier Coefficients B
! Loop over all Elements
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  IF (aElem%nNodes.NE. 8) THEN
    WRITE(*,*)' WARNING: checkJacobian only for Hexas...',aElem%nnodes
    aElem=>aElem%nextElem
    CYCLE
  END IF ! Only Elements with nNodes
  iElem=iElem+1 ! Element Counter

  CALL PackGeo(N,aElem,Xgeo)
  Jac(:,:,:) = EvalJac(Xgeo)
  maxJac=MAXVAL(ABS(Jac))
  minJac=MINVAL(Jac)
   
  scaledJac(iElem)=minJac/maxJac
  ! Check that Nodal Jacobians are positive
  !DO k=0,nAnalyze; DO j=0,nAnalyze; DO i=0,nAnalyze
  !  IF (Jac(i,j,k)/maxJac.LE.0.01) THEN
  !     WRITE(*,*)'Negative Jacobian, iELem=',iElem,'Jac=',Jac(i,j,k),'maxJac',maxJac
  !  END IF
  !END DO; END DO; END DO !i,j,k=0,nAnalyze

  aElem=>aElem%nextElem
END DO
IF (iElem.NE.nMeshElems) THEN
  WRITE(*,*) 'iElem /= nMeshElems',iElem,nMeshElems
END IF

! Error Section
IF(ANY(scaledJac.LE.jacobianTolerance))THEN
  OPEN(UNIT=100,FILE='Jacobian_Error.out',STATUS='UNKNOWN',ACTION='WRITE')
  WRITE(100,*) 'Corrupt Elems Found, i.e. the Jacobian is negative, '
  WRITE(100,*)
  WRITE(100,'(4(A12,4X))') 'iElem','xBary','yBary','zBary'
  WRITE(100,*)
  iElem=0
  aElem=>firstElem
  DO WHILE(ASSOCIATED(aElem))
    iElem=iElem+1
    IF(scaledJac(iElem).LT.jacobianTolerance) THEN
      xBary(:)=0.
      DO iNode=1,aElem%nNodes
        xBary(:)=xBary(:)+aElem%node(inode)%np%x
      END DO
      xBary(:)=xBary(:)/REAL(aElem%nNodes)
      WRITE(100,'(I12,3(4X,F12.8))') iElem, xBary(:)
      negativeJacobians=negativeJacobians+1
    END IF
    aElem=>aElem%nextElem
  END DO
  CLOSE(UNIT=100)
  WRITE(*,*) 'WARNING!!! Elements with negative Jacobian found, see Jacobian_Error.out for details.'
END IF
scaledJacStat(:)=0

DO iElem=1,nMeshElems
  i=CEILING(MAX(0.,scaledJac(iElem)*10))
  scaledJacStat(i)=scaledJacStat(i)+1 
END DO
WRITE(Unit_StdOut,'(A)') ' Number of element with scaled Jacobians ranging between:'
WRITE(Unit_StdOut,'(A)') '   <  0.0  <  0.1  <  0.2  <  0.3  <  0.4  <  0.5  <  0.6  <  0.7  <  0.8  <  0.9  <  1.0 '
DO i=0,10
  WRITE(Unit_StdOut,'(I6,X,A1)',ADVANCE='NO')scaledJacStat(i),'|'
END DO
WRITE(Unit_StdOut,'(A1)')' '
 

CALL Timer(.FALSE.)
WRITE(UNIT_StdOut,'(132("="))')
WRITE(UNIT_StdOut,*)

CONTAINS

  FUNCTION EvalJac(Xgeo_in) RESULT(Jac_out)
    USE MOD_ChangeBasis,ONLY:ChangeBasis3D
    !uses N,nAnalyze, VdmGLToAP and D_EqToGL from main routine
    IMPLICIT NONE
    !-------------------------------------------------------------
    !INPUT/ OUTPUT VARIABLES
    REAL,INTENT(IN)    :: Xgeo_in(3,0:N,0:N,0:N)
    REAL               :: Jac_out(0:nAnalyze,0:nAnalyze,0:nAnalyze)
    !-------------------------------------------------------------
    !LOCAL VARIABLES
    REAL    :: dXdxiGL  (3,0:N,0:N,0:N)
    REAL    :: dXdetaGL (3,0:N,0:N,0:N)
    REAL    :: dXdzetaGL(3,0:N,0:N,0:N)
    REAL    :: dXdxiAP  (3,0:nAnalyze,0:nAnalyze,0:nAnalyze)
    REAL    :: dXdetaAP (3,0:nAnalyze,0:nAnalyze,0:nAnalyze)
    REAL    :: dXdzetaAP(3,0:nAnalyze,0:nAnalyze,0:nAnalyze)
    INTEGER :: i,j,k,l
    !-------------------------------------------------------------
    dXdxiGL  =0.
    dXdetaGL =0.
    dXdzetaGL=0.
    DO k=0,N; DO j=0,N; DO i=0,N
      DO l=0,N
        dXdxiGL  (:,i,j,k) = dXdxiGL  (:,i,j,k) + D_EqToGL(i,l)*Xgeo_in(:,l,j,k)
        dXdetaGL (:,i,j,k) = dXdetaGL (:,i,j,k) + D_EqToGL(j,l)*Xgeo_in(:,i,l,k)
        dXdzetaGL(:,i,j,k) = dXdzetaGL(:,i,j,k) + D_EqToGL(k,l)*Xgeo_in(:,i,j,l)
      END DO !l=0,N
    END DO; END DO; END DO !i,j,k=0,N
    CALL ChangeBasis3D(3,N,nAnalyze,VdmGLToAP,dXdxiGL  ,dXdxiAP  )
    CALL ChangeBasis3D(3,N,nAnalyze,VdmGLToAP,dXdetaGL ,dXdetaAP )
    CALL ChangeBasis3D(3,N,nAnalyze,VdmGLToAP,dXdzetaGL,dXdzetaAP)
    DO k=0,nAnalyze; DO j=0,nAnalyze; DO i=0,nAnalyze
      Jac_out(i,j,k) = SUM(dXdxiAP(:,i,j,k)*CROSS(dXdetaAP(:,i,j,k),dXdzetaAP(:,i,j,k)))
    END DO; END DO; END DO !i,j,k=0,nAnalyze
  END FUNCTION EvalJac

END SUBROUTINE CheckJacobians

END MODULE MOD_Mesh_Jacobians
