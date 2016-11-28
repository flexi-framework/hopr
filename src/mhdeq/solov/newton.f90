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
MODULE MOD_Newton
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE NewtonMin1D 
  MODULE PROCEDURE NewtonMin1D 
END INTERFACE

INTERFACE NewtonRoot1D 
  MODULE PROCEDURE NewtonRoot1D 
END INTERFACE

INTERFACE NewtonMin2D 
  MODULE PROCEDURE NewtonMin2D 
END INTERFACE

ABSTRACT INTERFACE 
  FUNCTION i_f1x1(x) RESULT (y1x1)
    IMPLICIT NONE
    REAL :: x
    REAL :: y1x1
  END FUNCTION i_f1x1

  FUNCTION i_f1x2(x) RESULT (y1x2)
    IMPLICIT NONE
    REAL :: x(2)
    REAL :: y1x2
  END FUNCTION i_f1x2

  FUNCTION i_f2x2(x) RESULT (y2x2)
    IMPLICIT NONE
    REAL :: x(2)
    REAL :: y2x2(2)
  END FUNCTION i_f2x2

  FUNCTION i_f22x2(x) RESULT (y22x2)
    IMPLICIT NONE
    REAL :: x(2)
    REAL :: y22x2(2,2)
  END FUNCTION i_f22x2
END INTERFACE

!===================================================================================================================================

CONTAINS


FUNCTION NewtonMin1D(tol,a,b,x,FF,dFF,ddFF) RESULT (fmin)
!===================================================================================================================================
! Newton's iterative algorithm to find the minimimum of function f(x) in the interval [a,b], using df(x)=0 and the derivative 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: tol
REAL,INTENT(IN)    :: a,b
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: x    !initial guess on input, result on output
PROCEDURE(i_f1x1)  :: FF   !functional f(x) to minimize
PROCEDURE(i_f1x1)  :: dFF  !d/dx f(x)
PROCEDURE(i_f1x1)  :: ddFF !d^2/dx^2 f(x)
REAL               :: fmin !on output =f(x)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: x0
!===================================================================================================================================
x0=x
x=NewtonRoot1D(tol,a,b,x0,0.,dFF,ddFF)
fmin=FF(x)

END FUNCTION NewtonMin1D


FUNCTION NewtonRoot1D(tol,a,b,xin,F0,FR,dFR) RESULT (xout)
!===================================================================================================================================
! Newton's iterative algorithm to find the root of function FR(x(:)) in the interval [a(:),b(:)], using d/dx(:)F(x)=0 and the derivative 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL   ,INTENT(IN) :: tol
REAL   ,INTENT(IN) :: a,b
REAL   ,INTENT(IN) :: F0     ! function to find root is FR(x)-F0 
REAL   ,INTENT(IN) :: xin    !initial guess on input, result on output
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
PROCEDURE(i_f1x1)  :: FR      ! function to find root
PROCEDURE(i_f1x1)  :: dFR     ! multidimensional derivative d/dx f(x), size dim1
REAL               :: xout    !on output =f(x)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iter,maxiter
REAL                :: x,dx
REAL                :: tolx
LOGICAL             :: converged
!===================================================================================================================================

converged=.FALSE.
x=xin
tolx=tol*ABS(x)
maxiter=50
DO iter=1,maxiter
  dx=-(FR(x)-F0)/dFR(x)
  dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
  x = x+dx
  converged=(ABS(dx).LT.tolx).AND.(x.GT.a).AND.(x.LT.b)
  IF(converged) EXIT
END DO !iter
IF(.NOT.converged) STOP 'NewtonRoot1D not converged'
xout=x

END FUNCTION NewtonRoot1D


FUNCTION NewtonMin2D(tol,a,b,x,FF,dFF,ddFF) RESULT (fmin)
!===================================================================================================================================
! Newton's iterative algorithm to find the minimimum of function f(x,y) in the interval x(i)[a(i),b(i)],
! using grad(f(x)=0 and the derivative 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN):: tol
REAL,INTENT(IN):: a(2),b(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT) :: x(2) !initial guess on input, result on output
PROCEDURE(i_f1x2)  :: FF   !functional f(x,y) to minimize
PROCEDURE(i_f2x2)  :: dFF  !d/dx f(x,y),d/dyf(x,y)
PROCEDURE(i_f22x2) :: ddFF !d^2/dx^2 f(x)
REAL               :: fmin !on output =f(x,y)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iter,maxiter
REAL                :: dx(2)
REAL                :: det_Hess
REAL                :: gradF(2),Hess(2,2),HessInv(2,2)
REAL                :: tolx
LOGICAL             :: converged
!===================================================================================================================================
converged=.FALSE.
tolx=tol*SQRT(SUM(x*x))
maxiter=50
DO iter=1,maxiter
  Hess=ddFF(x)
  det_Hess = Hess(1,1)*Hess(2,2)-Hess(1,2)*Hess(2,1)
  IF(det_Hess.LT.1.0E-12) STOP 'det Hessian=0 in NewtonMin'
  HessInv(1,1)= Hess(2,2)
  HessInv(1,2)=-Hess(1,2)
  HessInv(2,1)=-Hess(2,1)
  HessInv(2,2)= Hess(1,1)
  HessInv=HessInv/det_Hess
  gradF=dFF(x) 
  dx=-MATMUL(HessInv,gradF)
  dx = MAX(-(x-a),MIN(b-x,dx)) !respect bounds
  x = x+dx
  converged=(SQRT(SUM(dx*dx)).LT.tolx).AND.ALL(x.GT.a).AND.ALL(x.LT.b)
  IF(converged) EXIT
END DO !iter
IF(.NOT.converged) STOP 'NewtonMin not converged'
fmin=FF(x)

END FUNCTION NewtonMin2D

END MODULE MOD_Newton
