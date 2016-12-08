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
MODULE MOD_CCInt
!===================================================================================================================================
! This module contains routines to integrate 1D functions recursively until a given tolerance, using clenshaw-curtis quadrature
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

INTERFACE CCint_Init
  MODULE PROCEDURE CCint_Init 
END INTERFACE

INTERFACE CCint 
  MODULE PROCEDURE CCint 
END INTERFACE

ABSTRACT INTERFACE 
  FUNCTION i_fxn(n,x) RESULT (yn)
    IMPLICIT NONE
    INTEGER :: n
    REAL :: x(n)
    REAL :: yn(n)
  END FUNCTION i_fxn

END INTERFACE

!===================================================================================================================================

CONTAINS


SUBROUTINE CCint_Init()
!===================================================================================================================================
! do a nested integration using clenshaw curtis recursive quadrature
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Basis1D, ONLY: ClenshawCurtisNodesAndWeights
USE MOD_CCInt_vars,ONLY: Imax,Rcc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!REAL, INTENT(IN)   :: a,b
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!REAL, INTENT(OUT)   :: fi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE tCC 
  INTEGER           :: N
  REAL,ALLOCATABLE  :: x(:)
  REAL,ALLOCATABLE  :: w(:)
END TYPE tCC
TYPE(tCC),ALLOCATABLE:: CC(:)
INTEGER             :: s,i,j,k,Nr
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'  INIT CLENSHAW-CURTIS ...'

Imax=7
ALLOCATE(CC(1:Imax))
Nr=2
DO i=1,Imax
  CC(i)%N=Nr
  ALLOCATE(CC(i)%x(0:Nr))
  ALLOCATE(CC(i)%w(0:Nr))
  CALL ClenshawCurtisNodesAndWeights(Nr,CC(i)%x(0:Nr),CC(i)%w(0:Nr))
  !inteval [0,1]
  CC(i)%x(0:Nr)=0.5*(1.+CC(i)%x(0:Nr))
  CC(i)%w(0:Nr)=0.5*CC(i)%w(0:Nr)
  Nr=Nr*2
END DO !i=0

ALLOCATE(Rcc(1:Imax))

! start with mid-point rule
s=1
Rcc(1)%N=2
Rcc(1)%np=3
ALLOCATE(Rcc(1)%x(3))
Rcc(1)%x(:)=CC(1)%x(:)
ALLOCATE(Rcc(1)%w(1:Imax,3))
DO j=1,Imax
  Rcc(1)%w(j,1)=CC(j)%w(0)
  Rcc(1)%w(j,2)=CC(j)%w(CC(j)%N/2)
  Rcc(1)%w(j,3)=CC(j)%w(CC(j)%N)
END DO
!recursive integration rules
DO s=2,Imax
  Nr=CC(s)%N
  Rcc(s)%N=Nr
  Rcc(s)%np=Nr/2
  ALLOCATE(Rcc(s)%x(Rcc(s)%np))
  ALLOCATE(Rcc(s)%w(s:Imax,Rcc(s)%np))
  DO k=1,Rcc(s)%np
    Rcc(s)%x(k)=CC(Imax)%x((2*k-1)*CC(Imax)%N/Nr)
    DO j=s,Imax
      Rcc(s)%w(j,k)=CC(j)%w((2*k-1)*CC(j)%N/Nr)
    END DO
  END DO !k
END DO 


!s=1
!  Nr=Rcc(s)%N
!  WRITE(*,'(A,I3,A,I3,A,I5)')'=== DEBUG,s=',s,' N= ', Rcc(s)%N, ' np= ',Rcc(s)%np
!    WRITE(*,'(A,I3,A,I5,A)')'DEBUG,x(k=',1,' )= xN(', 0,')'
!    DO j=s,Imax
!      WRITE(*,'(A,I3,A,I3,A,I3,A,I5,A)')'DEBUG,w(j=',j,',k= ',1,' )= w_',CC(j)%N,'(', 0,')'
!    END DO
!    WRITE(*,'(A,I3,A,I5,A)')'DEBUG,x(k=',1,' )= xN(', CC(Imax)%N/2,')'
!    DO j=s,Imax
!      WRITE(*,'(A,I3,A,I3,A,I3,A,I5,A)')'DEBUG,w(j=',j,',k= ',2,' )= w_',CC(j)%N,'(', CC(j)%N/2,')'
!    END DO
!    WRITE(*,'(A,I3,A,I5,A)')'DEBUG,x(k=',3,' )= xN(', CC(Imax)%N,')'
!    DO j=s,Imax
!      WRITE(*,'(A,I3,A,I3,A,I3,A,I5,A)')'DEBUG,w(j=',j,',k= ',3,' )= w_',CC(j)%N,'(', CC(j)%N,')'
!    END DO
!DO s=2,Imax
!  Nr=Rcc(s)%N
!  WRITE(*,'(A,I3,A,I3,A,I5)')'=== DEBUG,s=',s,' N= ', Rcc(s)%N, ' np= ',Rcc(s)%np
!  DO k=1,Rcc(s)%np
!    WRITE(*,'(A,I3,A,I5,A)')'DEBUG,x(k=',k,' )= xN(', (2*k-1)*CC(Imax)%N/Nr,')'
!    DO j=s,Imax
!      WRITE(*,'(A,I3,A,I3,A,I3,A,I5,A)')'DEBUG,w(j=',j,',k= ',k,' )= w_',CC(j)%N,'(', (2*k-1)*CC(j)%N/Nr,')'
!    END DO
!  END DO
!END DO

DEALLOCATE(CC)

CALL CCintTest()

WRITE(UNIT_stdOut,'(A)')'  DONE.'
END SUBROUTINE CCInt_Init


FUNCTION CCint(tol,FINT,converged) RESULT(res)
!===================================================================================================================================
! do a nested integration using clenshaw curtis recursive quadrature
! each stage s, the function is evaluated at the points of this stage (1st stage 3 points, i'th stage 2^(i-1) midpoints)
! and then multiplied with the weights of ALL stages! this way the function value does not need to be stored and 
! that for each stage, one only needs the values at the new mid-points.
! integration function has to be already transformed to have limits [0,1], since integration
! points and weights are defined in [0,1]
!===================================================================================================================================
! MODULES
USE MOD_CCInt_vars,ONLY: Imax,Rcc
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)   :: tol !tolerance
PROCEDURE(i_fxn)   :: FINT  ! function to be integrated FINT(1:n)=FINT(n,x)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: res  !integration result
LOGICAL, INTENT(OUT):: converged
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: s
REAL                :: ibuf(1:Imax)
!===================================================================================================================================

ibuf(1:Imax)=MATMUL(Rcc(1)%w(1:Imax,:),FINT(Rcc(1)%np,Rcc(1)%x(:))) !boundaries
DO s=2,Imax
  ibuf(s:Imax)=ibuf(s:Imax)+MATMUL(Rcc(s)%w(s:Imax,:),FINT(Rcc(s)%np,Rcc(s)%x(:)))
  !WRITE(*,*)'s= ',s ,'err',ABS(ibuf(s)-ibuf(s-1))
  converged=(ABS(ibuf(s)-ibuf(s-1))).LT.tol
  IF(converged)EXIT
END DO!s 
res=ibuf(s)
IF(.NOT.converged) THEN
  WRITE(*,'(A,I4)') 'WARNING RECURSIVE CC INT: HIGHEST LEVEL REACHED WITHOUT CONVERGENCE, N=',2**Imax
END IF

END FUNCTION CCInt

SUBROUTINE CCintTest()
!===================================================================================================================================
! do a nested integration using clenshaw curtis recursive quadrature
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: a,b
REAL                :: tol,res,res_ex
LOGICAL             :: converged
!===================================================================================================================================
a=0.1
b=1.9
tol=1.0E-10
res= (b-a)*CCint(tol/(b-a),FI,converged)

IF(.NOT.converged) THEN
  STOP 'TEST CCINT: Imax reached without convergence'
ELSE
  !FI=(b-a)*( (x*1.1-0.48)**13+(x*1.3+0.33)**14+(x-0.3)**3 )
  !res=(b-a)*( (x*1.1-0.48)**14/(14.*1.1)+ (x*1.3+0.33)**15/(15.*1.3)+(x-0.3)**4/4. )
  res_ex= ( (b*1.1-0.48)**14/(14.*1.1)+ (b*1.3+0.33)**15/(15.*1.3)+(b-0.3)**4/4. ) &
        -( (a*1.1-0.48)**14/(14.*1.1)+ (a*1.3+0.33)**15/(15.*1.3)+(a-0.3)**4/4. )
  IF(ABS(res-res_ex).LT.10*tol) THEN
    WRITE(*,*) 'TEST CCINT SUCCESSFUL: tol= ',tol, '|r-r_ex|= ',ABS(res-res_ex),' < 10*tol' 
  ELSE
    WRITE(*,'(A,E21.11)') ' CCINT NOT CORRECT', ABS(res-res_ex)
    STOP 
  END IF
END IF

CONTAINS 

  FUNCTION FI(np,x)
    IMPLICIT NONE
    INTEGER:: np
    REAL   :: x(1:np) !always [0,1]
    REAL   :: FI(1:np)
    !local
    REAL   :: xloc(1:np)
    xloc=a+(b-a)*x
    FI= (xloc*1.1-0.48)**13+(xloc*1.3+0.33)**14+(xloc-0.3)**3 
  END FUNCTION FI

END SUBROUTINE CCIntTest

END MODULE MOD_CCInt
