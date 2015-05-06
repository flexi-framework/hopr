#include "defines.f90"
MODULE MOD_Basis1D
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
INTERFACE Vandermonde1D
  MODULE PROCEDURE Vandermonde1D
END INTERFACE

INTERFACE GradVandermonde1D
  MODULE PROCEDURE GradVandermonde1D
END INTERFACE

INTERFACE JacobiP
  MODULE PROCEDURE JacobiP
END INTERFACE

INTERFACE GradJacobiP
  MODULE PROCEDURE GradJacobiP
END INTERFACE

INTERFACE ChebyGaussLobNodesAndWeights
   MODULE PROCEDURE ChebyGaussLobNodesAndWeights
END INTERFACE

INTERFACE PolynomialDerivativeMatrix
   MODULE PROCEDURE PolynomialDerivativeMatrix
END INTERFACE

INTERFACE BarycentricWeights 
   MODULE PROCEDURE BarycentricWeights 
END INTERFACE

INTERFACE InitializeVandermonde 
   MODULE PROCEDURE InitializeVandermonde
END INTERFACE

PUBLIC:: Vandermonde1D
PUBLIC:: GradVandermonde1D
PUBLIC:: JacobiP
PUBLIC:: GradJacobiP
PUBLIC:: ChebyGaussLobNodesAndWeights
PUBLIC:: PolynomialDerivativeMatrix
PUBLIC:: BarycentricWeights
PUBLIC:: InitializeVandermonde
!===================================================================================================================================

CONTAINS
SUBROUTINE  JacobiP(nNodes,x,alpha,beta,Deg,P)
!===================================================================================================================================
! evaluates the Nth Jacobi-polynomial at position xi, Algorithm in book of hesthaven and found in his matlab code
! The Jacobi Polynomials P_i^{(alpha,beta)}(x) are orthonormal with respect to the weighting function in the interval [-1,1]
! w(x)=(1-x)^alpha(1+x)^beta
! \int_{-1}^{1} P_i^{(alpha,beta)}(x) P_j^{(alpha,beta)}(x) w(x) dx = \delta_{ij}
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN):: nNodes  ! ?
REAL,INTENT(IN)    :: x(nNodes)  ! evaluation positions
INTEGER,INTENT(IN) :: alpha,beta ! coefficients of the weighting function w(x)=(1-x)^alpha(1+x)^beta
INTEGER,INTENT(IN) :: Deg          ! polynomial DEGREE of Jacobi Polynomial
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: P(nNodes)      ! value of Jacobi polynomial N at all positions x
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(nNodes) :: P_0,P_1  ! ?
REAL                   :: gamma0,gamma1,gammaf(1:alpha+beta+2)  ! ?
REAL                   :: aold,anew,bnew  ! ?
REAL                   :: ri,ralpha,rbeta       !temp
INTEGER                :: i,h1,h2  ! ?
!===================================================================================================================================
!fill gamma function, only for integer values, replace by real gamma function, if needed. Intrinsic gamma function only with GNU 
ralpha=REAL(alpha)
rbeta=REAL(beta)
gammaf(1:2)=1
DO i=3,alpha+beta+2
  gammaf(i)=(i-1)*gammaf(i-1)
END DO
gamma0=2.**(alpha+beta+1)/(ralpha+rbeta+1.)*gammaf(alpha+1)*gammaf(beta+1)/gammaf(alpha+beta+1)
P_0(:)=1./SQRT(gamma0)
IF(Deg.EQ.0) THEN
  P=P_0
  RETURN
END IF
gamma1=(ralpha+1.)*(rbeta+1.)/(ralpha+rbeta+3.)*gamma0
P_1(:)=0.5*((ralpha+rbeta+2.)*x(:) + (ralpha-rbeta))/SQRT(gamma1)
IF(Deg.EQ.1) THEN
  P(:)=P_1(:)
  RETURN
END IF
! a_i= 2/(2+aplha+beta)*sqrt( ( i*(i+alpha+beta)*(i+alpha)*(i+beta) ) / ( (2i+alpha+beta-1)(2i+alpha+beta+1) ) )
! b_i= (alpha**2-beta**2)/( (2i+alpha+beta)(2i+alpha+beta+2) )
! a_i for i=1
h1=alpha+beta
h2=beta*beta-alpha*alpha
aold= 2./REAL(2.+h1)*SQRT(REAL((1.+alpha)*(1+beta))/REAL(h1+3))
!start recurrence
DO i=2,Deg
  ri=REAL(i)
  h1= h1+2
  !a_i
  anew=2./REAL(h1+2)*SQRT(REAL(i*(i+alpha+beta)*(i+alpha)*(i+beta)) / REAL((h1+1)*(h1+3)) )
  !b_i
  bnew=REAL(h2)/REAL(h1*(h1+2))
  ! recurrence P(i)= ((x-b_i) P(i-1) - a_(i-1) P(i-2) )/a_i
  P(:)=((x(:)-bnew)*P_1(:)-aold*P_0(:))/anew
  P_0=P_1
  P_1=P
  aold=anew
END DO
END SUBROUTINE JacobiP


SUBROUTINE  GradJacobiP(nNodes,x,alpha,beta,Deg,GradP)
!===================================================================================================================================
! evaluates the first derivative of the Nth Jacobi-polynomial at position xi, 
! Algorithm in book of hesthaven and found in his matlab code
! The Jacobi Polynomials P_i^{(alpha,beta)}(x) are orthonormal with respect to the weighting function in the interval [-1,1]
! w(x)=(1-x)^alpha(1+x)^beta
! \int_{-1}^{1} P_i^{(alpha,beta)}(x) P_j^{(alpha,beta)}(x) w(x) dx = \delta_{ij}
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN):: nNodes  ! ?
REAL,INTENT(IN)    :: x(nNodes)  ! evaluation positions
INTEGER,INTENT(IN) :: alpha,beta ! coefficients of the weighting function w(x)=(1-x)^alpha(1+x)^beta
INTEGER,INTENT(IN) :: Deg        ! polynomial DEGREE of Jacobi Polynomial
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: GradP(nNodes)      ! value of the gradient of Jacobi polynomial N at all positions x
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(nNodes) :: P  ! ?
!===================================================================================================================================
IF(Deg.EQ.0)THEN
  gradP=0.
  RETURN
END IF
CALL JacobiP(nNodes,x,alpha+1,beta+1, Deg-1,P)
gradP=SQRT(REAL(Deg*(Deg+alpha+beta+1)))*P
END SUBROUTINE GradJacobiP


SUBROUTINE  JacobiP_all(nNodes,x,alpha,beta,Deg,P)
!===================================================================================================================================
! evaluates the Nth Jacobi-polynomial at position xi, Algorithm in book of hesthaven and found in his matlab code
! The Jacobi Polynomials P_i^{(alpha,beta)}(x) are orthonormal with respect to the weighting function in the interval [-1,1]
! w(x)=(1-x)^alpha(1+x)^beta
! \int_{-1}^{1} P_i^{(alpha,beta)}(x) P_j^{(alpha,beta)}(x) w(x) dx = \delta_{ij}
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN):: nNodes  ! ?
REAL,INTENT(IN)    :: x(nNodes)  ! evaluation positions
INTEGER,INTENT(IN) :: alpha,beta ! coefficients of the weighting function w(x)=(1-x)^alpha(1+x)^beta
INTEGER,INTENT(IN) :: Deg          ! polynomial DEGREE of Jacobi Polynomial
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: P(nNodes,0:Deg)      ! value of Jacobi polynomial N at all positions x
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL                   :: gamma0,gamma1,gammaf(1:alpha+beta+2)  ! ?
REAL                   :: aold,anew,bnew  ! ?
REAL                   :: ri,ralpha,rbeta       !temp
INTEGER                :: i,h1,h2  ! ?
!===================================================================================================================================
!fill gamma function, only for integer values, replace by real gamma function, if needed. Intrinsic gamma function only with GNU 
ralpha=REAL(alpha)
rbeta=REAL(beta)
gammaf(1:2)=1
DO i=3,alpha+beta+2
  gammaf(i)=(i-1)*gammaf(i-1)
END DO
gamma0=2.**(alpha+beta+1)/(ralpha+rbeta+1.)*gammaf(alpha+1)*gammaf(beta+1)/gammaf(alpha+beta+1)
P(:,0)=1./SQRT(gamma0)
IF(Deg.EQ.0) THEN
  RETURN
END IF
gamma1=(ralpha+1.)*(rbeta+1.)/(ralpha+rbeta+3.)*gamma0
P(:,1)=0.5*((ralpha+rbeta+2.)*x(:) + (ralpha-rbeta))/SQRT(gamma1)
IF(Deg.EQ.1) THEN
  RETURN
END IF
! a_i= 2/(2+aplha+beta)*sqrt( ( i*(i+alpha+beta)*(i+alpha)*(i+beta) ) / ( (2i+alpha+beta-1)(2i+alpha+beta+1) ) )
! b_i= (alpha**2-beta**2)/( (2i+alpha+beta)(2i+alpha+beta+2) )
! a_i for i=1
h1=alpha+beta
h2=beta*beta-alpha*alpha
aold= 2./REAL(2.+h1)*SQRT(REAL((1.+alpha)*(1+beta))/REAL(h1+3))
!start recurrence
DO i=2,Deg
  ri=REAL(i)
  h1= h1+2
  !a_i
  anew=2./REAL(h1+2)*SQRT(REAL(i*(i+alpha+beta)*(i+alpha)*(i+beta)) / REAL((h1+1)*(h1+3)) )
  !b_i
  bnew=REAL(h2)/REAL(h1*(h1+2))
  ! recurrence P(i)= ((x-b_i) P(i-1) - a_(i-1) P(i-2) )/a_i
  P(:,i)=((x(:)-bnew)*P(:,i-1)-aold*P(:,i-2))/anew
  aold=anew
END DO
END SUBROUTINE JacobiP_all


SUBROUTINE  GradJacobiP_all(nNodes,x,alpha,beta,Deg,GradP)
!===================================================================================================================================
! evaluates the first derivative of the Nth Jacobi-polynomial at position xi, 
! Algorithm in book of hesthaven and found in his matlab code
! The Jacobi Polynomials P_i^{(alpha,beta)}(x) are orthonormal with respect to the weighting function in the interval [-1,1]
! w(x)=(1-x)^alpha(1+x)^beta
! \int_{-1}^{1} P_i^{(alpha,beta)}(x) P_j^{(alpha,beta)}(x) w(x) dx = \delta_{ij}
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN):: nNodes  ! ?
REAL,INTENT(IN)    :: x(nNodes)  ! evaluation positions
INTEGER,INTENT(IN) :: alpha,beta ! coefficients of the weighting function w(x)=(1-x)^alpha(1+x)^beta
INTEGER,INTENT(IN) :: Deg        ! polynomial DEGREE of Jacobi Polynomial
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: GradP(nNodes,0:Deg)      ! value of the gradient of Jacobi polynomial N at all positions x
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
REAL,DIMENSION(nNodes,0:Deg) :: P  ! ?
INTEGER                      :: i  ! ?
!===================================================================================================================================
gradP=0.
IF(Deg.EQ.0)THEN
  RETURN
END IF
CALL JacobiP_all(nNodes,x,alpha+1,beta+1, Deg-1,P(:,1:Deg))
DO i=1,Deg
  gradP(:,i)=SQRT(REAL(i*(i+alpha+beta+1)))*P(:,i)
END DO
END SUBROUTINE GradJacobiP_all


SUBROUTINE Vandermonde1D(nNodes1D,Deg,r1D,VdM1D)
!===================================================================================================================================
! computes on a given set of nodes and a given polynomial degree
! the Vandermonde Matrix to 1D orthonormal Legendre polyonomials in reference space [-1,1]
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: nNodes1D,Deg  ! ?
REAL,INTENT(IN)              :: r1D(nNodes1D)   ! node positions in reference space [-1,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: VdM1D(0:nNodes1D-1,0:Deg)   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
CALL JacobiP_all(nNodes1D,r1D, 0, 0, Deg,Vdm1D(:,:))
END SUBROUTINE Vandermonde1D


SUBROUTINE GradVandermonde1D(nNodes1D,Deg,r1D,gradVdM1D)
!===================================================================================================================================
! computes on a given set of nodes and a given polynomial degree
! the Gradient Vandermonde Matrix to 1D orthonormal Legendre polyonomials in reference space [-1,1]
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)          :: nNodes1D,Deg  ! ?
REAL,INTENT(IN)              :: r1D(nNodes1D)   ! node positions in reference space [-1,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)             :: gradVdM1D(0:nNodes1D-1,0:Deg)   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
CALL GradJacobiP_all(nNodes1D,r1D, 0, 0,Deg,gradVdM1D(:,:))
END SUBROUTINE GradVandermonde1D

SUBROUTINE ChebyGaussLobNodesAndWeights(N_in,xGP,wGP)
!===================================================================================================================================
! algorithm 27, Kopriva
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in       ! polynomial degree, (N_in+1) CLpoints 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)          :: xGP(0:N_in)  ! Gausspoint positions for the reference interval [-1,1]
REAL,INTENT(OUT),OPTIONAL :: wGP(0:N_in)  ! Gausspoint weights
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: iGP  ! ?
!===================================================================================================================================
DO iGP=0,N_in
  xGP(iGP)=-COS(iGP/REAL(N_in)*ACOS(-1.))
END DO
IF(PRESENT(wGP))THEN
  DO iGP=0,N_in
    wGP(iGP)=ACOS(-1.)/REAL(N_in)
  END DO
  wGP(0)=wGP(0)*0.5
  wGP(N_in)=wGP(N_in)*0.5
END IF

! Debugging:

!SWRITE(*,*) "########### I am HERE: ChebyGaussLobNodesAndWeights"


END SUBROUTINE ChebyGaussLobNodesAndWeights


SUBROUTINE BarycentricWeights(N_in,xGP,wBary)
!===================================================================================================================================
! algorithm 30, Kopriva
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_in               ! polynomial degree 
REAL,INTENT(IN)    :: xGP(0:N_in)        ! Gausspoint positions for the reference interval [-1,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: wBary(0:N_in)      ! barycentric weights
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: iGP,jGP  ! ?
!===================================================================================================================================
wBary(:)=1.
DO iGP=1,N_in
  DO jGP=0,iGP-1
    wBary(jGP)=wBary(jGP)*(xGP(jGP)-xGP(iGP))
    wBary(iGP)=wBary(iGP)*(xGP(iGP)-xGP(jGP))
  END DO ! jGP
END DO ! iGP
wBary(:)=1./wBary(:)
END SUBROUTINE BarycentricWeights


SUBROUTINE PolynomialDerivativeMatrix(N_in,xGP,D)
!===================================================================================================================================
! algorithm 37, Kopriva
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_in              ! polynomial degree
REAL,INTENT(IN)    :: xGP(0:N_in)       ! Gausspoint positions for the reference interval [-1,1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: D(0:N_in,0:N_in)     ! differentiation Matrix
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: iGP,iLagrange  ! ?
REAL               :: wBary(0:N_in)   ! ?
!===================================================================================================================================
CALL BarycentricWeights(N_in,xGP,wBary)
D(:,:)=0.
DO iLagrange=0,N_in
  DO iGP=0,N_in
    IF(iLagrange.NE.iGP)THEN
      D(iGP,iLagrange)=wBary(iLagrange)/(wBary(iGP)*(xGP(iGP)-xGP(iLagrange)))
      D(iGP,iGP)=D(iGP,iGP)-D(iGP,iLagrange)
    END IF ! (iLagrange.NE.iGP)
  END DO ! iGP
END DO ! iLagrange
END SUBROUTINE PolynomialDerivativeMatrix


SUBROUTINE LagrangeInterpolationPolys(x,N_in,xGP,wBary,L)
!===================================================================================================================================
! Algorithm 34, Kopriva
! Computes all Lagrange functions evaluated at position x in [-1;1]
! Uses function ALMOSTEQUAL
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)   :: x          ! Coordinate
INTEGER,INTENT(IN) :: N_in          ! polynomial degree
REAL,INTENT(IN)    :: xGP(0:N_in)   ! Gausspoint positions for the reference interval [-1,1]
REAL,INTENT(IN)    :: wBary(0:N_in) ! Barycentric weights
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: L(0:N_in)     ! Lagrange basis functions evaluated at x
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                   :: iGP  ! ?
LOGICAL                   :: xEqualGP ! is x equal to a Gauss Point
REAL                      :: DummySum  ! ?
!===================================================================================================================================
xEqualGP=.FALSE.
DO iGP=0,N_in
  L(iGP)=0.
  IF(ALMOSTEQUAL(x,xGP(iGP))) THEN
    L(iGP)=1.
    xEqualGP=.TRUE.
  END IF ! (ALMOSTEQUAL(x,xGP(iGP)))
END DO ! iGP
! if x is equal to a Gauss point, L=(0,....,1,....0)
IF(xEqualGP) RETURN
DummySum=0.
DO iGP=0, N_in
  L(iGP)=wBary(iGP)/(x-xGP(iGP))
  DummySum=DummySum+L(iGP)
END DO

DO iGP=0,N_in
  L(iGP)=L(iGP)/DummySum
END DO

END SUBROUTINE LagrangeInterpolationPolys


SUBROUTINE InitializeVandermonde(N_In,N_Out,wBary_In,xi_In,xi_Out,Vdm)
!===================================================================================================================================
! build a 1D Vandermonde matrix using the lagrange basis functions of degree
! N_In, evaluated at the interpolation points xi_Out
!===================================================================================================================================
! MODULES
!USE MOD_Interpolation_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: N_In,N_Out  ! ?
REAL,INTENT(IN)    :: xi_In(0:N_In)  ! ?
REAL,INTENT(IN)    :: xi_Out(0:N_Out)  ! ?
REAL,INTENT(IN)    :: wBary_In(0:N_In)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: Vdm(0:N_Out,0:N_In)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER            :: iXi  ! ?
!===================================================================================================================================
DO iXi=0,N_Out
  CALL LagrangeInterpolationPolys(xi_Out(iXi),N_In,xi_In,wBary_In,Vdm(iXi,:))
!l(0:N_In)
END DO
END SUBROUTINE InitializeVandermonde

FUNCTION ALMOSTEQUAL(x,y)
!===================================================================================================================================
! Based on Algorithm 139, Kopriva
! Compares two real numbers
! Depends on PP_RealTolerance
! Takes into account that x,y is located in-between [-1;1]
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: x,y         ! 2 scalar real numbers
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL         :: AlmostEqual ! TRUE if |x-y| < 2*PP_RealTolerance
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================

AlmostEqual=.FALSE.
IF((x.EQ.0.).OR.(y.EQ.0.)) THEN
  IF(ABS(x-y).LE.2.*RealTolerance) AlmostEqual=.TRUE.
ELSE ! x, y not zero
  IF((ABS(x-y).LE.RealTolerance*ABS(x)).AND.((ABS(x-y).LE.RealTolerance*ABS(y)))) AlmostEqual=.TRUE.
END IF ! x,y zero

END FUNCTION ALMOSTEQUAL

END MODULE MOD_Basis1D
