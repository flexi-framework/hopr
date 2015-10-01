#include "defines.f90"
PROGRAM YPlusEstimator
!===================================================================================================================================
! tool to estimate first cell height from flow parameters
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Basis1D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! PROGRAM VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER          :: PP_N
REAL             :: dens,vel,mu,Re,Cf,L,height
REAL             :: yp,ypLG,ypLGL,ypCG,ypCGL,ypEqui
REAL,ALLOCATABLE :: xLG(:),xLGL(:),xCG(:),xCGL(:),w(:)
CHARACTER(LEN=1) :: yesno
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' HOPR Y+ ESTIMATOR'
WRITE(UNIT_stdOut,'(132("="))')


WRITE(UNIT_stdOut,'(A)') ' Enter the target polynomial degree of your simulation:'
READ(UNIT_stdIn,*) PP_N

WRITE(UNIT_stdOut,'(A)') ' Density:'
READ(UNIT_stdIn,*) dens
dens=ABS(dens)

WRITE(UNIT_stdOut,'(A)') ' Reference length:'
READ(UNIT_stdIn,*) L
L=ABS(L)

WRITE(UNIT_stdOut,'(A)') ' Cell height estimate based on Reynolds number (1) or flow parameters (2):'
READ(UNIT_stdIn,*) yesno

IF(yesno.EQ.'1')THEN
  WRITE(UNIT_stdOut,'(A)') ' Reynolds:'
  READ(UNIT_stdIn,*) Re
ELSE
  WRITE(UNIT_stdOut,'(A)') ' Flow velocity:'
  READ(UNIT_stdIn,*) vel
  vel=ABS(vel)

  WRITE(UNIT_stdOut,'(A)') ' Dynamic viscosity:'
  READ(UNIT_stdIn,*) mu
  mu =ABS(mu)

  Re=dens*vel*L/mu
END IF
Re=ABS(Re)


Cf=0.0260*Re**(-1.0/7.0)
Cf=0.0576*Re**(-1.0/5.0)
Cf=0.370*(LOG(Re)/LOG(10.))     **(-2.584)
Cf=(2.   *LOG(Re)/LOG(10.)-0.65)**(-2.3)
!tau_w = Cf*0.5*dens*vel*vel
!Uf = SQRT(tau_w/dens)
!yp=mu/(Uf*dens)

yp = L/Re*SQRT(2./Cf)

ALLOCATE(xLG(0:PP_N),xLGL(0:PP_N),xCG(0:PP_N),xCGL(0:PP_N),w(0:PP_N))
CALL LegendreGaussNodesAndWeights( PP_N,xLG, w)
CALL LegGaussLobNodesAndWeights(   PP_N,xLGL,w)
CALL ChebyshevGaussNodesAndWeights(PP_N,xCG, w)
CALL ChebyGaussLobNodesAndWeights( PP_N,xCGL,w)

ypEqui=yp*PP_N
ypLG  =yp/((xLG( 0)+1.)*0.5) 
ypLGL =yp/((xLGL(1)+1.)*0.5)
ypCG  =yp/((xCG( 0)+1.)*0.5)
ypCGL =yp/((xCGL(1)+1.)*0.5)
DEALLOCATE(xLG,xLGL,xCG,xCGL,w)

WRITE(UNIT_stdOut,'(A)')
WRITE(UNIT_stdOut,'(A)') ' CAUTION: The cell height estimates are valid for linear meshes only!'
WRITE(UNIT_stdOut,'(A)') '          Curved (especially internally stretched) grid cells and/or '
WRITE(UNIT_stdOut,'(A)') '          overintegration and filters will strongly influence the result !'
WRITE(UNIT_stdOut,'(A)')
WRITE(UNIT_stdOut,'(A)') ' Height of first cell for y+=1 for various node distributions:'
WRITE(UNIT_stdOut,'(A)')

WRITE(UNIT_stdOut,'(A,ES10.3)') ' First cell point        :  h=',yp
WRITE(UNIT_stdOut,'(A,ES10.3)') ' Equidistant             :  h=',ypEqui
WRITE(UNIT_stdOut,'(A,ES10.3)') ' Legendre-Gauss          :  h=',ypLG
WRITE(UNIT_stdOut,'(A,ES10.3)') ' Legendre-Gauss-Lobatto  :  h=',ypLGL
WRITE(UNIT_stdOut,'(A,ES10.3)') ' Chebyshev-Gauss         :  h=',ypCG
WRITE(UNIT_stdOut,'(A,ES10.3)') ' Chebyshev-Gauss-Lobatto :  h=',ypCGL

WRITE(UNIT_stdOut,'(A)')
WRITE(UNIT_stdOut,'(A)') ' If desired enter your grids cell height to compute Y+ for the grid:'
READ(UNIT_stdIn,*) height
WRITE(UNIT_stdOut,'(A)')

WRITE(UNIT_stdOut,'(A,ES10.3)') ' Equidistant             : y+=',height/ypEqui
WRITE(UNIT_stdOut,'(A,ES10.3)') ' Legendre-Gauss          : y+=',height/ypLG
WRITE(UNIT_stdOut,'(A,ES10.3)') ' Legendre-Gauss-Lobatto  : y+=',height/ypLGL
WRITE(UNIT_stdOut,'(A,ES10.3)') ' Chebyshev-Gauss         : y+=',height/ypCG
WRITE(UNIT_stdOut,'(A,ES10.3)') ' Chebyshev-Gauss-Lobatto : y+=',height/ypCGL

END PROGRAM YPlusEstimator

