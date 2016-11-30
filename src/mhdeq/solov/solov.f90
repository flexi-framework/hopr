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
MODULE MOD_Solov
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

INTERFACE InitSolov 
  MODULE PROCEDURE InitSolov 
END INTERFACE

! allow different dimensions of input/output arrays
!INTERFACE MapToSolov 
!  MODULE PROCEDURE MapToSolov 
!END INTERFACE

PUBLIC::InitSolov
PUBLIC::MapToSolov
!===================================================================================================================================

CONTAINS

SUBROUTINE InitSolov 
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals, ONLY:UNIT_StdOut
USE MOD_ReadInTools
USE MOD_Solov_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'  INIT SOLOV INPUT ...'

setup =GETINT("Solov_setup","0")
IF(setup.EQ.0)THEN
  P_R0    =GETREAL("Solov_R0")
  P_eps   =GETREAL("Solov_eps")
  P_kappa =GETREAL("Solov_kappa")
  P_delta =GETREAL("Solov_delta")
  P_A     =GETREAL("Solov_A")
  P_B0    =GETREAL("Solov_B0")
  P_qaxis =GETREAL("Solov_qaxis")
  P_paxis =GETREAL("Solov_paxis")
ELSE
  !preselected setups (defaults are set, can still be modified in inifile)
  SELECT CASE(setup)
  CASE(1) !circular
    WRITE(*,*)'circular equilibrium setup:'
    P_R0    =GETREAL("Solov_R0","10.")
    P_eps   =GETREAL("Solov_eps","0.4")
    P_kappa =GETREAL("Solov_kappa","1.")
    P_delta =GETREAL("Solov_delta","0.")
    P_A     =GETREAL("Solov_A","0.955")
    P_B0    =GETREAL("Solov_B0","1.0")
    P_qaxis =GETREAL("Solov_qaxis","1.89")
    P_paxis =GETREAL("Solov_paxis","2.0e-03")
  CASE(2) !parameter ITER-like
    WRITE(*,*)'Iter-like equilibrium setup:'
    P_R0    =GETREAL("Solov_R0","6.2")
    P_eps   =GETREAL("Solov_eps","0.32")
    P_kappa =GETREAL("Solov_kappa","1.7")
    P_delta =GETREAL("Solov_delta","0.33")
    P_A     =GETREAL("Solov_A","-0.155")
    P_B0    =GETREAL("Solov_B0","1.0")
    P_qaxis =GETREAL("Solov_qaxis","1.6")
    P_paxis =GETREAL("Solov_paxis","8.0e-02")
  CASE DEFAULT
    WRITE(*,*) "WARNING: This soloviev setup does not exist" ,setup
    STOP
  END SELECT
END IF !setup=0
IF(p_A.GT.1.) WRITE(*,*) 'WARNING, USING POSITIVE PRESSURE GRADIENT p_axis< PresEdge'
IF(p_eps.GT.0.98) WRITE(*,*) 'WARNING, EPSILON >0.98, makes no sense!'
IF(p_delta.GT.0.98) WRITE(*,*) 'WARNING, DELTA >0.98, makes no sense!'
  
asin_delta=ASIN(P_delta)

CALL InitSolovievEquilibrium()

nRhoCoefs=GETINT("nRhoCoefs","0")
IF(nRhoCoefs.GT.0)THEN
  ALLOCATE(RhoCoefs(nRhoCoefs))
  RhoCoefs=GETREALARRAY("RhoCoefs",nRhoCoefs)
END IF
WRITE(UNIT_stdOut,'(A)')'  ... DONE'
END SUBROUTINE InitSolov


SUBROUTINE InitSolovievEquilibrium 
!===================================================================================================================================
! Solve for psiCoefs depending on input data
!===================================================================================================================================
! MODULES
USE MOD_Solov_Vars
USE MOD_PsiEval, ONLY:EvalPsi
USE MOD_PsiEval, ONLY:EvaldPsi
USE MOD_Newton, ONLY:NewtonMin1D
USE MOD_Newton, ONLY:NewtonMin2D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: tol,a(2),b(2)
REAL    :: qaxis_fact,mu0,dp_dpsin,qedge
REAL    :: Faxis,deltaF2,Fedge
!===================================================================================================================================

CALL SolveForPsiCoefs()
WRITE(*,'(A,4(E22.15,1X)))')'   PsiCoefs(1:4)        : ', PsiCoefs(1:4)
WRITE(*,'(A,3(E22.15,1X)))')'   PsiCoefs(5:7)        : ', PsiCoefs(5:7)

!                tolerance,
tol=1.0e-15
a(1)=1-0.1*p_eps !bound x- 
b(1)=1+0.4*p_eps !bound x+
a(2)=-0.1*p_kappa*p_eps !bound y-
b(2)=-a(2)              !bound y+
             

!1D minimum find, since we know here y=0
xaxis(1)=1.
xaxis(2)=0.
psi_axis=NewtonMin1D(tol,a(1),b(1),xaxis(1),F1,dF1,ddF1 )
WRITE(*,'(A,2E22.15)')'   magentic axis x,y    : ',xaxis(:)

!2D minimum search
xaxis(1)=1.
xaxis(2)=0.
psi_axis=NewtonMin2D(tol,a,b,xaxis,F2,dF2,ddF2 )

psi_edge=EvalPsi(1.+p_eps,0.)
WRITE(*,'(A,2E22.15)')'   magentic axis x,y    : ',xaxis(:)
WRITE(*,'(A,E22.15)') '   psi at magentic axis : ',psi_axis
WRITE(*,'(A,E22.15)') '   psi at domain edge   : ',psi_edge

!find out scaling of psi function
Faxis=p_B0*p_R0
F2_axis=Faxis**2

! q_axis=F_axis/(R_axis*(d_rr(psireal)*d_zz(psireal))^1/2), psireal=psi*psi0, evaluated at axis...
!       =F_axis/psi_scale* 1/(R_axis*(d_rr(psi)*d_zz(psi))^1/2), from book Freidberg, ideal MHD, p.134
! R=x*R0,Z=y*R0
! d/dR=(d/dx)*1/R0, d/dZ=(d/dy)*1/R0

qaxis_fact=p_R0/(xaxis(1)*SQRT(EvaldPsi(1,2,xaxis(1),xaxis(2))*EvaldPsi(2,2,xaxis(1),xaxis(2))) )

psi_scale=qaxis_fact*Faxis/p_qaxis

WRITE(*,'(A,E22.15)') '   psi_scale            : ',psi_scale
WRITE(*,'(A,E22.15)') '   real psi on axis     : ',psi_scale*psi_axis

deltaF2 = 2*p_A*psi_scale**2*psi_axis/(p_R0**2)
F2_edge=F2_axis+deltaF2
Fedge=SQRT(F2_edge)


!gradient of the pressure over normalized flux, p=p0+dpdpsin*psin 
!mu0=1
mu0=1.
dp_dpsin=(1-p_A)*psi_scale**2*psi_axis/(mu0*p_R0**4)
PresEdge=p_paxis+dp_dpsin
WRITE(*,'(A,E22.15)') '   pressure on axis     : ',p_paxis
WRITE(*,'(A,E22.15)') '   pressure on edge     : ',PresEdge

IF(PresEdge.LT.0.) STOP 'Pressure negative in domain, increase paxis input'

CALL Eval_qedge(Fedge,qedge)

WRITE(*,'(A,E22.15)') '   q-factor on axis     : ',p_qaxis
WRITE(*,'(A,E22.15)') '   q-factor on edge     : ',qedge

CONTAINS
  !for 1d NewtonMin
  FUNCTION F1(x)
    IMPLICIT NONE
    REAL:: x
    REAL:: F1
    F1=EvalPsi(x,0.)
  END FUNCTION F1

  FUNCTION dF1(x)
    IMPLICIT NONE
    REAL:: x
    REAL:: dF1
    dF1=EvaldPsi(1,1,x,0.)
  END FUNCTION dF1

  FUNCTION ddF1(x)
    IMPLICIT NONE
    REAL:: x
    REAL:: ddF1
    ddF1=EvaldPsi(1,2,x,0.)
  END FUNCTION ddF1

  !for 2d NewtonMin
  FUNCTION F2(x)
    IMPLICIT NONE
    REAL:: x(2)
    REAL:: F2
    F2=EvalPsi(x(1),x(2))
  END FUNCTION F2

  FUNCTION dF2(x)
    IMPLICIT NONE
    REAL:: x(2)
    REAL:: dF2(2)
    dF2(1)=EvaldPsi(1,1,x(1),x(2))
    dF2(2)=EvaldPsi(2,1,x(1),x(2))
  END FUNCTION dF2

  FUNCTION ddF2(x)
    IMPLICIT NONE
    REAL:: x(2)
    REAL:: ddF2(2,2)
 
    ddF2(1,1)=EvaldPsi(1,2,x(1),x(2))
    ddF2(1,2)=EvaldPsi(3,1,x(1),x(2))
    ddF2(2,1)=ddF2(1,2)
    ddF2(2,2)=EvaldPsi(2,2,x(1),x(2))
  END FUNCTION ddF2

END SUBROUTINE InitSolovievEquilibrium 


SUBROUTINE Eval_qedge(Fedge,qedge)
!===================================================================================================================================
! compute qfactor on edge, needs an circular integration over a flux surface, here given by boundary curve  
!===================================================================================================================================
! MODULES
USE MOD_globals, ONLY:Pi
USE MOD_Solov_Vars,ONLY:asin_delta,p_eps,p_kappa,xaxis,p_R0,psi_scale
USE MOD_PsiEval,ONLY:EvalPsi,EvaldPsi
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Fedge
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT) :: qedge
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,nmax
REAL    :: x,y,rcoord,BthetaR
REAL    :: tau,dtau
REAL    :: psimax,psimin,psi_dx,psi_dy,intgr
!===================================================================================================================================
nmax=1000
dtau=2*pi/REAL(nmax) 
tau=0.
qedge=0.
DO i=1,nmax
  x = 1. + p_eps * COS(tau + asin_delta* SIN(tau))
  y = p_eps * p_kappa * SIN(tau)
  rcoord=SQRT((x-xaxis(1))**2+(y-xaxis(2))**2)
  psimin=MIN(psimin,EvalPsi(x,y))
  psimax=MAX(psimax,EvalPsi(x,y))
  psi_dx=EvaldPsi(1,1,x,y)
  psi_dy=EvaldPsi(2,1,x,y)
  !BR=-1/R*dpsireal_dZ, Bz=1/R*dpsireal_dR 
  ! Btheta=-BR*sin(theta)+BZ*cos(theta)
  ! Btheta*R= dpsireal_dZ*sin(theta)+dpsireal_dR*cos(theta), dpsiReal_dZ=psi_scale/R0*dpsi_dy,  sin=(y-y_a)/r, cos=(x-x_a)/r
  BthetaR=psi_scale/(p_R0*rcoord)*(psi_dy*(y-xaxis(2))+psi_dx*(x-xaxis(1))) !/R=(x*R0)
  
  intgr = rcoord/(x*BthetaR)  ! oint r/(R^2*Btheta) dtheta, r=rcoord*R0, R=x*R0
  qedge = qedge+intgr
  tau=tau+dtau
END DO !i
qedge=Fedge*qedge/REAL(nmax) !=Fedge/2*pi*qedge*dtau

END SUBROUTINE Eval_qedge


FUNCTION ApproxFluxMap(rr,beta) RESULT(xyCoords)
!===================================================================================================================================
! Approximately flux-aligned coordinates, rr[0,1] ~sqrt(psi), theta [0,2pi]
!===================================================================================================================================
! MODULES
USE MOD_Solov_Vars,ONLY:asin_delta,p_eps,p_kappa,xaxis
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: rr     !linear radial coordinate [0,1] ~ sqrt(psi)
REAL, INTENT(IN) :: beta   ! poloidal angle (not atan(y,x)!) 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL      :: xyCoords(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL      :: sinbeta
!===================================================================================================================================
sinbeta=SIN(beta)
xyCoords(1) =1.+(xaxis(1)-1)*(1-rr**2)+p_eps * rr*COS(beta + asin_delta*rr * sinbeta ) !edge_X-1, mit p_delta=rr*p_delta
xyCoords(2) = rr*p_eps * p_kappa * sinbeta
         
END FUNCTION ApproxFluxMap


FUNCTION ApproxFluxMap_dr(rr,beta) RESULT(xyCoords_dr)
!===================================================================================================================================
! Approximately flux-aligned coordinates, rr[0,1] ~sqrt(psi), theta [0,2pi]
!===================================================================================================================================
! MODULES
USE MOD_Solov_Vars,ONLY:asin_delta,p_eps,p_kappa,xaxis
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: rr     !linear radial coordinate [0,1] ~ sqrt(psi)
REAL, INTENT(IN) :: beta   ! poloidal angle (not atan(y,x)!) 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL      :: xyCoords_dr(2)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL      :: sinbeta
!===================================================================================================================================
sinbeta=SIN(beta)

xyCoords_dr(1) =-2.*rr*(xaxis(1)-1)+p_eps * (COS(beta + asin_delta*rr * sinbeta ) &
                                             -rr*(SIN(beta + asin_delta*rr * sinbeta ))*asin_delta*sinbeta)
xyCoords_dr(2) = p_eps * p_kappa * sinbeta
         
END FUNCTION ApproxFluxMap_dr


FUNCTION PsiToPsiNorm(psi) RESULT(PsiNorm)
!===================================================================================================================================
! evaluate normalized flux [0,1], 0 on axis, 1 on edge
!===================================================================================================================================
! MODULES
USE MOD_Solov_Vars,ONLY:psi_axis,psi_edge
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: psi
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: PsiNorm
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
PsiNorm=(Psi-psi_axis)/(psi_edge-psi_axis)
END FUNCTION PsiToPsiNorm


FUNCTION PsiNormToPsi(psiNorm) RESULT(Psi)
!===================================================================================================================================
! evaluate normalized flux [0,1], 0 on axis, 1 on edge
!===================================================================================================================================
! MODULES
USE MOD_Solov_Vars,ONLY:psi_axis,psi_edge
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN) :: psiNorm
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL :: Psi
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
Psi=PsiNorm*(psi_edge-psi_axis)+psi_axis
END FUNCTION PsiNormToPsi


SUBROUTINE SolveForPsiCoefs()
!===================================================================================================================================
! Builds up the linear system and solves for psiCoefs(0:7) 
!===================================================================================================================================
! MODULES
USE MOD_Basis1D,ONLY:SOLVE
USE MOD_Solov_Vars
USE MOD_PsiEval, ONLY:EvalPsiVec,EvaldPsidxVec,Evald2PsidxVec
USE MOD_PsiEval, ONLY:EvaldPsidyVec,Evald2PsidyVec
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
REAL                :: SysMat(7,7),RHS(7)
REAL                :: xm,xc,xp,yp
REAL                :: N1,N2,N3
REAL,DIMENSION(0:7) :: Psi_xmy0,dPsidx_xmy0,d2Psidy_xmy0
REAL,DIMENSION(0:7) :: Psi_xpy0,dPsidx_xpy0,d2Psidy_xpy0
REAL,DIMENSION(0:7) :: Psi_xcyp,dPsidy_xcyp,d2Psidx_xcyp,dPsidx_xcyp
!===================================================================================================================================
PsiCoefs(0)=1.

xm=1.-p_eps
xc=1.-p_delta*p_eps
xp=1.+p_eps
yp=p_kappa*p_eps
N1=  - (1 + asin_delta)**2 / (p_eps * p_kappa**2)
N2 = + (1 - asin_delta)**2 / (p_eps * p_kappa**2)
N3 = - p_kappa / (p_eps * COS(asin_delta)**2)

    Psi_xmy0=EvalPsiVec(xm,0.)
 dPsidx_xmy0=EvaldPsidxVec(xm,0.)
d2Psidy_xmy0=Evald2PsidyVec(xm,0.)

    Psi_xpy0=EvalPsiVec(xp,0.)
 dPsidx_xpy0=EvaldPsidxVec(xp,0.)
d2Psidy_xpy0=Evald2PsidyVec(xp,0.)

    Psi_xcyp=EvalPsiVec(xc,yp)
 dPsidx_xcyp=EvaldPsidxVec(xc,yp)
 dPsidy_xcyp=EvaldPsidyVec(xc,yp)
d2Psidx_xcyp=Evald2PsidxVec(xc,yp)

RHS(1) = -Psi_xpy0(0)
RHS(2) = -Psi_xmy0(0)
RHS(3) = -Psi_xcyp(0)
RHS(4) = -dPsidx_xcyp(0)
RHS(5) = -(d2Psidy_xpy0(0)+ N1*dPsidx_xpy0(0))
RHS(6) = -(d2Psidy_xmy0(0)+ N2*dPsidx_xmy0(0))
RHS(7) = -(d2Psidx_xcyp(0)+ N3*dPsidy_xcyp(0))

DO i=1,7
  sysmat(1,i) = Psi_xpy0(i)
  sysmat(2,i) = Psi_xmy0(i)
  sysmat(3,i) = Psi_xcyp(i)
  sysmat(4,i) = dPsidx_xcyp(i)
  sysmat(5,i) = d2Psidy_xpy0(i)+ N1*dPsidx_xpy0(i)
  sysmat(6,i) = d2Psidy_xmy0(i)+ N2*dPsidx_xmy0(i)
  sysmat(7,i) = d2Psidx_xcyp(i)+ N3*dPsidy_xcyp(i)
END DO !i=1,7

PsiCoefs(1:7)=SOLVE(sysmat,RHS)

END SUBROUTINE SolveForPsiCoefs


SUBROUTINE MapToSolov(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
!===================================================================================================================================
! Maps a cylinder (r,z,phi) to a toroidal closed flux surface configuration derived from VMEC data. 
! Surfaces with constant r become flux surfaces. z [0;1] is mapped to [0;2*pi] 
! for a fixed theta, uses a newton method to find the exact location  in r of psi(x(r),y(r))-Psi0=0, psi0=psi(psinorm=r_p^2)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MHDEQ_Vars,  ONLY:nVarMHDEQ
USE MOD_MHDEQ_Tools, ONLY: Eval1DPoly
USE MOD_Newton,      ONLY:NewtonRoot1D
USE MOD_Solov_Vars,  ONLY:p_R0,p_kappa,p_paxis,PresEdge,nRhoCoefs,RhoCoefs,psi_scale,F2_axis,F2_edge
USE MOD_PsiEval,     ONLY:EvalPsi,EvaldPsi
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: nTotal         ! total number of points
REAL, INTENT(IN)   :: x_in(3,nTotal) ! input coordinates represent a cylinder: 
INTEGER, INTENT(IN):: InputCoordSys  ! =0: x_in(1:3) are (x,y,z) coordinates in a cylinder of size r=[0;1], z=[0;1]
                                     ! =1: x_in(1:3) are (r,z,phi) coordinates r= [0;1], z= [0;1], phi=[0;1]
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)   :: x_out(3,nTotal) ! mapped x,y,z coordinates with vmec data
REAL,INTENT(OUT)   :: MHDEQdata(nVarMHDEQ,nTotal) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: r_p    ! raduis in cylindrical coordinate system
REAL    :: theta  ! poloidal angle [0,2pi]
REAL    :: zeta ! toroidal angle [0,2pi]
INTEGER :: iNode
INTEGER :: percent
REAL    :: coszeta,sinzeta
REAL    :: psiVal,psiNorm 
REAL    :: tol,rNewton,xNewton(2)
REAL    :: Density
REAL    :: R
REAL    :: BR,BZ,Bphi,Bcart(3) 
REAL    :: AR,AZ,Aphi,Acart(3)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,I8,A,A,A)')'  MAP ', nTotal,' NODES TO SOLOVIEV EQUILIBRIUM'
percent=0
tol=1.0E-12
DO iNode=1,nTotal
  ! output of progress in %
  IF((nTotal.GT.10000).AND.(MOD(iNode,(nTotal/100)).EQ.0)) THEN
    percent=percent+1
    WRITE(0,'(I4,A23,A1)',ADVANCE='NO')percent, ' % of nodes evaluated...',ACHAR(13)
  END IF
  SELECT CASE(InputCoordSys)
  CASE(0)!x_in(1:3) = x,y,z of cylinder with r<1 and z=[0;1]
    r_p   = SQRT(x_in(1,iNode)**2+x_in(2,iNode)**2) 
    Theta = ATAN2(x_in(2,iNode),x_in(1,iNode))
    zeta  = 2.*Pi*x_in(3,iNode) 
  CASE(1) !x_in(1:3) = (r,z,phi) with r= [0;1], z= [0;1], phi=[0;1] 
    r_p =  x_in(1,iNode) !=r
    Theta = 2.*Pi*x_in(3,iNode) !=2*pi*phi
    zeta  = 2.*Pi*x_in(2,iNode) !=2*pi*z
  END SELECT 
  coszeta=COS(zeta)
  sinzeta=SIN(zeta)
  psiNorm=r_p**2 !r ~ sqrt(psiNorm), r^2 ~ psiNorm 

  theta = theta -0.1*(p_kappa-1.)*SIN(2*theta) !slight angle correction with ellipticity
  IF(r_p.LT.1.0E-08) THEN
   rNewton=r_p
  ELSE
    psiVal = PsiNormToPsi(psiNorm)  
    rNewton= NewtonRoot1D(tol,0.,1.1,r_p,psiVal,FR1,dFR1)
  END IF
  xNewton=ApproxFluxMap(rNewton,theta)
  
  R=p_R0*xNewton(1)
  !Z=p_R0*xNewton(2)
  x_out(1,iNode)= R*COS(zeta)
  x_out(2,iNode)=-R*SIN(zeta)
  x_out(3,iNode)= (p_R0*xNewton(2))

  Density=Eval1DPoly(nRhoCoefs,RhoCoefs,psiNorm) 
  !BR=-1/R*dpsiReal_dZ = -1/(x*R0)*psi_scale*dpsi_dy*dy/dZ = -psiscale/(x*R0^2)*dpsi_dy
  !BZ= 1/R*dpsiReal_dR =  1/(x*R0)*psi_scale*dpsi_dx*dx/dR =  psiscale/(x*R0^2)*dpsi_dx
  !Bphi=F/R
  BR=-psi_scale/(p_R0*R)*EvaldPsi(2,1,xNewton(1),xNewton(2))
  BZ= psi_scale/(p_R0*R)*EvaldPsi(1,1,xNewton(1),xNewton(2))
  Bphi=SQRT(F2_axis+(F2_edge-F2_axis)*psiNorm)/R


  Bcart(1)= BR*coszeta-Bphi*sinzeta
  Bcart(2)= BR*sinzeta+Bphi*coszeta
  Bcart(3)= BZ
  
  !aphi= psireal/R = psi_scale*psi/(x*R0)

  Aphi=psi_scale/R*EvalPsi(xNewton(1),xNewton(2))
  AR  = 0.
  AZ  = 0.

  Acart(1)= AR*coszeta-Aphi*sinzeta
  Acart(2)= AR*sinzeta+Aphi*coszeta
  Acart(3)= AZ
  MHDEQdata(:,iNode)=0.
  MHDEQdata(  1,iNode)=Density
  MHDEQdata(  2,iNode)=p_paxis + (PresEdge-p_paxis)*psiNorm !linear pressure profile
  MHDEQdata( 3:5,iNode)=Bcart(:)
  MHDEQdata(   6,iNode)=psiNorm
  MHDEQdata(   7,iNode)=-1.
  MHDEQdata(8:10,iNode)=Acart(:)
END DO !iNode


WRITE(UNIT_stdOut,'(A)')'  ...DONE.                             '

!for newton search
CONTAINS
  FUNCTION FR1(r)
    USE MOD_PsiEval,ONLY:EvalPsi
    IMPLICIT NONE
    REAL     :: r
    REAL     :: FR1
    !local
    REAL     :: xloc(2)
    xloc=ApproxFluxMap(r,Theta) !theta from subroutine!
    FR1=EvalPsi(xloc(1),xloc(2))
  END FUNCTION FR1

  FUNCTION dFR1(r)
    USE MOD_PsiEval,ONLY:EvaldPsi
    IMPLICIT NONE
    REAL     :: r
    REAL     :: dFR1
    !local
    REAL     :: xloc(2),dxloc_dr(2)
    xloc=ApproxFluxMap(r,Theta) !theta from subroutine!
    dxloc_dr=ApproxFluxMap_dr(r,Theta) !theta from subroutine!
    !d/dr(psi(x(r),y(r))) = dpsi_dx(r)*dx_dr(r) +dpsi_dy*dy_dr(r)
    dFR1= EvaldPsi(1,1,xloc(1),xloc(2))*dxloc_dr(1) &
          +EvaldPsi(2,1,xloc(1),xloc(2))*dxloc_dr(2)
  END FUNCTION dFR1

END SUBROUTINE MapToSolov 


END MODULE MOD_Solov
