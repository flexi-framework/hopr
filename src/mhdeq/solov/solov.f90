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
    P_qaxis =GETREAL("Solov_qaxis")
    P_paxis =GETREAL("Solov_paxis")
  ELSE
    !preselected setups (defaults are set, can still be modified in inifile)
    SELECT CASE(setup)
    CASE(1) !circular
      P_R0    =GETREAL("Solov_R0","10.")
      P_eps   =GETREAL("Solov_eps","1.")
      P_kappa =GETREAL("Solov_kappa","1.")
      P_delta =GETREAL("Solov_delta","0.")
      P_A     =GETREAL("Solov_A","0.955")
      P_qaxis =GETREAL("Solov_qaxis","1.89")
      P_paxis =GETREAL("Solov_paxis","2.0e-03")
    CASE(2) !parameter ITER-like
      P_R0    =GETREAL("Solov_R0","6.2")
      P_eps   =GETREAL("Solov_eps","0.32")
      P_kappa =GETREAL("Solov_kappa","1.7")
      P_delta =GETREAL("Solov_delta","0.33")
      P_A     =GETREAL("Solov_A","-0.155")
      P_qaxis =GETREAL("Solov_qaxis","1.6")
      P_paxis =GETREAL("Solov_paxis","8.0e-02")
    CASE DEFAULT
      WRITE(*,*) "this soloviev setup does not exist" ,setup
      STOP
    END SELECT
  END IF !setup=0
  IF(p_A.GT.1.) WRITE(*,*) 'WARNING, USING POSITIVE PRESSURE GRADIENT p_axis< p_edge'
    
  CALL InitSolovievEquilibrium()

  nRhoCoefs=GETINT("Solov_nRhoCoefs","0")
  IF(nRhoCoefs.GT.0)THEN
    ALLOCATE(RhoCoefs(nRhoCoefs))
    RhoCoefs=GETREALARRAY("Solov_RhoCoefs",nRhoCoefs)
  END IF
  WRITE(UNIT_stdOut,'(A)')'  ... DONE'
END SUBROUTINE InitSolov


SUBROUTINE InitSolovievEquilibrium 
!===================================================================================================================================
! Solve for psiCoefs depending on input data
!===================================================================================================================================
! MODULES
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
END SUBROUTINE InitSolovievEquilibrium 


SUBROUTINE SolveForPsiCoefs()
!===================================================================================================================================
! Builds up the linear system and solves for psiCoefs 
!===================================================================================================================================
! MODULES
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
REAL                :: SysMat(7,7),InvSysMat(7,7),RHS(7)
REAL                :: xm,xc,xp,yp
REAL                :: alpha,N1,N2,N3
REAL,DIMENSION(0:7) :: Psi_xmy0,Psi_xcyp,Psi_xpy0
REAL,DIMENSION(0:7) :: dPsidx_xmy0,dPsidx_xcyp,dPsidx_xpy0
REAL,DIMENSION(0:7) :: d2Psidx_xcyp,d2Psidy_xpy0,d2Psidy_xmy0
!===================================================================================================================================
PsiCoefs(0)=1.

xm=1.-p_eps
xc=1.-p_delta*p_eps
xp=1.+p_eps
yp=p_kappa*p_eps
alpha = ASIN(p_delta)
N1=  - (1 + alpha)**2 / (p_eps * p_kappa**2)
N2 = + (1 - alpha)**2 / (p_eps * p_kappa**2)
N3 = - p_kappa / (p_eps * COS(alpha)**2)

Psi_xmy0=EvalPsiVec(xm,0.)
Psi_xpy0=EvalPsiVec(xp,0.)
Psi_xcyp=EvalPsiVec(xc,yp)


END SUBROUTINE SolveForPsiCoefs


SUBROUTINE MapToSolov(nTotal,x_in,InputCoordSys,x_out,MHDEQdata)
!===================================================================================================================================
! Maps a cylinder (r,z,phi) to a toroidal closed flux surface configuration derived from VMEC data. 
! Surfaces with constant r become flux surfaces. z [0;1] is mapped to [0;2*pi] 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_MHDEQ_Vars, ONLY:nVarMHDEQ
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
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,I8,A,A,A)')'  MAP ', nTotal,' NODES TO SOLOVIEV EQUILIBRIUM'
percent=0
DO iNode=1,nTotal
  ! output of progress in %
  IF((nTotal.GT.10000).AND.(MOD(iNode,(nTotal/100)).EQ.0)) THEN
    percent=percent+1
    WRITE(0,'(I4,A23,A1)',ADVANCE='NO')percent, ' % of nodes evaluated...',ACHAR(13)
  END IF
  SELECT CASE(InputCoordSys)
  CASE(0)!x_in(1:3) = x,y,z of cylinder with r<1 and z=[0;1]
    r_p   = SQRT(x_in(1,iNode)**2+x_in(2,iNode)**2) 
    theta = ATAN2(x_in(2,iNode),x_in(1,iNode))
    zeta  = -2.*Pi*x_in(3,iNode) 
  CASE(1) !x_in(1:3) = (r,z,phi) with r= [0;1], z= [0;1], phi=[0;1] 
    r_p =  x_in(1,iNode) !=r
    theta = 2.*Pi*x_in(3,iNode) !=2*pi*phi
    zeta  = -2.*Pi*x_in(2,iNode) !=2*pi*z
  END SELECT 
  !dummy
  x_out(:,iNode)=x_in(:,iNode)
  MHDEQdata(:,iNode)=0.
END DO !iNode
WRITE(UNIT_stdOut,'(A)')'  ...DONE.                             '
END SUBROUTINE MapToSolov 


END MODULE MOD_Solov
