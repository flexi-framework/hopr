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
#include "defines.f90"
MODULE MOD_VMEC
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

INTERFACE InitVMEC 
  MODULE PROCEDURE InitVMEC 
END INTERFACE

INTERFACE MapToVMEC 
  MODULE PROCEDURE MapToVMEC 
END INTERFACE

PUBLIC::InitVMEC
PUBLIC::MapToVMEC
!===================================================================================================================================

CONTAINS
SUBROUTINE InitVMEC 
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut
USE MOD_ReadInTools
USE MOD_VMEC_Mappings
USE MOD_VMEC_Vars
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
CHARACTER(LEN = 256) :: dataFile
INTEGER              :: ioError
INTEGER              :: iMode,iEven,iOdd
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'INIT VMEC INPUT ...'
useVMEC      = GETLOGICAL('useVMEC','.FALSE.')   ! Use / reconstruct spline boundaries
IF(useVMEC)THEN

  me_rank = 0
  me_size = 1
  !VMEC "wout*.nc"  file
  dataFile=GETSTR("VMECwoutfile")
  ! use internal remapping of VMEC output 
  corVMEC =GETLOGICAL("corVMEC",".FALSE.") 
  ! grid to get fourier modes (default: 128, 128)
  intPointsU=GETINT("VMEC_intPointsU","128")
  intPointsV=GETINT("VMEC_intPointsV","128")
  ! use scaled |B|*sqrtG instead of |B| (default: .FALSE.)
  useScaledB=GETLOGICAL("VMEC_useScaledB",".FALSE.")
  ! use last value for extrapolation or linear extrapolation (default: .FALSE.)
  useLastVal=GETLOGICAL("VMEC_useLastVal",".FALSE.,")
  ! use error bars for R, z smoothing (default: .TRUE.)
  useVar =GETLOGICAL("VMEC_useVar", ".TRUE.")
  ! smoothing weight for R, z smoothing (default: 1.D-3)
  weight =GETREAL("VMEC_weight","1.0-3")
  ! error intervals (default: 1.D-1, 9.5D-1)
  errInt= GETREALARRAY("VMEC_errInt",2,"1.0-1, 9.5-1")
  ! error bars for intervals (default: 1.D-1, 1.D-3, 1.D-2)
  errVal=GETREALARRAY("VMEC_errVal",3,"1.0-1,1.0-3,1.0-2")
  ! debug output (files/messages) (default: .FALSE.)
  debug = .FALSE.

 
  !! read VMEC 2000 output (netcdf)
  CALL ReadVmecOutput(dataFile)

  ALLOCATE(phinorm(nFluxVMEC))
  phinorm=(phi-phi(1))/phi(nFluxVMEC-1)

  !find even and odd m-modes, to seperately evalute them
  mn_mEven=0
  DO iMode=1,mn_mode
    IF(MOD(xm(iMode),2.).EQ.0.) mn_mEven=mn_mEven+1
  END DO ! i=1,mn_mode

  mn_mOdd=mn_mode-mn_mEven
  ALLOCATE(mn_mapOdd(mn_mOdd),mn_mapEven(mn_mEven))
  iEven=0
  iOdd=0
  DO iMode=1,mn_mode
    IF(MOD(xm(iMode),2.).EQ.0.) THEN
      iEven=iEven+1
      mn_mapEven(iEven)=iMode
    ELSE
      iOdd=iOdd+1
      mn_mapOdd(iOdd)=iMode
    END IF !even
  END DO ! i=1,mn_mode
  WRITE(UNIT_stdOut,*)'   Total Number of mn-modes:',mn_mode
  WRITE(UNIT_stdOut,*)'   Number of even(m) and odd(m) mn-modes:',mn_mEven,mn_mOdd
END IF !useVMEC

WRITE(UNIT_stdOut,'(A)')'... DONE'
END SUBROUTINE InitVMEC


FUNCTION MapToVMEC(xcyl) RESULT(xvmec)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_VMEC_Mappings, ONLY: mn_mode,xm,xn,Rmnc, Zmns,lUmnc,lVmnc
USE MOD_VMEC_Mappings, ONLY: nFluxVMEC,phi,phipf,iotas
USE MOD_VMEC_Vars,     ONLY: phinorm,mn_mEven,mn_mOdd,mn_mapOdd,mn_mapEven
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)  :: xcyl(3) ! x,y,z coordinates in a cylinder of size r=0,1, z=0,1
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL              :: xvmec(3) ! mapped x,y,z coordinates with vmec data
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: s1,s2
REAL    :: CosMN(mn_mode),SinMN(mn_mode)
REAL    :: phi_p  ! flux coordinate [0,1] (use radial distance of point position)
REAL    :: theta  ! poloidal angle [0,2pi]
REAL    :: zeta ! toroidal angle [0,2pi]
REAL    :: R,Z   
REAL    :: f1,f2,w1,w2 
REAL    :: phimin,phimax   
REAL    :: phipf_int,iotas_int
REAL    :: Btheta,Bzeta
!===================================================================================================================================
phi_p  = SQRT(xcyl(1)**2+xcyl(2)**2) 
theta   = ATAN2(xcyl(2),xcyl(1))
zeta   = 2*Pi*xcyl(3) 

CosMN  = COS(xm(:) * theta - xn(:) * zeta)
SinMN  = SIN(xm(:) * theta - xn(:) * zeta)

phi_p=phi_p**2 ! use scaling of radius to phi evaluation variable

!! look for the nearest supporting point
s1=1
DO s1 = 1, nFluxVMEC
  IF (phi_p < phinorm(s1)) EXIT
END DO
s2=MIN(s1+1,nFluxVMEC)

IF(s1.NE.s2)THEN
  !interpolation factor
  f2=(phi_p-phinorm(s1))/(phinorm(s2)-phinorm(s1))
  f1=1.-f2 
  ! interpolation of odd modes with weighting sqrt(phi) * ((1-frac) * r1/sqrt(phi_1) + frac* r2/sqrt(phi_2))
  w1=SQRT(phi_p/phinorm(s1))
  w2=SQRT(phi_p/phinorm(s2))

  R=InterpolateData(f1,f2,w1,w2,CosMN,Rmnc(:,s1),Rmnc(:,s2))

  Z=InterpolateData(f1,f2,w1,w2,SinMN,Zmns(:,s1),Zmns(:,s2))

  !compute magnetic field: 
  ! B^s     = 0 
  ! B^theta = dphi/ds*(iota-dlambda/dzeta)  : phipf*(iotas -  (lVmnc) ) (iotas & lmns are overwritten to full mesh)
  ! B^zeta = dphi/ds*(1+dlambda/dtheta)     : phipf*(1+ (lUmnc) )
  phipf_int = f1*phipf(s1)+f2*phipf(s2)
  iotas_int = f1*iotas(s1)+f2*iotas(s2) 

  Btheta = phipf_int*(iotas_int - InterpolateData(f1,f2,w1,w2,CosMN,lVmnc(:,s1),lVmnc(:,s2)))

  Bzeta  = phipf_int*(1.        + InterpolateData(f1,f2,w1,w2,CosMN,lUmnc(:,s1),lUmnc(:,s2)))
ELSE
  !weighting with sqrt(s) cancels, evaluate all modes at s2.
  R = SUM(Rmnc(:, s2)*CosMN(:))
  Z = SUM(Zmns(:, s2)*SinMN(:))
  Btheta = phipf(s2)*(iotas(s2) - SUM(lVmnc(:, s2)*CosMN(:)))
  Bzeta  = phipf(s2)*(1.        + SUM(lUmnc(:, s2)*CosMN(:)))
END IF !s1/=s2

!test circular torus
!R=xcyl(1)+10.
!Z=xcyl(2)

!test deformed torus
!R=(0.5+0.2*SIN(zeta))*xcyl(1)
!Z=-SIN(3*zeta)*R+COS(3*zeta)*xcyl(2)
!R=10+COS(3*zeta)*R+SIN(3*zeta)*xcyl(2)

xvmec(1)= R*COS(zeta)
xvmec(2)=-R*SIN(zeta)
xvmec(3)= Z

END FUNCTION MapToVmec 


FUNCTION InterpolateData(f1,f2,w1,w2,trig_mn,Xmn_s1,Xmn_s2)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_VMEC_Mappings, ONLY: xm,xn,mn_mode
USE MOD_VMEC_Vars, ONLY: mn_mEven,mn_mOdd,mn_mapOdd,mn_mapEven
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)  :: f1,f2,w1,w2
REAL, INTENT(IN)  :: trig_mn(mn_mode) ! sin or cos evaluated at theta,zeta
REAL, INTENT(IN)  :: Xmn_s1( mn_mode) ! variable in fourier space to interpolate
REAL, INTENT(IN)  :: Xmn_s2( mn_mode) ! variable in fourier space to interpolate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL              :: InterpolateData
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: val1e,val2e,val1o,val2o
!===================================================================================================================================
  !evaluate only modes with m= even
  val1e = SUM(trig_mn(mn_mapEven)*xmn_s1(mn_mapEven))
  val2e = SUM(trig_mn(mn_mapEven)*xmn_s2(mn_mapEven))
  !evaluate only modes with m= odd
  val1o = SUM(trig_mn(mn_mapOdd)*xmn_s1(mn_mapOdd))
  val2o = SUM(trig_mn(mn_mapOdd)*xmn_s2(mn_mapOdd))

  InterpolateData=f1*(val1e+w1*val1o) + f2*(val2e+w2*val2o)
END FUNCTION InterpolateData

END MODULE MOD_VMEC
