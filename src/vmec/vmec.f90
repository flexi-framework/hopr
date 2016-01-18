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

  NAMELIST / input / outDir, intPointsU, intPointsV, mapEuterpe, useArcTan,&
       useRho, gyroperp, rGridPoints, zGridPoints, nPhiCuts, startCut, endCut,&
       constFac, nFlux, nTheta, lastExtS, nExtraFlux, nInTheta, inSVal,&
       nPointsU, useConstTheta, nThetaInt, nPhiInt, nu, nv, dataFile, lowBR,&
       highBR, lowBZ, highBZ, corVMEC, useSquared, useScaledB, useLastVal,&
       useVar, weight, errInt, errVal, useSecant, newtonMax, maxRec, relax,&
       debug
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'INIT VMEC INPUT ...'
useVMEC      = GETLOGICAL('useVMEC','.FALSE.')   ! Use / reconstruct spline boundaries
IF(useVMEC)THEN

  me_rank = 0
  me_size = 1
  OPEN(1, ACTION = "READ", FILE = "vmec_input", FORM = "FORMATTED",&
       IOSTAT = ioError, POSITION = "REWIND", STATUS = "OLD")
  IF (ioError == 0) THEN
    READ(1, NML = input, IOSTAT = ioError)
    IF (me_rank == 0) THEN
      IF (ioError == 0) THEN
        PRINT "(A)", '  Using "input" for input parameters!'
      ELSE
        PRINT "(A)", '  Error while reading input parameters!'
        CALL EXIT(21)
      END IF
    END IF
    CLOSE(1)
  ELSE
    PRINT *, '  Missing input file!'
    CALL EXIT(22)
  END IF

  !! set default values
  IF (startCut < 0) startCut = 1
  IF ((endCut < 0) .OR. (endCut > nPhiCuts)) endCut = nPhiCuts

  !! read VMEC 2000 output (netcdf)
  CALL ReadVmecOutput(dataFile)

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
USE MOD_VMEC_Mappings, ONLY: Rmnc, Zmns,phi,xm,xn,nFluxVMEC,mn_mode
USE MOD_VMEC_Mappings, ONLY: CosTransFullMesh,SinTransFullMesh 
USE MOD_VMEC_Vars, ONLY: mn_mEven,mn_mOdd,mn_mapOdd,mn_mapEven
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
REAL    :: phi_p  ! flux coordinate [0,1] (use radial distance of point position)
REAL    :: theta  ! poloidal angle [0,2pi]
REAL    :: zeta ! toroidal angle [0,2pi]
REAL    :: R,Z   
REAL    :: r1e,r2e,z1e,z2e
REAL    :: r1o,r2o,z1o,z2o,frac,w1,w2 
REAL    :: phimin,phimax   
REAL    :: phinorm(nFluxVMEC)   
INTEGER :: s1,s2
!===================================================================================================================================
phi_p  = SQRT(xcyl(1)**2+xcyl(2)**2)
theta   = ATAN2(xcyl(2),xcyl(1))
zeta   = 2*Pi*xcyl(3) 


phinorm=(phi-phi(1))/phi(nFluxVMEC-1)

phi_p=phi_p**2

!! look for the nearest supporting point
s1=1
DO s1 = 1, nFluxVMEC
  IF (phi_p < phinorm(s1)) EXIT
END DO
s2=MIN(s1+1,nFluxVMEC)

IF(s1.NE.s2)THEN
  !evaluate only modes with m= even
  r1e = CosTransFullMesh(mn_mEven, Rmnc(mn_mapEven, s1), xm(mn_mapEven), xn(mn_mapEven), theta, zeta)
  r2e = CosTransFullMesh(mn_mEven, Rmnc(mn_mapEven, s2), xm(mn_mapEven), xn(mn_mapEven), theta, zeta)
  z1e = SinTransFullMesh(mn_mEven, Zmns(mn_mapEven, s1), xm(mn_mapEven), xn(mn_mapEven), theta, zeta)
  z2e = SinTransFullMesh(mn_mEven, Zmns(mn_mapEven, s2), xm(mn_mapEven), xn(mn_mapEven), theta, zeta)
  !evaluate only modes with m= odd
  r1o = CosTransFullMesh(mn_mOdd , Rmnc(mn_mapOdd , s1), xm(mn_mapOdd ), xn(mn_mapOdd ), theta, zeta)
  r2o = CosTransFullMesh(mn_mOdd , Rmnc(mn_mapOdd , s2), xm(mn_mapOdd ), xn(mn_mapOdd ), theta, zeta)
  z1o = SinTransFullMesh(mn_mOdd , Zmns(mn_mapOdd , s1), xm(mn_mapOdd ), xn(mn_mapOdd ), theta, zeta)
  z2o = SinTransFullMesh(mn_mOdd , Zmns(mn_mapOdd , s2), xm(mn_mapOdd ), xn(mn_mapOdd ), theta, zeta)
  !interpolation factor
  frac=(phi_p-phinorm(s1))/(phinorm(s2)-phinorm(s1))
  ! interpolation of odd modes with weighting sqrt(phi) * ((1-frac) * r1/sqrt(phi_1) + frac* r2/sqrt(phi_2))
  w1=SQRT(phi_p/phinorm(s1))
  w2=SQRT(phi_p/phinorm(s2))
  R=(1.-frac)*(r1e+w1*r1o) + frac*(r2e+w2*r2o)
  Z=(1.-frac)*(z1e+w1*z1o) + frac*(z2e+w2*z2o)
ELSE
  !weighting with sqrt(s) cancels, evaluate all modes at s2.
  R = CosTransFullMesh(mn_mode, Rmnc(:, s2), xm(:), xn(:), theta, zeta)
  Z = SinTransFullMesh(mn_mode, Zmns(:, s2), xm(:), xn(:), theta, zeta)
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
xvmec(3)=Z

END FUNCTION MapToVmec 


END MODULE MOD_VMEC
