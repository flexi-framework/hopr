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

! allow different dimensions of input/output arrays
!INTERFACE MapToVMEC 
!  MODULE PROCEDURE MapToVMEC 
!END INTERFACE

PUBLIC::InitVMEC
PUBLIC::MapToVMEC
!===================================================================================================================================

CONTAINS
SUBROUTINE InitVMEC 
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:UNIT_stdOut,abort
USE MOD_ReadInTools
USE MOD_VMEC_Vars
USE SPLINE1_MOD,       ONLY:SPLINE1_FIT 
USE MOD_VMEC_Mappings !ReadVMECoutput,VMEC variables
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
INTEGER              :: iMode
LOGICAL              :: useFilter
REAL,ALLOCATABLE     :: lmns_half(:,:)
REAL,ALLOCATABLE     :: gmnc_half_nyq(:,:)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'INIT VMEC INPUT ...'
useVMEC      = GETLOGICAL('useVMEC','.FALSE.')   ! Use / reconstruct spline boundaries
IF(useVMEC)THEN

  me_rank = 0
  me_size = 1
  !VMEC "wout*.nc"  file
  VMECdataFile=GETSTR("VMECwoutfile")
  ! use internal remapping of VMEC output 
  corVMEC =GETLOGICAL("corVMEC",".FALSE.") 
  ! grid for evaluation of fourier modes (default: 4,4) , should be set to (128, 128) if smoothing is needed
  intPointsU=GETINT("VMEC_intPointsU","4") 
  intPointsV=GETINT("VMEC_intPointsV","4")
  ! use scaled |B|*sqrtG instead of |B| (default: .FALSE.)
  useScaledB=GETLOGICAL("VMEC_useScaledB",".FALSE.")
  ! use last value for extrapolation or linear extrapolation (default: .FALSE.)
  useLastVal=GETLOGICAL("VMEC_useLastVal",".FALSE.")
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

  !density coefficients of the polynomial coefficients: rho_1+rho_2*x + rho_3*x^2 ...
  nRhoCoefs=GETINT("nRhoCoefs","0")
  IF(nRhoCoefs.GT.0)THEN
    RhoFluxVar=GETINT("RhoFluxVar") ! dependant variable: =0: phinorm (tor. flux), =1:chinorm (pol. flux)
    ALLOCATE(RhoCoefs(nRhoCoefs))
    RhoCoefs=GETREALARRAY("RhoCoefs",nRhoCoefs)
  END IF

 
  !! read VMEC 2000 output (netcdf)
  CALL ReadVmecOutput(VMECdataFile)

  !data on half mesh (copy before its overwritten in precalcdata)
  ALLOCATE(lmns_half(1:nFluxVMEC,mn_mode)) 
  ALLOCATE(gmnc_half_nyq(1:nFluxVMEC,mn_mode_nyq))

  lmns_half         = lmns
  gmnc_half_nyq     = gmnc
  ! half data is stored from 2:nFluxVMEC

  !normalized toroidal flux (=flux variable s [0;1] in VMEC)
  ALLOCATE(phinorm(nFluxVMEC))
  phinorm=(phi-phi(1))/phi(nFluxVMEC)
  WRITE(UNIT_stdOut,'(A,3F10.4)')'   normalized flux of first three flux surfaces',phinorm(2:4)
  !normalized poloidal flux (=flux variable in Grad-Shafranov equation)
  ALLOCATE(chinorm(nFluxVMEC))
  chinorm=(chi-chi(1))/(chi(nFluxVMEC)-chi(1))

  WRITE(UNIT_stdOut,*)'   Total Number of mn-modes:',mn_mode
  WRITE(UNIT_stdOut,*)'   Max Mode m,n: ',MAXVAL(xm),MAXVAL(xn)
  WRITE(UNIT_stdOut,*)'   Total Number of mn-modes (Nyquist):',mn_mode_nyq
  WRITE(UNIT_stdOut,*)'   Max Mode m,n: ',MAXVAL(xm_nyq),MAXVAL(xn_nyq)

  useFilter=.TRUE. !GETLOGICAL('VMECuseFilter','.TRUE.') !SHOULD BE ALWAYS TRUE...

  ALLOCATE(xmabs(mn_mode))
  DO iMode=1,mn_mode
    xmabs(iMode)=ABS(NINT(xm(iMode)))
    IF(useFilter)THEN
      IF(xmabs(iMode) > 3) THEN !Filtering for |m| > 3
        IF(MOD(xmabs(iMode),2) == 0) THEN
          xmabs(iMode)=2 !Even mode, remove rho**2
        ELSE
          xmabs(iMode)=3 !Odd mode, remove rho**3
        END IF
      END IF
    END IF !usefilter
  END DO !iMode=1,mn_mode

  ALLOCATE(xmabs_nyq(mn_mode_nyq))
  DO iMode=1,mn_mode_nyq
    xmabs_nyq(iMode)=ABS(NINT(xm_nyq(iMode)))
    IF(useFilter)THEN
      IF(xmabs_nyq(iMode) > 3) THEN !Filtering for |m| > 3
        IF(MOD(xmabs_nyq(iMode),2) == 0) THEN
          xmabs_nyq(iMode)=2 !Even mode, remove rho**2
        ELSE
          xmabs_nyq(iMode)=3 !Odd mode, remove rho**3
        END IF
      END IF
    END IF !usefilter
  END DO !iMode=1,mn_mode

  !prepare Spline interpolation
  ALLOCATE(rho(1:nFluxVMEC))
  rho(:)=SQRT(phinorm(:))
  

  ALLOCATE(Rmnc_Spl(4,1:nFluxVMEC,mn_mode)) !first dim is for spline interpolation
  ALLOCATE(Zmns_Spl(4,1:nFluxVMEC,mn_mode))
  ALLOCATE(lmns_Spl(4,1:nFluxVMEC,mn_mode))
  ALLOCATE(gmnc_nyq_Spl(4,1:nFluxVMEC,mn_mode_nyq))

  CALL FitSpline(mn_mode,xmAbs,Rmnc,Rmnc_Spl)
  CALL FitSpline(mn_mode,xmAbs,Zmns,Zmns_Spl)
  CALL FitSplineHalf(mn_mode,xmAbs,lmns_half,lmns_Spl)
  CALL FitSplineHalf(mn_mode_nyq,xmAbs_nyq,gmnc_half_nyq,gmnc_nyq_Spl)
  ALLOCATE(pres_spl(4,1:nFluxVMEC))
  pres_spl(1,:)=presf(:)
  CALL SPLINE1_FIT(nFluxVMEC,rho,pres_Spl(:,:), K_BC1=3, K_BCN=0)
  ALLOCATE(phipf_spl(4,1:nFluxVMEC))
  phipf_spl(1,:)=phipf(:)
  CALL SPLINE1_FIT(nFluxVMEC,rho,phipf_Spl(:,:), K_BC1=3, K_BCN=0)
  ALLOCATE(iota_spl(4,1:nFluxVMEC))
  iota_spl(1,:)=iotaf(:)
  CALL SPLINE1_FIT(nFluxVMEC,rho,iota_Spl(:,:), K_BC1=3, K_BCN=0)
  ALLOCATE(chinorm_spl(4,1:nFluxVMEC))
  chinorm_spl(1,:)=chinorm(:)
  CALL SPLINE1_FIT(nFluxVMEC,rho,chinorm_Spl(:,:), K_BC1=3, K_BCN=0)


  !OUTPUT
  nVarVMEC=14
  ALLOCATE(VMECVarNames(nVarVMEC))
  VMECvarnames( 1)='VMEC-torfluxnorm'
  VMECvarnames( 2)='VMEC-polfluxnorm'
  VMECvarnames( 3)='VMEC-A_X'
  VMECvarnames( 4)='VMEC-A_Y'
  VMECvarnames( 5)='VMEC-A_Z'
  VMECvarnames( 6)='VMEC-B_X'
  VMECvarnames( 7)='VMEC-B_Y'
  VMECvarnames( 8)='VMEC-B_Z'
  VMECvarnames( 9)='VMEC-Pressure'
  VMECvarnames(10)='VMEC-Density'
  VMECvarnames(11)='VMEC-11'
  VMECvarnames(12)='VMEC-12'
  VMECvarnames(13)='VMEC-13'
  VMECvarnames(14)='VMEC-14'
  nVarOutVMEC=10
  ALLOCATE(VMECoutVarMap(1:nVarOutVMEC))
  !output variables for HDF5 file: Density,Pressure,Bx,By,Bz,normalized poloidal & toroidal flux, magnetic potential
  VMECoutVarMap(:)=(/10,9,6,7,8,2,1,3,4,5/) 
END IF !useVMEC

WRITE(UNIT_stdOut,'(A)')'... DONE'
END SUBROUTINE InitVMEC


SUBROUTINE FitSpline(modes,mAbs,Xmn,Xmn_Spl)
!===================================================================================================================================
! Fit disrete data along flux surfaces as spline for each fourier mode
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_VMEC_Mappings, ONLY: nFluxVMEC
USE MOD_VMEC_Vars,     ONLY: rho 
USE SPLINE1_MOD,       ONLY:SPLINE1_FIT 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: modes
INTEGER, INTENT(IN) :: mabs(modes)
REAL, INTENT(IN)    :: Xmn(modes,nFluxVMEC)  ! fourier coefficients at all flux surfaces 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)  :: Xmn_Spl(4,nFluxVMEC,modes)  ! spline fitted fourier coefficients 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iMode,iFlux
!===================================================================================================================================
Xmn_Spl=0.
DO iMode=1,modes
  !scaling with rho^|m|
  DO iFlux=2,nFluxVMEC
    IF(mabs(iMode).EQ.0)THEN
      Xmn_Spl(1,iFlux,iMode)=Xmn(iMode,iFlux)
    ELSE
      Xmn_Spl(1,iFlux,iMode)=Xmn(iMode,iFlux) /(rho(iFlux)**mabs(iMode))
    END IF
  END DO !i
  !Parabolic extrapolation to axis with dx'(rho=0)=0.
  Xmn_Spl(1,1,iMode)=(Xmn_Spl(1,2,iMode)*rho(3)**2-Xmn_Spl(1,3,iMode)*rho(2)**2) /(rho(3)**2-rho(2)**2)
!  !Quadratic extrapolation to axis with dx'(rho=0)=0.
!  r1=rho(2)**2*rho(3)**4-rho(2)**4*rho(3)**2
!  r2=rho(2)**2*rho(4)**4-rho(2)**4*rho(4)**2
!  Xmn_Spl(1,1,iMode)= ( r1*(Xmn_Spl(1,2,iMode)*rho(4)**4-Xmn_Spl(1,4,iMode)*rho(2)**4) &
!                       -r2*(Xmn_Spl(1,2,iMode)*rho(3)**4-Xmn_Spl(1,3,iMode)*rho(2)**4)) &
!                     /( r1*(rho(4)**4-rho(2)**4)-r2*(rho(3)**4-rho(2)**4))
  CALL SPLINE1_FIT(nFluxVMEC,rho,Xmn_Spl(:,:,iMode), K_BC1=3, K_BCN=0)
END DO !iMode 

END SUBROUTINE FitSpline


SUBROUTINE FitSplineHalf(modes,mabs,Xmn_half,Xmn_Spl)
!===================================================================================================================================
! Fit disrete data along flux surfaces as spline for each fourier mode
! input is given on the half mesh 2:nFluxVMEC 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_VMEC_Mappings, ONLY: nFluxVMEC
USE MOD_VMEC_Vars,     ONLY: rho,phinorm 
USE SPLINE1_MOD,       ONLY:SPLINE1_FIT 
USE SPLINE1_MOD,       ONLY:SPLINE1_INTERP 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: modes
INTEGER, INTENT(IN) :: mabs(modes)
REAL, INTENT(IN)    :: Xmn_half(modes,nFluxVMEC)  ! fourier coefficients at all flux surfaces 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL, INTENT(OUT)  :: Xmn_Spl(4,nFluxVMEC,modes)  ! spline fitted fourier coefficients 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: iMode,iFlux
REAL              :: Xmn_half_Spl(4,nFluxVMEC+1)  ! spline fitted fourier coefficients 
REAL              :: rho_half(1:nFluxVMEC+1)
INTEGER           :: iFlag
CHARACTER(len=100):: message
!===================================================================================================================================
DO iFlux=1,nFluxVMEC-1
  rho_half(iFlux+1)=SQRT(0.5*(phinorm(iFlux+1)+phinorm(iFlux))) !0.5*(rho(iFlux)+rho(iFlux+1))
END DO
!add end points
rho_half(1)=0.
rho_half(nFluxVMEC+1)=1.

DO iMode=1,modes
  !scaling with rho^|m|
  DO iFlux=2,nFluxVMEC
    IF(mabs(iMode).EQ.0)THEN
      Xmn_half_Spl(1,iFlux)=Xmn_half(iMode,iFlux)
    ELSE
      Xmn_half_Spl(1,iFlux)=Xmn_half(iMode,iFlux) /(rho_half(iFlux)**mabs(iMode))
    END IF
  END DO !i
  !Parabolic extrapolation to axis with dx'(rho=0)=0.
  Xmn_Half_Spl(1,1)=(Xmn_Half_Spl(1,2)*rho_half(3)**2-Xmn_Half_Spl(1,3)*rho_half(2)**2) /(rho_half(3)**2-rho_half(2)**2)
  !Extrapolate to Edge 
  Xmn_Half_Spl(1,nFluxVMEC+1)= ( Xmn_half_Spl(1,nFluxVMEC  )*(rho_half(nFluxVMEC+1)-rho_half(nFluxVMEC-1))     &
                                -Xmn_half_Spl(1,nFluxVMEC-1)*(rho_half(nFluxVMEC+1)-rho_half(nFluxVMEC  )) )   &
                                   /(  rho_half(nFluxVMEC)   -rho_half(nFluxVMEC-1) )
  CALL SPLINE1_FIT(nFluxVMEC+1,rho_half,Xmn_half_Spl(:,:), K_BC1=3, K_BCN=0)
  iflag=0
  message=''
  CALL SPLINE1_INTERP((/1,0,0/),nFluxVMEC+1,rho_half,Xmn_half_Spl, &
                                nFluxVMEC  ,rho     ,Xmn_Spl(:,:,iMode),       &
                          iflag,message, K_BC1=3,K_BCN=0)
  !respline
  Xmn_Spl(2:4,:,iMode)=0.
  Xmn_Spl(1,1,iMode)  =(Xmn_Spl(1,2,iMode)*rho(3)**2-Xmn_Spl(1,3,iMode)*rho(2)**2) /(rho(3)**2-rho(2)**2)
  CALL SPLINE1_FIT(nFluxVMEC,rho,Xmn_Spl(:,:,iMode), K_BC1=3, K_BCN=0)
END DO !iMode 

END SUBROUTINE FitSplineHalf



SUBROUTINE MapToVMEC(nTotal,x_in,InputCoordSys,xvmec,vmecData)
!===================================================================================================================================
! Maps a cylinder (r,z,phi) to a toroidal closed flux surface configuration derived from VMEC data. 
! Surfaces with constant r become flux surfaces. z [0;1] is mapped to [0;2*pi] 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_VMEC_Mappings, ONLY: mn_mode,xm,xn
USE MOD_VMEC_Mappings, ONLY: mn_mode_nyq,xm_nyq,xn_nyq
USE MOD_VMEC_Mappings, ONLY: nFluxVMEC
USE MOD_VMEC_Mappings, ONLY: mu0
USE MOD_VMEC_Mappings, ONLY: chi,phi
USE MOD_VMEC_Vars
USE SPLINE1_MOD, ONLY: SPLINE1_EVAL
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
REAL,INTENT(OUT)   :: xvmec(3,nTotal) ! mapped x,y,z coordinates with vmec data
REAL,INTENT(OUT)   :: vmecData(nVarVMEC,nTotal) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iNode,percent
INTEGER :: iMode
INTEGER :: iGuess=1
REAL    :: CosMN(mn_mode)
REAL    :: SinMN(mn_mode)
REAL    :: CosMN_nyq(mn_mode_nyq)
REAL    :: r_p    ! raduis in cylindrical coordinate system
REAL    :: phi_p  ! normalized toroidal flux (=flux coordinate s [0,1]) (use radial distance of point position)
REAL    :: chi_p  ! normalized poloidal flux [0,1] 
REAL    :: theta  ! poloidal angle [0,2pi]
REAL    :: zeta ! toroidal angle [0,2pi]
REAL    :: coszeta,sinzeta
REAL    :: R,Z   
REAL    :: sqrtG   
REAL    :: dRdrho,dRdtheta,dRdzeta      !derivatives
REAL    :: dZdrho,dZdtheta,dZdzeta 
REAL    :: phi_int,chi_int,phipf_int,iota_int
REAL    :: rho_p,rhom,drhom,splOut(3) !for interpolation
REAL    :: Density,Pressure
REAL    :: lam,dldtheta,dldzeta,dldrho
REAL    :: Btheta,Bzeta          ! also called B^u, B^v, contravariant components of the magnetic field
REAL    :: Br,Bz,Bphi            !mangetic field components in (R,Z,phi) system (phi=zeta)
                                 ! Br=dRdtheta*Btheta+dRdzeta*Bszeta 
                                 ! Bz=dZdtheta*Btheta+dZdzeta*Bzeta
                                 ! Bphi=R*Bzeta
REAL    :: Bcart(3)              !magnetic field components in (X,Y,Z) system 
                                 ! Bx=Br*cos(phi) - Bphi*sin(phi)
                                 ! By=Br*sin(phi) + Bphi*cos(phi)
                                 ! Bz=Bz
REAL    :: Arho,Atheta,Azeta
REAL    :: Ar,Az,Aphi,Acart(3)   !magnetic potential
REAL    :: ssqrtG
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,I8,A,A,A)')'  MAP ', nTotal,' NODES TO VMEC DATA FROM ',TRIM(VMECdataFile),' ...'
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
  
  CosMN(:)      = COS(    xm(:) * theta -     xn(:) * zeta)
  SinMN(:)      = SIN(    xm(:) * theta -     xn(:) * zeta) 
  CosMN_nyq(:)  = COS(xm_nyq(:) * theta - xn_nyq(:) * zeta)
  
  phi_p=r_p**2 ! use scaling of radius to phi evaluation variable
  
  phi_p=MIN(1.,MAX(phi_p,1.0E-08))
  rho_p=SQRT(phi_p)

  R         =0.
  Z         =0.
  dRdtheta  =0.
  dRdzeta   =0.
  dRdrho    =0.
  dZdtheta  =0.
  dZdzeta   =0.
  dZdrho    =0.
  lam       =0.
  dldtheta  =0.
  dldzeta   =0.
  dldrho    =0.
  DO iMode=1,mn_mode
    SELECT CASE(xmabs(iMode))
    CASE(0)
      rhom=1.
      drhom=0.
    CASE(1)
      rhom=rho_p
      drhom=1.
    CASE(2)
      rhom=rho_p*rho_p
      drhom=2*rho_p
    CASE DEFAULT
      rhom=rho_p**xmabs(iMode)
      drhom=xmabs(iMode)*rho_p**(xmabs(iMode)-1)
    END SELECT
    CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Rmnc_Spl(:,:,iMode),iGuess,splout) 
    R    = R    + rhom*splout(1)*CosMN(iMode)
    !dR/dtheta (cos(m*theta-n*zeta))=-m*sin(m*theta-n*zeta)
    dRdtheta = dRdtheta - rhom*splout(1)*SinMN(iMode)*xm(iMode) 
    !dR/dzeta  (cos(m*theta-n*zeta))=n*sin(m*theta-n*zeta)
    dRdzeta  = dRdzeta  + rhom*splout(1)*SinMN(iMode)*xn(iMode)
    !dR/drho = sum_mn [ rho**m * dR_mn/drho + R_mn *d(rho**m)/drho ] * cos(m*theta-n*zeta)
    dRdrho   = dRdrho   + (rhom*splout(2)+splout(1)*drhom)*CosMN(iMode)
    CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,Zmns_Spl(:,:,iMode),iGuess,splout) 
    Z    = Z    + rhom*splout(1)*SinMN(iMode)
    !dZ/dtheta (sin(m*theta-n*zeta))=m*cos(m*theta-n*zeta)
    dZdtheta = dZdtheta + rhom*splout(1)*CosMN(iMode)*xm(iMode) 
    !dZ/dzeta  (sin(m*theta-n*zeta))=-n*cos(m*theta-n*zeta)
    dZdzeta  = dZdzeta  - rhom*splout(1)*CosMN(iMode)*xn(iMode)
    !dZ/drho = sum_mn [ rho**m * dZ_mn/drho + Z_mn *d(rho**m)/drho ] * sin(m*theta-n*zeta)
    dZdrho   = dZdrho   + (rhom*splout(2)+splout(1)*drhom)*SinMN(iMode)
    !lambda
    CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho_p,rho,lmns_Spl(:,:,iMode),iGuess,splout) 
    lam      = lam      + rhom*splout(1)*SinMN(iMode)
    !dL/dtheta (sin(m*theta-n*zeta))=m*cos(m*theta-n*zeta)
    dldtheta = dldtheta + rhom*splout(1)*CosMN(iMode)*xm(iMode) 
    !dl/dzeta  (sin(m*theta-n*zeta))=-n*cos(m*theta-n*zeta)
    dldzeta  = dldzeta  - rhom*splout(1)*CosMN(iMode)*xn(iMode)
    !dl/drho = sum_mn [ rho**m * dl_mn/drho + l_mn *d(rho**m)/drho ] * sin(m*theta-n*zeta)
    dldrho   = dldrho   + (rhom*splout(2)+splout(1)*drhom)*SinMN(iMode)
  END DO !iMode=1,mn_mode 

  sqrtG     =0.
  DO iMode=1,mn_mode_nyq
    IF(xmabs_nyq(iMode).EQ.0)THEN
      rhom=1.
      drhom=0.
    ELSEIF(xmabs_nyq(iMode).EQ.1)THEN
      rhom=rho_p
      drhom=1.
    ELSE
      rhom=rho_p**xmabs_nyq(iMode)
      drhom=xmabs_nyq(iMode)*rho_p**(xmabs_nyq(iMode)-1)
    END IF
    CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_p,rho,gmnc_nyq_Spl(:,:,iMode),iGuess,splout) 
    sqrtG    = sqrtG    + rhom*splout(1)*CosMN_nyq(iMode)
  END DO !iMode=1,mn_mode_nyq 

  CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_p,rho,pres_Spl(:,:),iGuess,splout) 
  pressure=splout(1)
  CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_p,rho,chinorm_Spl(:,:),iGuess,splout) 
  chi_p=splout(1) !normalized chi (poloidal flux)
  CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_p,rho,phipf_Spl(:,:),iGuess,splout) 
  phipf_int=splout(1)
  CALL SPLINE1_EVAL((/1,0,0/), nFluxVMEC,rho_p,rho,iota_Spl(:,:),iGuess,splout) 
  iota_int=splout(1)
  

  !  !compute magnetic field, following Michael Kraus formulas: 
  !  ! B^s     = 0 
  !  ! B^theta = dphi/ds*(iota-dlambda/dzeta)  : phipf*(iotaf -  (lVmnc) ) 
  !  ! B^zeta = dphi/ds*(1+dlambda/dtheta)     : phipf*(1+ (lUmnc) )
  !  !   ...( lmns is overwritten to full mesh and then d/dtheta d/dzeta is applied )
  

  !do not use sqrtG
  !sqrtG=0.5/(rho_p)*R*(dRdtheta*dZdrho-dRdrho*dZdtheta) ! =drho/ds *R*(dRdtheta*dZdrho-dRdrho*dZdtheta)

  !contravariant components of B (B^s=0) !!! sqrtG=Jacobian: R*(dRdtheta*dZds-dRds*dZdtheta)
!  Btheta = phipf_int/sqrtG*(iota_int - dldzeta)
!  Bzeta  = phipf_int/sqrtG*(1.  + dldtheta)
!
!  Br   =  dRdtheta*Btheta+dRdzeta*Bzeta
!  Bz   =  dZdtheta*Btheta+dZdzeta*Bzeta
!  Bphi =  R*Bzeta

  !cylindrical components of B
!  Br   =  phipf_int*(dRdtheta*(iota_int-dldzeta)+dRdzeta*(1.  + dldtheta))*2*rho_p/(R*(dRdtheta*dZdrho-dRdrho*dZdtheta))
!  Bz   =  phipf_int*(dZdtheta*(iota_int-dldzeta)+dZdzeta*(1.  + dldtheta))*2*rho_p/(R*(dRdtheta*dZdrho-dRdrho*dZdtheta))
!  Bphi =  phipf_int*(1.  + dldtheta)*2*rho_p/(dRdtheta*dZdrho-dRdrho*dZdtheta)

  ! directly cylindrical components of B
  Br   =  phipf_int*(dRdtheta*(iota_int-dldzeta)+dRdzeta*(1.  + dldtheta))/sqrtG
  Bz   =  phipf_int*(dZdtheta*(iota_int-dldzeta)+dZdzeta*(1.  + dldtheta))/sqrtG
  Bphi =  phipf_int*(1.  + dldtheta)*R/sqrtG
 
  phi_int= phi(1)+(phi(nFluxVMEC)-phi(1))*phi_p 
  chi_int= chi(1)+(chi(nFluxVMEC)-chi(1))*chi_p 

!  !covariant components of A!!!
!  Arho   = 0.5/rho_p*phi_int*dldrho
!  Atheta = phi_int*(1+dldtheta)
!  Azeta  = phi_int*dldzeta-chi_int
!
!  !cylindrical components of A
!  Ar   = (-dZdtheta*Arho + 0.5/rho_p*dZdrho*Atheta )*R/sqrtG
!  Az   = ( dRdtheta*Arho - 0.5/rho_p*dRdrho*Atheta )*R/sqrtG
!  Aphi = -Azeta/R

  !directly cylindrical components of A
  Ar   = phi_int*(-dZdtheta*dldrho + dZdrho*(1+dldtheta) )/(dRdtheta*dZdrho-dRdrho*dZdtheta)
  Az   = phi_int*( dRdtheta*dldrho - dRdrho*(1+dldtheta) )/(dRdtheta*dZdrho-dRdrho*dZdtheta)
  Aphi = (chi_int-phi_int*dldzeta)/R
  
  !cylindrical components of A
!  Ar   = phi_int*(-dZdtheta*(0.5*dldrho/rho_p) + (0.5*dZdrho/rho_p)*(1+dldtheta) )*R/sqrtG
!  Az   = phi_int*( dRdtheta*(0.5*dldrho/rho_p) - (0.5*dRdrho/rho_p)*(1+dldtheta) )*R/sqrtG
!  Aphi = (chi_int-phi_int*dldzeta)/R



  coszeta=COS(zeta)
  sinzeta=SIN(zeta)
  
  xvmec(1,iNode)= R*coszeta
  xvmec(2,iNode)= R*sinzeta
  xvmec(3,iNode)= Z
  
  Bcart(1)= Br*coszeta-Bphi*sinzeta
  Bcart(2)= Br*sinzeta+Bphi*coszeta
  Bcart(3)= Bz

  Acart(1)= Ar*coszeta-Aphi*sinzeta
  Acart(2)= Ar*sinzeta+Aphi*coszeta
  Acart(3)= Az
  
  Density=EvalPoly(nRhoCoefs,RhoCoefs,MERGE(phi_p,chi_p,RhoFluxVar.EQ.0)) 

  vmecData(  1,iNode)=phi_p
  vmecData(  2,iNode)=chi_p
  vmecData(3:5,iNode)=Acart(:)
  vmecData(6:8,iNode)=Bcart(:)
  vmecData(  9,iNode)=pressure*mu0 !pressure transformed to mu0=1
  vmecData( 10,iNode)=Density
  vmecData( 11,iNode)=lam
  vmecData( 12,iNode)=dRdrho
  vmecData( 13,iNode)=dldrho
  vmecData( 14,iNode)=dZdrho
END DO !iNode=1,nTotal

WRITE(UNIT_stdOut,'(A)')'  ...DONE.                             '
END SUBROUTINE MapToVmec 
 

FUNCTION EvalPoly(nCoefs,Coefs,x)
!===================================================================================================================================
! evalute monomial polynomial c_1+c_2*x+c_3*x^2 ...
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)  :: nCoefs                   !number of coefficients 
REAL, INTENT(IN)     :: Coefs(nCoefs)            !coefficients
REAL, INTENT(IN)     :: x                        !evaluation position
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL              :: EvalPoly
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
!===================================================================================================================================
EvalPoly=0.
DO i=nCoefs,1,-1
  EvalPoly=EvalPoly*x+Coefs(i)
END DO

END FUNCTION EvalPoly


END MODULE MOD_VMEC
