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
REAL                 :: maxMode_m,maxMode_n
LOGICAL              :: killmode
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
  ! grid for evaluation of fourier modes (default: 4,4) , should be set to (128, 128) if smoothing is needed
  intPointsU=GETINT("VMEC_intPointsU","4") 
  intPointsV=GETINT("VMEC_intPointsV","4")
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

!!  llocate memory for new values
  CALL AllocArrays

  !! calculate data for mapping (include recalculation of lambda and magnetic
  !! field after smoothing r, z)
  CALL PreCalcData

  ALLOCATE(phinorm(nFluxVMEC))
  phinorm=(phi-phi(1))/phi(nFluxVMEC)

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
  maxmode_m=MAXVAL(xm)
  maxmode_n=MAXVAL(xn)
  WRITE(UNIT_stdOut,*)'   Total Number of mn-modes:',mn_mode
  WRITE(UNIT_stdOut,*)'   Max Mode m,n: ',maxmode_m,maxmode_n
  WRITE(UNIT_stdOut,*)'   Number of even(m) and odd(m) mn-modes:',mn_mEven,mn_mOdd
  !find even and odd m-modes, to seperately evalute them
  mn_mEven_nyq=0
  DO iMode=1,mn_mode_nyq
    IF(MOD(xm_nyq(iMode),2.).EQ.0.) mn_mEven_nyq=mn_mEven_nyq+1
  END DO ! i=1,mn_mode_nyq

  mn_mOdd_nyq=mn_mode_nyq-mn_mEven_nyq
  ALLOCATE(mn_mapOdd_nyq(mn_mOdd_nyq),mn_mapEven_nyq(mn_mEven_nyq))
  iEven=0
  iOdd=0
  DO iMode=1,mn_mode_nyq
    IF(MOD(xm_nyq(iMode),2.).EQ.0.) THEN
      iEven=iEven+1
      mn_mapEven_nyq(iEven)=iMode
    ELSE
      iOdd=iOdd+1
      mn_mapOdd_nyq(iOdd)=iMode
    END IF !even
  END DO ! i=1,mn_mode_nyq
  WRITE(UNIT_stdOut,*)'   Total Number of mn-modes (Nyquist):',mn_mode_nyq
  WRITE(UNIT_stdOut,*)'   Max Mode m,n: ',MAXVAL(xm_nyq),MAXVAL(xn_nyq)
  WRITE(UNIT_stdOut,*)'   Number of even(m) and odd(m) mn-modes (Nyquist):',mn_mEven_nyq,mn_mOdd_nyq

!  useFilter=GETLOGICAL('VMECuseFilter','.FALSE.')
!  !find even and odd m-modes, to seperately evalute them AND FILTER HIGH MODES
!  IF(.NOT.useFilter)THEN
!    ALLOCATE(filtMap(mn_mode_nyq))
!    mn_mode_filt=mn_mode_nyq
!    DO iMode=1,mn_mode_nyq
!        filtMap(iMode)=iMode
!    END DO ! i=1,mn_mode_nyq
!  ELSE  
!    mn_mode_filt=0
!    DO iMode=1,mn_mode_nyq
!      IF((xm_nyq(iMode).LE.maxmode_m).AND.(ABS(xn_nyq(iMode)).LE.maxmode_n)) THEN
!        mn_mode_filt=mn_mode_filt+1
!      END IF
!    END DO ! i=1,mn_mode_nyq
!    ALLOCATE(filtMap(mn_mode_filt))
!    iFilt=0
!    DO iMode=1,mn_mode_nyq
!      IF((xm_nyq(iMode).LE.maxmode_m).AND.(ABS(xn_nyq(iMode)).LE.maxmode_n)) THEN
!        iFilt=iFilt+1
!        filtMap(iFilt)=iMode
!      END IF
!    END DO ! i=1,mn_mode_nyq
!  END IF!useFilter
!
!  mn_mEven_filt=0
!  DO iFilt=1,mn_mode_filt
!    iMode=FiltMap(iFilt)
!    IF(MOD(xm_nyq(iMode),2.).EQ.0.) mn_mEven_filt=mn_mEven_filt+1
!  END DO ! i=1,mn_mode_filt
!  
!  mn_mOdd_filt=mn_mode_filt-mn_mEven_filt
!  ALLOCATE(mn_mapOdd_filt(mn_mOdd_filt),mn_mapEven_filt(mn_mEven_filt))
!  iEven=0
!  iOdd=0
!  DO iFilt=1,mn_mode_filt
!    iMode=FiltMap(iFilt)
!    IF(MOD(xm_nyq(iMode),2.).EQ.0.) THEN
!      iEven=iEven+1
!      mn_mapEven_filt(iEven)=iMode
!    ELSE
!      iOdd=iOdd+1
!      mn_mapOdd_filt(iOdd)=iMode
!    END IF !even
!  END DO ! i=1,mn_mode_nyq
!  WRITE(UNIT_stdOut,*)'   Total Number of mn-modes (filtered):',mn_mode_filt
!  WRITE(UNIT_stdOut,*)'   Max Mode m,n: ',maxmode_m,maxmode_n
!  WRITE(UNIT_stdOut,*)'   Number of even(m) and odd(m) mn-modes (filtered):',mn_mEven_filt,mn_mOdd_filt

  !OUTPUT
  VMECvarnames(1)='VMEC-phinorm'
  VMECvarnames(2)='VMEC-Br'
  VMECvarnames(3)='VMEC-Bz'
  VMECvarnames(4)='VMEC-Bphi'
  VMECvarnames(5)='VMEC-absB'
  VMECvarnames(6)='VMEC-B_X'
  VMECvarnames(7)='VMEC-B_Y'
  VMECvarnames(8)='VMEC-B_Z'
  VMECvarnames(9)='VMEC-Pressure'
  VMECvarnames(10)='VMEC-10'
  VMECvarnames(11)='VMEC-11'
  VMECvarnames(12)='VMEC-12'
  VMECvarnames(13)='VMEC-13'
  VMECvarnames(14)='VMEC-14'
  VMECvarnames(15)='VMEC-15'
END IF !useVMEC

WRITE(UNIT_stdOut,'(A)')'... DONE'
END SUBROUTINE InitVMEC


SUBROUTINE MapToVMEC(xcyl,xvmec,vmecData)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_VMEC_Mappings, ONLY: mn_mode,xm,xn,Rmnc, Zmns,dRdUmns,dRdVmns,dZdUmnc,dZdVmnc,lUmnc,lVmnc
USE MOD_VMEC_Mappings, ONLY: mn_mode_nyq,xm_nyq,xn_nyq
USE MOD_VMEC_Mappings, ONLY: Bsupumnc_nyq,Bsupvmnc_nyq,Bmnc_nyq
USE MOD_VMEC_Mappings, ONLY: nFluxVMEC,phi,presf
USE MOD_VMEC_Mappings, ONLY: phipf,iotas,lUmnc,lVmnc
USE MOD_VMEC_Vars,     ONLY: phinorm,mn_mEven,mn_mOdd,mn_mapOdd,mn_mapEven
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)  :: xcyl(3) ! x,y,z coordinates in a cylinder of size r=0,1, z=0,1
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)  :: xvmec(3) ! mapped x,y,z coordinates with vmec data
REAL,INTENT(OUT)  :: vmecData(PP_nVarVMEC) 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,s1,s2
REAL    :: CosMN(mn_mode),SinMN(mn_mode)
REAL    :: CosMN_nyq(mn_mode_nyq)
!REAL    :: CosMN_nyq(mn_mode_nyq)
!REAL    :: SinMN_nyq(mn_mode_nyq)
REAL    :: phi_p  ! flux coordinate [0,1] (use radial distance of point position)
REAL    :: theta  ! poloidal angle [0,2pi]
REAL    :: zeta ! toroidal angle [0,2pi]
REAL    :: coszeta,sinzeta
REAL    :: R,Z   
REAL    :: f1,f2,w1,w2 
REAL    :: Bsupu,Bsupv           !covariant magnetic field components, B^s=0, B^u=B^theta, B^v=B^zeta
REAL    :: dRdU,dRdV,dZdU,dZdV   !derivatives of R and Z
REAL    :: Bnorm                 !|B| from VMEC
REAL    :: Br,Bz,Bphi            !mangetic field components in (R,Z,phi) system (phi=zeta)
                                 ! Br=dRdu*Bsupu+dRdv*Bsupv, Bz=dZdu*Bsupu+dZdv*Bsupv, Bphi=R*Bsupv
REAL    :: Bcart(3)              !magnetic field components in (X,Y,Z) system 
                                 ! Bx=Br*cos(phi) - Bphi*sin(phi)
                                 ! By=Br*sin(phi) + Bphi*cos(phi)
                                 ! Bz=Bz
REAL    :: Pressure
REAL    :: phipf_int,iotas_int,Btheta,Bzeta
!===================================================================================================================================
phi_p  = SQRT(xcyl(1)**2+xcyl(2)**2) 
theta   = ATAN2(xcyl(2),xcyl(1))
zeta   = -2*Pi*xcyl(3) ! which coordinate system is not clear, this gives correct geometry and fields

CosMN(:)  = COS(xm(:) * theta - xn(:) * zeta)
SinMN(:)  = SIN(xm(:) * theta - xn(:) * zeta)

CosMN_nyq(:)  = COS(xm_nyq(:) * theta - xn_nyq(:) * zeta)
!SinMN_nyq(:)  = SIN(xm_nyq(:) * theta - xn_nyq(:) * zeta) !not yet needed

phi_p=phi_p**2 ! use scaling of radius to phi evaluation variable

!! look for the nearest supporting point phinorm(s1) < phi_p < phinorm(s1) 
IF(phi_p.LT.1.0E-08)THEN
  s1=1
  s2=1
ELSE
  s1=nFluxVMEC
  DO i = 1, nFluxVMEC
    IF (phi_p .LT. phinorm(i)) THEN
     s1=i
     EXIT
    END IF
  END DO
  s2=MIN(s1+1,nFluxVMEC)
END IF
!! exception: extrapolate for phi_p<phinorm(2) 
!!            -> use s1=2 instead of s1=1, because 1/sqrt(phinorm(s1))=1/0. is not defined.
!s1=nFluxVMEC
!DO i = 2, nFluxVMEC
!  IF (phi_p .LT. phinorm(i)) THEN
!   s1=i
!   EXIT
!  END IF
!END DO
!s2=MIN(s1+1,nFluxVMEC)

!WRITE(*,*)'DEBUG,s1,s2,phi1,phi2,phi',s1,s2,phinorm(s1),phinorm(s2),phi_p

IF(s1.NE.s2)THEN
  !interpolation factor
  f2=(phi_p-phinorm(s1))/(phinorm(s2)-phinorm(s1))
  f1=1.-f2 
  ! interpolation of odd modes with weighting sqrt(phi) * ((1-frac) * r1/sqrt(phi_1) + frac* r2/sqrt(phi_2))
  w1=MERGE(0.,SQRT(phi_p/phinorm(s1)),s1.EQ.1) !avoid division by 0 
  w2=SQRT(phi_p/phinorm(s2))

  R=InterpolateData(f1,f2,w1,w2,CosMN,Rmnc(:,s1),Rmnc(:,s2))

  Z=InterpolateData(f1,f2,w1,w2,SinMN,Zmns(:,s1),Zmns(:,s2))


  Bsupu  = InterpolateData_nyq(f1,f2,w1,w2,CosMN_nyq,bsupumnc_nyq(:,s1),bsupumnc_nyq(:,s2))
  Bsupv  = InterpolateData_nyq(f1,f2,w1,w2,CosMN_nyq,bsupvmnc_nyq(:,s1),bsupvmnc_nyq(:,s2))
  Bnorm  = InterpolateData_nyq(f1,f2,w1,w2,CosMN_nyq,Bmnc_nyq(:,s1),Bmnc_nyq(:,s2))


  dRdu =InterpolateData(f1,f2,w1,w2,SinMN,dRdUmns(:,s1),dRdUmns(:,s2))
  dRdv =InterpolateData(f1,f2,w1,w2,SinMN,dRdVmns(:,s1),dRdVmns(:,s2))

  dZdu =InterpolateData(f1,f2,w1,w2,cosMN,dZdUmnc(:,s1),dZdUmnc(:,s2))
  dZdv =InterpolateData(f1,f2,w1,w2,cosMN,dZdVmnc(:,s1),dZdVmnc(:,s2))

  pressure = f1*presf(s1)+f2*presf(s2)
!  !compute magnetic field, following Michael Kraus formulas: 
!  ! B^s     = 0 
!  ! B^theta = dphi/ds*(iota-dlambda/dzeta)  : phipf*(iotaf -  (lVmnc) ) 
!  ! B^zeta = dphi/ds*(1+dlambda/dtheta)     : phipf*(1+ (lUmnc) )
!  !   ...( lmns is overwritten to full mesh and then d/du d/dv is applied )


!  phipf_int = (f1*phipf(s1)+f2*phipf(s2))
!  iotas_int = f1*iotas(s1)+f2*iotas(s2) !iotas overwritten to full mesh 
!
!  Btheta = phipf_int*(iotas_int - InterpolateData(f1,f2,w1,w2,CosMN,lVmnc(:,s1),lVmnc(:,s2)))
!  Bzeta  = phipf_int*(1.        + InterpolateData(f1,f2,w1,w2,CosMN,lUmnc(:,s1),lUmnc(:,s2)))
!  ! ==> DO NOT MATCH WITH Bsupu, Bsupv !!!
ELSE
  !weighting with sqrt(s) cancels, evaluate all modes at s2.
  R      = SUM(CosMN(:)*Rmnc(:, s2))
  Z      = SUM(SinMN(:)*Zmns(:, s2))

  Bsupu  = SUM(CosMN_nyq(:)*bsupumnc_nyq(:,s2))
  Bsupv  = SUM(CosMN_nyq(:)*bsupvmnc_nyq(:,s2))
  Bnorm  = SUM(CosMN_nyq(:)*Bmnc_nyq(:,s2))


  dRdu   = SUM(SinMN(:)*dRdUmns(:,s2))
  dRdv   = SUM(SinMN(:)*dRdVmns(:,s2))

  dZdu   = SUM(cosMN(:)*dZdUmnc(:,s2))
  dZdv   = SUM(cosMN(:)*dZdVmnc(:,s2))

  pressure = presf(s2)

!  Btheta = phipf(s2)*(iotas(s2) - SUM(CosMN*lVmnc(:,s2)))
!  Bzeta  = phipf(s2)*(1.        + SUM(CosMN*lUmnc(:,s2)))
END IF !s1/=s2

Br   =dRdu*Bsupu+dRdv*Bsupv
Bz   =dZdu*Bsupu+dZdv*Bsupv
Bphi =R*Bsupv

!Br   =dRdu*Btheta+dRdv*Bzeta
!Bz   =dZdu*Btheta+dZdv*Bzeta
!Bphi =R*Bzeta

!IF(s1.EQ.2)THEN
!  WRITE(*,*)'s1,s2',s1,s2
!  WRITE(*,*)'f1,f2',f1,f2
!  WRITE(*,*)'w1,w2',w1,w2
!  WRITE(*,*)'R,Z',R,Z
!  WRITE(*,'(A,3E15.5,A,3E15.5)')' (s,theta/Pi,zeta/Pi) ',phi_p,theta/Pi,zeta/Pi,  &
!                           ' compare: |B|',Bnorm,SQRT(Br*Br+Bz*Bz+Bphi*Bphi), Bnorm-SQRT(Br*Br+Bz*Bz+Bphi*Bphi)
!END IF

!WRITE(*,*)'s1,s2',s1,s2,phinorm(s1),phinorm(s2)
!WRITE(*,'(A,3E14.5)') ' (s,theta/Pi,zeta/Pi) ',phi_p,theta/Pi,zeta/Pi
!WRITE(*,'(2(A,2E14.5))') ' Btheta,Bsupu ',Btheta,Bsupu,' Bzeta,Bsupv ', Bzeta,Bsupv


coszeta=COS(zeta)
sinzeta=SIN(zeta)

xvmec(1)= R*coszeta
xvmec(2)= R*sinzeta
xvmec(3)= Z

Bcart(1)= Br*coszeta-Bphi*sinzeta
Bcart(2)= Br*sinzeta+Bphi*coszeta
Bcart(3)= Bz

vmecData(1)=phi_p
vmecData(2)=Br
vmecData(3)=Bz
vmecData(4)=Bphi
vmecData(5)=Bnorm !VMEC |B|
vmecData(6:8)=Bcart(:)
vmecData(9)=pressure
vmecData(10)=Bsupu
vmecData(11)=Bsupv
vmecData(12)=dRdu
vmecData(13)=dZdu
vmecData(14)=dRdv
vmecData(15)=dZdv

END SUBROUTINE MapToVmec 


FUNCTION InterpolateData(f1,f2,w1,w2,trig_mn,Xmn_s1,Xmn_s2)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_VMEC_Mappings, ONLY: mn_mode
USE MOD_VMEC_Vars, ONLY: mn_mEven,mn_mOdd,mn_mapOdd,mn_mapEven
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)  :: f1,f2            ! linear interpolation weights (f1+f2=1)
REAL, INTENT(IN)  :: w1,w2            ! extra-weights for odd modes (phi/sqrt(phi_1),phi/sqrt(phi_2))
REAL, INTENT(IN)  :: trig_mn(mn_mode) ! sin or cos evaluated at  theta (m),zeta (n)
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

FUNCTION InterpolateData_nyq(f1,f2,w1,w2,trig_mn_nyq,Xmn_nyq_s1,Xmn_nyq_s2)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_VMEC_Mappings, ONLY: mn_mode_nyq
USE MOD_VMEC_Vars, ONLY: mn_mEven_nyq,mn_mOdd_nyq,mn_mapOdd_nyq,mn_mapEven_nyq
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL, INTENT(IN)  :: f1,f2                    ! linear interpolation weights (f1+f2=1)
REAL, INTENT(IN)  :: w1,w2                    ! extra-weights for odd modes (phi/sqrt(phi_1),phi/sqrt(phi_2))
REAL, INTENT(IN)  :: trig_mn_nyq(mn_mode_nyq) ! sin or cos evaluated at theta (m),zeta (n)
REAL, INTENT(IN)  :: Xmn_nyq_s1( mn_mode_nyq) ! variable in fourier space to interpolate
REAL, INTENT(IN)  :: Xmn_nyq_s2( mn_mode_nyq) ! variable in fourier space to interpolate
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL              :: InterpolateData_nyq
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: val1e,val2e,val1o,val2o
!===================================================================================================================================
  !evaluate only modes with m= even
  val1e = SUM(trig_mn_nyq(mn_mapEven_nyq)*xmn_nyq_s1(mn_mapEven_nyq))
  val2e = SUM(trig_mn_nyq(mn_mapEven_nyq)*xmn_nyq_s2(mn_mapEven_nyq))
  !evaluate only modes with m= odd
  val1o = SUM(trig_mn_nyq(mn_mapOdd_nyq)*xmn_nyq_s1(mn_mapOdd_nyq))
  val2o = SUM(trig_mn_nyq(mn_mapOdd_nyq)*xmn_nyq_s2(mn_mapOdd_nyq))

  InterpolateData_nyq=f1*(val1e+w1*val1o) + f2*(val2e+w2*val2o)
END FUNCTION InterpolateData_nyq


END MODULE MOD_VMEC
