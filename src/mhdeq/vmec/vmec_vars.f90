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
! Copyright (C) 2017  Florian Hindenlang <hindenlang@gmail.com>
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
MODULE MOD_VMEC_Vars
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
LOGICAL             :: useVMEC
CHARACTER(LEN = 256):: VMECdataFile
INTEGER,ALLOCATABLE :: xmAbs(:)                  ! abs |xm(iMode)|
REAL,ALLOCATABLE    :: Psi_prof(:)               ! TOROIDAL flux profile (called phi_pf in VMEC)
REAL,ALLOCATABLE    :: Psinorm_prof(:)           ! normalized TOROIDAL flux profile 
REAL,ALLOCATABLE    :: chi_prof(:)               ! POLOIDAL flux profile ( called chi_pf in VMEC)

REAL,ALLOCATABLE    :: rho(:)                    ! := sqrt(psinorm) at all flux surface 
REAL,ALLOCATABLE    :: pres_Spl(:,:)             ! Spline coefficients in (rho) for Pressure, iota 
REAL,ALLOCATABLE    :: Psi_Spl(:,:)       
REAL,ALLOCATABLE    :: chi_Spl(:,:)       
REAL,ALLOCATABLE    :: Rmnc_Spl(:,:,:)           ! modified spline coefficients of Rmnc
REAL,ALLOCATABLE    :: Zmns_Spl(:,:,:)           !
REAL,ALLOCATABLE    :: lmns_Spl(:,:,:)           !
INTEGER             :: nRhoCoefs                 ! number of density coefficients 
INTEGER             :: RhoFluxVar                ! =0: rho(psinorm) Normalized toroidal flux variable, =1: rho(chinorm) 
REAL,ALLOCATABLE    :: RhoCoefs(:)               ! density coefficients of the polynomial coefficients:
                                                 ! rho_1+rho_2*x + rho_3*x^2 ...
!not used anymore
!INTEGER,ALLOCATABLE :: xmAbs_nyq(:)              ! abs |xm(iMode)|
!REAL,ALLOCATABLE    :: gmnc_nyq_Spl(:,:,:)       !
!REAL,ALLOCATABLE    :: iota_Spl(:,:)       
!REAL,ALLOCATABLE    :: dPsi_ds_Spl(:,:)          ! two modes in vmec: if toroidal flux is used for profiles, 
!                                                 ! then s=Psinorm [0,1], Psi(s)=Psi(1)+(Psi(n)-Psi(1))*s,dPsi/ds=Psi(n)-Psi(1). 
!                                                 ! VMEC can also use the normalized polodial fluxes for s (RFP=.TRUE. option),
!                                                 ! then dPsi/ds = dPsi/dchinorm = (chi(n)-chi(1))/iota
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
! Vandermande and D matrices
!===================================================================================================================================
END MODULE MOD_VMEC_Vars

