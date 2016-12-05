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
MODULE MOD_Solov_Vars
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER             :: setup       ! use given parameter setup: 10: circular, 20: iterlike
REAL                :: p_R0        !major radius
REAL                :: p_eps       !minor/major radius ratio
REAL                :: p_kappa     ! elliticity
REAL                :: p_delta     ! triangularity
REAL                :: asin_delta  ! ASIN(p_delta)
REAL                :: p_A         ! pA <1 (else deltap>0)
REAL                :: p_B0        ! toroidal magnetic field strength at magn. axis
REAL                :: p_qaxis     ! q-factor on axis
REAL                :: p_paxis     ! pressure on axis
INTEGER             :: nRhoCoefs   ! number of density coefficients 
REAL,ALLOCATABLE    :: RhoCoefs(:) !density coefficients of the polynomial coefficients:
!equilibrium related
REAL                :: xaxis(2)    !xy-position o of magnetic axis
REAL                :: psi_scale   !psi scaling between Soloviev and physical psiReal=psi_scale*psi
REAL                :: psi_axis    !psi value at the magnetic axis
REAL                :: psi_edge
REAL                :: presEdge      !pressure on edge, soloviev pressure profile linear in psinorm
REAL                :: F_axis     ! F^2 on axis, f-profile, dF^2/dpsi=const for soloviev 
REAL                :: deltaF2    ! F=SQRT(F_axis^2+deltaF2*psinorm)
REAL                :: psiCoefs(0:7) !coefficients for representing flux variable psi

!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
END MODULE MOD_Solov_Vars

