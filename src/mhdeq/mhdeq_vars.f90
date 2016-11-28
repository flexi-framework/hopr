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
MODULE MOD_MHDEQ_Vars
!===================================================================================================================================
! Variables for the MHD equilibrium mapping / data
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
INTEGER             :: WhichEquilibrium                ! 0: do nothing (default), 1: use VMEC data, 2: use Solov'ev equilibrium
LOGICAL             :: useMHDEQ                        ! =(whichEquilibrium>0)
INTEGER             :: nVarMHDEQ=10
REAL,ALLOCATABLE    :: MHDEQoutdataGL(:,:,:,:,:)       ! MHD equilibrium data to be written to hdf5 file, on Gauss-Lobatto nodes
REAL,ALLOCATABLE    :: MHDEQdataEq(:,:,:,:,:)          ! VMEC data on equidistant nodes (forvisualization) 
CHARACTER(LEN=255),DIMENSION(10),PARAMETER :: MHDEQvarNames(10)=(/ CHARACTER(LEN=255) :: &
                      'MHDEQ-Density'     & ! 1 
                     ,'MHDEQ-Pressure'    & ! 2
                     ,'MHDEQ-B_X'         & ! 3
                     ,'MHDEQ-B_Y'         & ! 4
                     ,'MHDEQ-B_Z'         & ! 5
                     ,'MHDEQ-polfluxnorm' & ! 6
                     ,'MHDEQ-torfluxnorm' & ! 7
                     ,'MHDEQ-A_X'         & ! 8
                     ,'MHDEQ-A_Y'         & ! 9
                     ,'MHDEQ-A_Z'         & !10    
                                         /)

!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
END MODULE MOD_MHDEQ_Vars

