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
MODULE MOD_MHDEQ_Tools
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

INTERFACE Eval1DPoly 
  MODULE PROCEDURE Eval1DPoly
END INTERFACE

PUBLIC::Eval1DPoly
!===================================================================================================================================

CONTAINS
 

FUNCTION Eval1DPoly(nCoefs,Coefs,x)
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
REAL              :: Eval1DPoly
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i
!===================================================================================================================================
Eval1DPoly=0.
DO i=nCoefs,1,-1
  Eval1DPoly=Eval1DPoly*x+Coefs(i)
END DO

END FUNCTION Eval1DPoly

END MODULE MOD_MHDEQ_Tools
