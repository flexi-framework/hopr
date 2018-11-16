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
! Copyright (C) 2017 Claus-Dieter Munz <munz@iag.uni-stuttgart.de>
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
PROGRAM HOPR
!===================================================================================================================================
! control program of HOPR (high order preprocessor)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Basis,       ONLY:InitBasis
USE MOD_Mesh,        ONLY:InitMesh,FillMesh
USE MOD_Mesh_Vars,   ONLY:negativeJacobians,jacobianTolerance
USE MOD_Output,      ONLY:InitOutput
USE MOD_ReadInTools, ONLY:IgnoredStrings
USE MOD_Search,      ONLY:InitSearch
#ifdef _OPENMP
USE omp_lib
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! PROGRAM VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
!===================================================================================================================================
CALL SetStackSizeUnlimited()
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') &
"           _______     _______    ___________________    ___________________   ___________________                         xX    "
WRITE(UNIT_stdOut,'(A)') &
"          /      /)   /      /)  /                  /)  /                  /) /                  /)     .xXXXXXXXXx.       X     "
WRITE(UNIT_stdOut,'(A)') &
"         /      //   /      //  /      _____       //  /      _____       // /      _____       //    .XXXXXXXXXXXXXXXx  .XXXXx  "
WRITE(UNIT_stdOut,'(A)') &
"        /      //   /      //  /      /)___/      //  /      /)___/      // /      /)___/      //   .XXXXXXXXXXXXXXXXXXXXXXXXXXx "
WRITE(UNIT_stdOut,'(A)') &
"       /      //___/      //  /      //   /      //  /      //___/      // /      //___/      //  .XXXXXXXXXXXXXXXXXXXXXXXX´     "
WRITE(UNIT_stdOut,'(A)') &
"      /                  //  /      //   /      //  /                  // /                  //  .XX``XXXXXXXXXXXXXXXXXX´        "
WRITE(UNIT_stdOut,'(A)') &
"     /      _____       //  /      //   /      //  /      ____________// /      __      ____//   XX`  `XXXXX`      .X´           "
WRITE(UNIT_stdOut,'(A)') &
"    /      /)___/      //  /      //   /      //  /      /)___________) /      /)_|    |____)   XX     XXX`       .´             "
WRITE(UNIT_stdOut,'(A)') &
"   /      //   /      //  /      //___/      //  /      //             /      //  |    |__     ,X`    XXX´                       "
WRITE(UNIT_stdOut,'(A)') &
"  /      //   /      //  /                  //  /      //             /      //   |      /)   ,X`   .XX´                         "
WRITE(UNIT_stdOut,'(A)') &
" /______//   /______//  /__________________//  /______//             /______//    |_____//   ,X`   XX´                           "
WRITE(UNIT_stdOut,'(A)') &
" )______)    )______)   )__________________)   )______)              )______)     )_____)   xX    XXx                            "
WRITE(UNIT_stdOut,'(A)')
WRITE(UNIT_stdOut,'(132("="))')

! Initialization
! hm, didn't know were to put this...
PI  = ACOS(-1.)
rPI = SQRT(PI)
! read .ini file, build mesh, basis, initial condition,...
CALL InitOutput()
CALL InitMesh()
CALL InitBasis()
CALL InitSearch()
CALL IgnoredStrings()
! Now build mesh!
CALL FillMesh()
WRITE(*,'(132("="))')
IF(negativeJacobians.GT.0) THEN
  WRITE(UNIT_stdOut,'(A,A,A)')' HOPR finished: Mesh "',TRIM(ProjectName)//'_mesh.h5','" written to HDF5 file.'
  WRITE(UNIT_stdOut,'(A,I8,A,E11.3,A)')' WARNING: ',negativeJacobians, &
                                       ' ELEMENT(S) WITH SCALED JACOBIAN < ',jacobianTolerance,' FOUND!!!'
ELSE
  WRITE(UNIT_stdOut,'(A,A,A)')' HOPR successfully finished: Mesh "',TRIM(ProjectName)//'_mesh.h5','" written to HDF5 file.'
END IF
WRITE(*,'(132("="))')
END PROGRAM HOPR
