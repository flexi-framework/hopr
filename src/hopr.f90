#include "defines.f90"
PROGRAM HOPR
!===================================================================================================================================
! control program of HOPR (high order preprocessor)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Basis,       ONLY:InitBasis
USE MOD_Mesh,        ONLY:InitMesh,FillMesh
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
WRITE(UNIT_stdOut,'(a,a,a)')' HOPR successfully finished: Mesh "',TRIM(ProjectName)//'_mesh.h5','" written to HDF5 file.'
WRITE(*,'(132("="))')
END PROGRAM HOPR
