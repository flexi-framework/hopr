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
!===================================================================================================================================
! Here, preprocessor variables for different equation systems and abbrevbiations for specific expressions are defined
!===================================================================================================================================

! Abbreviations
!-----------------------------------------------------------------------------------------------------------------------------------
#ifdef SUN
#define __DATE__ '__TIME__ and __DATE__ not'
#define __TIME__ 'available for SUN COMPILER'
#define IEEE_ISNAN
#elif PGI
#define NO_ISNAN
#endif
#ifndef __FILENAME__ 
#define __FILENAME__ __FILE__
#endif
#define __STAMP__ __FILENAME__,__LINE__,__DATE__,__TIME__

#ifdef GNU
#  define IEEE_IS_NAN ISNAN
#endif

#define SIZEOF_F(x) STORAGE_SIZE(x)/8

#ifdef GNU
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  k).OR.x<-HUGE(1_  k))       CALL ABORT(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ k).OR.x<-HUGE(1._ k))       CALL ABORT(__STAMP__,'Real conversion failed: out of range!')
#else
#define CHECKSAFEINT(x,k)  IF(x>HUGE(1_  ## k).OR.x<-HUGE(1_  ## k)) CALL ABORT(__STAMP__,'Integer conversion failed: out of range!')
#define CHECKSAFEREAL(x,k) IF(x>HUGE(1._ ## k).OR.x<-HUGE(1._ ## k)) CALL ABORT(__STAMP__,'Real conversion failed: out of range!')
#endif

#define SDEALLOCATE(A) IF(ASSOCIATED(A)) DEALLOCATE(A)
#define ERRWRITE(a,b) WRITE(UNIT_errFile,b)
#define LOGWRITE(a,b) IF(Logging) WRITE(UNIT_logOut,b)

! General tolerance for mesh node comparisons and sfc (1.E.-13 should be ok too)
#define PP_MeshTolerance 1.E-12
! General tolerance for Real comparisons (double precision)
#define PP_RealTolerance 1.E-16

! Settings for different compilers
!-----------------------------------------------------------------------------------------------------------------------------------
! Colored output
#if defined(GNU)||defined(PGI)||defined(BLUEGENE)
! No colors available (still)
#define __FGCOLORRED__    ''
#define __FGCOLORGREEN__  ''
#define __FGCOLORYELLOW__ ''
#define __FGCOLORBLUE__   ''
#define __FGCOLOREND__    ''
#else
! Use default colors
#define __FGCOLORRED__    '\033[31m'
#define __FGCOLORGREEN__  '\033[32m'
#define __FGCOLORYELLOW__ '\033[33m'
#define __FGCOLORBLUE__   '\033[34m'
#define __FGCOLOREND__    '\033[0m'
#endif

#if(PP_CGNS_INT==64)
#  define PP_CGNS_INT_TYPE INTEGER(KIND=8)
#endif
#if(PP_CGNS_INT==32)
#  define PP_CGNS_INT_TYPE INTEGER
#endif

