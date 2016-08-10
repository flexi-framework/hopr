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
MODULE MOD_Readin_GMSH_Vars
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
PUBLIC

INTEGER                            :: bOrd
INTEGER,ALLOCATABLE                :: tetMapGMSH(:,:,:),  pyrMapGMSH(:,:,:),  priMapGMSH(:,:,:),  hexMapGMSH(:,:,:)
INTEGER                            :: tetMapCGNSToGMSH(4),pyrMapCGNSToGMSH(5),priMapCGNSToGMSH(6),hexMapCGNSToGMSH(8)
INTEGER                            :: GMSH_TYPES(6,131)
INTEGER                            :: nBCs_GMSH
INTEGER,ALLOCATABLE                :: MapBC(:),MapBCInd(:)

CONTAINS
SUBROUTINE buildTypes()
!===================================================================================================================================
! Build list of GMSH Types
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                   :: tmp(UBOUND(GMSH_TYPES,1)*UBOUND(GMSH_TYPES,2))  ! ?
!===================================================================================================================================
! BaseElement (GMSH), ElemType(our), nNodes, Degree, Complete, Dim 
tmp=&
(/2 ,  2   , 2   ,  1 ,  1, 2,& ! MSH_LIN_2   1  
  3 ,  3   , 3   ,  1 ,  1, 2,& ! MSH_TRI_3   2  
  4 ,  4   , 4   ,  1 ,  1, 2,& ! MSH_QUA_4   3  
  5 ,  104 , 4   ,  1 ,  1, 3,& ! MSH_TET_4   4  
  8 ,  108 , 8   ,  1 ,  1, 3,& ! MSH_HEX_8   5  
  7 ,  106 , 6   ,  1 ,  1, 3,& ! MSH_PRI_6   6  
  6 ,  105 , 5   ,  1 ,  1, 3,& ! MSH_PYR_5   7  
  2 ,  202 , 3   ,  2 ,  1, 2,& ! MSH_LIN_3   8  
  3 ,  6   , 6   ,  2 ,  1, 2,& ! MSH_TRI_6   9  
  4 ,  7   , 9   ,  2 ,  1, 2,& ! MSH_QUA_9   10 
  5 ,  204 , 10  ,  2 ,  1, 3,& ! MSH_TET_10  11 
  8 ,  208 , 27  ,  2 ,  1, 3,& ! MSH_HEX_27  12 
  7 ,  206 , 18  ,  2 ,  1, 3,& ! MSH_PRI_18  13 
  6 ,  205 , 14  ,  2 ,  1, 3,& ! MSH_PYR_14  14 
  1 ,  1   , 1   ,  1 ,  1, 1,& ! MSH_PNT     15 
  4 ,  7   , 8   ,  2 ,  0, 2,& ! MSH_QUA_8   16 
  8 ,  208 , 20  ,  2 ,  0, 3,& ! MSH_HEX_20  17 
  7 ,  206 , 15  ,  2 ,  0, 3,& ! MSH_PRI_15  18 
  6 ,  205 , 13  ,  2 ,  0, 3,& ! MSH_PYR_13  19 
  3 ,  5   , 9   ,  3 ,  0, 2,& ! MSH_TRI_9   20 
  3 ,  5   , 10  ,  3 ,  1, 2,& ! MSH_TRI_10  21 
  3 ,  5   , 12  ,  4 ,  0, 2,& ! MSH_TRI_12  22 
  3 ,  5   , 15  ,  4 ,  1, 2,& ! MSH_TRI_15  23 
  3 ,  5   , 15  ,  5 ,  0, 2,& ! MSH_TRI_15I 24 
  3 ,  5   , 21  ,  5 ,  1, 2,& ! MSH_TRI_21  25 
  2 ,  202 , 4   ,  3 ,  1, 2,& ! MSH_LIN_4   26 
  2 ,  202 , 5   ,  3 ,  1, 2,& ! MSH_LIN_5   27 
  2 ,  202 , 6   ,  3 ,  1, 2,& ! MSH_LIN_6   28 
  4 ,  204 , 20  ,  3 ,  1, 3,& ! MSH_TET_20  29 
  4 ,  204 , 35  ,  4 ,  1, 3,& ! MSH_TET_35  30 
  4 ,  204 , 56  ,  5 ,  1, 3,& ! MSH_TET_56  31 
  4 ,  204 , 34  ,  4 ,  0, 3,& ! MSH_TET_34  32 
  4 ,  204 , 52  ,  5 ,  0, 3,& ! MSH_TET_52  33 
  9 , -1   ,-1   , -1 , -1,-2,& ! MSH_POLYG_  34 
  10, -1   ,-1   , -1 , -1,-3,& ! MSH_POLYH_  35 
  4 ,  7   , 16  ,  3 ,  1, 2,& ! MSH_QUA_16  36 
  4 ,  7   , 25  ,  4 ,  1, 2,& ! MSH_QUA_25  37 
  4 ,  7   , 36  ,  5 ,  1, 2,& ! MSH_QUA_36  38 
  4 ,  7   , 12  ,  3 ,  0, 2,& ! MSH_QUA_12  39 
  4 ,  7   , 16  ,  4 ,  0, 2,& ! MSH_QUA_16I 40 
  4 ,  7   , 20  ,  5 ,  0, 2,& ! MSH_QUA_20  41 
  3 ,  5   , 28  ,  6 ,  1, 2,& ! MSH_TRI_28  42 
  3 ,  5   , 36  ,  7 ,  1, 2,& ! MSH_TRI_36  43 
  3 ,  5   , 45  ,  8 ,  1, 2,& ! MSH_TRI_45  44 
  3 ,  5   , 55  ,  9 ,  1, 2,& ! MSH_TRI_55  45 
  3 ,  5   , 66  ,  10,  1, 2,& ! MSH_TRI_66  46 
  4 ,  7   , 49  ,  6 ,  1, 2,& ! MSH_QUA_49  47 
  4 ,  7   , 64  ,  7 ,  1, 2,& ! MSH_QUA_64  48 
  4 ,  7   , 81  ,  8 ,  1, 2,& ! MSH_QUA_81  49 
  4 ,  7   , 100 ,  9 ,  1, 2,& ! MSH_QUA_100 50 
  4 ,  7   , 121 ,  10,  1, 2,& ! MSH_QUA_121 51 
  3 ,  5   , 18  ,  6 ,  0, 2,& ! MSH_TRI_18  52 
  3 ,  5   , 21  ,  7 ,  0, 2,& ! MSH_TRI_21I 53 
  3 ,  5   , 24  ,  8 ,  0, 2,& ! MSH_TRI_24  54 
  3 ,  5   , 27  ,  9 ,  0, 2,& ! MSH_TRI_27  55 
  3 ,  5   , 30  ,  10,  0, 2,& ! MSH_TRI_30  56 
  4 ,  7   , 24  ,  6 ,  0, 2,& ! MSH_QUA_24  57 
  4 ,  7   , 28  ,  7 ,  0, 2,& ! MSH_QUA_28  58 
  4 ,  7   , 32  ,  8 ,  0, 2,& ! MSH_QUA_32  59 
  4 ,  7   , 36  ,  9 ,  0, 2,& ! MSH_QUA_36I 60 
  4 ,  7   , 40  ,  10,  0, 2,& ! MSH_QUA_40  61 
  2 ,  202 , 7   ,  6 ,  1, 2,& ! MSH_LIN_7   62 
  2 ,  202 , 8   ,  7 ,  1, 2,& ! MSH_LIN_8   63 
  2 ,  202 , 9   ,  8 ,  1, 2,& ! MSH_LIN_9   64 
  2 ,  202 , 10  ,  9 ,  1, 2,& ! MSH_LIN_10  65 
  2 ,  202 , 11  ,  10,  1, 2,& ! MSH_LIN_11  66 
  2 , -1   ,-1   , -1 , -1,-2,& ! MSH_LIN_B   67 
  3 , -1   ,-1   , -1 , -1,-2,& ! MSH_TRI_B   68 
  9 , -1   ,-1   , -1 , -1,-2,& ! MSH_POLYG_B 69 
  2 , -1   ,-1   , -1 , -1,-2,& ! MSH_LIN_C   70 
  5 ,  204 , 84  ,  6 ,  1, 3,& ! MSH_TET_84  71 
  5 ,  204 , 120 ,  7 ,  1, 3,& ! MSH_TET_120 72 
  5 ,  204 , 165 ,  8 ,  1, 3,& ! MSH_TET_165 73 
  5 ,  204 , 220 ,  9 ,  1, 3,& ! MSH_TET_220 74 
  5 ,  204 , 286 ,  10,  1, 3,& ! MSH_TET_286 75 
  1 , -1   ,-1   , -1 , -1,-9,& ! dummy       76
  1 , -1   ,-1   , -1 , -1,-9,& ! dummy       77
  1 , -1   ,-1   , -1 , -1,-9,& ! dummy       78
  5 ,  204 , 74  ,  6 ,  0, 3,& ! MSH_TET_74  79 
  5 ,  204 , 100 ,  7 ,  0, 3,& ! MSH_TET_100 80 
  5 ,  204 , 130 ,  8 ,  0, 3,& ! MSH_TET_130 81 
  5 ,  204 , 164 ,  9 ,  0, 3,& ! MSH_TET_164 82 
  5 ,  204 , 202 ,  10,  0, 3,& ! MSH_TET_202 83 
  2 , -1   ,-1   , -1 , -1,-2,& ! MSH_LIN_1   84 
  3 , -1   ,-1   , -1 , -1,-2,& ! MSH_TRI_1   85 
  4 , -1   ,-1   , -1 , -1,-2,& ! MSH_QUA_1   86 
  5 , -1   ,-1   , -1 , -1,-3,& ! MSH_TET_1   87 
  8 , -1   ,-1   , -1 , -1,-3,& ! MSH_HEX_1   88 
  7 , -1   ,-1   , -1 , -1,-3,& ! MSH_PRI_1   89 
  7 ,  206 , 40  ,  3 ,  1, 3,& ! MSH_PRI_40  90 
  7 ,  206 , 75  ,  4 ,  1, 3,& ! MSH_PRI_75  91 
  8 ,  208 , 64  ,  3 ,  1, 3,& ! MSH_HEX_64  92 
  8 ,  208 , 125 ,  4 ,  1, 3,& ! MSH_HEX_125 93 
  8 ,  208 , 216 ,  5 ,  1, 3,& ! MSH_HEX_216 94 
  8 ,  208 , 343 ,  6 ,  1, 3,& ! MSH_HEX_343 95 
  8 ,  208 , 512 ,  7 ,  1, 3,& ! MSH_HEX_512 96 
  8 ,  208 , 729 ,  8 ,  1, 3,& ! MSH_HEX_729 97 
  8 ,  208 , 1000,  9 ,  1, 3,& ! MSH_HEX_100 98 
  8 ,  208 , 56  ,  3 ,  0, 3,& ! MSH_HEX_56  99 
  8 ,  208 , 98  ,  4 ,  0, 3,& ! MSH_HEX_98  100
  8 ,  208 , 152 ,  5 ,  0, 3,& ! MSH_HEX_152 101
  8 ,  208 , 222 ,  6 ,  0, 3,& ! MSH_HEX_222 102
  8 ,  208 , 296 ,  7 ,  0, 3,& ! MSH_HEX_296 103
  8 ,  208 , 386 ,  8 ,  0, 3,& ! MSH_HEX_386 104
  8 ,  208 , 488 ,  9 ,  0, 3,& ! MSH_HEX_488 105
  7 ,  206 , 126 ,  5 ,  1, 3,& ! MSH_PRI_126 106
  7 ,  206 , 196 ,  6 ,  1, 3,& ! MSH_PRI_196 107
  7 ,  206 , 288 ,  7 ,  1, 3,& ! MSH_PRI_288 108
  7 ,  206 , 405 ,  8 ,  1, 3,& ! MSH_PRI_405 109
  7 ,  206 , 550 ,  9 ,  1, 3,& ! MSH_PRI_550 110
  7 ,  206 , 38  ,  3 ,  0, 3,& ! MSH_PRI_38  111
  7 ,  206 , 66  ,  4 ,  0, 3,& ! MSH_PRI_66  112
  7 ,  206 , 102 ,  5 ,  0, 3,& ! MSH_PRI_102 113
  7 ,  206 , 146 ,  6 ,  0, 3,& ! MSH_PRI_146 114
  7 ,  206 , 198 ,  7 ,  0, 3,& ! MSH_PRI_198 115
  7 ,  206 , 258 ,  8 ,  0, 3,& ! MSH_PRI_258 116
  7 ,  206 , 326 ,  9 ,  0, 3,& ! MSH_PRI_326 117
  6 ,  205 , 30  ,  3 ,  1, 3,& ! MSH_PYR_30  118
  6 ,  205 , 55  ,  4 ,  1, 3,& ! MSH_PYR_55  119
  6 ,  205 , 91  ,  5 ,  1, 3,& ! MSH_PYR_91  120
  6 ,  205 , 140 ,  6 ,  1, 3,& ! MSH_PYR_140 121
  6 ,  205 , 204 ,  7 ,  1, 3,& ! MSH_PYR_204 122
  6 ,  205 , 285 ,  8 ,  1, 3,& ! MSH_PYR_285 123
  6 ,  205 , 385 ,  9 ,  1, 3,& ! MSH_PYR_385 124
  6 ,  205 , 29  ,  3 ,  0, 3,& ! MSH_PYR_29  125
  6 ,  205 , 50  ,  4 ,  0, 3,& ! MSH_PYR_50  126
  6 ,  205 , 77  ,  5 ,  0, 3,& ! MSH_PYR_77  127
  6 ,  205 , 110 ,  6 ,  0, 3,& ! MSH_PYR_110 128
  6 ,  205 , 149 ,  7 ,  0, 3,& ! MSH_PYR_149 129
  6 ,  205 , 194 ,  8 ,  0, 3,& ! MSH_PYR_194 130
  6 ,  205 , 245 ,  9 ,  0, 3/) ! MSH_PYR_245 131
GMSH_TYPES=RESHAPE(tmp,(/UBOUND(GMSH_TYPES,1),UBOUND(GMSH_TYPES,2)/))

END SUBROUTINE buildTypes

SUBROUTINE getGMSHVolumeMapping()
!===================================================================================================================================
! Define mapping from GMSH .msh structure to tensorproduct (i,j,k) structure
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
ALLOCATE(tetMapGMSH(0:bOrd-1,0:bOrd-1,0:bOrd-1),pyrMapGMSH(0:bOrd-1,0:bOrd-1,0:bOrd-1))
ALLOCATE(priMapGMSH(0:bOrd-1,0:bOrd-1,0:bOrd-1),hexMapGMSH(0:bOrd-1,0:bOrd-1,0:bOrd-1))
tetMapCGNSToGMSH = (/1,4,2,3/)
pyrMapCGNSToGMSH = (/1,2,3,4,5/)
priMapCGNSToGMSH = (/1,2,3,4,5,6/)
hexMapCGNSToGMSH = (/4,3,7,8,1,2,6,5/)

SELECT CASE(bOrd)
CASE(2)
  tetMapGMSH(0,0,0)= 1
  tetMapGMSH(1,0,0)= 4
  tetMapGMSH(0,1,0)= 2
  tetMapGMSH(0,0,1)= 3
  !
  pyrMapGMSH(0,0,0)= 1
  pyrMapGMSH(1,0,0)= 2
  pyrMapGMSH(0,1,0)= 4
  pyrMapGMSH(1,1,0)= 3
  pyrMapGMSH(0,0,1)= 5

  !priMapGMSH(0,0,0)= 1
  !priMapGMSH(1,0,0)= 2
  !priMapGMSH(0,1,0)= 3
  !priMapGMSH(0,0,1)= 4
  !priMapGMSH(1,0,1)= 5
  !priMapGMSH(0,1,1)= 6
  !
  hexMapGMSH(0,0,0)= 4
  hexMapGMSH(1,0,0)= 3
  hexMapGMSH(0,1,0)= 8
  hexMapGMSH(1,1,0)= 7
  hexMapGMSH(0,0,1)= 1
  hexMapGMSH(1,0,1)= 2
  hexMapGMSH(0,1,1)= 5
  hexMapGMSH(1,1,1)= 6

CASE(3)
  tetMapGMSH(0,0,0)=  1
  tetMapGMSH(1,0,0)=  8
  tetMapGMSH(2,0,0)=  4
  tetMapGMSH(0,1,0)=  5
  tetMapGMSH(1,1,0)= 10
  tetMapGMSH(0,2,0)=  2
  tetMapGMSH(0,0,1)=  7
  tetMapGMSH(1,0,1)=  9
  tetMapGMSH(0,1,1)=  6
  tetMapGMSH(0,0,2)=  3

  !pyrMapGMSH(0,0,0)=  1
  !pyrMapGMSH(2,0,0)=  2
  !pyrMapGMSH(2,2,0)=  3
  !pyrMapGMSH(0,2,0)=  4
  !pyrMapGMSH(0,0,2)=  5
  !pyrMapGMSH(1,0,0)=  6
  !pyrMapGMSH(0,1,0)=  7
  !pyrMapGMSH(0,0,1)=  8
  !pyrMapGMSH(2,1,0)=  9
  !pyrMapGMSH(1,0,1)= 10
  !pyrMapGMSH(1,2,0)= 11
  !pyrMapGMSH(1,1,1)= 12
  !pyrMapGMSH(0,1,1)= 13
  !pyrMapGMSH(1,1,0)= 14

  !priMapGMSH(0,0,0)=  1 
  !priMapGMSH(2,0,0)=  2 
  !priMapGMSH(0,2,0)=  3 
  !priMapGMSH(0,0,2)=  4 
  !priMapGMSH(2,0,2)=  5 
  !priMapGMSH(0,2,2)=  6 
  !priMapGMSH(1,0,0)=  7 
  !priMapGMSH(0,1,0)=  8 
  !priMapGMSH(0,0,1)=  9 
  !priMapGMSH(1,1,0)= 10
  !priMapGMSH(2,0,1)= 11
  !priMapGMSH(0,2,1)= 12
  !priMapGMSH(1,0,2)= 13
  !priMapGMSH(0,1,2)= 14
  !priMapGMSH(1,1,2)= 15
  !priMapGMSH(1,0,1)= 16
  !priMapGMSH(0,1,1)= 17
  !priMapGMSH(1,1,1)= 18

  hexMapGMSH(0,0,0)=  4
  hexMapGMSH(1,0,0)= 14
  hexMapGMSH(2,0,0)=  3
  hexMapGMSH(0,1,0)= 16
  hexMapGMSH(1,1,0)= 25
  hexMapGMSH(2,1,0)= 15
  hexMapGMSH(0,2,0)=  8
  hexMapGMSH(1,2,0)= 20
  hexMapGMSH(2,2,0)=  7
  hexMapGMSH(0,0,1)= 10
  hexMapGMSH(1,0,1)= 21
  hexMapGMSH(2,0,1)= 12
  hexMapGMSH(0,1,1)= 23
  hexMapGMSH(1,1,1)= 27
  hexMapGMSH(2,1,1)= 24
  hexMapGMSH(0,2,1)= 18
  hexMapGMSH(1,2,1)= 26
  hexMapGMSH(2,2,1)= 19
  hexMapGMSH(0,0,2)=  1
  hexMapGMSH(1,0,2)=  9
  hexMapGMSH(2,0,2)=  2
  hexMapGMSH(0,1,2)= 11
  hexMapGMSH(1,1,2)= 22
  hexMapGMSH(2,1,2)= 13
  hexMapGMSH(0,2,2)=  5
  hexMapGMSH(1,2,2)= 17
  hexMapGMSH(2,2,2)=  6
CASE(4)
  tetMapGMSH(0,0,0)=  1 
  tetMapGMSH(1,0,0)= 12 
  tetMapGMSH(2,0,0)= 11 
  tetMapGMSH(3,0,0)=  4 
  tetMapGMSH(0,1,0)=  5 
  tetMapGMSH(1,1,0)= 18 
  tetMapGMSH(2,1,0)= 15 
  tetMapGMSH(0,2,0)=  6 
  tetMapGMSH(1,2,0)= 16 
  tetMapGMSH(0,3,0)=  2
  tetMapGMSH(0,0,1)= 10
  tetMapGMSH(1,0,1)= 19
  tetMapGMSH(2,0,1)= 13
  tetMapGMSH(0,1,1)= 17
  tetMapGMSH(1,1,1)= 20
  tetMapGMSH(0,2,1)=  7
  tetMapGMSH(0,0,2)=  9
  tetMapGMSH(1,0,2)= 14
  tetMapGMSH(0,1,2)=  8
  tetMapGMSH(0,0,3)=  3

  hexMapGMSH(0,0,0)=  4
  hexMapGMSH(1,0,0)= 20
  hexMapGMSH(2,0,0)= 19
  hexMapGMSH(3,0,0)=  3
  hexMapGMSH(0,1,0)= 23
  hexMapGMSH(1,1,0)= 50
  hexMapGMSH(2,1,0)= 49
  hexMapGMSH(3,1,0)= 21
  hexMapGMSH(0,2,0)= 24
  hexMapGMSH(1,2,0)= 51
  hexMapGMSH(2,2,0)= 52
  hexMapGMSH(3,2,0)= 22
  hexMapGMSH(0,3,0)=  8
  hexMapGMSH(1,3,0)= 32
  hexMapGMSH(2,3,0)= 31
  hexMapGMSH(3,3,0)=  7
  hexMapGMSH(0,0,1)= 12
  hexMapGMSH(1,0,1)= 34
  hexMapGMSH(2,0,1)= 35
  hexMapGMSH(3,0,1)= 16
  hexMapGMSH(0,1,1)= 44
  hexMapGMSH(1,1,1)= 60
  hexMapGMSH(2,1,1)= 59
  hexMapGMSH(3,1,1)= 46
  hexMapGMSH(0,2,1)= 43
  hexMapGMSH(1,2,1)= 64
  hexMapGMSH(2,2,1)= 63
  hexMapGMSH(3,2,1)= 47
  hexMapGMSH(0,3,1)= 28
  hexMapGMSH(1,3,1)= 56
  hexMapGMSH(2,3,1)= 55
  hexMapGMSH(3,3,1)= 30
  hexMapGMSH(0,0,2)= 11
  hexMapGMSH(1,0,2)= 33
  hexMapGMSH(2,0,2)= 36
  hexMapGMSH(3,0,2)= 15
  hexMapGMSH(0,1,2)= 41
  hexMapGMSH(1,1,2)= 57
  hexMapGMSH(2,1,2)= 58
  hexMapGMSH(3,1,2)= 45
  hexMapGMSH(0,2,2)= 42
  hexMapGMSH(1,2,2)= 61
  hexMapGMSH(2,2,2)= 62
  hexMapGMSH(3,2,2)= 48
  hexMapGMSH(0,3,2)= 27
  hexMapGMSH(1,3,2)= 53
  hexMapGMSH(2,3,2)= 54
  hexMapGMSH(3,3,2)= 29
  hexMapGMSH(0,0,3)=  1
  hexMapGMSH(1,0,3)=  9
  hexMapGMSH(2,0,3)= 10
  hexMapGMSH(3,0,3)=  2
  hexMapGMSH(0,1,3)= 13
  hexMapGMSH(1,1,3)= 37
  hexMapGMSH(2,1,3)= 38
  hexMapGMSH(3,1,3)= 17
  hexMapGMSH(0,2,3)= 14
  hexMapGMSH(1,2,3)= 40
  hexMapGMSH(2,2,3)= 39
  hexMapGMSH(3,2,3)= 18
  hexMapGMSH(0,3,3)=  5
  hexMapGMSH(1,3,3)= 25
  hexMapGMSH(2,3,3)= 26
  hexMapGMSH(3,3,3)=  6

CASE DEFAULT
  CALL abort(__STAMP__,&
             'Elements of specified or higher order are not implemented yet. Order: ',bOrd)
END SELECT
END SUBROUTINE getGMSHVolumeMapping


END MODULE MOD_Readin_GMSH_Vars
