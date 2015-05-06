#include "defines.f90"
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
INTEGER,ALLOCATABLE                :: tetMapCGNSToGMSH(:),pyrMapCGNSToGMSH(:),priMapCGNSToGMSH(:),hexMapCGNSToGMSH(:)
INTEGER                            :: GMSH_TYPES(6,131)

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
  7 ,  107 , 6   ,  1 ,  1, 3,& ! MSH_PRI_6   6  
  6 ,  106 , 5   ,  1 ,  1, 3,& ! MSH_PYR_5   7  
  2 ,  202 , 3   ,  2 ,  1, 2,& ! MSH_LIN_3   8  
  3 ,  6   , 6   ,  2 ,  1, 2,& ! MSH_TRI_6   9  
  4 ,  7   , 9   ,  2 ,  1, 2,& ! MSH_QUA_9   10 
  5 ,  204 , 10  ,  2 ,  1, 3,& ! MSH_TET_10  11 
  8 ,  208 , 27  ,  2 ,  1, 3,& ! MSH_HEX_27  12 
  7 ,  207 , 18  ,  2 ,  1, 3,& ! MSH_PRI_18  13 
  6 ,  206 , 14  ,  2 ,  1, 3,& ! MSH_PYR_14  14 
  1 ,  1   , 1   ,  1 ,  1, 1,& ! MSH_PNT     15 
  4 ,  7   , 8   ,  2 ,  0, 2,& ! MSH_QUA_8   16 
  8 ,  208 , 20  ,  2 ,  0, 3,& ! MSH_HEX_20  17 
  7 ,  207 , 15  ,  2 ,  0, 3,& ! MSH_PRI_15  18 
  6 ,  206 , 13  ,  2 ,  0, 3,& ! MSH_PYR_13  19 
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
  7 ,  207 , 40  ,  3 ,  1, 3,& ! MSH_PRI_40  90 
  7 ,  207 , 75  ,  4 ,  1, 3,& ! MSH_PRI_75  91 
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
  7 ,  207 , 126 ,  5 ,  1, 3,& ! MSH_PRI_126 106
  7 ,  207 , 196 ,  6 ,  1, 3,& ! MSH_PRI_196 107
  7 ,  207 , 288 ,  7 ,  1, 3,& ! MSH_PRI_288 108
  7 ,  207 , 405 ,  8 ,  1, 3,& ! MSH_PRI_405 109
  7 ,  207 , 550 ,  9 ,  1, 3,& ! MSH_PRI_550 110
  7 ,  207 , 38  ,  3 ,  0, 3,& ! MSH_PRI_38  111
  7 ,  207 , 66  ,  4 ,  0, 3,& ! MSH_PRI_66  112
  7 ,  207 , 102 ,  5 ,  0, 3,& ! MSH_PRI_102 113
  7 ,  207 , 146 ,  6 ,  0, 3,& ! MSH_PRI_146 114
  7 ,  207 , 198 ,  7 ,  0, 3,& ! MSH_PRI_198 115
  7 ,  207 , 258 ,  8 ,  0, 3,& ! MSH_PRI_258 116
  7 ,  207 , 326 ,  9 ,  0, 3,& ! MSH_PRI_326 117
  6 ,  206 , 30  ,  3 ,  1, 3,& ! MSH_PYR_30  118
  6 ,  206 , 55  ,  4 ,  1, 3,& ! MSH_PYR_55  119
  6 ,  206 , 91  ,  5 ,  1, 3,& ! MSH_PYR_91  120
  6 ,  206 , 140 ,  6 ,  1, 3,& ! MSH_PYR_140 121
  6 ,  206 , 204 ,  7 ,  1, 3,& ! MSH_PYR_204 122
  6 ,  206 , 285 ,  8 ,  1, 3,& ! MSH_PYR_285 123
  6 ,  206 , 385 ,  9 ,  1, 3,& ! MSH_PYR_385 124
  6 ,  206 , 29  ,  3 ,  0, 3,& ! MSH_PYR_29  125
  6 ,  206 , 50  ,  4 ,  0, 3,& ! MSH_PYR_50  126
  6 ,  206 , 77  ,  5 ,  0, 3,& ! MSH_PYR_77  127
  6 ,  206 , 110 ,  6 ,  0, 3,& ! MSH_PYR_110 128
  6 ,  206 , 149 ,  7 ,  0, 3,& ! MSH_PYR_149 129
  6 ,  206 , 194 ,  8 ,  0, 3,& ! MSH_PYR_194 130
  6 ,  206 , 245 ,  9 ,  0, 3/) ! MSH_PYR_245 131
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
ALLOCATE(hexMapCGNSToGMSH(bOrd**3))
SELECT CASE(bOrd)
CASE(2)
  tetMapGMSH(0,0,0)= 1
  tetMapGMSH(1,0,0)= 2
  tetMapGMSH(0,1,0)= 3
  tetMapGMSH(0,0,1)= 4
  
  pyrMapGMSH(0,0,0)= 1
  pyrMapGMSH(1,0,0)= 2
  pyrMapGMSH(1,1,0)= 3
  pyrMapGMSH(0,1,0)= 4
  pyrMapGMSH(0,0,1)= 5

  priMapGMSH(0,0,0)= 1
  priMapGMSH(1,0,0)= 2
  priMapGMSH(0,1,0)= 3
  priMapGMSH(0,0,1)= 4
  priMapGMSH(1,0,1)= 5
  priMapGMSH(0,1,1)= 6
 ! hexaeder auf cgns angepasst, andere stehen noch aus 
  hexMapGMSH(0,0,0)= 1
  hexMapGMSH(0,0,1)= 2
  hexMapGMSH(0,1,1)= 3
  hexMapGMSH(0,1,0)= 4
  hexMapGMSH(1,0,0)= 5
  hexMapGMSH(1,0,1)= 6
  hexMapGMSH(1,1,1)= 7
  hexMapGMSH(1,1,0)= 8

  hexMapCGNSToGMSH = (/1,5,8,4,2,6,7,3/)
CASE(3)
  tetMapGMSH(0,0,0)=  1
  tetMapGMSH(2,0,0)=  2
  tetMapGMSH(0,2,0)=  3
  tetMapGMSH(0,0,2)=  4
  tetMapGMSH(1,0,0)=  5
  tetMapGMSH(1,1,0)=  6
  tetMapGMSH(0,1,0)=  7
  tetMapGMSH(0,0,1)=  8
  tetMapGMSH(0,1,1)=  9
  tetMapGMSH(1,0,1)= 10

  pyrMapGMSH(0,0,0)=  1
  pyrMapGMSH(2,0,0)=  2
  pyrMapGMSH(2,2,0)=  3
  pyrMapGMSH(0,2,0)=  4
  pyrMapGMSH(0,0,2)=  5
  pyrMapGMSH(1,0,0)=  6
  pyrMapGMSH(0,1,0)=  7
  pyrMapGMSH(0,0,1)=  8
  pyrMapGMSH(2,1,0)=  9
  pyrMapGMSH(1,0,1)= 10
  pyrMapGMSH(1,2,0)= 11
  pyrMapGMSH(1,1,1)= 12
  pyrMapGMSH(0,1,1)= 13
  pyrMapGMSH(1,1,0)= 14

  priMapGMSH(0,0,0)=  1 
  priMapGMSH(2,0,0)=  2 
  priMapGMSH(0,2,0)=  3 
  priMapGMSH(0,0,2)=  4 
  priMapGMSH(2,0,2)=  5 
  priMapGMSH(0,2,2)=  6 
  priMapGMSH(1,0,0)=  7 
  priMapGMSH(0,1,0)=  8 
  priMapGMSH(0,0,1)=  9 
  priMapGMSH(1,1,0)= 10
  priMapGMSH(2,0,1)= 11
  priMapGMSH(0,2,1)= 12
  priMapGMSH(1,0,2)= 13
  priMapGMSH(0,1,2)= 14
  priMapGMSH(1,1,2)= 15
  priMapGMSH(1,0,1)= 16
  priMapGMSH(0,1,1)= 17
  priMapGMSH(1,1,1)= 18

  hexMapGMSH(0,0,0)=  1 
  hexMapGMSH(2,0,0)=  2 
  hexMapGMSH(2,2,0)=  3 
  hexMapGMSH(0,2,0)=  4 
  hexMapGMSH(0,0,2)=  5 
  hexMapGMSH(2,0,2)=  6 
  hexMapGMSH(2,2,2)=  7 
  hexMapGMSH(0,2,2)=  8 
  hexMapGMSH(1,0,0)=  9 
  hexMapGMSH(0,1,0)= 10
  hexMapGMSH(0,0,1)= 11
  hexMapGMSH(2,1,0)= 12
  hexMapGMSH(2,0,1)= 13
  hexMapGMSH(1,2,0)= 14
  hexMapGMSH(2,2,1)= 15
  hexMapGMSH(0,2,1)= 16
  hexMapGMSH(1,0,2)= 17
  hexMapGMSH(0,1,2)= 18
  hexMapGMSH(2,1,2)= 19
  hexMapGMSH(1,2,2)= 20
  hexMapGMSH(1,1,0)= 21
  hexMapGMSH(1,0,1)= 22
  hexMapGMSH(0,1,1)= 23
  hexMapGMSH(2,1,1)= 24
  hexMapGMSH(1,2,1)= 25
  hexMapGMSH(1,1,2)= 26
  hexMapGMSH(1,1,1)= 27
CASE(4)
  tetMapGMSH(0,0,0)=  1 
  tetMapGMSH(3,0,0)=  2 
  tetMapGMSH(0,3,0)=  3 
  tetMapGMSH(0,0,3)=  4 
  tetMapGMSH(1,0,0)=  5 
  tetMapGMSH(2,0,0)=  6 
  tetMapGMSH(2,1,0)=  7 
  tetMapGMSH(1,2,0)=  8 
  tetMapGMSH(0,2,0)=  9 
  tetMapGMSH(0,1,0)= 10
  tetMapGMSH(0,0,2)= 11
  tetMapGMSH(0,0,1)= 12
  tetMapGMSH(0,1,2)= 13
  tetMapGMSH(0,2,1)= 14
  tetMapGMSH(1,0,2)= 15
  tetMapGMSH(2,0,1)= 16
  tetMapGMSH(1,1,0)= 17
  tetMapGMSH(1,0,1)= 18
  tetMapGMSH(0,1,1)= 19
  tetMapGMSH(1,1,1)= 20

  hexMapGMSH(0,0,0)=  1 
  hexMapGMSH(3,0,0)=  2 
  hexMapGMSH(3,3,0)=  3 
  hexMapGMSH(0,3,0)=  4 
  hexMapGMSH(0,0,3)=  5 
  hexMapGMSH(3,0,3)=  6 
  hexMapGMSH(3,3,3)=  7 
  hexMapGMSH(0,3,3)=  8 
  hexMapGMSH(1,0,0)=  9 
  hexMapGMSH(2,0,0)= 10
  hexMapGMSH(0,1,0)= 11
  hexMapGMSH(0,2,0)= 12
  hexMapGMSH(0,0,1)= 13
  hexMapGMSH(0,0,2)= 14
  hexMapGMSH(3,1,0)= 15
  hexMapGMSH(3,2,0)= 16
  hexMapGMSH(3,0,1)= 17
  hexMapGMSH(3,0,2)= 18
  hexMapGMSH(2,3,0)= 19
  hexMapGMSH(1,3,0)= 20
  hexMapGMSH(3,3,1)= 21
  hexMapGMSH(3,3,2)= 22
  hexMapGMSH(0,3,1)= 23
  hexMapGMSH(0,3,2)= 24
  hexMapGMSH(1,0,3)= 25
  hexMapGMSH(2,0,3)= 26
  hexMapGMSH(0,1,3)= 27
  hexMapGMSH(0,2,3)= 28
  hexMapGMSH(3,1,3)= 29
  hexMapGMSH(3,2,3)= 30
  hexMapGMSH(2,3,3)= 31
  hexMapGMSH(1,3,3)= 32
  hexMapGMSH(1,1,0)= 33
  hexMapGMSH(1,2,0)= 34
  hexMapGMSH(2,2,0)= 35
  hexMapGMSH(2,1,0)= 36
  hexMapGMSH(1,0,1)= 37
  hexMapGMSH(2,0,1)= 38
  hexMapGMSH(2,0,2)= 39
  hexMapGMSH(1,0,2)= 40
  hexMapGMSH(0,1,1)= 41
  hexMapGMSH(0,1,2)= 42
  hexMapGMSH(0,2,2)= 43
  hexMapGMSH(0,2,1)= 44
  hexMapGMSH(3,1,1)= 45
  hexMapGMSH(3,2,1)= 46
  hexMapGMSH(3,2,2)= 47
  hexMapGMSH(3,1,2)= 48
  hexMapGMSH(2,3,1)= 49
  hexMapGMSH(1,3,1)= 50
  hexMapGMSH(1,3,2)= 51
  hexMapGMSH(2,3,2)= 52
  hexMapGMSH(1,1,3)= 53
  hexMapGMSH(2,1,3)= 54
  hexMapGMSH(2,2,3)= 55
  hexMapGMSH(1,2,3)= 56
  hexMapGMSH(1,1,1)= 57
  hexMapGMSH(2,1,1)= 58
  hexMapGMSH(2,2,1)= 59
  hexMapGMSH(1,2,1)= 60
  hexMapGMSH(1,1,2)= 61
  hexMapGMSH(2,1,2)= 62
  hexMapGMSH(2,2,2)= 63
  hexMapGMSH(1,2,2)= 64
CASE(5)
  tetMapGMSH(0,0,0)=  1 
  tetMapGMSH(4,0,0)=  2 
  tetMapGMSH(0,4,0)=  3 
  tetMapGMSH(0,0,4)=  4 
  tetMapGMSH(1,0,0)=  5 
  tetMapGMSH(2,0,0)=  6 
  tetMapGMSH(3,0,0)=  7 
  tetMapGMSH(3,1,0)=  8 
  tetMapGMSH(2,2,0)=  9 
  tetMapGMSH(1,3,0)= 10
  tetMapGMSH(0,3,0)= 11
  tetMapGMSH(0,2,0)= 12
  tetMapGMSH(0,1,0)= 13
  tetMapGMSH(0,0,3)= 14
  tetMapGMSH(0,0,2)= 15
  tetMapGMSH(0,0,1)= 16
  tetMapGMSH(0,1,3)= 17
  tetMapGMSH(0,2,2)= 18
  tetMapGMSH(0,3,1)= 19
  tetMapGMSH(1,0,3)= 20
  tetMapGMSH(2,0,2)= 21
  tetMapGMSH(3,0,1)= 22
  tetMapGMSH(1,1,0)= 23
  tetMapGMSH(1,2,0)= 24
  tetMapGMSH(2,1,0)= 25
  tetMapGMSH(1,0,1)= 26
  tetMapGMSH(2,0,1)= 27
  tetMapGMSH(1,0,2)= 28
  tetMapGMSH(0,1,1)= 29
  tetMapGMSH(0,1,2)= 30
  tetMapGMSH(0,2,1)= 31
  tetMapGMSH(1,1,2)= 32
  tetMapGMSH(2,1,1)= 33
  tetMapGMSH(1,2,1)= 34
  tetMapGMSH(1,1,1)= 35

  hexMapGMSH(0,0,0)=   1  
  hexMapGMSH(4,0,0)=   2  
  hexMapGMSH(4,4,0)=   3  
  hexMapGMSH(0,4,0)=   4  
  hexMapGMSH(0,0,4)=   5  
  hexMapGMSH(4,0,4)=   6  
  hexMapGMSH(4,4,4)=   7  
  hexMapGMSH(0,4,4)=   8  
  hexMapGMSH(1,0,0)=   9  
  hexMapGMSH(2,0,0)=  10 
  hexMapGMSH(3,0,0)=  11 
  hexMapGMSH(0,1,0)=  12 
  hexMapGMSH(0,2,0)=  13 
  hexMapGMSH(0,3,0)=  14 
  hexMapGMSH(0,0,1)=  15 
  hexMapGMSH(0,0,2)=  16 
  hexMapGMSH(0,0,3)=  17 
  hexMapGMSH(4,1,0)=  18 
  hexMapGMSH(4,2,0)=  19 
  hexMapGMSH(4,3,0)=  20 
  hexMapGMSH(4,0,1)=  21 
  hexMapGMSH(4,0,2)=  22 
  hexMapGMSH(4,0,3)=  23 
  hexMapGMSH(3,4,0)=  24 
  hexMapGMSH(2,4,0)=  25 
  hexMapGMSH(1,4,0)=  26 
  hexMapGMSH(4,4,1)=  27 
  hexMapGMSH(4,4,2)=  28 
  hexMapGMSH(4,4,3)=  29 
  hexMapGMSH(0,4,1)=  30 
  hexMapGMSH(0,4,2)=  31 
  hexMapGMSH(0,4,3)=  32 
  hexMapGMSH(1,0,4)=  33 
  hexMapGMSH(2,0,4)=  34 
  hexMapGMSH(3,0,4)=  35 
  hexMapGMSH(0,1,4)=  36 
  hexMapGMSH(0,2,4)=  37 
  hexMapGMSH(0,3,4)=  38 
  hexMapGMSH(4,1,4)=  39 
  hexMapGMSH(4,2,4)=  40 
  hexMapGMSH(4,3,4)=  41 
  hexMapGMSH(3,4,4)=  42 
  hexMapGMSH(2,4,4)=  43 
  hexMapGMSH(1,4,4)=  44 
  hexMapGMSH(1,1,0)=  45 
  hexMapGMSH(1,3,0)=  46 
  hexMapGMSH(3,3,0)=  47 
  hexMapGMSH(3,1,0)=  48 
  hexMapGMSH(1,2,0)=  49 
  hexMapGMSH(2,3,0)=  50 
  hexMapGMSH(3,2,0)=  51 
  hexMapGMSH(2,1,0)=  52 
  hexMapGMSH(2,2,0)=  53 
  hexMapGMSH(1,0,1)=  54 
  hexMapGMSH(3,0,1)=  55 
  hexMapGMSH(3,0,3)=  56 
  hexMapGMSH(1,0,3)=  57 
  hexMapGMSH(2,0,1)=  58 
  hexMapGMSH(3,0,2)=  59 
  hexMapGMSH(2,0,3)=  60 
  hexMapGMSH(1,0,2)=  61 
  hexMapGMSH(2,0,2)=  62 
  hexMapGMSH(0,1,1)=  63 
  hexMapGMSH(0,1,3)=  64 
  hexMapGMSH(0,3,3)=  65 
  hexMapGMSH(0,3,1)=  66 
  hexMapGMSH(0,1,2)=  67 
  hexMapGMSH(0,2,3)=  68 
  hexMapGMSH(0,3,2)=  69 
  hexMapGMSH(0,2,1)=  70 
  hexMapGMSH(0,2,2)=  71 
  hexMapGMSH(4,1,1)=  72 
  hexMapGMSH(4,3,1)=  73 
  hexMapGMSH(4,3,3)=  74 
  hexMapGMSH(4,1,3)=  75 
  hexMapGMSH(4,2,1)=  76 
  hexMapGMSH(4,3,2)=  77 
  hexMapGMSH(4,2,3)=  78 
  hexMapGMSH(4,1,2)=  79 
  hexMapGMSH(4,2,2)=  80 
  hexMapGMSH(3,4,1)=  81 
  hexMapGMSH(1,4,1)=  82 
  hexMapGMSH(1,4,3)=  83 
  hexMapGMSH(3,4,3)=  84 
  hexMapGMSH(2,4,1)=  85 
  hexMapGMSH(1,4,2)=  86 
  hexMapGMSH(2,4,3)=  87 
  hexMapGMSH(3,4,2)=  88 
  hexMapGMSH(2,4,2)=  89 
  hexMapGMSH(1,1,4)=  90 
  hexMapGMSH(3,1,4)=  91 
  hexMapGMSH(3,3,4)=  92 
  hexMapGMSH(1,3,4)=  93 
  hexMapGMSH(2,1,4)=  94 
  hexMapGMSH(3,2,4)=  95 
  hexMapGMSH(2,3,4)=  96 
  hexMapGMSH(1,2,4)=  97 
  hexMapGMSH(2,2,4)=  98 
  hexMapGMSH(1,1,1)=  99 
  hexMapGMSH(3,1,1)= 100
  hexMapGMSH(3,3,1)= 101
  hexMapGMSH(1,3,1)= 102
  hexMapGMSH(1,1,3)= 103
  hexMapGMSH(3,1,3)= 104
  hexMapGMSH(3,3,3)= 105
  hexMapGMSH(1,3,3)= 106
  hexMapGMSH(2,1,1)= 107
  hexMapGMSH(1,2,1)= 108
  hexMapGMSH(1,1,2)= 109
  hexMapGMSH(3,2,1)= 110
  hexMapGMSH(3,1,2)= 111
  hexMapGMSH(2,3,1)= 112
  hexMapGMSH(3,3,2)= 113
  hexMapGMSH(1,3,2)= 114
  hexMapGMSH(2,1,3)= 115
  hexMapGMSH(1,2,3)= 116
  hexMapGMSH(3,2,3)= 117
  hexMapGMSH(2,3,3)= 118
  hexMapGMSH(2,2,1)= 119
  hexMapGMSH(2,1,2)= 120
  hexMapGMSH(1,2,2)= 121
  hexMapGMSH(3,2,2)= 122
  hexMapGMSH(2,3,2)= 123
  hexMapGMSH(2,2,3)= 124
  hexMapGMSH(2,2,2)= 125
CASE(6)
  tetMapGMSH(0,0,0)=  1 
  tetMapGMSH(5,0,0)=  2 
  tetMapGMSH(0,5,0)=  3 
  tetMapGMSH(0,0,5)=  4 
  tetMapGMSH(1,0,0)=  5 
  tetMapGMSH(2,0,0)=  6 
  tetMapGMSH(3,0,0)=  7 
  tetMapGMSH(4,0,0)=  8 
  tetMapGMSH(4,1,0)=  9 
  tetMapGMSH(3,2,0)= 10
  tetMapGMSH(2,3,0)= 11
  tetMapGMSH(1,4,0)= 12
  tetMapGMSH(0,4,0)= 13
  tetMapGMSH(0,3,0)= 14
  tetMapGMSH(0,2,0)= 15
  tetMapGMSH(0,1,0)= 16
  tetMapGMSH(0,0,4)= 17
  tetMapGMSH(0,0,3)= 18
  tetMapGMSH(0,0,2)= 19
  tetMapGMSH(0,0,1)= 20
  tetMapGMSH(0,1,4)= 21
  tetMapGMSH(0,2,3)= 22
  tetMapGMSH(0,3,2)= 23
  tetMapGMSH(0,4,1)= 24
  tetMapGMSH(1,0,4)= 25
  tetMapGMSH(2,0,3)= 26
  tetMapGMSH(3,0,2)= 27
  tetMapGMSH(4,0,1)= 28
  tetMapGMSH(1,1,0)= 29
  tetMapGMSH(1,3,0)= 30
  tetMapGMSH(3,1,0)= 31
  tetMapGMSH(1,2,0)= 32
  tetMapGMSH(2,2,0)= 33
  tetMapGMSH(2,1,0)= 34
  tetMapGMSH(1,0,1)= 35
  tetMapGMSH(3,0,1)= 36
  tetMapGMSH(1,0,3)= 37
  tetMapGMSH(2,0,1)= 38
  tetMapGMSH(2,0,2)= 39
  tetMapGMSH(1,0,2)= 40
  tetMapGMSH(0,1,1)= 41
  tetMapGMSH(0,1,3)= 42
  tetMapGMSH(0,3,1)= 43
  tetMapGMSH(0,1,2)= 44
  tetMapGMSH(0,2,2)= 45
  tetMapGMSH(0,2,1)= 46
  tetMapGMSH(1,1,3)= 47
  tetMapGMSH(3,1,1)= 48
  tetMapGMSH(1,3,1)= 49
  tetMapGMSH(2,1,2)= 50
  tetMapGMSH(2,2,1)= 51
  tetMapGMSH(1,2,2)= 52
  tetMapGMSH(1,1,1)= 53
  tetMapGMSH(2,1,1)= 54
  tetMapGMSH(1,2,1)= 55
  tetMapGMSH(1,1,2)= 56

  hexMapGMSH(0,0,0)=   1  
  hexMapGMSH(5,0,0)=   2  
  hexMapGMSH(5,5,0)=   3  
  hexMapGMSH(0,5,0)=   4  
  hexMapGMSH(0,0,5)=   5  
  hexMapGMSH(5,0,5)=   6  
  hexMapGMSH(5,5,5)=   7  
  hexMapGMSH(0,5,5)=   8  
  hexMapGMSH(1,0,0)=   9  
  hexMapGMSH(2,0,0)=  10 
  hexMapGMSH(3,0,0)=  11 
  hexMapGMSH(4,0,0)=  12 
  hexMapGMSH(0,1,0)=  13 
  hexMapGMSH(0,2,0)=  14 
  hexMapGMSH(0,3,0)=  15 
  hexMapGMSH(0,4,0)=  16 
  hexMapGMSH(0,0,1)=  17 
  hexMapGMSH(0,0,2)=  18 
  hexMapGMSH(0,0,3)=  19 
  hexMapGMSH(0,0,4)=  20 
  hexMapGMSH(5,1,0)=  21 
  hexMapGMSH(5,2,0)=  22 
  hexMapGMSH(5,3,0)=  23 
  hexMapGMSH(5,4,0)=  24 
  hexMapGMSH(5,0,1)=  25 
  hexMapGMSH(5,0,2)=  26 
  hexMapGMSH(5,0,3)=  27 
  hexMapGMSH(5,0,4)=  28 
  hexMapGMSH(4,5,0)=  29 
  hexMapGMSH(3,5,0)=  30 
  hexMapGMSH(2,5,0)=  31 
  hexMapGMSH(1,5,0)=  32 
  hexMapGMSH(5,5,1)=  33 
  hexMapGMSH(5,5,2)=  34 
  hexMapGMSH(5,5,3)=  35 
  hexMapGMSH(5,5,4)=  36 
  hexMapGMSH(0,5,1)=  37 
  hexMapGMSH(0,5,2)=  38 
  hexMapGMSH(0,5,3)=  39 
  hexMapGMSH(0,5,4)=  40 
  hexMapGMSH(1,0,5)=  41 
  hexMapGMSH(2,0,5)=  42 
  hexMapGMSH(3,0,5)=  43 
  hexMapGMSH(4,0,5)=  44 
  hexMapGMSH(0,1,5)=  45 
  hexMapGMSH(0,2,5)=  46 
  hexMapGMSH(0,3,5)=  47 
  hexMapGMSH(0,4,5)=  48 
  hexMapGMSH(5,1,5)=  49 
  hexMapGMSH(5,2,5)=  50 
  hexMapGMSH(5,3,5)=  51 
  hexMapGMSH(5,4,5)=  52 
  hexMapGMSH(4,5,5)=  53 
  hexMapGMSH(3,5,5)=  54 
  hexMapGMSH(2,5,5)=  55 
  hexMapGMSH(1,5,5)=  56 
  hexMapGMSH(1,1,0)=  57 
  hexMapGMSH(1,4,0)=  58 
  hexMapGMSH(4,4,0)=  59 
  hexMapGMSH(4,1,0)=  60 
  hexMapGMSH(1,2,0)=  61 
  hexMapGMSH(1,3,0)=  62 
  hexMapGMSH(2,4,0)=  63 
  hexMapGMSH(3,4,0)=  64 
  hexMapGMSH(4,3,0)=  65 
  hexMapGMSH(4,2,0)=  66 
  hexMapGMSH(3,1,0)=  67 
  hexMapGMSH(2,1,0)=  68 
  hexMapGMSH(2,2,0)=  69 
  hexMapGMSH(2,3,0)=  70 
  hexMapGMSH(3,3,0)=  71 
  hexMapGMSH(3,2,0)=  72 
  hexMapGMSH(1,0,1)=  73 
  hexMapGMSH(4,0,1)=  74 
  hexMapGMSH(4,0,4)=  75 
  hexMapGMSH(1,0,4)=  76 
  hexMapGMSH(2,0,1)=  77 
  hexMapGMSH(3,0,1)=  78 
  hexMapGMSH(4,0,2)=  79 
  hexMapGMSH(4,0,3)=  80 
  hexMapGMSH(3,0,4)=  81 
  hexMapGMSH(2,0,4)=  82 
  hexMapGMSH(1,0,3)=  83 
  hexMapGMSH(1,0,2)=  84 
  hexMapGMSH(2,0,2)=  85 
  hexMapGMSH(3,0,2)=  86 
  hexMapGMSH(3,0,3)=  87 
  hexMapGMSH(2,0,3)=  88 
  hexMapGMSH(0,1,1)=  89 
  hexMapGMSH(0,1,4)=  90 
  hexMapGMSH(0,4,4)=  91 
  hexMapGMSH(0,4,1)=  92 
  hexMapGMSH(0,1,2)=  93 
  hexMapGMSH(0,1,3)=  94 
  hexMapGMSH(0,2,4)=  95 
  hexMapGMSH(0,3,4)=  96 
  hexMapGMSH(0,4,3)=  97 
  hexMapGMSH(0,4,2)=  98 
  hexMapGMSH(0,3,1)=  99 
  hexMapGMSH(0,2,1)= 100
  hexMapGMSH(0,2,2)= 101
  hexMapGMSH(0,2,3)= 102
  hexMapGMSH(0,3,3)= 103
  hexMapGMSH(0,3,2)= 104
  hexMapGMSH(5,1,1)= 105
  hexMapGMSH(5,4,1)= 106
  hexMapGMSH(5,4,4)= 107
  hexMapGMSH(5,1,4)= 108
  hexMapGMSH(5,2,1)= 109
  hexMapGMSH(5,3,1)= 110
  hexMapGMSH(5,4,2)= 111
  hexMapGMSH(5,4,3)= 112
  hexMapGMSH(5,3,4)= 113
  hexMapGMSH(5,2,4)= 114
  hexMapGMSH(5,1,3)= 115
  hexMapGMSH(5,1,2)= 116
  hexMapGMSH(5,2,2)= 117
  hexMapGMSH(5,3,2)= 118
  hexMapGMSH(5,3,3)= 119
  hexMapGMSH(5,2,3)= 120
  hexMapGMSH(4,5,1)= 121
  hexMapGMSH(1,5,1)= 122
  hexMapGMSH(1,5,4)= 123
  hexMapGMSH(4,5,4)= 124
  hexMapGMSH(3,5,1)= 125
  hexMapGMSH(2,5,1)= 126
  hexMapGMSH(1,5,2)= 127
  hexMapGMSH(1,5,3)= 128
  hexMapGMSH(2,5,4)= 129
  hexMapGMSH(3,5,4)= 130
  hexMapGMSH(4,5,3)= 131
  hexMapGMSH(4,5,2)= 132
  hexMapGMSH(3,5,2)= 133
  hexMapGMSH(2,5,2)= 134
  hexMapGMSH(2,5,3)= 135
  hexMapGMSH(3,5,3)= 136
  hexMapGMSH(1,1,5)= 137
  hexMapGMSH(4,1,5)= 138
  hexMapGMSH(4,4,5)= 139
  hexMapGMSH(1,4,5)= 140
  hexMapGMSH(2,1,5)= 141
  hexMapGMSH(3,1,5)= 142
  hexMapGMSH(4,2,5)= 143
  hexMapGMSH(4,3,5)= 144
  hexMapGMSH(3,4,5)= 145
  hexMapGMSH(2,4,5)= 146
  hexMapGMSH(1,3,5)= 147
  hexMapGMSH(1,2,5)= 148
  hexMapGMSH(2,2,5)= 149
  hexMapGMSH(3,2,5)= 150
  hexMapGMSH(3,3,5)= 151
  hexMapGMSH(2,3,5)= 152
  hexMapGMSH(1,1,1)= 153
  hexMapGMSH(4,1,1)= 154
  hexMapGMSH(4,4,1)= 155
  hexMapGMSH(1,4,1)= 156
  hexMapGMSH(1,1,4)= 157
  hexMapGMSH(4,1,4)= 158
  hexMapGMSH(4,4,4)= 159
  hexMapGMSH(1,4,4)= 160
  hexMapGMSH(2,1,1)= 161
  hexMapGMSH(3,1,1)= 162
  hexMapGMSH(1,2,1)= 163
  hexMapGMSH(1,3,1)= 164
  hexMapGMSH(1,1,2)= 165
  hexMapGMSH(1,1,3)= 166
  hexMapGMSH(4,2,1)= 167
  hexMapGMSH(4,3,1)= 168
  hexMapGMSH(4,1,2)= 169
  hexMapGMSH(4,1,3)= 170
  hexMapGMSH(3,4,1)= 171
  hexMapGMSH(2,4,1)= 172
  hexMapGMSH(4,4,2)= 173
  hexMapGMSH(4,4,3)= 174
  hexMapGMSH(1,4,2)= 175
  hexMapGMSH(1,4,3)= 176
  hexMapGMSH(2,1,4)= 177
  hexMapGMSH(3,1,4)= 178
  hexMapGMSH(1,2,4)= 179
  hexMapGMSH(1,3,4)= 180
  hexMapGMSH(4,2,4)= 181
  hexMapGMSH(4,3,4)= 182
  hexMapGMSH(3,4,4)= 183
  hexMapGMSH(2,4,4)= 184
  hexMapGMSH(2,2,1)= 185
  hexMapGMSH(2,3,1)= 186
  hexMapGMSH(3,3,1)= 187
  hexMapGMSH(3,2,1)= 188
  hexMapGMSH(2,1,2)= 189
  hexMapGMSH(3,1,2)= 190
  hexMapGMSH(3,1,3)= 191
  hexMapGMSH(2,1,3)= 192
  hexMapGMSH(1,2,2)= 193
  hexMapGMSH(1,2,3)= 194
  hexMapGMSH(1,3,3)= 195
  hexMapGMSH(1,3,2)= 196
  hexMapGMSH(4,2,2)= 197
  hexMapGMSH(4,3,2)= 198
  hexMapGMSH(4,3,3)= 199
  hexMapGMSH(4,2,3)= 200
  hexMapGMSH(3,4,2)= 201
  hexMapGMSH(2,4,2)= 202
  hexMapGMSH(2,4,3)= 203
  hexMapGMSH(3,4,3)= 204
  hexMapGMSH(2,2,4)= 205
  hexMapGMSH(3,2,4)= 206
  hexMapGMSH(3,3,4)= 207
  hexMapGMSH(2,3,4)= 208
  hexMapGMSH(2,2,2)= 209
  hexMapGMSH(3,2,2)= 210
  hexMapGMSH(3,3,2)= 211
  hexMapGMSH(2,3,2)= 212
  hexMapGMSH(2,2,3)= 213
  hexMapGMSH(3,2,3)= 214
  hexMapGMSH(3,3,3)= 215
  hexMapGMSH(2,3,3)= 216
CASE DEFAULT
  CALL abort(__STAMP__,&
             'Elements of specified or higher order are not implemented yet. Order: ',bOrd,999.)
END SELECT
END SUBROUTINE getGMSHVolumeMapping


END MODULE MOD_Readin_GMSH_Vars
