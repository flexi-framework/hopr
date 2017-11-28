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
#define SWRITE WRITE
MODULE MOD_ReadInTools
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE ISO_VARYING_STRING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE

PUBLIC::TRYREAD
PUBLIC::GETSTR
PUBLIC::CNTSTR
PUBLIC::GETINT
PUBLIC::GETREAL
PUBLIC::GETLOGICAL
PUBLIC::GETINTARRAY
PUBLIC::GETREALARRAY

PUBLIC::IgnoredStrings

!===================================================================================================================================

INTERFACE TRYREAD
  MODULE PROCEDURE TRYREAD
END INTERFACE

INTERFACE GETSTR
  MODULE PROCEDURE GETSTR
END INTERFACE

INTERFACE CNTSTR
  MODULE PROCEDURE CNTSTR
END INTERFACE

INTERFACE GETINT
  MODULE PROCEDURE GETINT
END INTERFACE

INTERFACE GETREAL
  MODULE PROCEDURE GETREAL
END INTERFACE

INTERFACE GETLOGICAL
  MODULE PROCEDURE GETLOGICAL
END INTERFACE

INTERFACE GETINTARRAY
  MODULE PROCEDURE GETINTARRAY
END INTERFACE

INTERFACE GETREALARRAY
  MODULE PROCEDURE GETREALARRAY
END INTERFACE

INTERFACE IgnoredStrings
  MODULE PROCEDURE IgnoredStrings
END INTERFACE

INTERFACE FillStrings
  MODULE PROCEDURE FillStrings
END INTERFACE

INTERFACE FindStr
  MODULE PROCEDURE FindStr
END INTERFACE

INTERFACE LowCase
  MODULE PROCEDURE LowCase
END INTERFACE

INTERFACE GetNewString
  MODULE PROCEDURE GetNewString
END INTERFACE

INTERFACE DeleteString
  MODULE PROCEDURE DeleteString
END INTERFACE

TYPE tString
  TYPE(Varying_String)::Str
  TYPE(tString),POINTER::NextStr,PrevStr
END TYPE tString

LOGICAL,PUBLIC::ReadInDone=.FALSE.
TYPE(tString),POINTER::FirstString

CONTAINS

FUNCTION TRYREAD(UnitLoc,Key,abortOpt)
!===================================================================================================================================
! Read string from specified unit
!===================================================================================================================================
! MODULES
USE MOD_Globals,ONLY:abort
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
INTEGER,INTENT(IN)                   :: unitLoc
LOGICAL,INTENT(IN),OPTIONAL          :: abortOpt
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: TRYREAD
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                              :: stat
CHARACTER(LEN=255)                   :: tmp
LOGICAL                              :: abortLoc=.TRUE.
!===================================================================================================================================
IF(PRESENT(abortOpt)) abortLoc=abortOpt
TRYREAD=.TRUE.
READ(unitLoc,*,IOSTAT=stat) tmp
IF(stat.NE.0)              TRYREAD=.FALSE.
IF(TRIM(Key).NE.TRIM(tmp)) TRYREAD=.FALSE.

IF(.NOT.TRYREAD.AND.abortLoc)&
  CALL abort(__STAMP__,&
             'Keyword '//TRIM(Key)//' not found in file.')
END FUNCTION TRYREAD


FUNCTION GETSTR(Key,Proposal)
!===================================================================================================================================
! Read string named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=255)                   :: GetStr   ! String read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=8)                     :: DefMsg  ! ?
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,GetStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,GetStr,DefMsg)
END IF
SWRITE(UNIT_StdOut,'(a3,a30,a3,a33,a3,a7,a3)')' | ',TRIM(Key),' | ', TRIM(GetStr),' | ',TRIM(DefMsg),' | '
END FUNCTION GETSTR

FUNCTION CNTSTR(Key,Proposal)
!===================================================================================================================================
! Counts all occurances of string named "key" from inifile and store in "GETSTR". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Inifile was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                              :: CntStr   ! Number of parameters named "Key" in inifile
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: IntProposal  ! ?
CHARACTER(LEN=LEN(Key))              :: TmpKey  ! ?
TYPE(tString),POINTER                :: Str1  ! ?
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings

CntStr=0
CALL LowCase(Key,TmpKey)
! Remove blanks
TmpKey=REPLACE(TmpKey," ","",Every=.TRUE.)

! Search
Str1=>FirstString
DO WHILE (ASSOCIATED(Str1))
  IF (INDEX(CHAR(Str1%Str),TRIM(TmpKey)//'=').EQ.1) CntStr=CntStr+1
  ! Next string in list
  Str1=>Str1%NextStr
END DO
IF (CntStr.EQ.0) THEN
  IF (PRESENT(Proposal)) THEN
    READ(Proposal,'(I3)')IntProposal
    CntStr=IntProposal
  ELSE
    SWRITE(UNIT_StdOut,*) 'Inifile missing necessary keyword item : ',TRIM(TmpKey)
    CALL abort(__STAMP__, &
         'Code stopped during inifile parsing!')
  END IF
END IF
END FUNCTION CNTSTR

FUNCTION GETINT(Key,Proposal)
!===================================================================================================================================
! Read integer named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                              :: GetInt  ! Integer read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)                   :: HelpStr  ! ?
CHARACTER(LEN=8)                     :: DefMsg  ! ?
INTEGER                              :: ioerr
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*,IOSTAT=ioerr)GetInt
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (integer):'
  WRITE(*,*) TRIM(key),' = ',TRIM(helpStr)
  STOP     
END IF
SWRITE(UNIT_StdOut,'(a3,a30,a3,i33,a3,a7,a3)')' | ',TRIM(Key),' | ', GetInt,' | ',TRIM(DefMsg),' | '
END FUNCTION GETINT


FUNCTION GETREAL(Key,Proposal)
!===================================================================================================================================
! Read real named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key      ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                                 :: GetReal  ! Real read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=500)                   :: HelpStr  ! ?
CHARACTER(LEN=8)                     :: DefMsg  ! ?
INTEGER                              :: ioerr
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
! Find values of pi in the string
CALL getPImultiplies(helpstr)
READ(HelpStr,*,IOSTAT=ioerr)GetReal
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (real):'
  WRITE(*,*) TRIM(key),' = ',TRIM(helpStr)
  STOP     
END IF
SWRITE(UNIT_StdOut,'(a3,a30,a3,e33.5,a3,a7,a3)')' | ',TRIM(Key),' | ', GetReal,' | ',TRIM(DefMsg),' | '
END FUNCTION GETREAL


FUNCTION GETLOGICAL(Key,Proposal)
!===================================================================================================================================
! Read logical named "key" from setup file and store in "GETINT". If keyword "Key" is not found in ini file,
! the default value "Proposal" is used for "GETINT" (error if "Proposal" not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key        ! Search for this keyword in ini file
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal   ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                              :: GetLogical ! Logical read from setup file or initialized with default value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)                   :: HelpStr  ! ?
CHARACTER(LEN=8)                     :: DefMsg  ! ?
INTEGER                              :: ioerr
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*,IOSTAT=ioerr)GetLogical
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (logical):'
  WRITE(*,*) TRIM(key),' = ',TRIM(helpStr)
  STOP     
END IF
SWRITE(UNIT_StdOut,'(a3,a30,a3,l33,a3,a7,a3)')' | ',TRIM(Key),' | ', GetLogical,' | ',TRIM(DefMsg),' | '
END FUNCTION GETLOGICAL


FUNCTION GETINTARRAY(Key,nIntegers,Proposal)
!===================================================================================================================================
! Read array of "nIntegers" integer values named "Key" from ini file. If keyword "Key" is not found in setup file, the default
! values "Proposal" are used to create the array (error if "Proposal" not given). Setup file was read in before and is stored as
! list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
    IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key              ! Search for this keyword in ini file
INTEGER,INTENT(IN)                   :: nIntegers        ! Number of values in array
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal         ! Default values as character string (as in setup file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                   :: GetIntArray(nIntegers)      ! Integer array read from setup file or initialized with default values
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)        :: HelpStr  ! ?
CHARACTER(LEN=8)          :: DefMsg  ! ?
INTEGER                   :: iInteger  ! ?
INTEGER                   :: ioerr
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
READ(HelpStr,*,IOSTAT=ioerr)GetIntArray
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (integer array):'
  WRITE(*,*) TRIM(key),' = ',TRIM(helpStr)
  STOP     
END IF
SWRITE(UNIT_stdOut,'(a3,a30,a3,a28,i4,a4,a7,a3)',ADVANCE='NO') ' | ',TRIM(Key),' | ',&
                                                               'Integer array of size (',nIntegers,') | ',TRIM(DefMsg),' | '
DO iInteger=0,nIntegers-1
  IF ((iInteger.GT.0) .AND. (MOD(iInteger,8).EQ.0)) THEN
    SWRITE(UNIT_stdOut,*)
    SWRITE(UNIT_stdOut,'(a80,a3)',ADVANCE='NO')'',' | '
  END IF
  SWRITE(UNIT_stdOut,'(i5)',ADVANCE='NO')GetIntArray(iInteger+1)
END DO
SWRITE(UNIT_stdOut,*)
END FUNCTION GETINTARRAY


FUNCTION GETREALARRAY(Key,nReals,Proposal)
!===================================================================================================================================
! Read array of "nReals" real values named "Key" from ini file. If keyword "Key" is not found in setup file, the default
! values "Proposal" are used to create the array (error if "Proposal" not given). Setup file was read in before and is stored as
! list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key              ! Search for this keyword in ini file
INTEGER,INTENT(IN)                   :: nReals           ! Number of values in array
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal         ! Default values as character string (as in setup file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                      :: GetRealArray(nReals)        ! Real array read from setup file or initialized with default values
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255+nReals*50) :: HelpStr  ! ?
CHARACTER(LEN=8)          :: DefMsg  ! ?
INTEGER                   :: iReal  ! ?
INTEGER                   :: ioerr  ! ?
!===================================================================================================================================
! Read-in ini file if not done already
CALL FillStrings

IF (PRESENT(Proposal)) THEN
  CALL FindStr(Key,HelpStr,DefMsg,Proposal)
ELSE
  CALL FindStr(Key,HelpStr,DefMsg)
END IF
CALL getPImultiplies(helpstr)
READ(HelpStr,*,IOSTAT=ioerr)GetRealArray
IF(ioerr.NE.0)THEN
  WRITE(*,*)'PROBLEM IN READIN OF LINE (RealArray):'
  WRITE(*,*) TRIM(key),' = ',TRIM(helpStr)
  STOP     
END IF
SWRITE(UNIT_stdOut,'(a3,a30,a3,a28,i4,a4,a7,a3)',ADVANCE='NO') ' | ',TRIM(Key),' | ',&
                                                               'Real array of size (',nReals,') | ',TRIM(DefMsg),' | '
DO iReal=0,nReals-1
  IF ((iReal.GT.0) .AND. (MOD(iReal,8).EQ.0)) THEN
    SWRITE(UNIT_stdOut,*)
    SWRITE(UNIT_stdOut,'(a80,a3)',ADVANCE='NO')'',' | '
  END IF
  SWRITE(UNIT_stdOut,'(f5.2)',ADVANCE='NO')GetRealArray(iReal+1)
END DO
SWRITE(UNIT_stdOut,*)
END FUNCTION GETREALARRAY


SUBROUTINE IgnoredStrings()
!===================================================================================================================================
! Prints out remaining strings in list after read-in is complete
!===================================================================================================================================
! MODULES
USE ISO_VARYING_STRING
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tString),POINTER                  :: Str1  ! ?
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)')" THE FOLLOWING INI-FILE PARAMETERS WERE IGNORED:"
Str1=>FirstString
DO WHILE(ASSOCIATED(Str1))
  SWRITE(UNIT_stdOut,'(A4,A)')" |- ",TRIM(CHAR(Str1%Str))
  Str1=>Str1%NextStr
END DO
SWRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE IgnoredStrings


SUBROUTINE FillStrings(IniFile)
!===================================================================================================================================
! Read ini file and put each line in a string object. All string objects are connected to a list of string objects starting
! with "firstString"
!===================================================================================================================================
! MODULES
USE ISO_VARYING_STRING
USE,INTRINSIC :: ISO_FORTRAN_ENV,ONLY:IOSTAT_END
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN),OPTIONAL   :: IniFile                    ! Name of ini file to be read in
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tString),POINTER                  :: Str1=>NULL(),Str2=>NULL()  ! ?
CHARACTER(LEN=255)                     :: HelpStr,Str  ! ?
CHARACTER(LEN=300)                     :: File  ! ?
TYPE(Varying_String)                   :: aStr,bStr,Separator  ! ?
INTEGER                                :: EOF  ! ?
!===================================================================================================================================
! Check if we have read in ini file already
IF (ReadInDone) RETURN
! Get name of ini file
IF (PRESENT(IniFile)) THEN
  File = TRIM(IniFile)
ELSE
  IF(COMMAND_ARGUMENT_COUNT().LT.1) STOP 'Parameter file not specified! Usage: "hopr <parameter.ini>"'
  CALL GET_COMMAND_ARGUMENT(1,File)
END IF
SWRITE(UNIT_StdOut,*)'| Reading from file "',TRIM(File),'":'

OPEN(UNIT   = 103,        &
     FILE   = File,       &
     STATUS = 'OLD',      &
     ACTION = 'READ',     &
     ACCESS = 'SEQUENTIAL')
EOF=0

NULLIFY(Str1,Str2)
DO WHILE(EOF.NE.IOSTAT_END)
  IF(.NOT.ASSOCIATED(Str1)) CALL GetNewString(Str1)
    ! Read line from file
    CALL Get(103,aStr,iostat=EOF)
    Str=aStr
!IPWRITE(*,*)'Reading: ',Str,EOF
    IF (EOF.NE.IOSTAT_END) THEN
      ! Remove comments with "!"
      CALL Split(aStr,Str1%Str,"!")
      ! Remove comments with "#"
      CALL Split(Str1%Str,bStr,"#")
      Str1%Str=bStr
      ! Remove "%" sign from old ini files, i.e. mesh% disc% etc.
      CALL Split(Str1%Str,bStr,"%",Separator,Back=.false.)
      ! If we have a newtype ini file, take the other part
      IF(LEN(CHAR(Separator)).EQ.0) Str1%Str=bStr
      ! Remove blanks
      Str1%Str=Replace(Str1%Str," ","",Every=.true.)
      ! Replace brackets
      Str1%Str=Replace(Str1%Str,"(/"," ",Every=.true.)
      Str1%Str=Replace(Str1%Str,"/)"," ",Every=.true.)
      ! Replace commas
      Str1%Str=Replace(Str1%Str,","," ",Every=.true.)
      ! Lower case
      CALL LowCase(CHAR(Str1%Str),HelpStr)
      ! If we have a remainder (no comment only)
      IF(LEN_TRIM(HelpStr).GT.2) THEN
        Str1%Str=Var_Str(HelpStr)
        IF(.NOT.ASSOCIATED(Str2)) THEN
          FirstString=>Str1
        ELSE
          Str2%NextStr=>Str1
          Str1%PrevStr=>Str2
        END IF
        Str2=>Str1
        CALL GetNewString(Str1)
      END IF
    END IF
END DO
CLOSE(103)
ReadInDone=.TRUE.

CALL UserDefinedVars()


END SUBROUTINE FillStrings

SUBROUTINE UserDefinedVars()
!===================================================================================================================================
! Get the user defined variables 
!===================================================================================================================================
! MODULES
USE iso_varying_string
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                          :: i,j,nDefVars
TYPE(Varying_String),ALLOCATABLE :: DefVar(:,:)
TYPE(tString),POINTER            :: Str1  ! ?
TYPE(Varying_String)             :: vStr,vStr1,vStr2
LOGICAL                          :: found
!===================================================================================================================================
nDefVars=CNTSTR('defvar','0')
IF(nDefVars.EQ.0) RETURN
SWRITE(UNIT_StdOut,'(A,I4,A)')' | Found ',nDefVars,' UserDefined variables: '
ALLOCATE(DefVar(2,nDefVars))
DO i=1,nDefVars
  CALL GetDefVar(DefVar(:,i))
  !check if part of the variable name was used before
  DO j=1,i-1
    IF (INDEX(TRIM(CHAR(DefVar(1,i))),TRIM(CHAR(DefVar(1,j)))).NE.0) THEN
      SWRITE(UNIT_StdOut,*) '!! WARNING !! Problem with DEFVAR ', TRIM(CHAR(DefVar(1,i)))
      SWRITE(UNIT_StdOut,*) '  a part of this variable name was already used in DEFVAR ' ,TRIM(CHAR(DefVar(1,j)))
      CALL abort(__STAMP__, &
         'DEFVAR: do not reuse same strings for variable names! Code stopped during inifile parsing!')
    END IF
  END DO
  Str1=>FirstString
  DO WHILE(ASSOCIATED(Str1))
    vStr=Str1%Str
    vStr2=Str1%Str
    CALL Split(vStr2,vStr1,"=",back=.FALSE.) 
    found=.FALSE.
    IF (INDEX(TRIM(CHAR(vStr2)),TRIM(CHAR(DefVar(1,i)))).NE.0) THEN
      found=.TRUE.
      vStr2=replace(vStr2,TRIM(CHAR(DefVar(1,i))),TRIM(CHAR(DefVar(2,i))),Every=.TRUE.)
    END IF
    IF(Found)THEN
      !SWRITE(UNIT_StdOut,*)'DEBUG, ',TRIM(CHAR(Str1%str))
      Str1%Str=CHAR(vStr1)//'='//CHAR(vStr2)
      !SWRITE(UNIT_StdOut,*)' >>>>>>',TRIM(CHAR(Str1%str))
    END IF !found
    ! Next string in list
    Str1=>Str1%NextStr
  END DO !WHILE Str1 associated
END DO !i=1,nDefVars

END SUBROUTINE UserDefinedVars 


SUBROUTINE GetDefVar(DefVar)
!===================================================================================================================================
! Get the user defined variables 
!===================================================================================================================================
! MODULES
USE iso_varying_string
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(Varying_String):: DefVar(2)   !Name, Value
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=255)  :: HelpStr  ! ?
CHARACTER(LEN=255)  :: aStr  ! ?
CHARACTER(LEN=8)    :: DefMsg  ! ?
TYPE(Varying_String):: vStr,vStr1,vStr2,vStrTmp,vStr_narr
INTEGER             :: DefVarInt
REAL                :: DefVarReal
INTEGER,ALLOCATABLE :: DefVarIntArr(:)
REAL,ALLOCATABLE    :: DefVarRealArr(:)
INTEGER             :: nArr
LOGICAL             :: DefVarIsInt
LOGICAL             :: DefVarIsIntarray
LOGICAL             :: DefVarIsReal
LOGICAL             :: DefVarIsRealarray
!===================================================================================================================================
CALL FindStr('DEFVAR',HelpStr,DefMsg)
CALL LowCase(HelpStr,aStr)
vStr=aStr

CALL Split(vStr,vStr1,":",back=.FALSE.) !split after first occurence in vStr, and put first part in vStr1, second part in vStr
CALL Split(vStr,vStr2,":",back=.TRUE.)  !split after last occurence in Vstr and put second part in vStr1 and first part in Vstr

vStrTmp=vStr1
CALL Split(vStrtmp,vStr1,"~",back=.FALSE.) !first part in vStr1, second part in VStrtmp
CALL Split(vStrtmp,vStr_narr,"~",back=.TRUE.) !second part in VStr_narr, first part in VStrtmp
CALL Split(vStr_narr,vStrTmp,")",back=.TRUE.) !first part VStr_narr

DefVarIsInt      =(CHAR(vStr1).EQ.'(int)') 
DefVarIsIntarray =(CHAR(vStr1).EQ.'(int') !array  must be of format (int~n)
DefVarIsReal     =(CHAR(vStr1).EQ.'(real)') 
DefVarIsRealarray=(CHAR(vStr1).EQ.'(real') 


IF(.NOT.((DefVarIsInt).OR.(DefVarIsIntArray).OR.(DefVarIsReal).OR.(defVarIsRealarray) ))THEN
  SWRITE(UNIT_StdOut,*) 'DEFVAR not correctly defined: ',TRIM(HelpStr)
    CALL abort(__STAMP__, &
         'Code stopped during inifile parsing!')
END IF

IF(DefVarIsIntArray.OR.DefVarIsRealArray)THEN
  aStr=CHAR(vstr_narr)
  READ(aStr,*) nArr
END IF

!now take the second part of the definition ( nvar = xxx)
vStr=VStr2
CALL Split(vStr,vStr1,"=",back=.FALSE.) 
CALL Split(vStr,vStr2,"=",back=.TRUE.) 
DefVar(1)=vStr1
DefVar(2)=vStr2
aStr=CHAR(DefVar(2))
IF(DefVarIsInt) THEN
  READ(aStr,*)DefVarInt
  SWRITE(UNIT_StdOut,'(A3,A30,A3,I33,A3,A7,A3)')  ' | ',TRIM(CHAR(DefVar(1))),' | ', DefVarInt      ,' | ','=>INT  ',' | '
ELSE IF(DefVarIsReal)THEN
  READ(aStr,*)DefVarReal
  SWRITE(UNIT_StdOut,'(A3,A30,A3,E33.5,A3,A7,A3)')' | ',TRIM(CHAR(DefVar(1))),' | ', DefVarReal     ,' | ','=>REAL ',' | '
ELSE IF(DefVarIsIntArray)THEN
  ALLOCATE(DefVarIntArr(nArr))
  READ(aStr,*)DefVarIntArr(:)
  SWRITE(UNIT_StdOut,'(A3,A30,A31,I4,A4,A7,A3,'//TRIM(CHAR(vStr_narr))//'(X,I4))')' | ',TRIM(CHAR(DefVar(1))), &
        ' |      Integer array of size (', nArr,') | ','=>INT  ',' | ',DefVarIntArr(1:narr)
  WRITE(aStr,'('//TRIM(CHAR(vStr_narr))//'(X,I8))') DefVarIntArr !overwrite
  DEALLOCATE(DefVarIntArr)
  DefVar(2)=aStr
ELSE IF(DefVarIsRealArray)THEN
  ALLOCATE(DefVarRealArr(nArr))
  READ(aStr,*)DefVarRealArr(:)
  SWRITE(UNIT_StdOut,'(A3,A30,A31,I4,A4,A7,A3,'//TRIM(CHAR(vStr_narr))//'(X,E8.1))')' | ',TRIM(CHAR(DefVar(1))), &
        ' |         Real array of size (', nArr,') | ','=>REAL ',' | ',DefVarRealArr(1:narr)
  WRITE(aStr,'('//TRIM(CHAR(vStr_narr))//'(X,E23.15))') DefVarRealArr !overwrite
  DEALLOCATE(DefVarRealArr)
  DefVar(2)=aStr
END IF

END SUBROUTINE GetDefVar 


SUBROUTINE GetNewString(Str)
!===================================================================================================================================
! Create and initialize new string object.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
TYPE(tString),POINTER,INTENT(INOUT) :: Str ! New string
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
NULLIFY(Str)
ALLOCATE(Str)
NULLIFY(Str%NextStr,Str%PrevStr)
END SUBROUTINE GetNewString


SUBROUTINE DeleteString(Str)
!===================================================================================================================================
! Remove string "Str" from list of strings witFirstString,h first element "DirstString" and delete string.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tString),POINTER,INTENT(INOUT) :: Str         ! String to delete
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF (ASSOCIATED(Str%NextStr)) Str%NextStr%PrevStr=>Str%PrevStr
IF (ASSOCIATED(Str,FirstString)) THEN
  FirstString=>Str%NextStr
ELSE
  Str%PrevStr%NextStr=>Str%NextStr
END IF
DEALLOCATE(Str)
NULLIFY(Str)
END SUBROUTINE DeleteString


SUBROUTINE FindStr(Key,Str,DefMsg,Proposal)
!===================================================================================================================================
! Find parameter string containing keyword "Key" in list of strings starting with "FirstString" and return string "Str" without
! keyword. If keyword is not found in list of strings, return default values "Proposal" (error if not given).
! Ini file was read in before and is stored as list of character strings starting with "FirstString".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)          :: Key         ! Search for this keyword in ini file
CHARACTER(LEN=8),INTENT(INOUT)       :: DefMsg      ! Default message = keyword not found, return default parameters (if available)
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proposal    ! Default values as character string (as in ini file)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(OUT)         :: Str         ! Parameter string without keyword
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=LEN(Key))              :: TmpKey   ! ?
TYPE(tString),POINTER                :: Str1  ! ?
LOGICAL                              :: Found  ! ?
!===================================================================================================================================
DefMsg='*CUSTOM'
! Convert to lower case
CALL LowCase(Key,TmpKey)
! Remove blanks
TmpKey=REPLACE(TmpKey," ","",Every=.TRUE.)
Found=.FALSE.
Str1=>FirstString
DO WHILE(.NOT.Found)
  IF (.NOT.ASSOCIATED(Str1)) THEN
    IF (.NOT.PRESENT(Proposal)) THEN
      SWRITE(UNIT_StdOut,*) 'Inifile missing necessary keyword item : ',TRIM(TmpKey)
      CALL abort(__STAMP__, &
           'Code stopped during inifile parsing!')
    ELSE ! Return default value
!      CALL LowCase(TRIM(Proposal),Str)
      Str=TRIM(Proposal)
      IF (Str(1:1).NE.'@') THEN
        DefMsg='DEFAULT'
      END IF
      RETURN
    END IF ! (.NOT.PRESENT(Proposal))
  END IF ! (.NOT.ASSOCIATED(Str1))

  IF (INDEX(CHAR(Str1%Str),TRIM(TmpKey)//'=').EQ.1) THEN
    Found=.TRUE.
    Str1%Str=replace(Str1%Str,TRIM(TmpKey)//'=',"",Every=.TRUE.)
    Str=TRIM(CHAR(Str1%Str))
    ! Remove string from list
    CALL DeleteString(Str1)
  ELSE
    ! Next string in list
    Str1=>Str1%NextStr
  END IF

END DO
END SUBROUTINE FindStr


SUBROUTINE LowCase(Str1,Str2)
!===================================================================================================================================
! Transform upper case letters in "Str1" into lower case letters, result is "Str2"
!==================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: Str1 ! Input string 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(OUT) :: Str2 ! Output string, lower case letters only
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                      :: iLen,nLen,Upper  ! ?
CHARACTER(LEN=*),PARAMETER   :: lc='abcdefghijklmnopqrstuvwxyz'  ! ?
CHARACTER(LEN=*),PARAMETER   :: UC='ABCDEFGHIJKLMNOPQRSTUVWXYZ'  ! ?
LOGICAL                      :: HasEq  ! ?
!===================================================================================================================================
HasEq=.FALSE.
Str2=Str1
nLen=LEN_TRIM(Str1)
DO iLen=1,nLen
  ! Transformation stops at "="
  IF(Str1(iLen:iLen).EQ.'=') HasEq=.TRUE.
  Upper=INDEX(UC,Str1(iLen:iLen))
  IF ((Upper > 0).AND. .NOT. HasEq) Str2(iLen:iLen) = lc(Upper:Upper)
END DO
END SUBROUTINE LowCase

SUBROUTINE getPImultiplies(helpstr)
!===================================================================================================================================
! Searches for the occurence of 'PI','Pi','pi' and 'pI' in a helpstr and replaces
! it with the value of pi=3.1415... etc. and oes a multiplication.
!===================================================================================================================================
! MODULES
USE iso_varying_string
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(INOUT) :: helpstr   ! Input character string
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(varying_string)      :: separator  ! ?
TYPE(varying_string)      :: astr,bstr,cstr,dstr  ! ?
CHARACTER(LEN=1000)       :: dummystr  ! ?
REAL                      :: PI  ! ?
REAL                      :: adummy  ! ?
LOGICAL                   :: finished  ! ?
!===================================================================================================================================
!Initialiazation
dstr=var_str("")
PI=ACOS(-1.)
finished=.false.
!Replace all occurences of pi in the string by one symbol
helpstr=replace(helpstr,"PI","@",every=.true.)  ! Insert value of Pi
helpstr=replace(helpstr,"pi","@",every=.true.)  ! Insert value of Pi
helpstr=replace(helpstr,"pI","@",every=.true.)  ! Insert value of Pi
helpstr=replace(helpstr,"Pi","@",every=.true.)  ! Insert value of Pi
astr=var_str(helpstr)
! loop over string
DO WHILE(.NOT. finished)
  !split sting at "@-occurences"
  CALL split(astr,bstr,"@",separator,back=.false.) !bStr is string in front of @
  IF(len(char(separator)) .NE. 0)THEN 
    ! we have found something, bnow get the factor in front of @
    CALL split(bstr,cstr," ",separator,back=.true.)
    IF(LEN(char(cstr)) .EQ. 0)THEN
      !no factor
      adummy=1
    ELSE
      !extract factor 
      dummystr=trim(char(cstr))
      READ(dummystr,*)adummy
    ENDIF
    !do the multiplication and recombine the string into "dstr"
    adummy=PI*adummy
    WRITE(dummystr,'(2a,1X,E23.15)')trim(char(dstr)),trim(char(bstr)),adummy
    dstr=var_str(dummystr)
  ELSE
    ! we did not find anything now recombine the remaining string into "dstr"
    WRITE(dummystr,'(2a)')trim(char(dstr)),trim(char(bstr))
    dstr=var_str(dummystr)
    finished=.true.
  END IF
END DO
helpstr=trim(char(dstr))
END SUBROUTINE getPImultiplies

END MODULE MOD_ReadInTools
