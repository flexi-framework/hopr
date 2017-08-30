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
MODULE MOD_Output
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOd_Globals
USE MOD_Output_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitOutput
  MODULE PROCEDURE InitOutput
END INTERFACE

!INTERFACE Visualize
!  MODULE PROCEDURE Visualize
!END INTERFACE


PUBLIC::InitOutput
PUBLIC::Visualize
!===================================================================================================================================

CONTAINS
SUBROUTINE InitOutput()
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_ReadInTools
USE MOD_Output_vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: ioStatus  ! ?
CHARACTER(LEN=255)  :: strFileName  ! ?
CHARACTER(LEN=8)    :: StrDate  ! ?
CHARACTER(LEN=10)   :: StrTime  ! ?
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT OUTPUT...'
ProjectName=GETSTR('projectname')
! Open file for error output
WRITE(strFileName,'(A)')TRIM(ProjectName)//'_ERRORS.out'
OPEN(UNIT=UNIT_errFile, &
     FILE=strFileName,  &
     STATUS='REPLACE',  &
     ACTION='WRITE',    &
     IOSTAT=ioStatus)
Logging=GETLOGICAL('Logging','.FALSE.')
! Open file for logging
IF(Logging)THEN
  WRITE(strFileName,'(A)')TRIM(ProjectName)//'.log'
  OPEN(UNIT=UNIT_logOut,  &
       FILE=strFileName,  &
       STATUS='UNKNOWN',  &
       ACTION='WRITE',    &
       POSITION='APPEND', &
       IOSTAT=ioStatus)
  CALL DATE_AND_TIME(StrDate,StrTime)
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,'(132("#"))')
  WRITE(UNIT_logOut,*)
  WRITE(UNIT_logOut,*)'STARTED LOGGING ON ',StrDate(7:8),'.',StrDate(5:6),'.',StrDate(1:4),' | ',&
                      StrTime(1:2),':',StrTime(3:4),':',StrTime(5:10)
END IF  ! Logging
DebugVisu   =GETLOGICAL('DebugVisu','.FALSE.')
OutputFormat=GETINT('OutputFormat','0')  
IF(DebugVisu) THEN
  DebugVisuLevel=GETINT('DebugVisuLevel','0')  
  IF(DebugVisuLevel.GE.2) &
    Visu_sJ_limit  =GETREAL('Visu_sJ_limit','1.')  
END IF

sfc_type  =GETSTR('sfc_type','hilbert')
dosortIJK=GETLOGICAL('doSortIJK','.FALSE.')
useSpaceFillingCurve=GETLOGICAL('useSpaceFillingCurve','.TRUE.')
sfc_boundbox=GETINT('sfc_boundbox','2')

OutputInitDone=.TRUE.
WRITE(UNIT_stdOut,'(A)')' INIT OUTPUT DONE!'
WRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitOutput


SUBROUTINE Visualize(dim1,nVal,Nplot,nElems,VarNames,Coord,Values,FileName)
!===================================================================================================================================
! Call of software specific output routines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Output_vars   ,ONLY:OutputFormat
USE MOD_Output_VTK    ,ONLY:WriteDataToVTK
USE MOD_Output_tecplot,ONLY:WriteDataToASCIITecplot
USE MOD_Output_CGNS   ,ONLY:WriteDataToCGNS
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: dim1                    ! dimension of the data (either 2 or 3)
INTEGER,INTENT(IN)            :: nVal                    ! Number of nodal output variables
INTEGER,INTENT(IN)            :: NPlot                   ! Number of output points per element: (Nplot+1)**dim1
INTEGER,INTENT(IN)            :: nElems                  ! Number of output elements
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVal)          ! Names of all variables that will be written out
REAL,INTENT(IN)                :: Coord(3,1:(Nplot+1)**dim1,1:nElems)  ! ?
REAL,INTENT(IN)                :: Values(1:nVal,1:(Nplot+1)**dim1,1:nElems)  ! ?
CHARACTER(LEN=255),INTENT(IN)  :: FileName  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: strOutputFile  ! ?
!===================================================================================================================================
SELECT CASE(OutputFormat)
CASE(0) ! VTK Output Format
  strOutputFile=TRIM(FileName)//'.vtu'
  CALL WriteDataToVTK(dim1,3,nVal,Nplot,nElems,VarNames,Coord,Values,strOutputFile)
CASE(1) ! Tecplot ASCII Output
  strOutputFile=TRIM(Filename)//'.dat'
  CALL WriteDataToASCIITecplot(dim1,nVal,Nplot,nElems,VarNames,Coord,Values,strOutputFile)
CASE(2) ! CGNS output format
  strOutputFile=TRIM(FileName)//'.cgns'
  CALL WriteDataToCGNS(dim1,nVal,Nplot,nElems,VarNames,Coord,Values,TRIM(strOutputFile))
END SELECT
END SUBROUTINE Visualize


END MODULE MOD_Output
