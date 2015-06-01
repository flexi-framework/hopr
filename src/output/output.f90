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
OPEN(UNIT=UNIT_errOut,  &
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
DebugVisu  =GETLOGICAL('DebugVisu','.FALSE.')
IF(DebugVisu) THEN
  DebugVisuLevel=GETINT('DebugVisuLevel','0')  
  OutputFormat  =GETINT('OutputFormat','0')  
END IF

sfc_type  =GETSTR('sfc_type','hilbert')
dosortIJK=GETLOGICAL('doSortIJK','.FALSE.')
useSpaceFillingCurve=GETLOGICAL('useSpaceFillingCurve','.TRUE.')

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
  CALL WriteDataToVTK(dim1,nVal,Nplot,nElems,VarNames,Coord,Values,strOutputFile)
CASE(1) ! Tecplot ASCII Output
  strOutputFile=TRIM(Filename)//'.dat'
  CALL WriteDataToASCIITecplot(dim1,nVal,Nplot,nElems,VarNames,Coord,Values,strOutputFile)
CASE(2) ! CGNS output format
  strOutputFile=TRIM(FileName)//'.cgns'
  CALL WriteDataToCGNS(dim1,nVal,Nplot,nElems,VarNames,Coord,Values,TRIM(strOutputFile))
END SELECT
END SUBROUTINE Visualize


END MODULE MOD_Output
