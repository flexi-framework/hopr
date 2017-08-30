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
MODULE MOD_Output_CGNS
!===================================================================================================================================
! Module for generic data output in CGNS format
!===================================================================================================================================
! MODULES
USE CGNS
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=5)::ProgramName='HOPR'

INTERFACE WriteDataToCGNS
  MODULE PROCEDURE WriteDataToCGNS
END INTERFACE

PUBLIC::WriteDataToCGNS
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteDataToCGNS(dim1,nVal,NPlot,nElems,VarNames,Coord,Values,FileString)
!===================================================================================================================================
! Subroutine to write 3D point data to CGNS format
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: dim1                    ! dimension of the data (either 2 or 3)
INTEGER,INTENT(IN)            :: nVal                    ! Number of nodal output variables
INTEGER,INTENT(IN)            :: NPlot                   ! Number of output points
INTEGER,INTENT(IN)            :: nElems                  ! Number of output elements
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVal)          ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Coord(3,1:(NPlot+1)**dim1,1:nElems) ! CoordsVector 
REAL,INTENT(IN)               :: Values(1:nVal,1:(NPlot+1)**dim1,1:nElems) ! Statevector
CHARACTER(LEN=*),INTENT(IN)   :: FileString              ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
PP_CGNS_INT_TYPE   :: CGNSFile,CGNSBase,CGNSZone,CGNSCoords,CGNSsection,CGNSFlowSol,CGNSFieldInd ! ?
PP_CGNS_INT_TYPE   :: iSize(1,1:3)  ! ?
PP_CGNS_INT_TYPE   :: iErr  ! ?
INTEGER            :: iVal,i,j,k,iElem  ! ?
INTEGER            :: NPlot_p1,NPlot_p1_2,NPlot_p1_3,NPlot_2  ! ?
INTEGER            :: NodeIDElem  ! ?
PP_CGNS_INT_TYPE   :: ElemConn(1:2**dim1,1:NPlot**dim1,1:nElems)  ! ?
PP_CGNS_INT_TYPE   :: one,zero   ! ?
CHARACTER(LEN=255) :: CGName  ! ?
CHARACTER(LEN=32)  :: VarNames32(nVal)                   ! CGNS uses only 32 characters for (variable-)names
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE DATA TO CGNS FILE... "//TRIM(FileString)
one=1
zero=0

! Open the cgns file
CALL cg_open_f(TRIM(FileString), MODE_WRITE, CGNSFile, iErr)
IF (iErr .NE. CG_OK) CALL my_cg_error_exit('Error Opening CGNS File.',CGNSFile)
! Create base
CALL cg_base_write_f(CGNSFile,ProgramName//'Base',dim1,3,CGNSBase,iErr)
IF (iErr .NE. CG_OK) CALL my_cg_error_exit('Error Creating CGNS Base.',CGNSFile)

! Compute number of elements and nodes 
isize(1,1) = ((NPlot+1)**dim1)*nElems
isize(1,2) = (NPlot**dim1)*nElems
isize(1,3) = 0

! Create new zone in file
CGName=ProgramName//'VisuData'
CALL cg_zone_write_f(CGNSfile,CGNSBase,TRIM(CGname),isize,Unstructured,CGNSZone,iErr)
IF (iErr .NE. CG_OK) CALL my_cg_error_exit('Error Creating CGNS Zone.',CGNSFile)

! Write x-coordinates
CALL cg_coord_write_f(CGNSFile,CGNSBase,CGNSZone,RealDouble,'CoordinateX',             &
                      Coord(1,:,:),                                                    &
                      CGNSCoords, iErr)
IF (iErr .NE. CG_OK) CALL my_cg_error_exit('Error Writing x-Coordinates.',CGNSFile)

! Write y-coordinates
CALL cg_coord_write_f(CGNSFile,CGNSBase,CGNSZone,RealDouble,'CoordinateY',             &
                      Coord(2,:,:),                                                    &
                      CGNSCoords, iErr)
IF (iErr .NE. CG_OK) CALL my_cg_error_exit('Error Writing y-Coordinates.',CGNSFile)

! Write z-coordinates
CALL cg_coord_write_f(CGNSFile,CGNSBase,CGNSZone,RealDouble,'CoordinateZ',             &
                      Coord(3,:,:),                                                    &
                      CGNSCoords, iErr)
IF (iErr .NE. CG_OK) CALL my_cg_error_exit('Error Writing z-Coordinates.',CGNSFile)

! Build element connectivity
NPlot_p1  =(NPlot+1) 
NPlot_p1_2=(NPlot+1)**2
SELECT CASE(dim1)
CASE(2)
  NodeIDElem=0
  DO iElem=1,nElems
    DO j=1,NPlot
      DO i=1,NPlot
        ElemConn(1,i+Nplot*(j-1),iElem) = NodeIDElem+i+  (j-1)*Nplot_p1  !P1(CGNS=tecplot standard)
        ElemConn(2,i+Nplot*(j-1),iElem) = NodeIDElem+i+1+(j-1)*Nplot_p1  !P2
        ElemConn(3,i+Nplot*(j-1),iElem) = NodeIDElem+i+1+ j   *Nplot_p1  !P3     
        ElemConn(4,i+Nplot*(j-1),iElem) = NodeIDElem+i+   j   *Nplot_p1  !P4
      END DO 
    END DO 
    NodeIDElem=NodeIDElem+NPlot_p1_2
  END DO
  ! Write Element Connectivity 
  CALL cg_section_write_f(CGNSFile,CGNSBase,CGNSZone,'Elements',QUAD_4,one,isize(1,2),zero, &
                          ElemConn,CGNSsection,iErr)
CASE(3)
  NPlot_p1_3=(NPlot+1)**3
  NPlot_2=NPlot**2
  
  NodeIDElem=0
  DO iElem=1,nElems
    DO k=1,NPlot
      DO j=1,NPlot
        DO i=1,NPlot
          !visuHexaElem  (CGNS=tecplot standard)
          ElemConn(1,i+Nplot*(j-1)+Nplot_2*(k-1),iElem) = NodeIDElem+i+  (j-1)*NPlot_p1+(k-1)*NPlot_p1_2       !P1
          ElemConn(2,i+Nplot*(j-1)+Nplot_2*(k-1),iElem) = NodeIDElem+i+1+(j-1)*NPlot_p1+(k-1)*NPlot_p1_2       !P2
          ElemConn(3,i+Nplot*(j-1)+Nplot_2*(k-1),iElem) = NodeIDElem+i+1+ j   *NPlot_p1+(k-1)*NPlot_p1_2       !P3     
          ElemConn(4,i+Nplot*(j-1)+Nplot_2*(k-1),iElem) = NodeIDElem+i+   j   *NPlot_p1+(k-1)*NPlot_p1_2       !P4
          ElemConn(5,i+Nplot*(j-1)+Nplot_2*(k-1),iElem) = NodeIDElem+i+  (j-1)*NPlot_p1+ k   *NPlot_p1_2       !P5
          ElemConn(6,i+Nplot*(j-1)+Nplot_2*(k-1),iElem) = NodeIDElem+i+1+(j-1)*NPlot_p1+ k   *NPlot_p1_2       !P6
          ElemConn(7,i+Nplot*(j-1)+Nplot_2*(k-1),iElem) = NodeIDElem+i+1+ j   *NPlot_p1+ k   *NPlot_p1_2       !P7     
          ElemConn(8,i+Nplot*(j-1)+Nplot_2*(k-1),iElem) = NodeIDElem+i+   j   *NPlot_p1+ k   *NPlot_p1_2       !P8
        END DO 
      END DO 
    END DO 
    NodeIDElem=NodeIDElem+NPlot_p1_3
  END DO
  ! Write Element Connectivity 
  CALL cg_section_write_f(CGNSFile,CGNSBase,CGNSZone,'Elements',HEXA_8,one,isize(1,2),zero, &
                          ElemConn,CGNSsection,iErr)
END SELECT


IF (iErr .NE. CG_OK) CALL my_cg_error_exit('Error Writing Connectivity.',CGNSFile)
! Write out point data
CGname='FlowSolution'
CALL cg_sol_write_f(CGNSFile,CGNSBase,CGNSZone,TRIM(CGname),Vertex,CGNSFlowSol,iErr)
IF (iErr .NE. CG_OK) CALL my_cg_error_exit('Error Creating CGNS Flow Solution Node.',CGNSFile)
DO iVal=1,nVal
  VarNames32(iVal)=VarNames(iVal)(1:32)
  CALL cg_field_write_f(CGNSFile,CGNSBase,CGNSZone,CGNSFlowSol,RealDouble,TRIM(VarNames(iVal)),     &
                        Values(iVal,:,:),                                                            &
                        CGNSFieldInd,iErr)
  IF (iErr .NE. CG_OK) CALL my_cg_error_exit('Error Writing CGNS Variable '//TRIM(VarNames(iVal))//'.',CGNSFile)
END DO


! Close CGNS file
CALL cg_close_f(CGNSfile,iErr)
IF (iErr .NE. CG_OK) CALL my_cg_error_exit('Error closing CGNS File.',CGNSFile)
WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"  DONE"
END SUBROUTINE WriteDataToCGNS




SUBROUTINE my_cg_error_exit(ErrorMessage,CGNSFile)
!===================================================================================================================================
! Provides a controlled code abort in case of an CGNS error
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)              :: ErrorMessage ! name of the file that caused an error
PP_CGNS_INT_TYPE ,INTENT(IN)  :: CGNSFile  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)  :: CGNSmessage ! CGNS error message
INTEGER             :: iErr  ! ?
!===================================================================================================================================
! Retrieve error message
CALL cg_get_error_f(CGNSmessage)      
! Exit calculation
CALL cg_close_f(CGNSFile,iErr)
WRITE(Unit_StdOut,'(A)') ErrorMessage
WRITE(Unit_StdOut,'(A)') 'CGNS Error Output: '
WRITE(Unit_StdOut,'(A)') CGNSmessage
CALL abort(__STAMP__, &
          'CGNS error!')
END SUBROUTINE my_cg_error_exit

END MODULE MOD_Output_CGNS
