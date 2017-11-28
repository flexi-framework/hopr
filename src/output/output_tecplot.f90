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

MODULE MOD_Output_Tecplot
!===================================================================================================================================
! Module for generic data output in Tecplot fromat
!
! WARNING: WriteDataToTecplot works only for POSTPROCESSING
!
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

INTERFACE WriteDataToASCIITecplot
  MODULE PROCEDURE WriteDataToASCIITecplot
END INTERFACE


PUBLIC::WriteDataToASCIITecplot
!==================================================================================================================================

CONTAINS

SUBROUTINE WriteDataToASCIITecplot(dim1,nVal,NPlot,nElems,VarNames,Coord,Values,FileString)
!===================================================================================================================================
! Subroutine to write 3D point data to Tecplot format
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: dim1                    ! dimension of the data (either 2 or 3)
INTEGER,INTENT(IN)            :: nVal                    ! Number of nodal output variables
INTEGER,INTENT(IN)            :: NPlot                   ! Number of output points .EQ. NAnalyze
INTEGER,INTENT(IN)            :: nElems                  ! Number of output elements
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVal)          ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Coord(3,1:(NPlot+1)**dim1,nElems)      ! CoordsVector 
REAL,INTENT(IN)               :: Values(nVal,1:(Nplot+1)**dim1,nElems)   ! Statevector 
CHARACTER(LEN=*),INTENT(IN)   :: FileString              ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iElem,iVal,Offset  ! ?
CHARACTER(LEN=255) :: Format_nVal  ! ?
CHARACTER(LEN=500) :: Format_Title  ! ?
CHARACTER(LEN=255) :: VarString  ! ?
INTEGER            :: NodeIDElem,NPlot_p1,nPlot_p1_2,nPlot_p1_3  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"     WRITE DATA TO TECPLOT ASCII FILE... "//TRIM(FileString)

!assemble format strings
WRITE(Format_nVal,'(A1,I4,A19)')'(',3+nVal-1,'(E21.14,1X),E21.14)'
Format_Title(1:51)='VARIABLES="CoordinateX","CoordinateY","CoordinateZ"'
Offset = 52
DO iVal=1,nVal
  WRITE(VarString,'(A2,A,A1)')',"',TRIM(VarNames(iVal)),'"'
  Format_Title(Offset:Offset+LEN(TRIM(VarString)))=TRIM(VarString)
  Offset=Offset+LEN(TRIM(VarString))
END DO

!connected 3D FEM data
NPlot_p1=(NPlot+1)
NPlot_p1_2=NPlot_p1*NPlot_p1
OPEN(44,FILE=TRIM(FileString),Status="REPLACE")
WRITE(44,'(A)')Format_Title(1:Offset-1)
SELECT CASE(dim1)
CASE(2)
  WRITE(44,'(A,I8,A,I8,A)')'ZONE T="",DATAPACKING=POINT, NODES=',nElems*nPlot_p1_2,', ELEMENTS=',    &
                           nElems*(NPlot)**2,',ZONETYPE=FEQUADRILATERAL'
  DO iElem=1,nElems
    DO i=1,NPlot_p1_2
        WRITE(44,Format_nVal) Coord(:,i,iElem),Values(:,i,iElem)
    END DO
  END DO 
  !element connectivity
  NodeIDElem=0
  DO iElem=1,nElems
    DO j=1,NPlot
      DO i=1,NPlot
        !visuQuadElem  
        WRITE(44,'(4(I8,1X))')&
          NodeIDElem+i+  (j-1)*(NPlot+1),      & !P1(CGNS=tecplot standard)
          NodeIDElem+i+1+(j-1)*(NPlot+1),      & !P2
          NodeIDElem+i+1+ j   *(NPlot+1),      & !P3     
          NodeIDElem+i+   j   *(NPlot+1)         !P4
      END DO 
    END DO 
    NodeIDElem=NodeIDElem+NPlot_p1_2
  END DO 
CASE(3)
  NPlot_p1_3=NPlot_p1_2*NPlot_p1
  WRITE(44,'(A,I8,A,I8,A)')'ZONE T="",DATAPACKING=POINT, NODES=',nElems*NPlot_p1_3,', ELEMENTS=',nElems*(NPlot)**3,&
                         ',ZONETYPE=FEBRICK'
  DO iElem=1,nElems
    DO i=1,NPlot_p1_3
      WRITE(44,Format_nVal) Coord(:,i,iElem),Values(:,i,iElem)
    END DO 
  END DO 
  !element connectivity
  NodeIDElem=0
  DO iElem=1,nElems
    DO k=1,NPlot
      DO j=1,NPlot
        DO i=1,NPlot
          !visuHexaElem  
          WRITE(44,'(8(I8,1X))')&
            NodeIDElem+i+  (j-1)*(NPlot+1)+(k-1)*NPlot_p1_2,      & !P1(CGNS=tecplot standard)
            NodeIDElem+i+1+(j-1)*(NPlot+1)+(k-1)*NPlot_p1_2,      & !P2
            NodeIDElem+i+1+ j   *(NPlot+1)+(k-1)*NPlot_p1_2,      & !P3     
            NodeIDElem+i+   j   *(NPlot+1)+(k-1)*NPlot_p1_2,      & !P4
            NodeIDElem+i+  (j-1)*(NPlot+1)+ k   *NPlot_p1_2,      & !P5
            NodeIDElem+i+1+(j-1)*(NPlot+1)+ k   *NPlot_p1_2,      & !P6
            NodeIDElem+i+1+ j   *(NPlot+1)+ k   *NPlot_p1_2,      & !P7     
            NodeIDElem+i+   j   *(NPlot+1)+ k   *NPlot_p1_2         !P8
        END DO 
      END DO 
    END DO 
    NodeIDElem=NodeIDElem+NPlot_p1_3
  END DO 
END SELECT

CLOSE(44)
WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"   DONE"
END SUBROUTINE WriteDataToASCIITecplot




END MODULE MOD_Output_Tecplot
