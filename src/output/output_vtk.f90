
MODULE MOD_Output_VTK
!===================================================================================================================================
! Module for generic data output in vtk xml fromat
!
! WARNING: WriteDataToVTK works only for POSTPROCESSING
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

INTERFACE WriteDataToVTK
  MODULE PROCEDURE WriteDataToVTK
END INTERFACE

!INTERFACE LinkVTKFiles
!  MODULE PROCEDURE LinkVTKFiles
!END INTERFACE

PUBLIC::WriteDataToVTK
!PUBLIC::LinkVTKFiles
!===================================================================================================================================

CONTAINS

SUBROUTINE WriteDataToVTK(dim1,nVal,NPlot,nElems,VarNames,Coord,Values,FileString)
!===================================================================================================================================
! Subroutine to write 3D point data to VTK format
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
INTEGER,INTENT(IN)            :: nElems               ! Number of output elements
REAL,INTENT(IN)               :: Coord(3,1:(Nplot+1)**dim1,nElems)      ! CoordsVector 
CHARACTER(LEN=*),INTENT(IN)   :: VarNames(nVal)          ! Names of all variables that will be written out
REAL,INTENT(IN)               :: Values(nVal,1:(Nplot+1)**dim1,nElems)   ! Statevector 
CHARACTER(LEN=*),INTENT(IN)   :: FileString              ! Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k,iVal,iElem,Offset,nBytes,nVTKElems,nVTKCells,ivtk=44  ! ?
INTEGER            :: INT  ! ?
INTEGER            :: Vertex(2**dim1,(NPlot+1)**dim1*nElems)  ! ?
INTEGER            :: NPlot_p1_3,NPlot_p1_2,NPlot_p1,NodeID,NodeIDElem,ElemType  ! ?
CHARACTER(LEN=35)  :: StrOffset,TempStr1,TempStr2  ! ?
CHARACTER(LEN=200) :: Buffer  ! ?
CHARACTER(LEN=1)   :: lf  ! ?
REAL(KIND=4)       :: Float  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)',ADVANCE='NO')"   WRITE DATA TO VTX XML BINARY (VTU) FILE... "//TRIM(FileString)

NPlot_p1  =(Nplot+1)
NPlot_p1_2=Nplot_p1*Nplot_p1
NPlot_p1_3=NPlot_p1_2*Nplot_p1

! Line feed character
lf = char(10)

! Write file
OPEN(UNIT=ivtk,FILE=TRIM(FileString),ACCESS='STREAM')
! Write header
Buffer='<?xml version="1.0"?>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify file type
nVTKElems=NPlot_p1**dim1*nElems
nVTKCells=NPlot**dim1*nElems
Buffer='  <UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
WRITE(TempStr1,'(I16)')nVTKElems
WRITE(TempStr2,'(I16)')nVTKCells
Buffer='    <Piece NumberOfPoints="'//TRIM(ADJUSTL(TempStr1))//'" &
NumberOfCells="'//TRIM(ADJUSTL(TempStr2))//'">'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify point data
Buffer='      <PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=0
WRITE(StrOffset,'(I16)')Offset
DO iVal=1,nVal
  Buffer='        <DataArray type="Float32" Name="'//TRIM(VarNames(iVal))//'" &
format="appended" offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
  Offset=Offset+SIZEOF(INT)+nVTKElems*SIZEOF(FLOAT)
  WRITE(StrOffset,'(I16)')Offset
END DO
Buffer='      </PointData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify cell data
Buffer='      <CellData> </CellData>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify coordinate data
Buffer='      <Points>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='        <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended" &
offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF(INT)+3*nVTKElems*SIZEOF(FLOAT)
WRITE(StrOffset,'(I16)')Offset
Buffer='      </Points>'//lf;WRITE(ivtk) TRIM(Buffer)
! Specify necessary cell data
Buffer='      <Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
! Connectivity
Buffer='        <DataArray type="Int32" Name="connectivity" format="appended" &
offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF(INT)+2**dim1*nVTKElems*SIZEOF(INT)
WRITE(StrOffset,'(I16)')Offset
! Offsets
Buffer='        <DataArray type="Int32" Name="offsets" format="appended" &
offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Offset=Offset+SIZEOF(INT)+nVTKElems*SIZEOF(INT)
WRITE(StrOffset,'(I16)')Offset
! Elem types
Buffer='        <DataArray type="Int32" Name="types" format="appended" &
offset="'//TRIM(ADJUSTL(StrOffset))//'"/>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='      </Cells>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='    </Piece>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='  </UnstructuredGrid>'//lf;WRITE(ivtk) TRIM(Buffer)
! Prepare append section
Buffer='  <AppendedData encoding="raw">'//lf;WRITE(ivtk) TRIM(Buffer)
! Write leading data underscore
Buffer='_';WRITE(ivtk) TRIM(Buffer)

! Write binary raw data into append section
! Point data
nBytes = nVTKElems*SIZEOF(FLOAT)
DO iVal=1,nVal
  WRITE(ivtk) nBytes,REAL(Values(iVal,:,:),4)
END DO
! Points
nBytes = nBytes * 3
WRITE(ivtk) nBytes
WRITE(ivtk) REAL(Coord(:,:,:),4)
! Connectivity
SELECT CASE(dim1)
CASE(2)
  NodeID = 0
  NodeIDElem = 0
  DO iElem=1,nElems
    DO j=1,NPlot
      DO i=1,NPlot
        NodeID = NodeID+1
        !visuQuadElem
        Vertex(:,NodeID) = (/                  &
          NodeIDElem+i+   j   *(NPlot+1)-1,    & !P4
          NodeIDElem+i+  (j-1)*(NPlot+1)-1,    & !P1(CGNS=tecplot standard)
          NodeIDElem+i+1+(j-1)*(NPlot+1)-1,    & !P2
          NodeIDElem+i+1+ j   *(NPlot+1)-1    /) !P3
      END DO
    END DO
    NodeIDElem=NodeIDElem+NPlot_p1_2
  END DO
CASE(3)
  NodeID=0
  NodeIDElem=0
  DO iElem=1,nElems
    DO k=1,NPlot
      DO j=1,NPlot
        DO i=1,NPlot
          NodeID=NodeID+1
          !
          Vertex(:,NodeID)=(/                                       &
            NodeIDElem+i+   j   *(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P4(CGNS=tecplot standard)
            NodeIDElem+i+  (j-1)*(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P1
            NodeIDElem+i+1+(j-1)*(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P2
            NodeIDElem+i+1+ j   *(NPlot+1)+(k-1)*NPlot_p1_2-1,      & !P3
            NodeIDElem+i+   j   *(NPlot+1)+ k   *NPlot_p1_2-1,      & !P8
            NodeIDElem+i+  (j-1)*(NPlot+1)+ k   *NPlot_p1_2-1,      & !P5
            NodeIDElem+i+1+(j-1)*(NPlot+1)+ k   *NPlot_p1_2-1,      & !P6
            NodeIDElem+i+1+ j   *(NPlot+1)+ k   *NPlot_p1_2-1      /) !P7
        END DO
      END DO
    END DO
    !
    NodeIDElem=NodeIDElem+NPlot_p1_3
  END DO
END SELECT
nBytes = 2**dim1*nVTKElems*SIZEOF(INT)
WRITE(ivtk) nBytes
WRITE(ivtk) Vertex(:,:)
! Offset
nBytes = nVTKElems*SIZEOF(INT)
WRITE(ivtk) nBytes
WRITE(ivtk) (Offset,Offset=2**dim1,2**dim1*nVTKElems,2**dim1)
! Elem type
ElemType =3+3*dim1 !9 VTK_QUAD 12  VTK_HEXAHEDRON
WRITE(ivtk) nBytes
WRITE(ivtk) (ElemType,iElem=1,nVTKElems)
! Write footer
Buffer=lf//'  </AppendedData>'//lf;WRITE(ivtk) TRIM(Buffer)
Buffer='</VTKFile>'//lf;WRITE(ivtk) TRIM(Buffer)
CLOSE(ivtk)
WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"   DONE"
END SUBROUTINE WriteDataToVTK
 
 

END MODULE MOD_Output_VTK
