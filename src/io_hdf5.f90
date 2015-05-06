#include "defines.f90"
MODULE MOD_IO_HDF5
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE HDF5
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tNodePtr,tSidePtr,tElemPtr
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
TYPE(tNodePtr),POINTER         :: Nodes(:)
TYPE(tSidePtr),POINTER         :: Sides(:)
TYPE(tElemPtr),POINTER         :: Elems(:)
INTEGER, ALLOCATABLE           :: BCType(:,:)
CHARACTER(LEN=255), ALLOCATABLE:: BCNames(:)
!--- output_vars
INTEGER(HID_T)                 :: File_ID
CHARACTER(LEN=255)             :: OutputFileName
INTEGER(HID_T)                 :: Plist_ID,info
INTEGER                        :: iError
INTEGER(SIZE_T)                :: SizeSet

INTEGER                        :: nDims
INTEGER(HSIZE_T),POINTER       :: HSize(:)

INTEGER,PARAMETER              :: ElemInfoSize=6        !number of entry in each line of ElemInfo
INTEGER,PARAMETER              :: ELEM_Type=1           !entry position in ElemInfo
INTEGER,PARAMETER              :: ELEM_Zone=2           
INTEGER,PARAMETER              :: ELEM_FirstSideInd=3
INTEGER,PARAMETER              :: ELEM_LastSideInd=4
INTEGER,PARAMETER              :: ELEM_FirstNodeInd=5
INTEGER,PARAMETER              :: ELEM_LastNodeInd=6

INTEGER,PARAMETER              :: SideInfoSize=5
INTEGER,PARAMETER              :: SIDE_Type=1           !entry position in SideInfo
INTEGER,PARAMETER              :: SIDE_ID=2
INTEGER,PARAMETER              :: SIDE_nbElemID=3
INTEGER,PARAMETER              :: SIDE_nbLocSide_Flip=4
INTEGER,PARAMETER              :: SIDE_BCID=5
INTEGER,PARAMETER              :: nSidesElem(3:8)= & ! number of sides of an element nSides=nSidesElem(Elem%nNodes)
                                        (/3,4,5,5,0,6/) !tria,quad/tetra,pyra,prism,hex
INTEGER,PARAMETER              :: LinMap(1:8,4:8)=RESHAPE(        & !CGNS -> IJK ordering for element corner nodes 
                                                                    !jNode=LinMap(iNode,Elem%nNodes)
                                       (/1,2,3,4,0,0,0,0,         & !Tet, one-to-one
                                         1,2,4,3,5,0,0,0,         & !Pyra 
                                         1,2,3,4,5,6,0,0,         & !Prism one-to-one
                                         0,0,0,0,0,0,0,0,         & !nothing
                                         1,2,4,3,5,6,8,7/),(/8,5/)) !Hex

INTEGER,ALLOCATABLE            :: ElemInfo(:,:),SideInfo(:,:)
REAL,ALLOCATABLE               :: ElemWeight(:)
REAL,ALLOCATABLE               :: ElemBarycenters(:,:)
INTEGER,ALLOCATABLE            :: GlobalNodeIDs(:)
REAL,ALLOCATABLE               :: NodeCoords(:,:)
INTEGER,ALLOCATABLE            :: Elem_IJK(:,:)
INTEGER                        :: nElems_IJK(3)
INTEGER                        :: nGlobalElems
INTEGER                        :: nElems,nSides,nNodes
INTEGER                        :: ElemCounter(2,11)
INTEGER                        :: nSideIDs,nNodeIDs
INTEGER                        :: nBCs
LOGICAL                        :: curvedfound
LOGICAL                        :: initMesh=.FALSE. 

INTERFACE INVMAP
  MODULE PROCEDURE INVMAP
END INTERFACE

INTERFACE OpenHDF5File
  MODULE PROCEDURE OpenHDF5File
END INTERFACE

INTERFACE CloseHDF5File
  MODULE PROCEDURE CloseHDF5File
END INTERFACE

!===================================================================================================================================

CONTAINS


FUNCTION INVMAP(ID,nIDs,ArrID)
!===================================================================================================================================
! find the inverse Mapping p.e. NodeID-> entry in NodeMap (a sorted array of unique NodeIDs), using bisection 
! if Index is not in the range, -1 will be returned, if it is in the range, but is not found, 0 will be returned!!
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: ID            ! ID to search for
INTEGER, INTENT(IN)                :: nIDs          ! size of ArrID
INTEGER, INTENT(IN)                :: ArrID(nIDs)   ! 1D array of IDs
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                            :: INVMAP               ! index of ID in NodeMap array
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i,maxSteps,low,up,mid  ! ?
!===================================================================================================================================
INVMAP=0
maxSteps=INT(LOG(REAL(nIDs))*1.4426950408889634556)+1    !1/LOG(2.)=1.4426950408889634556
low=1
up=nIDs
IF((ID.LT.ArrID(low)).OR.(ID.GT.ArrID(up))) THEN
  !WRITE(*,*)'WARNING, Node Index Not in local range -> set to -1'
  INVMAP=-1  ! not in the range!
  RETURN
END IF 
IF(ID.EQ.ArrID(low))THEN
  INVMAP=low
ELSEIF(ID.EQ.ArrID(up))THEN
  INVMAP=up
ELSE
  !bisection
  DO i=1,maxSteps
    mid=(up-low)/2+low
    IF(ID .EQ. ArrID(mid))THEN
      INVMAP=mid                     !index found!
      EXIT
    ELSEIF(ID .GT. ArrID(mid))THEN ! seek in upper half
      low=mid
    ELSE
      up=mid
    END IF
  END DO
END IF
END FUNCTION INVMAP 

! HFD5 STUFF
SUBROUTINE OpenHDF5File(FileString,create)
!===================================================================================================================================
! Open HDF5 file and groups
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)   :: FileString  ! ?
LOGICAL,INTENT(IN)            :: create   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: Plist_ID  ! ?
!===================================================================================================================================
LOGWRITE(UNIT_stdOut,'(A)')'  OPEN HDF5 FILE "',TRIM(FileString),'" ...'
! Initialize FORTRAN predefined datatypes
CALL H5OPEN_F(iError)
! Setup file access property list with parallel I/O access (MPI) or with default property list.
CALL H5PCREATE_F(H5P_FILE_ACCESS_F, Plist_ID, iError)
! Open the file collectively.
IF(create)THEN
  CALL H5FCREATE_F(TRIM(FileString), H5F_ACC_TRUNC_F, File_ID, iError, access_prp = Plist_ID)
ELSE
  CALL H5FOPEN_F(TRIM(FileString), H5F_ACC_RDONLY_F, File_ID, iError, access_prp = Plist_ID)
END IF
CALL H5PCLOSE_F(Plist_ID, iError)
LOGWRITE(UNIT_stdOut,*)'...DONE!'
END SUBROUTINE OpenHDF5File


SUBROUTINE CloseHDF5File()
!===================================================================================================================================
! Close HDF5 file and groups 
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
LOGWRITE(UNIT_stdOut,'(A)')'  CLOSE HDF5 FILE...'
! Close file
CALL H5FCLOSE_F(File_ID, iError)
! Close FORTRAN predefined datatypes.
CALL H5CLOSE_F(iError)
File_ID=0
LOGWRITE(UNIT_stdOut,*)'...DONE!'
END SUBROUTINE CloseHDF5File

END MODULE MOD_IO_HDF5