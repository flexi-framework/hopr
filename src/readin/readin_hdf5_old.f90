#include "defines.f90"
MODULE MOD_Readin_HDF5_OLD
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE HDF5
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tNodePtr,tSidePtr,tElemPtr
!USE MOD_IO_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
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
INTEGER,PARAMETER              :: ELEM_Type=1           !entry position, 
INTEGER,PARAMETER              :: ELEM_Zone=2           
INTEGER,PARAMETER              :: ELEM_FirstSideInd=3
INTEGER,PARAMETER              :: ELEM_LastSideInd=4
INTEGER,PARAMETER              :: ELEM_FirstNodeInd=5
INTEGER,PARAMETER              :: ELEM_LastNodeInd=6

INTEGER,PARAMETER              :: SideInfoSize=4
INTEGER,PARAMETER              :: SIDE_Type=1           !entry position
INTEGER,PARAMETER              :: SIDE_ID=2
INTEGER,PARAMETER              :: SIDE_nbElemID=3
INTEGER,PARAMETER              :: SIDE_BCID=4

INTEGER,ALLOCATABLE            :: ElemInfo(:,:),SideInfo(:,:),NodeInfo(:)
REAL,ALLOCATABLE               :: ElemWeight(:)
REAL,ALLOCATABLE               :: ElemBarycenters(:,:)
REAL,ALLOCATABLE               :: NodeCoords(:,:)
INTEGER,ALLOCATABLE            :: NodeMap(:)
INTEGER,ALLOCATABLE            :: Elem_IJK(:,:)
INTEGER                        :: nElems_IJK(3)
INTEGER                        :: nGlobalElems
INTEGER                        :: nElems,nSides,nNodes,locnSides,locnNodes
INTEGER                        :: ElemCounter(11,2)
INTEGER                        :: offsetElem,offsetSideID,offsetNodeID,offsetSide,offsetNode
INTEGER                        :: iElem,iSide,iNode,i,j,k
INTEGER                        :: nSideIDs,nNodeIDs
INTEGER                        :: nBCs,BoundaryOrder_mesh
LOGICAL                        :: curvedfound
LOGICAL                        :: initMesh=.FALSE. 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadMeshFromHDF5_OLD
  MODULE PROCEDURE ReadMeshFromHDF5_OLD
END INTERFACE

PUBLIC::ReadMeshFromHDF5_OLD
!===================================================================================================================================

CONTAINS
SUBROUTINE ReadMeshFromHDF5_OLD(FileString,ElemWeightFile_in)
!===================================================================================================================================
! Subroutine to read the mesh from HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:usecurveds,N
USE MOD_Mesh_Vars,ONLY:nMeshElems
USE MOD_Mesh_Vars,ONLY:nNodesElemSideMapping,ElemSideMapping
USE MOD_Mesh_Vars,ONLY:BoundaryType
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewNode,getNewBC,getNewSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)  :: FileString  ! ?
CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: ElemWeightFile_in  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER            :: Elem  ! ?
TYPE(tSide),POINTER            :: aSide,bSide  ! ?
INTEGER                        :: i1,j1,k1,nBCSides,nPeriodicSides  ! ?
INTEGER                        :: iNodeP,NodeID,SideID,ElemID,iBC  ! ?
INTEGER                        :: offset,first,last  ! ?
LOGICAL                        :: oriented  ! ?
CHARACTER(LEN=255)             :: ElemWeightFile  ! ?
LOGICAL                        :: fileExists  ! ?
LOGICAL                        :: doConnection  ! ?
!===================================================================================================================================
IF(initMesh) RETURN
INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
IF(.NOT.FileExists)  CALL abort(__STAMP__, &
        'readMesh from HDF5, file "'//TRIM(FileString)//'" does not exist',999,999.)

WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
  !----------------------------------------------------------------------------------------------------------------------------
  !                              WEIGHTS
  !----------------------------------------------------------------------------------------------------------------------------
IF(.NOT.PRESENT(ElemWeightFile_in))THEN
  ElemWeightFile=TRIM(FileString)
ELSE
  ElemWeightFile=TRIM(ElemWeightFile_in)
END IF
WRITE(UNIT_stdOut,'(A)')'READ WEIGHTS FROM HDF5 FILE "'//TRIM(ElemWeightFile)//'" ...'
! Open HDF5 file
CALL OpenHDF5File(ElemWeightFile,create=.FALSE.)

CALL GetHDF5DataSize(File_ID,'ElemWeight',nDims,HSize)
nGlobalElems=HSize(1) !global number of elements
DEALLOCATE(HSize)
nElems=nGlobalElems   !local number of Elements 

CALL CloseHDF5File()
CALL Timer(.FALSE.)

WRITE(UNIT_stdOut,'(A)')'READ MESH FROM HDF5 FILE "'//TRIM(FileString)//'" ...'
CALL Timer(.TRUE.)
! Open HDF5 file
CALL OpenHDF5File(FileString,create=.FALSE.)
CALL ReadArrayFromHDF5(File_ID,'ElemCounter',2,(/11,2/),0,IntegerArray=ElemCounter)
WRITE(UNIT_stdOut,'(A)')  'Number of Elements'
DO i=1,11
  WRITE(UNIT_stdOut,'(I4,1X,I8)')ElemCounter(i,:)
END DO
WRITE(UNIT_stdOut,'(A5,I8)') 'SUM: ',nGlobalElems

CALL readBCs()

CALL GetHDF5Attribute(File_ID,'BoundaryOrder',1,IntegerScalar=BoundaryOrder_mesh)
IF(usecurveds) THEN
  IF(N+1.NE.BoundaryOrder_mesh) THEN
    WRITE(*,*) 'new boundary order = ',N+1,&
                ' Does not correspond to Boundary order in Mesh file = ',BoundaryOrder_mesh
  END IF
END IF
!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from HDF5
ALLOCATE(ElemInfo(1:nElems,ElemInfoSize))
WRITE(*,*)'READ ELEMENTS'
CALL ReadArrayFromHDF5(File_ID,'ElemInfo',2,(/nElems,ElemInfoSize/),0,IntegerArray=ElemInfo)

ALLOCATE(Elems(1:nElems))

nMeshElems=nElems
DO iElem=1,nElems
  iSide=ElemInfo(iElem,ELEM_FirstSideInd) !first index -1 in Sideinfo
  iNode=ElemInfo(iElem,ELEM_FirstNodeInd) !first index -1 in NodeInfo
  CALL getNewElem(Elems(iElem)%ep)
  Elem=>Elems(iElem)%ep
  Elem%Ind    = iElem
  Elem%Type   = ElemInfo(iElem,ELEM_Type)
  Elem%Zone   = ElemInfo(iElem,ELEM_Zone)
  Elem%nNodes = GETNNODES(elem%Type)
  ALLOCATE(Elem%Node(Elem%nNodes))
END DO

!ALLOCATE(ElemBarycenters(nElems,3)) 
!WRITE(*,*)'READ ELEMENT BARYCENTERS'
!CALL ReadArrayFromHDF5(File_ID,'ElemBarycenters',2,(/nElems,3/),0,RealArray=ElemBarycenters)
!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------

!read local Node Info from HDF5 
offset=ElemInfo(1,ELEM_FirstNodeInd) ! hdf5 array starts at 0-> -1
nNodes=ElemInfo(nElems,ELEM_LastNodeInd)-ElemInfo(1,ELEM_FirstNodeInd)
first=offset+1
last =offset+nNodes
ALLOCATE(NodeInfo(first:last))
WRITE(*,*)'READ NODEINFO'
CALL ReadArrayFromHDF5(File_ID,'NodeInfo',1,(/nNodes/),offset,IntegerArray=NodeInfo)

CALL GetNodeMap() !get nNodeIDs and NodeMap from NodeInfo array
LOGWRITE(*,*)'MIN,MAX,SIZE of NodeMap',MINVAL(NodeMap),MAXVAL(NodeMap),SIZE(NodeMap,1)

ALLOCATE(Nodes(1:nNodeIDs)) ! pointer list, entry is known by INVMAP(i,nNodeIDs,NodeMap)
DO i=1,nNodeIDs
  NULLIFY(Nodes(i)%np)
END DO
!assign nodes 
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  iNode=ElemInfo(iElem,ELEM_FirstNodeInd) !first index -1 in NodeInfo
  locnNodes=ElemInfo(iElem,Elem_LastNodeInd)-ElemInfo(iElem,Elem_FirstNodeInd)
  SELECT CASE(Elem%nNodes)
  CASE(3,4,5) !tria, quad, tetra, pyra
    locnSides=Elem%nNodes
  CASE(6) ! prism
    locnSides=5
  CASE(8) ! hexa
    locnSides=6
  END SELECT
  
  IF(N.EQ.1)THEN
    IF(Elem%nNodes.NE.locnNodes-locnSides) & 
            CALL abort(__STAMP__, &
            ' Sanity check, number of element nodes do not fit for Ngeo=1')
  END IF
  DO i=1,Elem%nNodes
    iNode=iNode+1
    NodeID=ABS(NodeInfo(iNode))     !global, unique NodeID
    iNodeP=INVMAP(NodeID,nNodeIDs,NodeMap)  ! index in local Nodes pointer array
    IF(iNodeP.LE.0) STOP 'Problem in INVMAP' 
    IF(.NOT.ASSOCIATED(Nodes(iNodeP)%np))THEN
      CALL GetNewNode(Nodes(iNodeP)%np,0) 
      Nodes(iNodeP)%np%ind=NodeID 
    END IF
    Elem%Node(i)%np=>Nodes(iNodeP)%np
    Nodes(iNodeP)%np%refCount=Nodes(iNodeP)%np%refCount+1
  END DO !i=1,Elem%nNodes
  IF(N.GT.1)THEN !Curved
    SELECT CASE(Elem%nNodes)
    CASE(4) !Tet
      Elem%nCurvedNodes=(N+1)*(N+2)*(N+3)/6
    CASE(5) !pyra
      Elem%nCurvedNodes=(N+1)*(N+2)*(2*N+3)/6
    CASE(6) !prism
      Elem%nCurvedNodes=(N+1)**2*(N+2)/2
    CASE(8) !hex
      Elem%nCurvedNodes=(N+1)**3
    END SELECT
    IF(Elem%nCurvedNodes.NE.LocnNodes-Elem%nNodes-locnSides) &
            CALL abort(__STAMP__, &
            ' Sanity check, number of curved element nodes do not fit',LocnNodes,REAL(N))
    ALLOCATE(Elem%CurvedNode(Elem%nCurvedNodes))
    DO i=1,Elem%nCurvedNodes
      iNode=iNode+1
      NodeID=ABS(NodeInfo(iNode))     !global, unique NodeID
      iNodeP=INVMAP(NodeID,nNodeIDs,NodeMap)  ! index in local Nodes pointer array
      IF(iNodeP.LE.0) STOP 'Problem in INVMAP' 
      IF(.NOT.ASSOCIATED(Nodes(iNodeP)%np))THEN
        CALL GetNewNode(Nodes(iNodeP)%np,0) 
        Nodes(iNodeP)%np%ind=NodeID 
      END IF
      Elem%CurvedNode(i)%np=>Nodes(iNodeP)%np
      Nodes(iNodeP)%np%refCount=Nodes(iNodeP)%np%refCount+1
    END DO !i=1,Elem%nCurvedNodes
  END IF !linear/Curved
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              SIDES
!----------------------------------------------------------------------------------------------------------------------------


offset=ElemInfo(1,ELEM_FirstSideInd) ! hdf5 array starts at 0-> -1  
nSides=ElemInfo(nElems,ELEM_LastSideInd)-ElemInfo(1,ELEM_FirstSideInd)
!read local SideInfo from HDF5 
first=offset+1
last =offset+nSides
ALLOCATE(SideInfo(first:last,SideInfoSize))
CALL ReadArrayFromHDF5(File_ID,'SideInfo',2,(/nSides,SideinfoSize/),offset,IntegerArray=SideInfo)

!WRITE(*,*)'DEBUG, SideInfo'
!DO i=FirstSideInd,LastSideInd
!  WRITE(*,*)i,':',SideInfo(i,:)
!END DO

DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  iNode=ElemInfo(iElem,ELEM_FirstNodeInd) !first index -1 in NodeInfo
  iNode=iNode+Elem%nNodes+Elem%nCurvedNodes
  iSide=ElemInfo(iElem,ELEM_FirstSideInd) !first index -1 in Sideinfo
  !build up sides of the element using element Nodes and CGNS standard
  SELECT CASE(Elem%nNodes)
  CASE(3,4,5) !tria, quad, tetra, pyra
  locnSides=Elem%nNodes
  CASE(6) ! prism
  locnSides=5
  CASE(8) ! hexa
  locnSides=6
  END SELECT
  ALLOCATE(Sides(locnSides))
  DO i=1,locnSides
    locnNodes=nNodesElemSideMapping(Elem%nNodes,i)
    CALL getNewSide(aSide,locnNodes)
    aSide%LocSide=i
    DO j=1,aSide%nNodes
      aSide%Node(j)%np=>Elem%Node(ElemSideMapping(Elem%nNodes,i,j))%np
      aSide%Node(j)%np%RefCount=aSide%Node(j)%np%RefCount+1
    END DO !j=1,Side%nNodes
    Sides(i)%sp=>aSide
    NULLIFY(aSide)
  END DO

  IF(locnSides.NE.ElemInfo(iElem,ELEM_LastSideInd)-ElemInfo(iElem,ELEM_FirstSideInd))THEN
    ! we have hanging sides
    WRITE(*,*)'DEBUG, Hanging sides not implemented yet!!!'
    STOP 
  END IF
  ! assign oriented nodes
  DO i=1,locnSides
    aSide=>Sides(i)%sp
    iSide=iSide+1
    aSide%Elem=>Elem
    oriented=(Sideinfo(iSide,SIDE_ID).GT.0)
    
    aSide%Ind=ABS(SideInfo(iSide,SIDE_ID))
    iNode=iNode+1
    NodeID=NodeInfo(iNode) !first oriented corner node
    IF(oriented)THEN !oriented side
      DO j=1,aSide%nNodes
        aSide%OrientedNode(j)%np=>aSide%Node(j)%np
        aSide%OrientedNode(j)%np%RefCount=aSide%OrientedNode(j)%np%RefCount+1
      END DO !j=1,Side%nNodes
    ELSE !not oriented
      DO k=1,aSide%nNodes
        IF(aSide%Node(k)%np%ind.EQ.ABS(NodeID)) EXIT
      END DO
      IF(k.GT.aSide%nNodes) STOP 'NodeID doesnt belng to side'
      DO j=1,aSide%nNodes
        aSide%OrientedNode(j)%np=>aSide%Node(k)%np
        aSide%OrientedNode(j)%np%RefCount=aSide%OrientedNode(j)%np%RefCount+1
        k=k-1
        IF(k.EQ.0)k=aSide%nNodes
      END DO !j=1,Side%nNodes
    END IF
  END DO !i=1,locnSides
  !transform to side pointer list
  Elem%firstSide=>Sides(1)%sp
  aSide=>Elem%firstSide
  DO i=2,locnSides
    aSide%nextElemSide => Sides(i)%sp
    NULLIFY(Sides(i)%sp)
    aSide=>aSide%nextElemSide
  END DO !i=2,locnsides
  DEALLOCATE(Sides)
END DO !iElem


! build up side connection 
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  iSide=ElemInfo(iElem,ELEM_FirstSideInd) !first index -1 in Sideinfo
  aSide=>Elem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    iSide=iSide+1
    sideID = ABS(SideInfo(iSide,SIDE_ID))
    elemID = SideInfo(iSide,SIDE_nbElemID)
    iBC    = SideInfo(iSide,SIDE_BCID)
    doConnection=.TRUE. ! for periodic sides if BC is reassigned as non periodic
    IF(iBC.NE.0)THEN !BC
      CALL getNewBC(aSide%BC)
      aSide%BC%BCType     = BoundaryType(iBC,1)
      aSide%curveIndex    = BoundaryType(iBC,2)
      aSide%BC%BCState    = BoundaryType(iBC,3)
      aSide%BC%BCalphaInd = BoundaryType(iBC,4)
      aSide%BC%BCindex    = iBC
      IF(aSide%BC%BCType.NE.1)doConnection=.FALSE. 
    END IF
    IF(.NOT.ASSOCIATED(aSide%connection))THEN
      IF((elemID.NE.0).AND.doConnection)THEN !connection 
        IF((elemID.LE.nElems).AND.(elemID.GE.1))THEN !local connection
          bSide=>Elems(elemID)%ep%firstSide
          DO WHILE(ASSOCIATED(bSide))
            IF(bSide%ind.EQ.aSide%ind)THEN
              aSide%connection=>bSide
              bSide%connection=>aSide
              EXIT
            END IF
            bSide=>bSide%nextElemSide
          END DO
        ELSE !MPI connection
          CALL abort(__STAMP__, &
            ' elemID of neighbor not in global Elem list ',999,999.)
        END IF
      END IF
    END IF !connection associated
    aSide=>aSide%nextElemSide
  END DO 
END DO !iElem
!transform to element pointer list
firstElem=>Elems(1)%ep
Elem=>firstElem
DO iElem=2,nElems
  Elem%nextElem          => Elems(iElem)%ep
  Elem%nextElem%prevElem => Elem
  Elem=>Elem%nextElem
  NULLIFY(Elems(iElem)%ep)
END DO

DEALLOCATE(Elems,ElemInfo,SideInfo,NodeInfo)

! get physical coordinates
ALLOCATE(NodeCoords(nNodeIDs,3))

CALL ReadCoordsFromHDF5(File_ID,'NodeCoords',(/nNodeIDs,3/),NodeMap,NodeCoords)
!assign to pointer
DO i=1,nNodeIDs
  IF(ASSOCIATED(Nodes(i)%np))Nodes(i)%np%x=NodeCoords(i,:)
END DO

DEALLOCATE(NodeCoords)
DEALLOCATE(NodeMap)

LOGWRITE(*,*)'DEBUG,check connectivity'
! Check connectivity
i1=0
j1=0
k1=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  aSide=>Elem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    IF(aSide%LocSide .LE. 0) &
      CALL abort(__STAMP__, &
        'io_hdf5: Side%LocSide not set!',999,999.)
    IF(ASSOCIATED(aSide%BC))THEN
      IF((.NOT. ASSOCIATED(aSide%Connection)) .AND. (aSide%BC%BCType .EQ. 1))THEN
        i1=i1+1
      END IF
    ELSE
      k1=k1+1
      IF(.NOT. ASSOCIATED(aSide%Connection))THEN
        DO iNode=1,aSide%nNodes
          LOGWRITE(*,*)aSide%Node(iNode)%np%x
        END DO
        j1=j1+1
      END IF
    END IF
    aSide=>aSide%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO
IF(i1+j1 .GT. 0) THEN
  LOGWRITE(*,*)'..',k1
  CALL abort(__STAMP__, &
    'missing Connection of Side: with/without BC',i1,REAL(j1))
END IF

i1=0
j1=0
k1=0
nBCSides=0
nPeriodicSides=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  aSide=>Elem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    DO iNode=1,aSide%nNodes
      aSide%Node(iNode)%np%tmp=-1
    END DO
    aSide=>aSide%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  i1=i1+1
  aSide=>Elem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    j1=j1+1
    IF(ASSOCIATED(aSide%BC))THEN
      nBCSides=nBCSides+1
      IF(ASSOCIATED(aSide%connection)) nPeriodicSides=nPeriodicSides+1
    END IF
    DO iNode=1,aSide%nNodes
      IF(aSide%Node(iNode)%np%tmp .NE. 0)THEN
        k1=k1+1
        aSide%Node(iNode)%np%tmp=0
      END IF
    END DO
    aSide=>aSide%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO
LOGWRITE(*,*)'nElems:',i1
LOGWRITE(*,*)'nSides:',j1
LOGWRITE(*,*)'nNodes:',k1
LOGWRITE(*,*)'nBCSides:',nBCSides
LOGWRITE(*,*)'nPeriodicSides:',nPeriodicSides

CALL CloseHDF5File() 
initMesh=.TRUE.
CALL Timer(.FALSE.)
WRITE(UNIT_stdOut,'(132("~"))')
END SUBROUTINE ReadMeshFromHDF5_OLD


SUBROUTINE ReadBCs()
!===================================================================================================================================
! Subroutine to write Data to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:nUserDefinedBoundaries,BoundaryName,BoundaryType
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER, ALLOCATABLE           :: BCMapping(:)  ! ?
INTEGER                        :: iBC,iUserBC,nBCs  ! ?
INTEGER                        :: nDims  ! ?
INTEGER                        :: Offset=0 ! Every process reads all BCs
!===================================================================================================================================
! Read boundary names from HDF5 file
CALL GetHDF5DataSize(File_ID,'BCNames',nDims,HSize)
nBCs=HSize(1)
DEALLOCATE(HSize)
ALLOCATE(BCNames(nBCs), BCMapping(nBCs))
CALL ReadArrayFromHDF5(File_ID,'BCNames',1,(/nBCs/),Offset,StrArray=BCNames)  ! Type is a dummy type only
! User may have redefined boundaries in the ini file. So we have to create mappings for the boundaries.
IF(nUserDefinedBoundaries.EQ.0)THEN !no boundary mapping
  ALLOCATE(BoundaryName(nBCs))
  BoundaryName(:)=BCNames(:)
END IF
BCMapping=0
IF(nUserDefinedBoundaries .GT. 0)THEN
  DO iBC=1,nBCs
    DO iUserBC=1,nUserDefinedBoundaries
      IF(INDEX(TRIM(BCNames(iBC)),TRIM(BoundaryName(iUserBC))).NE.0) BCMapping(iBC)=iUserBC
    END DO
  END DO
END IF

! Read boundary types from HDF5 file
CALL GetHDF5DataSize(File_ID,'BCType',nDims,HSize)
IF(HSize(1).NE.nBCs) STOP 'Problem in readBC'
ERRWRITE(*,*)'BCType: ',nDims,HSize(:)
DEALLOCATE(HSize)
ALLOCATE(BCType(nBCs,4))
offset=0
CALL ReadArrayFromHDF5(File_ID,'BCType',2,(/nBCs,4/),Offset,IntegerArray=BCType)
IF(nUserDefinedBoundaries.EQ.0)THEN !no boundary mapping
  ALLOCATE(BoundaryType(nBCs,4))
  BoundaryType(:,:)=BCType(:,:)
END IF
! Now apply boundary mappings
IF(nUserDefinedBoundaries .GT. 0)THEN
  DO iBC=1,nBCs
    IF(BCMapping(iBC) .NE. 0)THEN
      WRITE(*,*)'   Boundary in HDF file found :',TRIM(BCNames(iBC))
      WRITE(*,*)'                    was       :',BCType(iBC,:)
      WRITE(*,*)'                    is set to :',BoundaryType(BCMapping(iBC),:)
      BCType(iBC,:)=BoundaryType(BCMapping(iBC),:)
    END IF
  END DO
END IF

IF(nUserDefinedBoundaries .EQ. 0) nUserDefinedBoundaries=nBCs
DEALLOCATE(BCNames,BCType,BCMapping)
END SUBROUTINE ReadBCs


FUNCTION GETNNODES(ElementType)
!===================================================================================================================================
! Get nNodes from Element Type 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: ElementType  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER                            :: GETNNODES  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
  SELECT CASE(ElementType)
  CASE(3)
    GETNNODES=3
  CASE(6)
    GETNNODES=BoundaryOrder_mesh*(BoundaryOrder_mesh+1)/2
  CASE(4,5)
    GETNNODES=4
  CASE(7)
    GETNNODES=BoundaryOrder_mesh*BoundaryOrder_mesh
  CASE(104,204)
    GETNNODES=4
  CASE(105,115,205)
    GETNNODES=5
  CASE(106,116,206)
    GETNNODES=6
  CASE(108,118,208)
    GETNNODES=8
  END SELECT
END FUNCTION GETNNODES



! HFD5 STUFF
SUBROUTINE GetHDF5DataSize(Loc_ID,DSetName,nDims,Size)
!===================================================================================================================================
! Subroutine to...
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)               :: DSetName  ! ?
INTEGER(HID_T),INTENT(IN)      :: Loc_ID  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)            :: nDims  ! ?
INTEGER(HSIZE_T),POINTER,INTENT(OUT) :: Size(:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: DSet_ID,FileSpace  ! ?
INTEGER(HSIZE_T), POINTER      :: SizeMax(:)  ! ?
!===================================================================================================================================
!WRITE(UNIT_stdOut,'(A,A,A)')'GET SIZE OF "',TRIM(DSetName),'" IN HDF5 FILE... '
! Initialize FORTRAN predefined datatypes

! Get size of array ----------------------------------------------------------------------------------------------------------------
! Open the dataset with default properties.
CALL H5DOPEN_F(Loc_ID, TRIM(DSetName) , DSet_ID, iError)
! Get the data space of the dataset.
CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
! Get number of dimensions of data space
CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(FileSpace, nDims, iError)
! Get size and max size of data space
ALLOCATE(Size(nDims),SizeMax(nDims))
CALL H5SGET_SIMPLE_EXTENT_DIMS_F(FileSpace, Size, SizeMax, iError)
CALL H5SCLOSE_F(FileSpace, iError)
CALL H5DCLOSE_F(DSet_ID, iError)

!WRITE(UNIT_stdOut,*)'...DONE!'
!WRITE(UNIT_stdOut,'(132("-"))')

END SUBROUTINE GetHDF5DataSize


SUBROUTINE ReadArrayFromHDF5(Loc_ID,ArrayName,Rank,nVal,Offset_in,RealArray,IntegerArray,StrArray)
!===================================================================================================================================
! Subroutine to read arrays of rank "Rank" with dimensions "Dimsf(1:Rank)".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: Rank ! ?
INTEGER, INTENT(IN)                :: Offset_in  ! ?
INTEGER, INTENT(IN)                            :: nVal(Rank)  ! ?
INTEGER(HID_T), INTENT(IN)         :: Loc_ID  ! ?
CHARACTER(LEN=*),INTENT(IN)        :: ArrayName  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL              ,DIMENSION(Rank),OPTIONAL,INTENT(OUT) :: RealArray  ! ?
INTEGER           ,DIMENSION(Rank),OPTIONAL,INTENT(OUT) :: IntegerArray  ! ?
CHARACTER(LEN=255),DIMENSION(Rank),OPTIONAL,INTENT(OUT) :: StrArray  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: DSet_ID, Type_ID, MemSpace, FileSpace, PList_ID  ! ?
INTEGER(HSIZE_T)               :: Offset(Rank),Dimsf(Rank)  ! ?
!===================================================================================================================================
LOGWRITE(UNIT_stdOut,'(A,I1.1,A,A,A)')'    READ ',Rank,'D ARRAY "',TRIM(ArrayName),'"'
! Read array -----------------------------------------------------------------------------------------------------------------------
Dimsf=nVal
LOGWRITE(*,*)'Dimsf,Offset=',Dimsf,Offset_in
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
CALL H5DOPEN_F(Loc_ID, TRIM(ArrayName) , DSet_ID, iError)
! Define and select the hyperslab to use for reading.
CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)
Offset(:)=0
Offset(1)=Offset_in
CALL H5SSELECT_HYPERSLAB_F(FileSpace, H5S_SELECT_SET_F, Offset, Dimsf, iError)
! Create property list
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
! Read the data
IF(PRESENT(RealArray))THEN
  CALL H5DREAD_F(DSet_ID,H5T_NATIVE_DOUBLE,&
                  RealArray    ,Dimsf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
END IF
IF(PRESENT(IntegerArray))THEN
  CALL H5DREAD_F(DSet_ID,H5T_NATIVE_INTEGER, &
                  IntegerArray ,Dimsf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
END IF
IF(PRESENT(StrArray))THEN
  ! Get datatype for the character string array
  CALL H5DGET_TYPE_F(DSet_ID, Type_ID, iError)
  CALL H5DREAD_F(DSet_ID,Type_ID,&
                  StrArray     ,Dimsf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
  CALL H5TCLOSE_F(Type_ID, iError)
END IF

! Close the property list
CALL H5PCLOSE_F(PList_ID,iError)
! Close the file dataspace
CALL H5SCLOSE_F(FileSpace,iError)
! Close the dataset
CALL H5DCLOSE_F(DSet_ID, iError)
! Close the memory dataspace
CALL H5SCLOSE_F(MemSpace,iError)
LOGWRITE(UNIT_stdOut,*)'...DONE!'
END SUBROUTINE ReadArrayFromHDF5


SUBROUTINE ReadCoordsFromHDF5(Loc_ID,ArrayName,nVal,ElementList,CoordArray)
!===================================================================================================================================
! Subroutine to read arrays of rank "Rank" with dimensions "Dimsf(1:Rank)".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN)                :: ElementList(:)  ! ?
INTEGER, INTENT(IN)                            :: nVal(2)  ! ?
INTEGER(HID_T), INTENT(IN)         :: Loc_ID  ! ?
CHARACTER(LEN=*),INTENT(IN)        :: ArrayName  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                   :: CoordArray(nVal(1),nVal(2))  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: rank,i,j  ! ?
INTEGER(HID_T)                 :: DSet_ID, MemSpace, FileSpace, PList_ID  ! ?
INTEGER(SIZE_T)                :: num_elements  ! ?
INTEGER(HSIZE_T)               :: Coords(2,nVal(1)*nVal(2)),Dimsf(2)  ! ?
!===================================================================================================================================
rank=2
LOGWRITE(UNIT_stdOut,'(A,I1.1,A,A,A)')'    READ ',Rank,'D ARRAY "',TRIM(ArrayName),'"'
! Read array -----------------------------------------------------------------------------------------------------------------------
Dimsf=nVal
CALL H5SCREATE_SIMPLE_F(Rank, Dimsf, MemSpace, iError)
CALL H5DOPEN_F(Loc_ID, TRIM(ArrayName) , DSet_ID, iError)
! Define and select the hyperslab to use for reading.
CALL H5DGET_SPACE_F(DSet_ID, FileSpace, iError)

!select elements from Element List
num_elements=0
DO j=1,Dimsf(2)
  DO i=1,Dimsf(1)
    num_elements=num_elements+1
    Coords(1,num_elements)=ElementList(i)
    Coords(2,num_elements)=j 
  END DO
END DO
CALL H5SSELECT_ELEMENTS_F(FileSpace, H5S_SELECT_SET_F, rank,num_elements,Coords,iError)

! Create property list
CALL H5PCREATE_F(H5P_DATASET_XFER_F, PList_ID, iError)
! Read the data
CALL H5DREAD_F(DSet_ID,H5T_NATIVE_DOUBLE,CoordArray,Dimsf,iError,mem_space_id=MemSpace,file_space_id=FileSpace,xfer_prp=PList_ID)
! Close the property list
CALL H5PCLOSE_F(PList_ID,iError)
! Close the file dataspace
CALL H5SCLOSE_F(FileSpace,iError)
! Close the dataset
CALL H5DCLOSE_F(DSet_ID, iError)
! Close the memory dataspace
CALL H5SCLOSE_F(MemSpace,iError)

LOGWRITE(UNIT_stdOut,*)'...DONE!'

END SUBROUTINE ReadCoordsFromHDF5


SUBROUTINE GetHDF5Attribute(Loc_ID_in,AttribName,nVal,DatasetName,RealScalar,IntegerScalar,StrScalar,LogicalScalar,&
                                                                  RealArray,IntegerArray)
!===================================================================================================================================
! Subroutine to read attributes from HDF5 file.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(HID_T), INTENT(IN)           :: Loc_ID_in  ! ?
INTEGER,INTENT(IN)                              :: nVal  ! ?
CHARACTER(LEN=*), INTENT(IN)         :: AttribName  ! ?
CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: DatasetName  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL              ,OPTIONAL,INTENT(OUT) :: RealArray(nVal)  ! ?
INTEGER           ,OPTIONAL,INTENT(OUT) :: IntegerArray(nVal)  ! ?
REAL              ,OPTIONAL,INTENT(OUT) :: RealScalar  ! ?
INTEGER           ,OPTIONAL,INTENT(OUT) :: IntegerScalar  ! ?
LOGICAL           ,OPTIONAL,INTENT(OUT) :: LogicalScalar  ! ?
CHARACTER(LEN=255),OPTIONAL,INTENT(OUT) :: StrScalar  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                 :: Attr_ID, Type_ID,Loc_ID  ! ?
INTEGER(HSIZE_T), DIMENSION(1) :: Dimsf  ! ?
INTEGER                        :: inttolog  ! ?
!===================================================================================================================================
LOGWRITE(UNIT_stdOut,*)' READ ATTRIBUTE "',TRIM(AttribName),'" FROM HDF5 FILE...'
Dimsf(1)=nVal
IF(PRESENT(DatasetName))THEN
 ! Open dataset
  CALL H5DOPEN_F(File_ID, TRIM(DatasetName),Loc_ID, iError)
ELSE
  Loc_ID=Loc_ID_in
END IF
! Create scalar data space for the attribute.
! Create the attribute for group Loc_ID.
CALL H5AOPEN_F(Loc_ID, TRIM(AttribName), Attr_ID, iError)
! Write the attribute data.
IF(PRESENT(RealArray))THEN
  CALL H5AREAD_F(Attr_ID, H5T_NATIVE_DOUBLE, RealArray, Dimsf, iError)
END IF
IF(PRESENT(RealScalar))THEN
  CALL H5AREAD_F(Attr_ID, H5T_NATIVE_DOUBLE, RealScalar, Dimsf, iError)
END IF
IF(PRESENT(IntegerArray))THEN
  CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER , IntegerArray, Dimsf, iError)
END IF
IF(PRESENT(IntegerScalar))THEN
  CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER , IntegerScalar, Dimsf, iError)
END IF
IF(PRESENT(LogicalScalar))THEN
  CALL H5AREAD_F(Attr_ID, H5T_NATIVE_INTEGER , inttolog, Dimsf, iError)
  LogicalScalar=(inttolog.EQ.1)
END IF
IF(PRESENT(StrScalar))THEN
  CALL H5AGET_TYPE_F(Attr_ID, Type_ID, iError)  ! Get HDF5 data type for character string
  CALL H5AREAD_F(Attr_ID, Type_ID, StrScalar, Dimsf, iError)
  CALL H5TCLOSE_F(Type_ID, iError)
  LOGWRITE(UNIT_stdOut,*)' SCALAR STRING READ "',TRIM(StrScalar)
END IF
! Close the attribute.
CALL H5ACLOSE_F(Attr_ID, iError)
IF(PRESENT(DataSetName))THEN
  ! Close the dataset and property list.
  CALL H5DCLOSE_F(Loc_ID, iError)
END IF
LOGWRITE(UNIT_stdOut,*)'...DONE!'
END SUBROUTINE GetHDF5Attribute


SUBROUTINE GetNodeMap()
!===================================================================================================================================
! take NodeInfo array, sort it, eliminate mulitple IDs and return the Mapping 1->NodeID1, 2->NodeID2, ... 
! this is useful if the NodeID list of the mesh are not contiguous, essentially occuring when using domain decomposition (MPI)
!===================================================================================================================================
! MODULES
USE MOD_SortingTools,ONLY:Qsort1Int
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: temp(nNodes+1),i,nullpos
!===================================================================================================================================
temp(1)=0
temp(2:nNodes+1)=NodeInfo
!sort
CALL Qsort1Int(temp)
nullpos=INVMAP(0,nNodes+1,temp)
IF(nullpos.LT.1) CALL abort(__STAMP__,'INVMAP failed.',nullpos)
!count unique entries
nNodeIDs=1
DO i=nullpos+2,nNodes+1
  IF(temp(i).NE.temp(i-1)) nNodeIDs = nNodeIDs+1
END DO
!associate unique entries
ALLOCATE(NodeMap(nNodeIDs))
nNodeIDs=1
NodeMap(1)=temp(nullpos+1)
DO i=nullpos+2,nNodes+1
  IF(temp(i).NE.temp(i-1)) THEN
    nNodeIDs = nNodeIDs+1
    NodeMap(nNodeIDs)=temp(i)
  END IF
END DO
END SUBROUTINE GetNodeMap 


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
END MODULE MOD_Readin_HDF5_OLD
