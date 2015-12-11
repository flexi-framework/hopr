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
#include "defines.f90"
MODULE MOD_Readin_HDF5
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE HDF5
USE MOD_IO_HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE ReadMeshFromHDF5
  MODULE PROCEDURE ReadMeshFromHDF5
END INTERFACE

INTERFACE DatasetExists
  MODULE PROCEDURE DatasetExists
END INTERFACE

PUBLIC::ReadMeshFromHDF5,DatasetExists
!===================================================================================================================================

CONTAINS
SUBROUTINE ReadMeshFromHDF5(FileString)
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
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER            :: Elem  ! ?
TYPE(tSide),POINTER            :: aSide,bSide  ! ?
INTEGER                        :: i1,j1,k1,nBCSides,nPeriodicSides  ! ?
INTEGER                        :: NodeID,SideID,ElemID,iBC  ! ?
INTEGER                        :: iNode,iSide,iElem,i,j,k
INTEGER                        :: offset,first,last  ! ?
INTEGER                        :: Ngeo 
INTEGER                        :: locnNodes 
INTEGER                        :: locnSides 
LOGICAL                        :: oriented  ! ?
LOGICAL                        :: fileExists  ! ?
LOGICAL                        :: doConnection  ! ?
!===================================================================================================================================
IF(initMesh) RETURN
INQUIRE (FILE=TRIM(FileString), EXIST=fileExists)
IF(.NOT.FileExists)  CALL abort(__STAMP__, &
        'readMesh from HDF5, file "'//TRIM(FileString)//'" does not exist',999,999.)

WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')'READ MESH FROM HDF5 FILE "'//TRIM(FileString)//'" ...'
! Open HDF5 file
CALL OpenHDF5File(FileString,create=.FALSE.)

CALL GetHDF5Attribute(File_ID,'nElems',1,IntegerScalar=nGlobalElems)
nElems=nGlobalElems   !local number of Elements 


CALL ReadArrayFromHDF5(File_ID,'ElemCounter',2,(/2,11/),0,IntegerArray=ElemCounter)
WRITE(UNIT_stdOut,'(A)')  'Number of Elements'
DO i=1,11
  WRITE(UNIT_stdOut,'(I4,1X,I8)')ElemCounter(:,i)
END DO
WRITE(UNIT_stdOut,'(A5,I8)') 'SUM: ',nGlobalElems

CALL readBCs()

CALL GetHDF5Attribute(File_ID,'Ngeo',1,IntegerScalar=Ngeo)
IF(usecurveds) THEN
  IF(N.NE.Ngeo) THEN
    WRITE(*,*) 'boundary order of inifile = ',N+1,&
               ' does not match the boundary order in the mesh file:', Ngeo+1
    STOP
  END IF
END IF
!----------------------------------------------------------------------------------------------------------------------------
!                              ELEMENTS
!----------------------------------------------------------------------------------------------------------------------------

!read local ElemInfo from HDF5
ALLOCATE(ElemInfo(ElemInfoSize,1:nElems))
CALL ReadArrayFromHDF5(File_ID,'ElemInfo',2,(/ElemInfoSize,nElems/),0,IntegerArray=ElemInfo)

ALLOCATE(Elems(1:nElems))

nMeshElems=nElems
DO iElem=1,nElems
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  iNode=ElemInfo(ELEM_FirstNodeInd,iElem) !first index -1 in NodeInfo
  CALL getNewElem(Elems(iElem)%ep)
  Elem=>Elems(iElem)%ep
  Elem%Ind    = iElem
  Elem%Type   = ElemInfo(ELEM_Type,iElem)
  Elem%Zone   = ElemInfo(ELEM_Zone,iElem)
  Elem%nNodes = MOD(elem%Type,10)
  ALLOCATE(Elem%Node(Elem%nNodes))
END DO

!----------------------------------------------------------------------------------------------------------------------------
!                              NODES
!----------------------------------------------------------------------------------------------------------------------------

!read local Node Info from HDF5 
offset=ElemInfo(ELEM_FirstNodeInd,1) ! hdf5 array starts at 0-> -1
nNodes=ElemInfo(ELEM_LastNodeInd,nElems)-ElemInfo(ELEM_FirstNodeInd,1)
first=offset+1
last =offset+nNodes

ALLOCATE(GlobalNodeIDs(first:last),NodeCoords(1:3,first:last))
CALL ReadArrayFromHDF5(File_ID,'GlobalNodeIDs',1,(/nNodes/),offset,IntegerArray=GlobalNodeIDs)
CALL ReadArrayFromHDF5(File_ID,'NodeCoords',2,(/3,nNodes/),offset,RealArray=NodeCoords)


CALL GetHDF5Attribute(File_ID,'nUniqueNodes',1,IntegerScalar=nNodeIDs)

ALLOCATE(Nodes(1:nNodeIDs)) ! pointer list of unique nodes
DO i=1,nNodeIDs
  NULLIFY(Nodes(i)%np)
END DO

!assign nodes 
DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  iNode=ElemInfo(ELEM_FirstNodeInd,iElem) !first index -1 in NodeCoords
  locnNodes=ElemInfo(Elem_LastNodeInd,iElem)-ElemInfo(Elem_FirstNodeInd,iElem)
  IF(N.EQ.1)THEN
    IF(Elem%nNodes.NE.locnNodes) CALL abort(__STAMP__, &
            ' Sanity check, number of element nodes do not fit for Ngeo=1')
    DO i=1,Elem%nNodes
      iNode=iNode+1
      NodeID=GlobalNodeIDs(iNode)     !global, unique NodeID
      IF(.NOT.ASSOCIATED(Nodes(NodeID)%np))THEN
        CALL GetNewNode(Nodes(NodeID)%np,0) 
        Nodes(NodeID)%np%ind=NodeID 
        Nodes(NodeID)%np%x(:) = NodeCoords(:,iNode)
      END IF
      Elem%Node(LinMap(i,Elem%nNodes))%np=>Nodes(NodeID)%np
      Nodes(NodeID)%np%refCount=Nodes(NodeID)%np%refCount+1
    END DO !i=1,Elem%nNodes
  ELSE !Curved
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
    IF(Elem%nCurvedNodes.NE.LocnNodes) CALL abort(__STAMP__, &
            ' Sanity check, number of curved element nodes do not fit')
    ALLOCATE(Elem%CurvedNode(Elem%nCurvedNodes))
    DO i=1,Elem%nCurvedNodes
      iNode=iNode+1
      NodeID=GlobalNodeIDs(iNode)     !global, unique NodeID
      IF(.NOT.ASSOCIATED(Nodes(NodeID)%np))THEN
        CALL GetNewNode(Nodes(NodeID)%np,0) 
        Nodes(NodeID)%np%ind=NodeID 
        Nodes(NodeID)%np%x(:) = NodeCoords(:,iNode)
      END IF
      Elem%CurvedNode(i)%np=>Nodes(NodeID)%np
      Nodes(NodeID)%np%refCount=Nodes(NodeID)%np%refCount+1
    END DO !i=1,Elem%nCurvedNodes
    Elem%Node(1)%np=>Elem%CurvedNode(1)%np
    Elem%Node(2)%np=>Elem%CurvedNode(N+1)%np
    SELECT CASE(Elem%nNodes)
    CASE(4) !Tet
      Elem%Node(3)%np=>Elem%CurvedNode((N+1)*(N+2)/2)%np
      Elem%Node(4)%np=>Elem%CurvedNode(Elem%nCurvedNodes)%np
    CASE(5) !pyra
      Elem%Node(3)%np=>Elem%CurvedNode((N+1)**2)%np
      Elem%Node(4)%np=>Elem%CurvedNode(N*(N+1)+1)%np
      Elem%Node(5)%np=>Elem%CurvedNode(Elem%nCurvedNodes)%np
    CASE(6) !prism
      Elem%Node(3)%np=>Elem%CurvedNode((N+1)*(N+2)/2)%np
      Elem%Node(4)%np=>Elem%CurvedNode(N*(N+1)*(N+2)/2+1)%np
      Elem%Node(5)%np=>Elem%CurvedNode(N*(N+1)*(N+2)/2+N+1)%np
      Elem%Node(6)%np=>Elem%CurvedNode(Elem%nCurvedNodes)%np
    CASE(8) !hex
      Elem%Node(3)%np=>Elem%CurvedNode((N+1)**2)%np
      Elem%Node(4)%np=>Elem%CurvedNode(N*(N+1)+1)%np
      Elem%Node(5)%np=>Elem%CurvedNode(N*(N+1)**2+1)%np
      Elem%Node(6)%np=>Elem%CurvedNode(N*(N+1)**2+N+1)%np
      Elem%Node(7)%np=>Elem%CurvedNode(Elem%nCurvedNodes)%np
      Elem%Node(8)%np=>Elem%CurvedNode(N*(N+1)*(N+2)+1)%np
    END SELECT
    DO i=1,Elem%nNodes
      NodeID=Elem%Node(i)%np%ind
      Nodes(NodeID)%np%refCount=Nodes(NodeID)%np%refCount+1
    END DO
  END IF !linear/Curved
END DO

DEALLOCATE(NodeCoords,GlobalNodeIDs)
!----------------------------------------------------------------------------------------------------------------------------
!                              SIDES
!----------------------------------------------------------------------------------------------------------------------------


offset=ElemInfo(ELEM_FirstSideInd,1) ! hdf5 array starts at 0-> -1  
nSides=ElemInfo(ELEM_LastSideInd,nElems)-ElemInfo(ELEM_FirstSideInd,1)
!read local SideInfo from HDF5 
first=offset+1
last =offset+nSides
ALLOCATE(SideInfo(SideInfoSize,first:last))
CALL ReadArrayFromHDF5(File_ID,'SideInfo',2,(/SideInfoSize,nSides/),offset,IntegerArray=SideInfo)


DO iElem=1,nElems
  Elem=>Elems(iElem)%ep
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  !build up sides of the element using element Nodes and CGNS standard
  locnSides=nSidesElem(Elem%nNodes)
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

  IF(locnSides.NE.ElemInfo(ELEM_LastSideInd,iElem)-ElemInfo(ELEM_FirstSideInd,iElem))THEN
    ! we have hanging sides
    WRITE(*,*)'DEBUG, more sides in element as standard, not implemented yet!!!'
    STOP 
  END IF
  ! assign oriented nodes
  DO i=1,locnSides
    aSide=>Sides(i)%sp
    iSide=iSide+1
    aSide%Elem=>Elem
    oriented=(Sideinfo(SIDE_ID,iSide).GT.0)
    
    aSide%Ind=ABS(SideInfo(SIDE_ID,iSide))
    aSide%flip=MOD(SideInfo(SIDE_nbLocSide_flip,iSide),10)  !this entry is 10*nbLocSide+Flip
    IF(oriented)THEN !oriented side
      DO j=1,aSide%nNodes
        aSide%OrientedNode(j)%np=>aSide%Node(j)%np
        aSide%OrientedNode(j)%np%RefCount=aSide%OrientedNode(j)%np%RefCount+1
      END DO !j=1,Side%nNodes
    ELSE !not oriented
      k=aSide%flip
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
  iSide=ElemInfo(ELEM_FirstSideInd,iElem) !first index -1 in Sideinfo
  aSide=>Elem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    iSide=iSide+1
    sideID = ABS(SideInfo(SIDE_ID,iSide))
    elemID = SideInfo(SIDE_nbElemID,iSide)
    iBC    = SideInfo(SIDE_BCID,iSide)
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

DEALLOCATE(Elems,ElemInfo,SideInfo)

!ALLOCATE(ElemBarycenters(3,nElems)) 
!CALL ReadArrayFromHDF5(File_ID,'ElemBarycenters',2,(/3,nElems/),0,RealArray=ElemBarycenters)

CALL CloseHDF5File() 
!======================== READ IN FINISHED =================================

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

initMesh=.TRUE.
CALL Timer(.FALSE.)
WRITE(UNIT_stdOut,'(132("~"))')
END SUBROUTINE ReadMeshFromHDF5


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
IF(HSize(2).NE.nBCs) STOP 'Problem in readBC'
ERRWRITE(*,*)'BCType: ',nDims,HSize(:)
DEALLOCATE(HSize)
ALLOCATE(BCType(4,nBCs))
offset=0
CALL ReadArrayFromHDF5(File_ID,'BCType',2,(/4,nBCs/),Offset,IntegerArray=BCType)
IF(nUserDefinedBoundaries.EQ.0)THEN !no boundary mapping
  ALLOCATE(BoundaryType(nBCs,4))
  BoundaryType(:,:)=TRANSPOSE(BCType(:,:))
END IF
! Now apply boundary mappings
IF(nUserDefinedBoundaries .GT. 0)THEN
  DO iBC=1,nBCs
    IF(BCMapping(iBC) .NE. 0)THEN
      WRITE(*,*)'   Boundary in HDF file found :',TRIM(BCNames(iBC))
      WRITE(*,*)'                    was       :',BCType(:,iBC)
      WRITE(*,*)'                    is set to :',BoundaryType(BCMapping(iBC),:)
      BCType(:,iBC)=BoundaryType(BCMapping(iBC),:)
    END IF
  END DO
END IF

IF(nUserDefinedBoundaries .EQ. 0) nUserDefinedBoundaries=nBCs
DEALLOCATE(BCNames,BCType,BCMapping)
END SUBROUTINE ReadBCs




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


SUBROUTINE DatasetExists(Loc_ID,DSetName,Exists,attrib)
!===================================================================================================================================
! Subroutine to check wheter a dataset on the hdf5 file exists
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*)                     :: DSetName
INTEGER(HID_T),INTENT(IN)            :: Loc_ID
LOGICAL,INTENT(IN),OPTIONAL          :: attrib
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)                  :: Exists
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(HID_T)                       :: DSet_ID
INTEGER                              :: hdferr
!===================================================================================================================================
! we have no "h5dexists_f", so we use the error given by h5dopen_f.
! this produces hdf5 error messages even if everything is ok, so we turn the error msgs off 
! during this operation.
! auto error messages off
CALL h5eset_auto_f(0, hdferr)
! Open the dataset with default properties.
IF(PRESENT(attrib).AND.attrib)THEN
  CALL H5AOPEN_F(Loc_ID, TRIM(DSetName), DSet_ID, iError)
  CALL H5ACLOSE_F(DSet_ID, iError)
ELSE
  CALL H5DOPEN_F(Loc_ID, TRIM(DSetName), DSet_ID, iError)
  CALL H5DCLOSE_F(DSet_ID, iError)
END IF
Exists=.TRUE.
IF(iError.LT.0) Exists=.FALSE.
! auto error messages on
CALL h5eset_auto_f(1, hdferr)
END SUBROUTINE DatasetExists


END MODULE MOD_Readin_HDF5
