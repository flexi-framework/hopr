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
#include "hopr.h"


MODULE MOD_Readin_CGNS
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE CGNS
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTERFACE ReadCGNSmesh
  MODULE PROCEDURE ReadCGNSmesh
END INTERFACE

INTERFACE ReadCGNSSurfaceMesh
  MODULE PROCEDURE ReadCGNSSurfaceMesh
END INTERFACE

INTERFACE openBase
  MODULE PROCEDURE openBase 
END INTERFACE

INTERFACE abortCGNS 
  MODULE PROCEDURE abortCGNS
END INTERFACE

PUBLIC::ReadCGNSmesh
PUBLIC::ReadCGNSSurfaceMesh
PUBLIC::openBase
PUBLIC::abortCGNS
!===================================================================================================================================

CONTAINS
SUBROUTINE ReadCGNSmesh()
!===================================================================================================================================
! This subroutine reads the restart data from the CGNS file and prepares the the calculation for continuation. It's pretty long!
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:nMeshFiles,MeshFileName
USE MOD_Mesh_Vars,ONLY:MeshDim
USE MOD_Mesh_Vars,ONLY:n2dNodes
USE MOD_Mesh_Vars,ONLY:nZones
USE MOD_Mesh_Vars,ONLY:FirstElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                    :: iFile  ! ?
PP_CGNS_INT_TYPE           :: iZone  ! ?
PP_CGNS_INT_TYPE           :: nBases,nCGNSZones  ! Number of bases / zones in CGNS file
PP_CGNS_INT_TYPE           :: nNodesGlob         ! Total number of nodes in the mesh file (in 2d: only first layer)
PP_CGNS_INT_TYPE           :: nZonesGlob         ! Total number of zones in the mesh.
PP_CGNS_INT_TYPE           :: CGNSFile,CGNSBase  ! CGNS file handles
PP_CGNS_INT_TYPE           :: CellDim, PhysDim   ! Dimesnion of elements,physical dimension
PP_CGNS_INT_TYPE           :: iError             ! Error flag
PP_CGNS_INT_TYPE           :: md  ! ?
PP_CGNS_INT_TYPE           :: file_type  ! ?
REAL(KIND=4)               :: version  ! ?

CHARACTER(LEN=32)          :: CGName             ! necessary data for CGNS
PP_CGNS_INT_TYPE           :: ZoneType  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')'Reading CGNS mesh...'
md=meshDim
nNodesGlob = 0
nZonesGlob = 0
! Let's fetz now
DO iFile=1,nMeshFiles
  ! Open CGNS file
  CALL OpenBase(TRIM(MeshFileName(iFile)),MODE_READ,md,md,CGNSFile,CGNSBase,.TRUE.)
  CALL CG_VERSION_F(CGNSFile, version, iError)
  WRITE(UNIT_stdOut,*)'CGNS version:',version
  CALL CG_IS_CGNS_F(TRIM(MeshFileName(iFile)), file_type, iError)

  ! Get number of bases in CGNS file
  CALL CG_NBASES_F(CGNSfile,nBases,iError)
  IF (iError .NE. CG_OK) &
    CALL abortCGNS(__STAMP__,CGNSFile)
  IF(nBases .GT. 1) WRITE(UNIT_stdOut,*)'WARNING - Found ',nBases,' bases in CGNS file. Reading only base 1!'
  CGNSBase=1  ! Multiple bases not supported at the moment
  
  ! Check dimensions of CGNS base
  CALL CG_BASE_READ_F(CGNSfile,CGNSBase,CGname,CellDim,PhysDim,iError)
  IF(iError .NE. CG_OK) &
    CALL abortCGNS(__STAMP__,CGNSFile)
  IF((INT(CellDim) .NE. MeshDim))THEN! .OR. (PhysDim .NE. MeshDim))THEN
    WRITE(UNIT_stdOut,*)'ERROR-Invalid dimensions in CGNS file: CellDim=',CellDim,', PhysDim=',PhysDim,'(MeshDim=',MeshDim,')'
    STOP
  END IF
  ! Get number of zones in CGNSBase
  CALL CG_NZONES_F(CGNSfile,CGNSBase,nCGNSZones,iError)
  IF(iError .NE. CG_OK) CALL abortCGNS(__STAMP__,CGNSFile)
  DO iZone=1,nCGNSZones
    nZonesGlob=nZonesGlob+1
    IF(nZonesGlob.GT.nZones)&
      WRITE(UNIT_stdOut,*)'ERROR: number of zones in inifile does not correspond to number of zones in meshfile(s)',nZones

    ! Check structured / unstructured
    CALL cg_zone_type_f(CGNSFile, CGNSBase, iZone, ZoneType, iError)
    IF (iError .NE. CG_OK) CALL cg_error_exit_f()
    IF (ZoneType.EQ.Structured)THEN
      CALL ReadCGNSMeshStruct(FirstElem,CGNSFile,CGNSBase,iZone,nZonesGlob,nNodesGlob)
    ELSEIF(ZoneType.EQ.Unstructured)THEN
      CALL ReadCGNSMeshUnstruct(FirstElem,CGNSFile,CGNSBase,iZone,nZonesGlob,nNodesGlob)
    ELSE
      STOP 'Wrong zone type specifier, should be structured or unstructured.'
    END IF


  END DO ! iZone
  CALL closeFile(CGNSFile)
END DO ! iFile=1,nMeshFiles
! Set total number of 2d nodes. This is needed to generate unique node indices in fill25DMesh.
IF(MeshDim .EQ. 2) n2dNodes=nNodesGlob
  
CALL Timer(.FALSE.)
END SUBROUTINE ReadCGNSmesh


SUBROUTINE ReadCGNSMeshUnstruct(FirstElem_in,CGNSFile,CGNSBase,iZone,nZonesGlob,nNodesGlob)
!===================================================================================================================================
! This subroutine reads unstructured 3D meshes from the CGNS file and prepares the element list. 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tElemPtr,tSide
USE MOD_Mesh_Vars,ONLY:MeshDim
USE MOD_Mesh_Vars,ONLY:BoundaryType
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewNode,getNewBC
USE MOD_Mesh_Vars,ONLY:BugFix_ANSA_CGNS
USE MOD_Mesh_Basis,ONLY:createSides,GetBoundaryIndex
USE MOD_SortingTools,ONLY:Qsort1Int
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(INOUT)                  :: FirstElem_in  ! ?
PP_CGNS_INT_TYPE ,INTENT(IN)         :: CGNSFile,CGNSBase  ! ?
PP_CGNS_INT_TYPE ,INTENT(IN)         :: iZone  ! ?
PP_CGNS_INT_TYPE ,INTENT(IN)         :: nZonesGlob  ! ?
PP_CGNS_INT_TYPE ,INTENT(INOUT)      :: nNodesGlob  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElemPtr),ALLOCATABLE   :: Elems(:)                            ! Pointer array to elements
TYPE(tElem),POINTER          :: ActualElem  ! ?
TYPE(tSide),POINTER          :: Side                                ! Local side pointer
REAL, ALLOCATABLE            :: NodeCoords(:,:)                     ! Coordinates of grid nodes (nDim,nNodes)
REAL                         :: NorVec(3)  ! ?
PP_CGNS_INT_TYPE ,ALLOCATABLE:: ElemMapping(:)                      ! Global element index -> surface or volume element index
PP_CGNS_INT_TYPE ,ALLOCATABLE:: ElemConnect(:,:), SurfElemConnect(:,:) ! Connectivity data for (surface) element (nNodes+1,nElems)
PP_CGNS_INT_TYPE ,ALLOCATABLE:: LocalConnect(:), ParentData(:,:)    ! Arrays for reading the (local=one section) connectivity data
PP_CGNS_INT_TYPE ,ALLOCATABLE:: BCPoints(:), BCElemList(:)          ! Boundary node /side indices
PP_CGNS_INT_TYPE             :: iNode,iElem,iBCElem,iElemGlob,iVolElem,iSurfElem
PP_CGNS_INT_TYPE             :: dm, iSect, DimVec(3),j  ! ?
PP_CGNS_INT_TYPE             :: IndMin, IndMax, iStart, iEnd        ! Start and end index
PP_CGNS_INT_TYPE             :: nSect                               ! Number of sections in CGNS file
PP_CGNS_INT_TYPE             :: nZoneElems, nSectElems, nSurfElems  ! Number of elems in zone (volume+surface)/section/on surface
PP_CGNS_INT_TYPE             :: nNodes,nElems                       ! Number of nodes / elements in a zone (volume only)
PP_CGNS_INT_TYPE             :: nELemNodes                          ! Number of nodes of an element
PP_CGNS_INT_TYPE             :: iBC, nCGNSBC, nBCElems, nBCPoints  ! ?
PP_CGNS_INT_TYPE             :: BCTypeIndex                         ! Index of boundary condition defined in parameter file
PP_CGNS_INT_TYPE             :: LocDim, LocType, nNodesLoc          ! Dimension, type, number of nodes of local=section elements
PP_CGNS_INT_TYPE             :: SizeZone(3)                         ! CGNS datastructure variables
PP_CGNS_INT_TYPE             :: SectionElemType                     ! Type of elements in CGNS file
PP_CGNS_INT_TYPE             :: ParentDataFlag                      ! 0=no parent data for elems available, 1=parent data available
PP_CGNS_INT_TYPE             :: PntSetType                          ! BC data format (points or surface elemnents)
PP_CGNS_INT_TYPE             :: NormalListFlag,NormalIndex(3)          ! CGNS datastructure variables
PP_CGNS_INT_TYPE             :: DataType                            ! CGNS datastructure variables
PP_CGNS_INT_TYPE             ::                         nDataSet    ! CGNS datastructure variables
PP_CGNS_INT_TYPE             :: CellDim, PhysDim                    ! Dimesnion of elements,physical dimension
PP_CGNS_INT_TYPE             :: iError                              ! Error flag
PP_CGNS_INT_TYPE             :: FirstElemInd  ! ?
CHARACTER(LEN=32)            :: CGName                              ! necessary data for CGNS
CHARACTER(LEN=30)            :: coordNameCGNS(3)                 ! List of CGNS names for the coordinates
LOGICAL, ALLOCATABLE         :: NodeIsBCNode(:)                     ! .TRUE. = node iNode is boundary node
LOGICAL                      :: SideIsBCSide                        ! .TRUE. = Side is a boundary side
PP_CGNS_INT_TYPE             :: skip  ! ?
PP_CGNS_INT_TYPE             :: one  ! ?
LOGICAL                      :: orient2D  ! ?
INTEGER,ALLOCATABLE          :: nBCNodes(:),BCInds(:,:)
INTEGER                      :: locInds(4),nUnique
LOGICAL,ALLOCATABLE          :: BCFound(:)
!===================================================================================================================================
coordNameCGNS(1) = 'CoordinateX'
coordNameCGNS(2) = 'CoordinateY'
coordNameCGNS(3) = 'CoordinateZ'
one=1
! Check dimensions of CGNS base
CALL CG_BASE_READ_F(CGNSfile,CGNSBase,CGname,CellDim,PhysDim,iError)
IF(iError .NE. CG_OK) &
  CALL abortCGNS(__STAMP__,CGNSFile)
IF((INT(CellDim) .NE. MeshDim) .OR. (INT(PhysDim) .NE. MeshDim))THEN
  WRITE(UNIT_stdOut,*)'ERROR-Invalid dimensions in CGNS file: CellDim=',CellDim,', PhysDim=',PhysDim,'(MeshDim=',MeshDim,')'
  STOP
END IF
! Start with reading zones: total number of Nodes and Elems
CALL CG_ZONE_READ_F(CGNSfile,CGNSBase,iZone,CGname,SizeZone,iError)
IF (iError .NE. CG_OK) &
  CALL abortCGNS(__STAMP__,CGNSFile)
WRITE(UNIT_stdOut,*)'Read Zone ',TRIM(CGname)

! Read node coordinates
nNodes=SizeZone(1)
ALLOCATE(NodeCoords(3,nNodes))
NodeCoords=0.
DO dm=1,MeshDim
  CGname=TRIM(CoordNameCGNS(dm))
  CALL CG_COORD_READ_F(CGNSfile,CGNSBase,iZone,CGName,RealDouble,one,nNodes,NodeCoords(dm,:),iError)
  IF (iError .NE. CG_OK)THEN
    WRITE(UNIT_stdOut,*)'ERROR - Could not read coordinate(',dm,'): ',TRIM(CoordNameCGNS(dm))
    CALL CG_NCOORDS_F(CGNSFile,CGNSBase,iZone,PhysDim,iError )  ! Here we use PhysDim as nCoords
    WRITE(UNIT_stdOut,*)'Number of coordinates in CGNS file: ',PhysDim
    WRITE(UNIT_stdOut,*)'Mesh dimension: ',MeshDim
    DO j=1,PhysDim
      CALL CG_COORD_INFO_F(CGNSFile,CGNSBase,iZone,j,LocType,CGName,iError)
      WRITE(UNIT_stdOut,*)'Coordname (',j,')=',TRIM(CGName)
      WRITE(UNIT_stdOut,*)'Datatype=',LocType
    END DO
    CALL abortCGNS(__STAMP__,CGNSFile)
  END IF
END DO
! Read in CGNS sections and count volume and boundary face elems
CALL CG_NSECTIONS_F(CGNSfile,CGNSBase,iZone,nSect,iError)
IF (iError .NE. CG_OK) CALL abortCGNS(__STAMP__,CGNSFile)
nZoneElems=0
DO iSect=1,nSect ! Vol. and boundary face elems
  CALL CG_SECTION_READ_F(CGNSfile,CGNSBase,iZone,iSect,CGname,SectionElemType,IndMin,IndMax,nBCElems,ParentDataFlag,iError)
  IF (iError .NE. CG_OK) CALL abortCGNS(__STAMP__,CGNSFile)
  nSectElems=1+IndMax-IndMin
  nZoneElems=nZoneElems+nSectElems
END DO
ALLOCATE(ElemMapping(nZoneElems))   ! Global element index -> volume / face element index
ElemMapping(:)=0
nElems=SizeZone(2)


! Read element connectivity
ALLOCATE(Elems(nElems))
ALLOCATE(ElemConnect(13,nElems)) ! max 8 + 1 schalter
nSurfElems=nZoneElems-nElems
ALLOCATE(SurfElemConnect(5,nSurfElems))

iVolElem =0
iSurfElem=0
DO iSect=1,nSect ! Vol. and Face elems
  ! Read in Elem indMin & indMax
  CALL CG_SECTION_READ_F(CGNSfile,CGNSBase,iZone,iSect,CGname,SectionElemType,IndMin,IndMax,nBCElems,ParentDataFlag,iError)
  WRITE(UNIT_StdOut,*)'   read section',TRIM(CGname)
  IF (iError .NE. CG_OK) CALL abortCGNS(__STAMP__,CGNSFile)
  IF(SectionElemType .LT. TRI_3) CYCLE !ignore additional sections with data <nDim-1
  CALL CG_ELEMENTDATASIZE_F(CGNSFile,CGNSBase,iZone,iSect,nSectElems,iError)  ! Get number of connectivity values
  ALLOCATE(LocalConnect(nSectElems))
  nSectElems=1+IndMax-IndMin ! Important for surface elements only 
                             ! (nSectElems, Parent1 | Parent2 | ParentSide1 | ParentSide2)...but we don't use it
  ALLOCATE(ParentData(nSectElems,4))
  ! Read in local connectivity data
  CALL CG_ELEMENTS_READ_F(CGNSfile,CGNSBase,iZone,iSect,LocalConnect,ParentData,iError)
  
  ! Check if 2D element is not oriented in z+, check only first element#
  IF(MeshDim .EQ. 2)THEN
    orient2D=.TRUE. !
    IF(SectionElemType .EQ. MIXED) THEN
      locType=LocalConnect(1)
      iStart=2
    ELSE
      locType=SectionElemType
      iStart=1
    END IF
    CALL CG_NPE_F(loctype,nNodesLoc,iError)
    IF(((nodeCoords(1,LocalConnect(iStart+1          ))-nodeCoords(1,LocalConnect(iStart)))* &
        (nodeCoords(2,LocalConnect(iStart+nNodesLoc-1))-nodeCoords(2,LocalConnect(iStart)))- &
        (nodeCoords(2,LocalConnect(iStart+1          ))-nodeCoords(2,LocalConnect(iStart)))* &
        (nodeCoords(1,LocalConnect(iStart+nNodesLoc-1))-nodeCoords(1,LocalConnect(iStart)))  ).LT. 0. ) THEN
        orient2D=.false. ! wrong orientation
    END IF
  END IF  !(MeshDim .EQ. 2)

  iStart=1
  DO iElem=1,nSectElems
    iEnd=iStart
    IF(SectionElemType .EQ. MIXED) THEN
      LocType=LocalConnect(iStart)  ! Mixed type elements, read elem type from array
      iStart =iStart+1              ! First value is elem type
    ELSE
      LocType=SectionElemType       ! Single type elements, elem type = SectionElemType
      iEnd   =iEnd-1                ! Only nElemNodes values
    END IF
    CALL CG_NPE_F(LocType,nNodesLoc,iError) ! Get number of nodes for iElem
    iEnd=iEnd+nNodesLoc

    LocDim=1
    IF(LocType .GT.  BAR_3) LocDim=2
    IF(LocType .GT. QUAD_9) LocDim=3
    IF(LocType .GT.  MIXED) LocDim=2 !NGON_n

    IF(LocDim .EQ. MeshDim) THEN ! volume element
      iVolElem=iVolElem+1
      IF(iVolElem.EQ.1) FirstElemInd=IndMin+iElem-1 !start of volume zone, only possible fro ONE VOLUME ZONE!
      IF(iVolElem .GT. nElems)THEN
        CALL closeFile(CGNSFile)
        CALL abort(__STAMP__,&
                       'Something wrrrrong with element numbers in CGNS File zone :',INT(iZone))
      END IF

      ElemConnect(1            ,iVolElem)=LocType
      ElemConnect(2:nNodesLoc+1,iVolElem)=LocalConnect(iStart:iEnd)
      IF(MeshDim .EQ. 2)THEN
        IF(.NOT. orient2D) THEN
          skip=ElemConnect(2,iVolElem)
          ElemConnect(2,iVolElem)=ElemConnect(1+nNodesLoc,iVolElem)
          ElemConnect(1+nNodesLoc,iVolElem)=skip
          IF(nNodesLoc .EQ. 4) THEN
            skip=ElemConnect(3,iVolElem)
            ElemConnect(3,iVolElem)=ElemConnect(4,iVolElem)
            ElemConnect(4,iVolElem)=skip
          END IF 
        END IF
      END IF  !(MeshDim .EQ. 2)
      ElemMapping(IndMin+iElem-1)        =iVolElem

    ELSEIF(LocDim .EQ. MeshDim-1) THEN ! surface element
      iSurfElem=iSurfElem+1
      IF(iSurfElem.GT.nSurfElems)THEN
        CALL closeFile(CGNSFile)
        CALL abort(__STAMP__,&
                       'Something wrrrrong with surf element numbers in CGNS File zone :',INT(iZone))
      END IF
      
      SurfElemConnect(1            ,iSurfElem)=LocType
      SurfElemConnect(2:nNodesLoc+1,iSurfElem)=LocalConnect(iStart:iEnd)
      ElemMapping(IndMin+iElem-1)             =iSurfElem

    END IF   ! LocDim .EQ. MeshDim
    iStart=iEnd+1
  END DO ! elements in section
  DEALLOCATE(LocalConnect,ParentData)
END DO !sections

! Rebuild the elements of zone iZone
DO iElem=1,nElems
  iVolElem=ElemMapping(FirstElemInd-1+iElem)
  CALL CG_NPE_F(ElemConnect(1,iVolElem),nElemNodes,iError)
  IF (iError .NE. CG_OK) CALL abortCGNS(__STAMP__,CGNSFile)
  CALL GetNewElem(Elems(iElem)%EP)
  Elems(iElem)%ep%Zone  =nZonesGlob
  Elems(iElem)%ep%Ind   =iVolElem
  Elems(iElem)%ep%nNodes=nElemNodes
  ALLOCATE(Elems(iElem)%EP%Node(nElemNodes))
  DO iNode=1,nElemNodes
    CALL GetNewNode(Elems(iElem)%EP%Node(iNode)%NP)
    Elems(iElem)%EP%Node(iNode)%NP%Ind     =ElemConnect(iNode+1,iVolElem)+nNodesGlob  ! Make node indices unique
    Elems(iElem)%EP%Node(iNode)%NP%x       =NodeCoords(:,ElemConnect(iNode+1,iVolElem))
    Elems(iElem)%EP%Node(iNode)%NP%RefCount=1
  END DO
  CALL CreateSides(Elems(iElem)%EP,.TRUE.)
END DO !nElems
DEALLOCATE(ElemConnect,NodeCoords)

! Now read in all boundary data
CALL CG_NBOCOS_F(CGNSfile,CGNSBase,iZone,nCGNSBC,iError)
IF (iError .NE. CG_OK) CALL abortCGNS(__STAMP__,CGNSFile)
ALLOCATE(NodeIsBCNode(nNodes))

DO iBC=1,nCGNSBC
  NodeIsBCNode(:)=.FALSE.    ! Reset
  CALL CG_BOCO_INFO_F(CGNSfile,CGNSBase,iZone,iBC,CGname,BCTypeIndex,PntSetType,nBCElems,NormalIndex,NormalListFlag, &
                                                                                         DataType,nDataSet,iError)
  IF (iError .NE. CG_OK) CALL  abortCGNS(__STAMP__,CGNSFile)
  IF (NormalListFlag .NE. 0)THEN
    CALL closeFile(CGNSFile)
    CALL abort(__STAMP__,&
              'incompatible Boundary definition NormalList')
  END IF

  IF(INDEX(TRIM(CGName),'_on_CURVES') .GT. 0)CYCLE  !ICEM fix, always writes this BC in addition
  
  BCTypeIndex=GetBoundaryIndex(CGName)
  IF(BCTypeIndex .EQ. -1)THEN
    WRITE(UNIT_stdOut,*)'ERROR - Could not find corresponding boundary definition of', CGName
    CYCLE
  END IF

  IF(Bugfix_ANSA_CGNS) PntSetType=ElementList

  ! Boundary is given as a list of boundary nodes
  IF((PntSetType .EQ. PointList).OR.(PntSetType .EQ. PointRange))THEN
    IF(PntSetType .EQ. PointRange) THEN
      CALL CG_BOCO_READ_F(CGNSfile,CGNSBase,iZone,iBC,DimVec,NorVec,iError)
      nBCPoints=1+DimVec(2)-DimVec(1)
      DO j=1,nBCPoints
        iNode=DimVec(1)-1+j  ! Global node index
        NodeIsBCNode(iNode)=.TRUE.
      END DO
    ELSE
      nBCPoints=nBCElems
      ALLOCATE(BCPoints(nBCPoints))
      CALL CG_BOCO_READ_F(CGNSfile,CGNSBase,iZone,iBC,BCPoints,NorVec,iError)
      DO j=1,nBCPoints
        NodeIsBCNode(BCPoints(j))=.TRUE.
      END DO
      DEALLOCATE(BCPoints)
    END IF

    ! I think we will finish the boundaries
    DO iElem=1,nElems
      Side=>Elems(iElem)%EP%firstSide
      DO WHILE(ASSOCIATED(Side))
        DO iNode=1,Side%nNodes
          SideIsBCSide=NodeIsBCNode(Side%Node(iNode)%NP%Ind-nNodesGlob) ! Check if node is a boundary node
          IF(.NOT. SideIsBCSide) EXIT                        ! All nodes of a side must be boundary nodes
        END DO
        IF(SideIsBCSide)THEN
          CALL getNewBC(Side%BC)
          Side%BC%BCType    =BoundaryType(BCTypeIndex,1)
          Side%CurveIndex   =BoundaryType(BCTypeIndex,2)
          Side%BC%BCstate   =BoundaryType(BCTypeIndex,3)
          Side%BC%BCalphaInd=BoundaryType(BCTypeIndex,4)
          Side%BC%BCIndex   =BCTypeIndex
        END IF
        Side=>Side%nextElemSide
      END DO  ! WHILE(ASSOCIATED(Side))
    END DO  ! iElem=1,nElems

  ! Boundary is given as a list of boundary elements
  ELSEIF ((PntSetType .EQ. ElementList).OR.(PntSetType .EQ. ElementRange)) THEN
    IF(PntSetType .EQ. ElementRange)THEN
      CALL CG_BOCO_READ_F(CGNSfile,CGNSBase,iZone,iBC,DimVec,NorVec,iError)
      nBCElems=1+DimVec(2)-DimVec(1)
      ALLOCATE(BCElemList(nBCElems))
      DO iElem=1,nBCElems
        BCElemList(iElem)=iElem-1+DimVec(1)
      END DO
    ELSE
      ALLOCATE(BCElemList(nBCElems))
      CALL CG_BOCO_READ_F(CGNSfile,CGNSBase,iZone,iBC,BCElemList,NorVec,iError)
    END IF

    ALLOCATE(nBCNodes(nBCElems))
    ALLOCATE(BCInds(4,nBCElems))
    ALLOCATE(BCFound(nBCElems))
    DO iElem=1,nBCElems
      iElemGlob=BCElemList(iElem)              ! Global element index
      iSurfElem=ElemMapping(iElemGlob)         ! Surface element index (boundary elements must be surface elements)
      LocType  =SurfElemConnect(1,iSurfElem)   ! Element type (important for mixed types mode)
      CALL CG_NPE_F(LocType,nNodesLoc,iError)  ! Get number of nodes of element
      nBCNodes(iElem)=nNodesLoc
      BCInds(1:nNodesLoc,iElem)=SurfElemConnect(2:nNodesLoc+1,iSurfElem)
      CALL Qsort1Int(BCInds(1:nNodesLoc,iElem))
      NodeIsBCNode(  BCInds(1:nNodesLoc,iElem))=.TRUE.
    END DO
    DEALLOCATE(BCElemList)

    ! I think we will finish the boundaries
    BCFound=.FALSE.
    DO iElem=1,nElems
      Side=>Elems(iElem)%EP%firstSide
      DO WHILE(ASSOCIATED(Side))
        DO iNode=1,Side%nNodes
          locInds(iNode)=Side%Node(iNode)%NP%Ind-nNodesGlob
        END DO
        IF(.NOT.ALL(NodeIsBCNode(locInds(1:side%nNodes))).OR.ASSOCIATED(Side%BC))THEN ! speed up search
          Side=>Side%nextElemSide
          CYCLE
        END IF
        SideIsBCSide=.FALSE.
        ! sort nodes and remove duplicates
        CALL Qsort1Int(locInds(1:Side%nNodes))
        nUnique=1
        DO iNode=2,Side%nNodes
          IF(locInds(nUnique).NE.locInds(iNode))THEN
            nUnique=nUnique+1
            locInds(nUnique)=locInds(iNode)
          END IF
        END DO
        DO iBCElem=1,nBCElems
          IF(BCFound(iBCElem)) CYCLE
          IF(nUnique.NE.nBCNodes(iBCElem)) CYCLE
          SideIsBCSide=(ALL(locInds(1:nUnique).EQ.BCInds(1:nUnique,iBCElem)))
          IF(SideIsBCSide)THEN
            BCFound(iBCElem)=.TRUE.
            EXIT
          END IF
        END DO
        IF(SideIsBCSide)THEN
          CALL getNewBC(Side%BC)
          Side%BC%BCType    =BoundaryType(BCTypeIndex,1)
          Side%CurveIndex   =BoundaryType(BCTypeIndex,2)
          Side%BC%BCstate   =BoundaryType(BCTypeIndex,3)
          Side%BC%BCalphaInd=BoundaryType(BCTypeIndex,4)
          Side%BC%BCIndex   =BCTypeIndex
        END IF
        Side=>Side%nextElemSide
      END DO  ! WHILE(ASSOCIATED(Side))
    END DO  ! iElem=1,nElems
    IF(.NOT.ALL(BCFound)) THEN
      PRINT*,'BC sides found/expected',COUNT(BCFound),nBCElems
      CALL abort(__STAMP__,'ERROR: not all BC sides could be associated')
    END IF
    DEALLOCATE(nBCNodes,BCInds,BCFound)
  ELSE
    CALL closeFile(CGNSFile)
    CALL abort(__STAMP__,&
               'unknown BC specification')
  END IF
END DO !nCGNSBC (boundaries done...)

DEALLOCATE(ElemMapping,SurfElemConnect,NodeIsBCNode)

DO iElem=1,nElems
  IF(.NOT. ASSOCIATED(FirstElem_in))THEN
    FirstElem_in   => Elems(iElem)%EP
  ELSE
    ActualElem                   => Elems(iElem)%EP
    ActualElem%nextElem          => FirstElem_in
    ActualElem%nextElem%prevElem => ActualElem
    FirstElem_in                    => ActualElem
  END IF
  NULLIFY(FirstElem_in%prevElem)
END DO !iElem=1,nElems
DEALLOCATE(Elems)
nNodesGlob = nNodesGlob+nNodes

END SUBROUTINE ReadCGNSMeshUnstruct


SUBROUTINE ReadCGNSMeshStruct(FirstElem_in,CGNSFile,CGNSBase,iZone,nZonesGlob,nNodesGlob)
!===================================================================================================================================
! This subroutine reads structured 2D and 3D meshes from the CGNS file and prepares the element list.
!===================================================================================================================================
! MODULES
USE MOD_CartMesh ,ONLY:GetNewHexahedron
USE MOD_Mesh_Vars,ONLY:tElem,tElemPtr,tSide,tNodePtr
USE MOD_Mesh_Vars,ONLY:DZ,nMeshElems,meshDim
USE MOD_Mesh_Vars,ONLY:BoundaryType,useCurveds,N,NBlock,MeshIsAlreadyCurved
USE MOD_Mesh_Vars,ONLY:nSkip,nSkipZ
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewNode,getNewBC,GETNEWQUAD,deleteNode
USE MOD_Mesh_Basis,ONLY:createSides,GetBoundaryIndex
USE MOD_Basis_Vars,ONLY:HexaMapInv
USE MOD_Basis     ,ONLY:GetVandermonde
USE MOD_ChangeBasis,ONLY:ChangeBasis2D,ChangeBasis3D
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(INOUT)                  :: FirstElem_in  ! ?
PP_CGNS_INT_TYPE ,INTENT(IN)         :: CGNSFile,CGNSBase  ! ?
PP_CGNS_INT_TYPE ,INTENT(IN)         :: iZone  ! ?
PP_CGNS_INT_TYPE ,INTENT(IN)         :: nZonesGlob  ! ?
PP_CGNS_INT_TYPE ,INTENT(INOUT)      :: nNodesGlob  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES

TYPE(tNodePtr),POINTER        :: Mnodes(:,:,:)  ! ?
TYPE(tNodePtr)                :: CornerNode(8)  ! temporary corner nodes of the element (hexaheadron)
TYPE(tElem),POINTER           :: aElem  ! ?
TYPE(tSide),POINTER           :: aSide  ! ?
PP_CGNS_INT_TYPE ,ALLOCATABLE :: DimVec(:,:)  ! ?
PP_CGNS_INT_TYPE              :: irmin(3),irmax(3),irmaxorg(3)  ! ?
PP_CGNS_INT_TYPE              :: nCGNSBC,iBC,iBCFace  ! Number of Boundary Conditions
PP_CGNS_INT_TYPE ,ALLOCATABLE :: isize(:,:)  ! ?
PP_CGNS_INT_TYPE              :: iSide,nBCElems,PntSetType  ! ?
PP_CGNS_INT_TYPE ,ALLOCATABLE :: BCIndex(:,:),BCTypeIndex(:),countBCs(:),nBCFaces(:)  ! ?
PP_CGNS_INT_TYPE              :: BCTypeI  ! ?
PP_CGNS_INT_TYPE              :: NormalListFlag,DataType,nDataSet  ! ?
PP_CGNS_INT_TYPE              :: iError                                 ! Error flag
PP_CGNS_INT_TYPE              :: i,j,k,l,m,step,kk,ll,mm  ! ?
PP_CGNS_INT_TYPE              :: k2,k3  ! ?
PP_CGNS_INT_TYPE              :: stepk,stepl,stepm  ! ?
PP_CGNS_INT_TYPE              :: whichdir  ! ?
PP_CGNS_INT_TYPE              :: nElems(3)
REAL                          :: dir(3,3),scalprod  ! ?
INTEGER                       :: N_loc  ! ?
PP_CGNS_INT_TYPE              :: SideMap(3,2)                        ! xi/eta/zeta plus/minus -> localSideInd
REAL,ALLOCATABLE              :: NodeCoords(:,:,:,:),NodeCoordsTmp(:,:,:,:) ! Coordinates of the Nodes in x,y,z direktion
REAL,ALLOCATABLE              :: Vdm_NBlock_N(:,:) ! Vdm from NBlock to N (if curved)
CHARACTER(LEN=32)             :: CGName, FamilyName,ZoneName  ! ?
LOGICAL                       :: onbnd  ! ?
PP_CGNS_INT_TYPE              :: one   ! ?
PP_CGNS_INT_TYPE              :: NormalIndex(3)  ! ?
PP_CGNS_INT_TYPE              :: NormalListSize  ! ?
PP_CGNS_INT_TYPE ,ALLOCATABLE :: BCElems(:,:)  ! ?
PP_CGNS_INT_TYPE ,ALLOCATABLE :: NormalList(:)  ! ?
LOGICAL                       :: zFit
!===================================================================================================================================
ALLOCATE(isize(meshDim,3))
ALLOCATE(DimVec(meshDim,2))
one=1
CALL CG_ZONE_READ_F(CGNSfile,CGNSBase,iZone,CGname,iSize,iError)
ZoneName=CGname
IF(useCurveds.AND.MeshIsAlreadyCurved)THEN 
  N_loc=N
  IF(NBlock.NE.-1) N_loc=NBlock
ELSE
  N_loc=1
END IF
Step=N_loc

irmin=1
IF(meshdim.EQ.3)THEN
  irmax=isize(:,1)
ELSE
  irmax(1:2)=isize(:,1)
  irmax(3)=(N_loc+1)*nSkip
END IF
irmaxorg=irmax
WRITE(UNIT_stdOut,'(A,A,A,3I8)')'Read Zone ',TRIM(CGname),', with block elem size:',irmax(:)-1

DO k=1,meshDim
  IF (MOD((irmax(k)-1),(step)).NE.0) THEN
    IF(useCurveds) THEN 
      WRITE(UNIT_StdOut,'(A)') 'WARNING: cannot read block, step=(order-1)*nSkip does not fit with block elem size.'
    ELSE
      WRITE(UNIT_StdOut,'(A)') 'WARNING: cannot read block, nSkip does not fit with block elem size.'
    END IF
    RETURN
  END IF
END DO

! Allocate list for node Coordinates
IF(meshdim.EQ.3) THEN
  ALLOCATE(NodeCoords(3,1:irmax(1),1:irmax(2),1:irmax(3)))
ELSEIF(meshdim.EQ.2) THEN
  ALLOCATE(nodeCoords(3,1:irmax(1),1:irmax(2),1))
ELSE
  STOP 'Incompatible meshDimension'
END IF

! Read Coordinates
CALL cg_coord_read_f(CGNSFile,CGNSBase,iZone,'CoordinateX',REALDOUBLE,irmin,irmax,NodeCoords(1,:,:,:),iError)
CALL cg_coord_read_f(CGNSFile,CGNSBase,iZone,'CoordinateY',REALDOUBLE,irmin,irmax,NodeCoords(2,:,:,:),iError)
CALL cg_coord_read_f(CGNSFile,CGNSBase,iZone,'CoordinateZ',REALDOUBLE,irmin,irmax,NodeCoords(3,:,:,:),iError)

! Apply skip only in z-DIrection
IF(nSkipZ.NE.1)THEN
  stepk=1
  stepl=1
  stepm=1
  !find z direction
  dir(1,:)=NORMALIZE(NodeCoords(:,2,1,1)-NodeCoords(:,1,1,1),3)
  dir(2,:)=NORMALIZE(NodeCoords(:,1,2,1)-NodeCoords(:,1,1,1),3)
  dir(3,:)=NORMALIZE(NodeCoords(:,1,1,2)-NodeCoords(:,1,1,1),3)
  DO whichdir=1,3
    scalprod=SUM(dir(whichdir,:)*(/0.,0.,1./))
    IF(ABS(scalprod).GT.0.95) EXIT
  END DO
  SELECT CASE(whichDir)
  CASE(1)
    stepk=nSkipZ
    zFit=.NOT.(MOD((irmax(1)-1),(stepk)).NE.0)
  CASE(2)
    stepl=nSkipZ
    zFit=.NOT.(MOD((irmax(2)-1),(stepl)).NE.0)
  CASE(3)
    stepm=nSkipZ
    zFit=.NOT.(MOD((irmax(3)-1),(stepm)).NE.0)
  END SELECT 
  IF(.NOT.zfit)THEN
    IF(useCurveds)THEN
      WRITE(UNIT_StdOut,'(A)') 'WARNING: cannot read block, step=(order-1)*nSkipZ does not fit with block elem size.'
    ELSE
      WRITE(UNIT_StdOut,'(A)') 'WARNING: cannot read block, nSkipZ does not fit with block elem size.'
    END IF
    RETURN
  END IF !zfit

  ! Now apply nSkip in z-dir
  ALLOCATE(NodeCoordsTmp(3,1:((irmax(1)-1)/stepk)+1,1:((irmax(2)-1)/stepl)+1,1:((irmax(3)-1)/stepm)+1))
  DO k=1,irmax(1),stepk; DO l=1,irmax(2),stepl; DO m=1,irmax(3),stepm
    NodeCoordsTmp(:,1+(k-1)/stepk,1+(l-1)/stepl,1+(m-1)/stepm)=NodeCoords(:,k,l,m)
  END DO; END DO; END DO !k,l,m
  DEALLOCATE(NodeCoords)
  irmax=(/SIZE(NodeCoordsTmp,2),SIZE(NodeCoordsTmp,3),SIZE(NodeCoordsTmp,4)/)
  ALLOCATE(NodeCoords(3,irmax(1),irmax(2),irmax(3)))
  NodeCoords=NodeCoordsTmp
  DEALLOCATE(NodeCoordsTmp)
END IF !nSkipZ NE 1

IF(useCurveds.AND.MeshIsAlreadyCurved.AND.N_loc.NE.N)THEN ! change basis to N
  nElems(1:3)=(irmax(1:3)-1)/NBlock
  ALLOCATE(Vdm_NBlock_N(0:N,0:NBlock))
  CALL GetVandermonde(N_loc,'VISU',N,'VISU',Vdm_NBlock_N,modal=.FALSE.)
  IF(meshdim.EQ.3) THEN
    ALLOCATE(NodeCoordsTmp(3,1:nElems(1)*N+1,1:nElems(2)*N+1,1:nElems(3)*N+1))
    DO k=0,nElems(3)-1; DO j=0,nElems(2)-1; DO i=0,nElems(1)-1
      CALL ChangeBasis3D(3,N_loc,N,Vdm_NBlock_N,                   &
                         NodeCoords(   :,i*N_loc+1:(i+1)*N_loc+1,  &
                                         j*N_loc+1:(j+1)*N_loc+1,  &
                                         k*N_loc+1:(k+1)*N_loc+1), &
                         NodeCoordsTmp(:,i*N    +1:(i+1)*N    +1,  &
                                         j*N    +1:(j+1)*N    +1,  &
                                         k*N    +1:(k+1)*N    +1))
    END DO; END DO; END DO
  ELSEIF(meshdim.EQ.2) THEN
    ALLOCATE(NodeCoordsTmp(3,1:nElems(1)*N+1,1:nElems(2)*N+1,1))
    DO j=0,nElems(2)-1; DO i=0,nElems(1)-1
      CALL ChangeBasis2D(3,N_loc,N,Vdm_NBlock_N,                     &
                         NodeCoords(   :,i*N_loc+1:(i+1)*N_loc+1,    &
                                         j*N_loc+1:(j+1)*N_loc+1, 1),&
                         NodeCoordsTmp(:,i*N    +1:(i+1)*N    +1,    &
                                         j*N    +1:(j+1)*N    +1, 1))
    END DO; END DO
  END IF
  DEALLOCATE(Vdm_NBlock_N)
  DEALLOCATE(NodeCoords)
  irmax=(/SIZE(NodeCoordsTmp,2),SIZE(NodeCoordsTmp,3),SIZE(NodeCoordsTmp,4)/)
  ALLOCATE(NodeCoords(3,irmax(1),irmax(2),irmax(3)))
  NodeCoords=NodeCoordsTmp
  DEALLOCATE(NodeCoordsTmp)
  N_loc=N
END IF


! Building temporary nodes
ALLOCATE(Mnodes(irmax(1),irmax(2),irmax(3)))
DO k=1,irmax(1)
  DO l=1,irmax(2)
    DO m=1,irmax(3)
      CALL GetNewNode(Mnodes(k,l,m)%np)
      IF(meshDim.EQ.3)THEN
        Mnodes(k,l,m)%np%x      =NodeCoords(:,k,l,m)      ! Node coordinates are assigned
      ELSE
        Mnodes(k,l,m)%np%x(1:2) =NodeCoords(1:2,k,l,1)    ! Node coordinates are assigned
        Mnodes(k,l,m)%np%x(3)   =REAL(m-1)*REAL(DZ)/REAL(step) ! Node coordinates are assigned       
      END IF
      nNodesGlob=nNodesGlob+1
      Mnodes(k,l,m)%np%ind=nNodesGlob ! Node index is assigned

      ! Marker for Boundary Conditions
      Mnodes(k,l,m)%np%tmp=0
      IF(m.EQ.1.AND.meshdim.EQ.3)        Mnodes(k,l,m)%np%tmp=Mnodes(k,l,m)%np%tmp+1      !zeta minus
      IF(l.EQ.1)                         Mnodes(k,l,m)%np%tmp=Mnodes(k,l,m)%np%tmp+20     !eta minus
      IF(k.EQ.irmax(1) )                 Mnodes(k,l,m)%np%tmp=Mnodes(k,l,m)%np%tmp+300    !xi plus
      IF(l.EQ.irmax(2) )                 Mnodes(k,l,m)%np%tmp=Mnodes(k,l,m)%np%tmp+4000   !eta plus
      IF(k.EQ.1)                         Mnodes(k,l,m)%np%tmp=Mnodes(k,l,m)%np%tmp+50000  !xi minus
      IF(m.EQ.irmax(3).AND.meshdim.EQ.3) Mnodes(k,l,m)%np%tmp=Mnodes(k,l,m)%np%tmp+600000 !zeta plus
    END DO
  END DO
END DO
DEALLOCATE(NodeCoords)

DO k=1,irmax(1)-N,N
  DO l=1,irmax(2)-N,N
    DO m=1,irmax(3)-N,N
      CornerNode(1)%np=>Mnodes(k  ,l  ,m  )%np
      CornerNode(2)%np=>Mnodes(k+N,l  ,m  )%np
      CornerNode(3)%np=>Mnodes(k+N,l+N,m  )%np
      CornerNode(4)%np=>Mnodes(k  ,l+N,m  )%np
      CornerNode(5)%np=>Mnodes(k  ,l  ,m+N)%np
      CornerNode(6)%np=>Mnodes(k+N,l  ,m+N)%np
      CornerNode(7)%np=>Mnodes(k+N,l+N,m+N)%np
      CornerNode(8)%np=>Mnodes(k  ,l+N,m+N)%np
      IF(meshdim.EQ.3)THEN
        CALL GetNewHexahedron(CornerNode)
        CALL CreateSides(FirstElem_in,.TRUE.)
      ELSE
        CALL GetNewQuad(FirstElem_in,CornerNode(1:4))
        CALL CreateSides(FirstElem_in,.TRUE.) ! Transforms 2D element to hexahedron using DZ
        DO kk=1,4
          CALL deleteNode(FirstElem_in%node(kk+4)%np)
          FirstElem_in%node(kk+4)%np => CornerNode(kk+4)%np ! overwrite remaining cornerNodes
        END DO
        CALL CreateSides(FirstElem_in,.FALSE.)
      END IF
      nMeshElems=nMeshElems+1
      FirstElem_in%ind =nMeshElems
      FirstElem_in%zone=nZonesGlob

      IF(useCurveds.AND.MeshIsAlreadyCurved)THEN !read in curvedNodes
        FirstElem_in%nCurvedNodes=(N_loc+1)**3
        ALLOCATE(FirstElem_in%curvedNode(FirstElem_in%nCurvedNodes))
        DO kk=0,N; DO ll=0,N; DO mm=0,N
          FirstElem_in%curvedNode(HexaMapInv(kk,ll,mm))%np=>Mnodes(k+kk,l+ll,m+mm)%np
        END DO; END DO; END DO
      END IF!useCurveds
    END DO !m
  END DO !l
END DO !k
DEALLOCATE(Mnodes)

!------------------ READ BCs ------------------------------------!
! Check for number of Boundary Conditions nCGNSBC 
CALL CG_NBOCOS_F(CGNSFile,CGNSBase,iZone,nCGNSBC,iError) !Number of Boundary conditions nCGNSBC
IF (iError .NE. CG_OK) CALL cg_error_exit_f()

IF (nCGNSBC.LT.1) RETURN ! exit if there are no boundary conditions

ALLOCATE(BCIndex(1:nCGNSBC,6),BCTypeIndex(1:nCGNSBC),countBCs(1:nCGNSBC),nBCFaces(1:nCGNSBC))
SideMap(1,1) = 5 ! xi minus
SideMap(1,2) = 3 ! xi plus
SideMap(2,1) = 2 ! eta minus 
SideMap(2,2) = 4 ! eta plus
SideMap(3,1) = 1 ! zeta minus
SideMap(3,2) = 6 ! zeta plus 

! Read in BC Data
BCIndex=-1
nBCFaces=0
BCTypeIndex=0
countBCs=0
DO iBC=1,nCGNSBC !Loop over all BCs
  CALL CG_BOCO_INFO_F(CGNSfile,CGNSBase,iZone,iBC,CGname,BCTypeI,PntSetType,nBCElems,NormalIndex, &
                      NormalListFlag, DataType,nDataSet,iError)
  IF (iError.NE.CG_OK) CALL  abortCGNS(__STAMP__,CGNSFile)
!  IF (NormalListFlag .NE. 0)THEN
!    CALL closeFile(CGNSFile)
!    CALL abort(__STAMP__,&
!               'Incompatible Boundary definition NormalList')
!  END IF
  IF (INDEX(TRIM(CGName),'_on_CURVES') .GT. 0) CYCLE !ICEM fix, always writes this BC in addition ??????

  CALL CG_GOTO_F(CGNSFile, CGNSBase, iError,'Zone_t',iZone,'ZoneBC_t',1,'BC_t',iBC,'end')
  CALL CG_FAMNAME_READ_F(FamilyName, iError) !Reading Family Name of BC
  IF(iError.NE.CG_OK)THEN  ! if family name is not available
    FamilyName=CGName
  END IF
  BCTypeIndex(iBC)=GetBoundaryIndex(FamilyName)
  IF (BCTypeIndex(iBC).EQ.-1) THEN  
    WRITE(UNIT_stdOut,*)'ERROR - Could not find corresponding boundary definition of ',FamilyName
    CYCLE
  END IF
  ALLOCATE(BCElems(MeshDim,nBCElems))
  NormalListSize=nBCElems*MeshDim
  ALLOCATE(NormalList(NormalListSize))
  CALL CG_BOCO_READ_F(CGNSfile,CGNSBase,iZone,iBC,BCElems,NormalList,iError) 
  IF(PntSetType.EQ.PointRange)THEN
    IF(ANY(BCElems.LE.0))THEN
      WRITE(UNIT_StdOut,'(A)') &
        'WARNING: corrupted pointrange found on Boundary '//TRIM(FamilyName)//' ( '//TRIM(CGName)//', '//TRIM(ZoneName)//' )'
    ELSE
      nBCFaces(iBC)  = 1
      IF(nBCElems.NE.2) STOP 'PointRange has only 2 entries!'
      DO k=1,meshDim
        IF(((BCElems(k,1).NE.irmaxorg(k)).AND.(BCElems(k,1).NE.1)).AND. &
           ((BCElems(k,2).NE.irmaxorg(k)).AND.(BCElems(k,2).NE.1))) THEN
          WRITE(UNIT_StdOut,*)'WARNING: Block face has multiple boundary faces, BoundaryName: ',TRIM(FamilyName)
        END IF
        IF(BCElems(k,1).EQ.BCElems(k,2))&
          BCIndex(iBC,1) = MERGE(SideMap(k,1),SideMap(k,2),BCElems(k,1).EQ.1) ! else irmax(k)
      END DO
      IF (BCindex(iBC,1).EQ.-1) STOP 'ERROR - pointrange does not allow association of BC'
    END IF
  END IF
  IF(PntSetType.EQ.PointList)THEN
    IF(nBCElems.EQ.1) THEN
       WRITE(UNIT_StdOut,*) 'Warning: Single point BC. Zone No,BC no, BCname, ',iZone,iBC,FamilyName
    ELSE
      ! loop through the point list and detect BC Faces
      ! ONLY WORKS WITH BLOCKS WITH irmax > 2 !!
      ! enables multiple block faces for one BC
      nBCFaces(iBC)  = 0
      DO k=1,meshDim
        k2=MAX(3-k,1)
        k3=MIN(3,5-k)  ! cyclic indices
        ! 1 side
        DO l=1,nBCElems
          IF(BCElems(k,l) .EQ. 1) THEN
            IF((BCElems(k2,l).GT. 1) .AND. (BCElems(k2,l).LT.irmaxorg(k2))) THEN
              IF((BCElems(k3,l).GT. 1) .AND. (BCElems(k3,l).LT.irmaxorg(k3))) THEN
                ! Inner Point (not on edge or corner) on Side defines BC
                nBCFaces(iBC)=nBCFaces(iBC)+1
                BCIndex(iBC,nBCFaces(iBC)) = SideMap(k,1)
                EXIT 
              END IF
            END IF
          END IF
        END DO  ! l 
        ! irmax side
        DO l=1,nBCElems
          IF(BCElems(k,l) .EQ. irmaxorg(k)) THEN
            IF((BCElems(k2,l).GT. 1) .AND. (BCElems(k2,l).LT.irmaxorg(k2))) THEN
              IF((BCElems(k3,l).GT. 1) .AND. (BCElems(k3,l).LT.irmaxorg(k3))) THEN
                ! Inner Point (not on edge or corner) on Side defines BC
                nBCFaces(iBC)=nBCFaces(iBC)+1
                BCIndex(iBC,nBCFaces(iBC)) = SideMap(k,2)
                EXIT 
              END IF
            END IF
          END IF
        END DO  ! l 
      END DO ! k

!      isdir=.TRUE.
!      DO k=1,meshDim
!        DO l=2,nBCElems
!          IF(BCElems(k,l).NE.BCElems(k,1)) THEN
!            isdir(k)=.FALSE.
!            EXIT
!          END IF
!        END DO !l
!      END DO ! k
!      IF(isdir(1).AND.(.NOT.isdir(2)) .AND. (.NOT.isdir(3)) ) THEN
!        k=1
!        BCIndex(iBC) = MERGE(SideMap(k,1),SideMap(k,2),BCElems(k,1).EQ.1) ! else irmax(k)
!      ELSEIF((.NOT.isdir(1)) .AND. isdir(2).AND. (.NOT.isdir(3)) ) THEN
!        k=2
!        BCIndex(iBC) = MERGE(SideMap(k,1),SideMap(k,2),BCElems(k,1).EQ.1) ! else irmax(k)
!      ELSEIF((.NOT.isdir(1)) .AND. (.NOT.isdir(2)).AND. isdir(3) ) THEN
!        k=3
!        BCIndex(iBC) = MERGE(SideMap(k,1),SideMap(k,2),BCElems(k,1).EQ.1) ! else irmax(k)
!      ELSEIF(.NOT.ANY(isdir(:))) THEN
!       WRITE(*,*)'Zone,BC,no BCs',ZoneName,CGName,nCGNSBC
!       STOP 'no BC direction found'
!      ELSE 
!       !last possible case: edge where two BC meet.
!        BCIndex(iBC)=-1
!       WRITE(*,*) 'Warning: Edge BC. Zone No,BC no, BCname, ',iZone,iBC,FamilyName
!      END IF
    END IF
  END IF
  DEALLOCATE(BCElems,NormalList)
END DO !iBC=1,nCGNSBC

aElem=>FirstElem_in
DO WHILE(ASSOCIATED(aElem))
  IF (aElem%Zone.NE.iZone)THEN
    aElem=>aElem%nextElem
    CYCLE
  END IF
  aSide=>aElem%FirstSide
  DO WHILE(ASSOCIATED(aSide))
    DO iBC=1,nCGNSBC                ! nCGNSBC Boundary Conditions 
      DO iBCFace=1,nBCFaces(iBC)    ! nBCFaces: number of block faces with this BC
        iSide=BCIndex(iBC,iBCFace) 
        IF(iSide.LT.0) CYCLE
        onBnd=.TRUE.
        DO l=1,4 !Loop over all Nodes of this side
          IF(INT(MOD(aSide%Node(l)%np%tmp,10**iSide)/(10**(iSide-1))).NE.iSide) THEN
            onBnd=.FALSE.
            EXIT !Loop
          END IF
        END DO !l=1,4
        IF(onBnd) THEN
          !WRITE(*,*) 'BCi',iBC
          IF(BoundaryType(BCTypeIndex(iBC),1).EQ.0) EXIT  !ignore BCType=0 (DEFAULT_SURFACES)
          CALL getNewBC(aSide%BC)
          aSide%BC%BCType    =BoundaryType(BCTypeIndex(iBC),1)
          aSide%CurveIndex   =BoundaryType(BCTypeIndex(iBC),2)
          aSide%BC%BCstate   =BoundaryType(BCTypeIndex(iBC),3)
          aSide%BC%BCalphaInd=BoundaryType(BCTypeIndex(iBC),4)
          aSide%BC%BCIndex   =BCTypeIndex(iBC)
          countBCs(iBC) = countBCs(iBC)+1
        END IF
      END DO !iBCFace=1,nBCFaces(iBC)
    END DO !iBC=1,nCGNSBC 
    aSide=>aSide%nextElemSide
  END DO !WHILE(ASSOCIATED(aSide))
  aElem=>aElem%nextElem
END DO !WHILE(ASSOCIATED(aElem))
DEALLOCATE(BCIndex,BCTypeIndex,countBCs,nBCFaces)
END SUBROUTINE ReadCGNSMeshStruct



SUBROUTINE ReadCGNSSurfaceMesh(FirstElem_in,FileName)
!===================================================================================================================================
! This subroutine reads the surface data from an unstructured CGNS file and prepares the element list.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tElemPtr,tSide
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewNode,getNewSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(INOUT)                  :: FirstElem_in
CHARACTER(LEN=255),INTENT(IN)                   :: FileName
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
PP_CGNS_INT_TYPE             :: iZone
PP_CGNS_INT_TYPE             :: nBases,nCGNSZones  ! Number of bases / zones in CGNS file
PP_CGNS_INT_TYPE             :: nNodesGlob         ! Total number of nodes in the mesh file
PP_CGNS_INT_TYPE             :: nElemsGlob         ! Total number of elements in mesh file
PP_CGNS_INT_TYPE             :: nZonesGlob         ! Total number of zones in the mesh.
PP_CGNS_INT_TYPE             :: CGNSFile,CGNSBase  ! CGNS file handles
PP_CGNS_INT_TYPE             :: md  ! ?
PP_CGNS_INT_TYPE             :: file_type  ! ?
REAL(KIND=4)                 :: version  ! ?
                             
CHARACTER(LEN=32)            :: CGName             ! necessary data for CGNS
PP_CGNS_INT_TYPE             :: ZoneType  ! ?
TYPE(tElemPtr),ALLOCATABLE   :: Elems(:)                            ! Pointer array to elements
TYPE(tElem),POINTER          :: ActualElem  ! ?
TYPE(tSide),POINTER          :: Side                                ! Local side pointer
REAL, ALLOCATABLE            :: NodeCoords(:,:)                     ! Coordinates of grid nodes (nDim,nNodes)
PP_CGNS_INT_TYPE ,ALLOCATABLE:: ElemConnect(:,:)                    ! Connectivity data for (surface) element (nNodes+1,nElems)
PP_CGNS_INT_TYPE ,ALLOCATABLE:: LocalConnect(:), ParentData(:,:)    ! Arrays for reading the (local=one section) connectivity data
PP_CGNS_INT_TYPE             :: iNode,iElem,iSecElem,iSide  ! ?
PP_CGNS_INT_TYPE             :: dm, iSect, j  ! ?
PP_CGNS_INT_TYPE             :: IndMin, IndMax, iStart, iEnd        ! Start and end index
PP_CGNS_INT_TYPE             :: nSect                               ! Number of sections in CGNS file
PP_CGNS_INT_TYPE             :: nNodes,nElems,nSectElems            ! Number of nodes / elements in a zone (volume only)
PP_CGNS_INT_TYPE             :: nElemNodes,nBCElems                 ! Number of nodes of an element
PP_CGNS_INT_TYPE             :: LocDim, LocType, nNodesLoc          ! Dimension, type, number of nodes of local=section elements
PP_CGNS_INT_TYPE             :: SizeZone(3)                         ! CGNS datastructure variables
PP_CGNS_INT_TYPE             :: SectionElemType                     ! Type of elements in CGNS file
PP_CGNS_INT_TYPE             :: ParentDataFlag                      ! 0=no parent data for elems available, 1=parent data available
PP_CGNS_INT_TYPE             :: CellDim, PhysDim                    ! Dimesnion of elements,physical dimension
PP_CGNS_INT_TYPE             :: iError                              ! Error flag
CHARACTER(LEN=30)            :: coordNameCGNS(3)                 ! List of CGNS names for the coordinates
PP_CGNS_INT_TYPE             :: one  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,*)'Read CGNS File: ',TRIM(FileName)
! Open CGNS file
CALL OpenBase(TRIM(FileName),MODE_READ,md,md,CGNSFile,CGNSBase,.TRUE.)
!CALL CG_OPEN_F(TRIM(MeshFileName(iFile)), CG_MODE_READ, CGNSFile, iError)
CALL CG_VERSION_F(CGNSFile, version, iError)
WRITE(UNIT_stdOut,*)'CGNS version:',version
CALL CG_IS_CGNS_F(TRIM(FileName), file_type, iError)

! Get number of bases in CGNS file
CALL CG_NBASES_F(CGNSfile,nBases,iError)
IF (iError .NE. CG_OK) &
  CALL abortCGNS(__STAMP__,CGNSFile)
IF(nBases .GT. 1) WRITE(UNIT_stdOut,*)'WARNING - Found ',nBases,' bases in CGNS file. Reading only base 1!'
CGNSBase=1  ! Multiple bases not supported at the moment

! Check dimensions of CGNS base
CALL CG_BASE_READ_F(CGNSfile,CGNSBase,CGname,CellDim,PhysDim,iError)
IF(iError .NE. CG_OK) &
  CALL abortCGNS(__STAMP__,CGNSFile)

! Get number of zones in CGNSBase
CALL CG_NZONES_F(CGNSfile,CGNSBase,nCGNSZones,iError)
IF(iError .NE. CG_OK) CALL abortCGNS(__STAMP__,CGNSFile)

nNodesGlob=0
nZonesGlob=0
nElemsGlob=0
nZonesGlob=0
DO iZone=1,nCGNSZones
  nZonesGlob=nZonesGlob+1
  ! Check structured / unstructured
  CALL cg_zone_type_f(CGNSFile, CGNSBase, iZone, ZoneType, iError)
  IF (iError .NE. CG_OK) CALL cg_error_exit_f()
  IF (ZoneType.EQ.Structured)THEN
    STOP 'no structured readin for surface data'
  END IF
  coordNameCGNS(1) = 'CoordinateX'
  coordNameCGNS(2) = 'CoordinateY'
  coordNameCGNS(3) = 'CoordinateZ'
  one=1
  ! Check dimensions of CGNS base
  CALL CG_BASE_READ_F(CGNSfile,CGNSBase,CGname,CellDim,PhysDim,iError)
  IF(iError .NE. CG_OK) &
    CALL abortCGNS(__STAMP__,CGNSFile)
  ! Start with reading zones: total number of Nodes and Elems
  CALL CG_ZONE_READ_F(CGNSfile,CGNSBase,iZone,CGname,SizeZone,iError)
  IF (iError .NE. CG_OK) &
    CALL abortCGNS(__STAMP__,CGNSFile)
  WRITE(UNIT_stdOut,*)'Read Zone ',TRIM(CGname)
  
  ! Read node coordinates
  nNodes=SizeZone(1)
  ALLOCATE(NodeCoords(3,nNodes))
  NodeCoords=0.
  DO dm=1,3
    CGname=TRIM(CoordNameCGNS(dm))
    CALL CG_COORD_READ_F(CGNSfile,CGNSBase,iZone,CGName,RealDouble,one,nNodes,NodeCoords(dm,:),iError)
    IF (iError .NE. CG_OK)THEN
      WRITE(UNIT_stdOut,*)'ERROR - Could not read coordinate(',dm,'): ',TRIM(CoordNameCGNS(dm))
      CALL CG_NCOORDS_F(CGNSFile,CGNSBase,iZone,PhysDim,iError )  ! Here we use PhysDim as nCoords
      WRITE(UNIT_stdOut,*)'Number of coordinates in CGNS file: ',PhysDim
      WRITE(UNIT_stdOut,*)'Mesh dimension: ',3
      DO j=1,PhysDim
        CALL CG_COORD_INFO_F(CGNSFile,CGNSBase,iZone,j,LocType,CGName,iError)
        WRITE(UNIT_stdOut,*)'Coordname (',j,')=',TRIM(CGName)
        WRITE(UNIT_stdOut,*)'Datatype=',LocType
      END DO
      CALL abortCGNS(__STAMP__,CGNSFile)
    END IF
  END DO
  ! Read in CGNS sections and count volume and boundary face elems
  CALL CG_NSECTIONS_F(CGNSfile,CGNSBase,iZone,nSect,iError)
  IF (iError .NE. CG_OK) CALL abortCGNS(__STAMP__,CGNSFile)
  nElems=SizeZone(2)
  
  
  ! Read element connectivity
  ALLOCATE(ElemConnect(5,nElems))
  ALLOCATE(Elems(nElems))
  
  iElem=0
  DO iSect=1,nSect ! Vol. and Face elems
    ! Read in Elem indMin & indMax
    CALL CG_SECTION_READ_F(CGNSfile,CGNSBase,iZone,iSect,CGname,SectionElemType,IndMin,IndMax,nBCElems,ParentDataFlag,iError)
    WRITE(UNIT_StdOut,*)'   read section',TRIM(CGname)
    IF (iError .NE. CG_OK) CALL abortCGNS(__STAMP__,CGNSFile)
    IF(SectionElemType .LT. TRI_3) CYCLE !ignore additional sections with data <nDim-1
    CALL CG_ELEMENTDATASIZE_F(CGNSFile,CGNSBase,iZone,iSect,nSectElems,iError)  ! Get number of connectivity values
    ALLOCATE(LocalConnect(nSectElems))
    ALLOCATE(ParentData(nSectElems,4))
    ! Read in local connectivity data
    CALL CG_ELEMENTS_READ_F(CGNSfile,CGNSBase,iZone,iSect,LocalConnect,ParentData,iError)
  
    nSectElems=1+IndMax-IndMin ! Important for surface elements only 
                               ! (nSectElems, Parent1 | Parent2 | ParentSide1 | ParentSide2)...but we don't use it
    iStart=1
    DO iSecElem=1,nSectElems
      iEnd=iStart
      IF(SectionElemType .EQ. MIXED) THEN
        LocType=LocalConnect(iStart)  ! Mixed type elements, read elem type from array
        iStart =iStart+1              ! First value is elem type
      ELSE
        LocType=SectionElemType       ! Single type elements, elem type = SectionElemType
        iEnd   =iEnd-1                ! Only nElemNodes values
      END IF
      CALL CG_NPE_F(LocType,nNodesLoc,iError) ! Get number of nodes for iElem
      iEnd=iEnd+nNodesLoc
  
      LocDim=1
      IF(LocType .GT.  BAR_3) LocDim=2
      IF(LocType .GT. QUAD_9) LocDim=3
      IF(LocType .GT.  MIXED) LocDim=2 !NGON_n
  
      IF(LocDim .EQ. 2) THEN ! surface element
        iElem=iElem+1
        IF(iElem.GT.nElems)THEN
          CALL closeFile(CGNSFile)
          CALL abort(__STAMP__,&
                         'Something wrrrrong with surf element numbers in CGNS File zone :',INT(iZone))
        END IF
        
        ElemConnect(1            ,iElem)=LocType
        ElemConnect(2:nNodesLoc+1,iElem)=LocalConnect(iStart:iEnd)
  
      END IF   
      iStart=iEnd+1
    END DO ! elements in section
    DEALLOCATE(LocalConnect,ParentData)
  END DO !sections
  
  ! Rebuild the elements of zone iZone
  DO iElem=1,nElems
    CALL CG_NPE_F(ElemConnect(1,iElem),nElemNodes,iError)
    IF (iError .NE. CG_OK) CALL abortCGNS(__STAMP__,CGNSFile)
    CALL GetNewElem(Elems(iElem)%EP)
    ActualElem=>Elems(iElem)%ep
    ActualElem%Zone  =nZonesGlob
    ActualElem%Ind   =iElem+nElemsGlob
    ActualElem%nNodes=nElemNodes
    ALLOCATE(ActualElem%Node(nElemNodes))
    DO iNode=1,nElemNodes
      CALL GetNewNode(ActualElem%Node(iNode)%NP)
      ActualElem%Node(iNode)%NP%Ind     =ElemConnect(iNode+1,iElem)+nNodesGlob  ! Make node indices unique
      ActualElem%Node(iNode)%NP%x       =NodeCoords(:,ElemConnect(iNode+1,iElem))
      ActualElem%Node(iNode)%NP%RefCount=1
    END DO
    ! buildSides
    DO iSide=1,nElemNodes
      IF (iSide .EQ. 1) THEN
        CALL getNewSide(ActualElem%firstSide,2)
        Side=>ActualElem%firstSide
      ELSE
        CALL getNewSide(Side%nextElemSide,2)  
        Side=>Side%nextElemSide
      END IF
      Side%LocSide=iSide
      Side%Elem=>ActualElem
      DO iNode=1,2
        Side%Node(iNode)%np=>ActualElem%Node(MOD(iSide+iNode-2,nElemNodes)+1)%np
        Side%Node(iNode)%np%refCount=Side%Node(iNode)%np%refCount+1
      END DO
    END DO
  END DO !nElems
  DEALLOCATE(ElemConnect,NodeCoords)
  
  DO iElem=1,nElems
    IF(.NOT. ASSOCIATED(FirstElem_in))THEN
      FirstElem_in   => Elems(iElem)%EP
    ELSE
      ActualElem                   => Elems(iElem)%EP
      ActualElem%nextElem          => FirstElem_in
      ActualElem%nextElem%prevElem => ActualElem
      FirstElem_in                    => ActualElem
    END IF
    NULLIFY(FirstElem_in%prevElem)
  END DO !iElem=1,nElems
  DEALLOCATE(Elems)
  nNodesGlob = nNodesGlob+nNodes
  nElemsGlob = nElemsGlob+nElems
END DO ! iZone

END SUBROUTINE ReadCGNSSurfaceMesh


SUBROUTINE openBase(ioName,mode,celldim,physdim,CGNSFile,CGNSBase,externBase_in)
!===================================================================================================================================
! Opens CGNS files and, if they do not yet exist, prepares their datastructure
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: ioName        ! name for/of the CGNS file
PP_CGNS_INT_TYPE,INTENT(IN)    :: mode          ! CGNS file mode: MODE_READ,MODE_MODIFY
PP_CGNS_INT_TYPE,INTENT(IN)    :: celldim       ! dimension of the cells: 3 for volume cells, 2 for surface cells
PP_CGNS_INT_TYPE,INTENT(IN)    :: physdim       ! physical dimension of the calculation
LOGICAL,OPTIONAL,INTENT(IN)    :: externBase_in ! base specified from external
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
PP_CGNS_INT_TYPE,INTENT(OUT)    :: CGNSFile      ! CGNS file handle
PP_CGNS_INT_TYPE,INTENT(OUT)    :: CGNSBase      ! CGNS base handle
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
PP_CGNS_INT_TYPE    :: nB,i !,j        ! loop/count variables
PP_CGNS_INT_TYPE    :: cd,pd           ! CGNS subroutine variables
PP_CGNS_INT_TYPE    :: ier             ! CGNS error variable
CHARACTER(LEN=100)  :: filename        ! necessary data for CGNS
CHARACTER(LEN=32)   :: CGname,basename ! necessary data for CGNS
LOGICAL             :: externBase      ! provide CGNS base from external
PP_CGNS_INT_TYPE    :: one   ! ?
!===================================================================================================================================
! read/write + force creation if not existing
one=1
externBase=.FALSE.
IF(PRESENT(externBase_in)) externBase=externBase_in
CGNSBase=-999
CGNSFile=-999
basename='NSDG_database'
IF(.NOT. externBase) THEN
  filename=TRIM(ioName)//'.cgns'
ELSE
  filename=TRIM(ioName)
END IF

! Open the file
CALL cg_open_f(filename,MODE_READ,CGNSfile,ier)
IF (ier .NE. CG_OK) THEN
  IF(mode .EQ. MODE_READ) RETURN
  IF(mode .EQ. MODE_MODIFY) THEN ! force file creation, since the file does not yet exist
    CALL cg_open_f(filename,MODE_WRITE,CGNSfile,ier)
    IF (ier .NE. CG_OK) &
      CALL abortCGNS(__STAMP__,CGNSFile)
    CALL closeFile(CGNSFile)
    CALL cg_open_f(filename,mode,CGNSfile,ier)
    IF (ier .NE. CG_OK) &
      CALL abortCGNS(__STAMP__,CGNSFile)
  END IF
END IF

! Find out how many bases we have (should be just one)
CALL cg_nbases_f(CGNSfile,nB,ier)
IF (ier .NE. CG_OK) &
  CALL abortCGNS(__STAMP__,CGNSFile)
IF(nB .GT. 1) THEN
  WRITE(UNIT_StdOut,*)'Warning, more than one bases in file : ',TRIM(fileName),' taking only first one'
END IF

IF(externBase) THEN
  CALL cg_base_read_f(CGNSfile,one,CGname,cd,pd,ier)
  IF (ier .NE. CG_OK) &
    CALL abortCGNS(__STAMP__,CGNSFile)
  CGNSbase=1
ELSE
  DO i=1,nB
    CALL cg_base_read_f(CGNSfile,i,CGname,cd,pd,ier)
    IF (ier .NE. CG_OK) &
      CALL abortCGNS(__STAMP__,CGNSFile)
    IF (CGname .EQ. basename) THEN
      CGNSbase=i
      IF( (cd .NE. cellDim) .OR. (pd .NE. physDim))THEN
        CALL closeFile(CGNSFile)
        CALL abort(__STAMP__,&
           'Wrong dimensionalities in database ' // TRIM(basename) // ' in file:'//TRIM(filename),INT(cellDim),REAL(physDim))
      END IF
    END IF
  END DO
END IF

! If no base present, the file is new, so we create the right one now
IF (CGNSBase .EQ. -999) THEN ! not found
  IF(mode .EQ. MODE_READ)THEN
     CALL closeFile(CGNSFile)
     CALL abort(__STAMP__,& 
                'Cannot find base '//TRIM(basename)//' in cgns file: '//TRIM(filename))
  END IF

  IF((mode .EQ. MODE_MODIFY) .OR. (mode.EQ.MODE_WRITE)) THEN ! create base anyway
    CALL cg_base_write_f(CGNSfile,basename,celldim,physdim,CGNSBase,ier)
    IF (ier .NE. CG_OK) &
      CALL abortCGNS(__STAMP__,CGNSFile)
  END IF
END IF
END SUBROUTINE openBase


SUBROUTINE closeFile(CGNSFile)
!===================================================================================================================================
! Aborts the current running calculation in a controlled way
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
PP_CGNS_INT_TYPE ,INTENT(IN)  :: CGNSFile   ! CGNS file handle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
PP_CGNS_INT_TYPE    :: ier  ! ?
!===================================================================================================================================
! close CGNS file first (if still open)
CALL cg_close_f(CGNSfile,ier)
IF (ier .NE. CG_OK) CALL cg_error_print_f()
END SUBROUTINE closeFile


SUBROUTINE abortCGNS(sourceFile,sourceLine,compDate,compTime,CGNSFile)
!===================================================================================================================================
! Provides a controlled code abort in case of an CGNS error
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: sourceFile ! name of the file that caused an error
INTEGER,INTENT(IN)             :: sourceLine ! number of line in file that caused an error
CHARACTER(LEN=*),INTENT(IN)    :: compDate   ! date of computation
CHARACTER(LEN=*),INTENT(IN)    :: compTime   ! computation time
PP_CGNS_INT_TYPE,INTENT(IN)    :: CGNSFile   ! CGNS file handle
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)  :: message    ! CGNS error message
!===================================================================================================================================
! Retrieve error message
CALL cg_get_error_f(message)      
! Exit calculation
CALL closeFile(CGNSFile)
CALL abort(sourceFile,sourceLine,compDate,compTime,message)
END SUBROUTINE abortCGNS

END MODULE MOD_Readin_CGNS

