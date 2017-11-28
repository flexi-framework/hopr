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
#include "hopr.h"
MODULE MOD_Search
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Search_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE getNewSearchMesh
   MODULE PROCEDURE getNewSearchMesh
END INTERFACE

INTERFACE getFirstToObject
   MODULE PROCEDURE getFirstToObject
END INTERFACE

INTERFACE GetNextToObject
   MODULE PROCEDURE getNextToObject
END INTERFACE

INTERFACE deleteSearchMesh
   MODULE PROCEDURE deleteSearchMesh
END INTERFACE

INTERFACE insertNode
   MODULE PROCEDURE insertNode
END INTERFACE

INTERFACE insertUniqueNode
   MODULE PROCEDURE insertUniqueNode
END INTERFACE

INTERFACE getUniqueNode
   MODULE PROCEDURE getUniqueNode
END INTERFACE

INTERFACE insertAllMeshSides
   MODULE PROCEDURE insertAllMeshSides
END INTERFACE

INTERFACE insertSide
   MODULE PROCEDURE insertSide
END INTERFACE

INTERFACE flushSearchMesh
   MODULE PROCEDURE flushSearchMesh
END INTERFACE
  
INTERFACE idxok
   MODULE PROCEDURE idxok
END INTERFACE

INTERFACE getIdx
   MODULE PROCEDURE getIdx
END INTERFACE

INTERFACE deleteToObjectList
   MODULE PROCEDURE deleteToObjectList
END INTERFACE

INTERFACE insertToObject
   MODULE PROCEDURE insertToObject
END INTERFACE

INTERFACE insertNextToObject
   MODULE PROCEDURE insertNextToObject
END INTERFACE

!===================================================================================================================================

CONTAINS
SUBROUTINE InitSearch()
!===================================================================================================================================
! Read search mesh parameters from ini file and initialize search mesh variables
!===================================================================================================================================
! MODULES
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT SEARCH...'
nElemsNodeSearch=GETINT('nElemsNodeSearch','25')
RefineSideSearch=GETREAL('RefineSideSearch','8.')
WRITE(UNIT_stdOut,'(A)')' INIT SEARCH DONE!'
WRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitSearch


SUBROUTINE getNewSearchMesh(searchMesh,safedx,xmin,xmax,dx)
!===================================================================================================================================
! Create new Cartesian search mesh and define search area around an arbitrary point
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:ConformConnect
USE MOD_Mesh_Vars,ONLY:maxDX,SpaceQuandt
USE MOD_Mesh_Vars,ONLY:useCurveds
USE MOD_Mesh_Tolerances,ONLY:SAMEPOINT
USE MOD_SortingTools,ONLY:QSort1Int1Pint
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
LOGICAL,INTENT(IN)                    :: safedx            ! .TRUE.: calculate dx from mesh geometry, .FALSE.: dx = (xmax-xmin)/50.
REAL,INTENT(IN),OPTIONAL              :: xmin(3)  ! Search mesh extents
REAL,INTENT(IN),OPTIONAL              :: xmax(3)  ! Search mesh extents
REAL,INTENT(IN),OPTIONAL              :: dx                ! specify dx manually 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tSearchMesh),POINTER,INTENT(OUT) :: searchMesh        ! New search mesh
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER       :: Elem                          ! Local element pointer
REAL                      :: deltaSafe(3)                    ! ?
INTEGER,ALLOCATABLE       :: dist(:)                         ! ?
INTEGER                   :: iNode,l,i,j,k                   ! ?
!===================================================================================================================================
IF(ASSOCIATED(searchMesh)) CALL abort(__STAMP__,&
                                     'Warning: search mesh already allocated.')
ALLOCATE(searchMesh)
searchMesh%actualInd   = -1
searchMesh%idx_ur      = -1
searchMesh%idxCounter  = -1
searchMesh%nMesh       = 1   ! Number of mesh points for search mesh
searchMesh%objectCount = 0
NULLIFY(searchMesh%actualToObject)
! searchMesh%nSearch defines the Spirale des Todes
! default values used for node search
searchMesh%nSearch         = 1  ! Initialize first independent of mesh dimension
IF(.NOT.ConformConnect)searchMesh%nSearch(1:3) = nElemsNodeSearch
! Get search mesh extents if not given 
IF(.NOT. PRESENT(xmin)) THEN
  ! Calculate from mesh geometry
  searchMesh%xmin = 1.e12
  searchMesh%xmax = -1.e12
  Elem=>firstElem
  DO WHILE(ASSOCIATED(Elem))
    DO iNode=1,Elem%nNodes
      searchMesh%xmin=MIN(Elem%Node(iNode)%np%x,searchMesh%xmin)
      searchMesh%xmax=MAX(Elem%Node(iNode)%np%x,searchMesh%xmax)
    END DO
    IF(useCurveds)THEN
      DO iNode=1,Elem%nCurvedNodes
        searchMesh%xmin=MIN(Elem%curvedNode(iNode)%np%x,searchMesh%xmin)
        searchMesh%xmax=MAX(Elem%curvedNode(iNode)%np%x,searchMesh%xmax)
      END DO
    END IF
    Elem=>Elem%nextElem
  END DO
ELSE
  ! Use given arguments
  searchMesh%xmin=xmin
  searchMesh%xmax=xmax
END IF

! Calculate dx for search mesh
IF (safedx) THEN !determine/estimate dx from mesh geometry, used for side search, i.e. rasterfahndung/connect 
  searchMesh%dx=maxDX
ELSE ! use nElemsNodeSearch to determine dx, used for node search, i.e. getuniquenode, at moment: homogenous in each direction
  searchMesh%dx=(searchMesh%xmax-searchMesh%xmin)/(2.*nElemsNodeSearch)
  IF(PRESENT(dx)) searchMesh%dx=dx
END IF
IF(ConformConnect)searchMesh%dx=MAX(searchMesh%dx,1.1*SpaceQuandt*PP_MeshTolerance) !  !

! Numerical tolerances: enlarge slightly for numerical security
deltaSafe        = 0.00001*(searchMesh%xmax-searchMesh%xmin)  
searchMesh%xmin  = searchMesh%xmin-deltaSafe
searchMesh%xmax  = searchMesh%xmax+deltaSafe
searchMesh%dx    = 1.00001*searchMesh%dx
searchMesh%Mid   = 0.5*(searchMesh%xmin+searchMesh%xmax)
searchMesh%Radius= 0.5*SQRT(SUM((searchMesh%xmax-searchMesh%xmin)**2))*1.01  ! 10% security
IF (safedx) THEN
  ! Introduce Refine Factor for side search mesh
  searchMesh%dx = searchMesh%dx/RefineSideSearch
  ! Determine the new size of the Spirale des Todes
  searchMesh%nSearch(1:3)= INT(maxDX/searchMesh%dx)+1 !+1 for security 
ENDIF 
! Determine the cartesian search mesh
searchMesh%nMesh(1:3) = INT(MAX(1.,(searchMesh%xmax-searchMesh%xmin)/searchMesh%dx))+1
! Number of Spirale des Todes links
searchMesh%n    = PRODUCT((2*searchMesh%nSearch+1))  
! Define Spirale des Todes. If Search is not successfull in the determined i,j,k grid cell, then 
!  proceed search in neighbour grid cells following the Spirale des Todes
ALLOCATE(searchMesh%mapping(searchMesh%n,3),dist(searchMesh%n),searchMesh%sortidx(searchMesh%n))
! This corresponds to the search area, not to the mesh
l=0
DO i=-searchMesh%nSearch(1),searchMesh%nSearch(1)   
  DO j=-searchMesh%nSearch(2),searchMesh%nSearch(2)  
    DO k=-searchMesh%nSearch(3),searchMesh%nSearch(3)  
      l=l+1
      ! Map search mesh indices to one-dimensional search index.
      searchMesh%mapping(l,:)=(/i,j,k/)
      dist(l)=SUM(searchMesh%mapping(l,:)**2)
      searchMesh%sortidx(l)=l
    END DO
  END DO
END DO
! Sort mapping by distance "dist" in ascending order
CALL Qsort1Int1PInt(dist,searchMesh%sortidx)
DEALLOCATE(dist)
! Reduced Spirale des Todes, covers only the DIRECT neighbour grid cells
searchMesh%nred=27 
! Allocate and initialize search mesh
ALLOCATE(searchMesh%sm(searchMesh%nMesh(1),searchMesh%nMesh(2),searchMesh%nMesh(3)))
DO i=1,searchMesh%nMesh(1)  ! Search mesh
  DO j=1,searchMesh%nMesh(2)  ! Search mesh
    DO k=1,searchMesh%nMesh(3)  ! Search mesh
      NULLIFY(searchMesh%sm(i,j,k)%ToObject)
    END DO
  END DO
END DO
END SUBROUTINE getNewSearchMesh


FUNCTION getfirstToObject(searchMesh,reduced,coords,idx)
!===================================================================================================================================
! Get first search mesh object around a point with coordinates "coords" or search mesh indices "idx"
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSearchMesh),POINTER,INTENT(INOUT) :: searchMesh   ! ?
LOGICAL,INTENT(IN)                   :: reduced       ! Use reduced search area (faster / less secure)
REAL,INTENT(IN),OPTIONAL             :: coords(3)     ! Point coordinates
INTEGER,INTENT(IN),OPTIONAL          :: idx(3)        ! Point indices
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tToObject),POINTER   :: getfirstToObject ! First search mesh object at point with "coord" / "idx"
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                   :: maxn   ! ?
!===================================================================================================================================
searchMesh%actualInd=1  ! Reset
IF(PRESENT(idx)) THEN
  searchMesh%actualInd=idx  ! Use given indices
ELSEIF(PRESENT(coords)) THEN
  CALL getIdx(searchmesh,coords,searchMesh%actualInd)  ! Calculate indices from coordinates
ELSE
  CALL abort(__STAMP__, & 
    'error calling getFirstToObject')
END IF
maxn=searchMesh%n  ! Number of points in search area
IF(Reduced) maxn=searchMesh%nRed  ! Use reduced number of points
NULLIFY(getFirstToObject)
IF(.NOT.idxok(searchMesh,searchMesh%actualInd)) THEN
  ! Point is not located within search mesh
  searchMesh%actualInd  = -1
  NULLIFY(searchMesh%actualToObject)
  searchMesh%idxCounter = -1
ELSE
  searchMesh%idxCounter = 0
  searchMesh%idx_ur     = searchMesh%actualInd  ! Current search index used as starting index
  DO WHILE((.NOT. ASSOCIATED(getFirstToObject)) .AND. (searchMesh%idxCounter .LT. maxn))
    searchMesh%idxCounter = searchMesh%idxCounter + 1  ! Next element in search area around starting point
    searchMesh%actualInd(1:3) = searchMesh%idx_ur(1:3)+searchMesh%mapping(searchMesh%sortidx(searchMesh%idxCounter),1:3)
    IF(idxok(searchmesh,searchMesh%actualInd))THEN
      getFirstToObject => searchMesh%sm(searchMesh%actualInd(1),searchMesh%actualInd(2),searchMesh%actualInd(3))%ToObject
    END IF
  END DO
  IF(ASSOCIATED(getFirstToObject)) THEN
    searchMesh%actualToObject=>getFirstToObject
  ELSE
    NULLIFY(searchMesh%actualToObject)
  END IF
END IF
END FUNCTION getfirstToObject


FUNCTION getNextToObject(searchMesh,reduced)
!===================================================================================================================================
! Get next search mesh object
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! searchMesh%idx_ur           : Starting point of search
! searchMesh%actualInd        : Indices of current search mesh point
! searchMesh%idxCounter       : Current point in search area around point "searchMesh%idx_ur"
! searchMesh%sortidx          : Search area index of current point in search area
! searchMesh%mapping          : Mapping from search area point index to search mesh point indices
! SearchMesh%actualToObject   : Current search mesh object
TYPE(tSearchMesh),POINTER,INTENT(IN) :: searchMesh
LOGICAL,INTENT(IN)        :: reduced         ! Reduced search area (faster / less secure)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tToObject),POINTER   :: getNextToObject ! Next object in search mesh
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                   :: idx(3),maxn   ! ?
!===================================================================================================================================
idx  = 1             ! Reset indices
maxn = searchMesh%n  ! Number of points in search area
IF(Reduced) maxn=searchMesh%nRed  ! Reduced search area
NULLIFY(getNextToObject)
IF(ASSOCIATED(SearchMesh%actualToObject)) THEN
  IF(ASSOCIATED(searchMesh%actualToObject%nextToObject)) THEN
    getNextToObject=>searchMesh%actualToObject%nextToObject  ! Simply use next search mesh object in list
  ELSE
    DO WHILE((.NOT. ASSOCIATED(getNextToObject)) .AND. (searchMesh%idxCounter .LT. maxn))
      searchMesh%idxCounter = searchMesh%idxCounter + 1  ! Next point in search area
      idx(1:3) = searchMesh%idx_ur(1:3) + searchMesh%mapping(searchMesh%sortidx(searchMesh%idxCounter),1:3)
      IF(idxok(searchMesh,idx)) THEN
        getNextToObject => searchMesh%sm(idx(1),idx(2),idx(3))%ToObject
        searchMesh%actualInd = idx
      END IF
    END DO
    IF(.NOT. ASSOCIATED(getNextToObject))THEN  ! No next object found
      searchMesh%idxCounter = -1
      searchMesh%actualInd  = -1
    END IF
  END IF
END IF
searchMesh%actualToObject=>getNextToObject

END FUNCTION getNextToObject


SUBROUTINE flushSearchMesh(searchMesh,withObjects)
!===================================================================================================================================
! Delete search mesh objects and associated mesh objects (if withObjects = .TRUE.)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSearchMesh),POINTER,INTENT(INOUT) :: searchMesh  ! ?
LOGICAL,INTENT(IN)                      :: withObjects ! Delete mesh objects (nodes, sides, elements) too
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER                   :: i,j,k   ! ?
!===================================================================================================================================
DO i=1,searchMesh%nMesh(1)
  DO j=1,searchMesh%nMesh(2)
    DO k=1,searchMesh%nMesh(3)
      CALL deleteToObjectList(searchMesh%sm(i,j,k)%ToObject,withObjects)
    END DO
  END DO
END DO
END SUBROUTINE flushSearchMesh


SUBROUTINE deleteSearchMesh(searchMesh,withObjects)
!===================================================================================================================================
! Delete "search mesh" and associated mesh objects (if withObjects =.TRUE.)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSearchMesh),POINTER,INTENT(INOUT) :: searchMesh    ! ?
LOGICAL,INTENT(IN),OPTIONAL             :: withObjects ! Delete mesh objects (nodes, sides, elements) too
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
! Delete objects
IF(PRESENT(withObjects)) THEN
  CALL flushSearchMesh(searchMesh,withObjects)
ELSE
  CALL flushSearchMesh(searchMesh,.FALSE.)
END IF
! Delete search mesh
DEALLOCATE(searchMesh%sortidx,searchMesh%mapping,searchMesh%sm)
DEALLOCATE(searchMesh)
END SUBROUTINE deleteSearchMesh


SUBROUTINE insertNode(searchMesh,Node,Elem)
!===================================================================================================================================
! Insert "Node" in search mesh
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSearchMesh),POINTER,INTENT(INOUT) :: searchMesh  ! ?
TYPE(tNode),POINTER,INTENT(IN)          :: Node   ! ?
TYPE(tElem),POINTER,INTENT(IN),OPTIONAL :: Elem   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tToObject),POINTER   :: ToObject   ! ?
INTEGER                   :: idx(3)   ! ?
!===================================================================================================================================
CALL getIdx(searchmesh,Node%x,idx)  ! Calculate search mesh indices of node
IF(.NOT.idxok(searchMesh,idx)) CALL abort(__STAMP__, &
                                          'search Mesh error')
CALL getNewToObject(ToObject)  ! Create new search mesh object
ToObject%Node=>Node
IF(PRESENT(Elem))ToObject%elem=>Elem
Node%refCount=Node%refCount+1
CALL insertToObject(searchMesh%sm(idx(1),idx(2),idx(3))%ToObject,ToObject)  ! Insert Toobject in search mesh
searchMesh%objectCount=searchMesh%objectCount+1
END SUBROUTINE insertNode


SUBROUTINE insertUniqueNode(node,searchMesh,redundantNodes,marker)
!===================================================================================================================================
! Insert unique node in uniqueSearchMesh, if already there in deleteSearchMesh,uses flag marker to check if node is already in mesh 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Tolerances,ONLY:SAMEPOINT
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNode),POINTER,INTENT(INOUT)       :: Node  ! ?
TYPE(tSearchMesh),POINTER,INTENT(INOUT) :: searchMesh,redundantNodes  ! ?
INTEGER,INTENT(IN)                      :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tToObject),POINTER   :: toObj   ! ?
TYPE(tNode),POINTER       :: uniqueNode   ! ?
!===================================================================================================================================
IF(.NOT.ASSOCIATED(node)) WRITE(*,*)'WARNING; NODE NOT ASSOCIATED'
IF(node%tmp.EQ.marker) RETURN ! node to be checked is already unique node, don't check twice

NULLIFY(uniqueNode)
toObj => getFirstToObject(searchMesh,.TRUE.,node%x)    ! Get first search mesh object
DO WHILE(ASSOCIATED(toObj))
  IF(SAMEPOINT(toObj%node%x,node%x,PP_MeshTolerance))THEN  ! Compare node coordinates
    uniqueNode=>toObj%node
    EXIT
  END IF
  toObj=>getNextToObject(searchMesh,.TRUE.)
END DO

IF(ASSOCIATED(uniqueNode))THEN
  CALL insertNode(redundantNodes,node)
  uniqueNode%ind=MIN(uniqueNode%ind,node%ind)
  node => uniqueNode
ELSE
  CALL insertNode(searchMesh,node)
  node%tmp=marker
END IF

END SUBROUTINE insertUniqueNode

FUNCTION getUniqueNode(searchMesh,coords,acceptOutSide,ind)
!===================================================================================================================================
! Check if node with coordinates "coords" already exists and create a new node if necessary
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:getNewNode
USE MOD_Mesh_Tolerances,ONLY:SAMEPOINT
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSearchMesh),POINTER,INTENT(INOUT) :: searchMesh  ! ?
REAL,INTENT(IN)                         :: coords(3)  ! Coordinates of node
LOGICAL,INTENT(IN)                      :: acceptOutSide ! Accept point outside search mesh
INTEGER,OPTIONAL,INTENT(IN)             :: ind           ! node ind
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tNode),POINTER                     :: getUniqueNode ! New / existing node
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tNode),POINTER       :: bNode      ! Local node pointer
TYPE(tToObject),POINTER   :: ToObject   ! Search mesh object
INTEGER                   :: idx_ur(3)  ! Indices of starting point
LOGICAL                   :: outside    ! Object is outside search mesh (possible for MPI case)
!===================================================================================================================================
CALL getIdx(searchmesh,coords,idx_ur)  ! Calculate search mesh indices from coordinates -> starting point
outside=.FALSE.
IF(.NOT.idxok(searchMesh,idx_ur)) THEN
  ! Coordinates are outside search mesh
  IF(acceptOutSide) THEN
    idx_ur  = MAX((/1,1,1/),MIN(idx_ur,searchMesh%nMesh))
    outside = .TRUE.
  ELSE
    ERRWRITE(*,*)'node coords',idx_ur,coords
    ERRWRITE(*,*)'searchmesh',searchMesh%xmin, searchMesh%xmax, searchMesh%nSearch
    CALL abort(__STAMP__, 'search Mesh error,node lies outside searchmesh')
  END IF
END IF

ToObject=>getFirstToObject(searchMesh,.TRUE.,idx=idx_ur)  ! Get pointer to first search mesh object at point idx_ur
DO WHILE(ASSOCIATED(ToObject) .AND. .NOT.outside)
  bNode=>ToObject%Node
  IF(ASSOCIATED(bNode))THEN
    IF (SAMEPOINT(bNode%x,coords,PP_MeshTolerance)) THEN  ! Compare node coordinates
      getUniqueNode=>bNode
      RETURN
    END IF
  END IF
  ToObject=>getNextToObject(searchMesh,.TRUE.)
END DO
! No existing node, create new one and insert in search mesh
CALL getNewNode(getUniqueNode,1)
CALL getNewToObject(ToObject)
ToObject%Node=>getUniqueNode
getUniqueNode%x=coords
IF(PRESENT(ind)) getUniqueNode%ind=ind

CALL insertToObject(searchMesh%sm(idx_ur(1),idx_ur(2),idx_ur(3))%ToObject,ToObject)

END FUNCTION getUniqueNode


SUBROUTINE insertSide(searchMesh,Side,mode,VVpolicy,success)
!===================================================================================================================================
! Insert "Side" into search mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:VV
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSearchMesh),POINTER,INTENT(INOUT) :: searchMesh  ! ?
TYPE(tSide),POINTER,INTENT(IN)          :: Side  ! ?
INTEGER,INTENT(IN)        :: mode       ! Reference point for insertion: 1=nodes of side 2,3=pseudo-barycenter 4=first node
INTEGER,INTENT(IN)        :: VVpolicy   ! Treatment of periodic boundary sides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL,INTENT(OUT)       :: success    ! Operation was successful
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tToObject),POINTER   :: newToObject,toObject   ! ?
REAL                      :: psudoBary(3),VV_loc(3)   ! ?
INTEGER                   :: idx(3),iNode !idx_ur(3),l,jNode,dm   ! ?
LOGICAL                   :: alreadyIn !Found,   ! ?
!===================================================================================================================================
success=.TRUE.
SELECT CASE(mode)
CASE(1) ! at each node ! fatal
  DO iNode=1,Side%nNodes
    CALL getIdx(searchmesh,Side%Node(iNode)%np%x,idx) ! Calculate search mesh indices of node iNode
    IF(idxok(searchMesh,idx)) THEN  ! Check if indices lie inside the limits of the search mesh
      alreadyIn=.FALSE.
      ToObject=>searchMesh%sm(idx(1),idx(2),idx(3))%ToObject
      DO WHILE(ASSOCIATED(ToObject))
        IF(ASSOCIATED(ToObject%side)) THEN
          ! Compare sides
          IF(ASSOCIATED(ToObject%Side,Side)) THEN
            alReadyIn=.TRUE.
          END IF
        END IF
        ToObject=>ToObject%nextToObject
      END DO
      IF(.NOT. alreadyIn) THEN
        ! Insert side into search mesh
        CALL getNewToObject(newToObject)
        newToObject%Side=>Side
        CALL insertToObject(searchMesh%sm(idx(1),idx(2),idx(3))%ToObject,newToObject) ! fill Local side search mesh
        searchMesh%objectCount=searchMesh%objectCount+1
      END IF
    ELSE
      CALL abort(__STAMP__,'search mesh error')
    END IF
  END DO
CASE(2,3) ! at psudoBary only one instance ! not fatal
  ! Calculate pseudo-barycenter
  psudobary=0.
  DO iNode=1,Side%nNodes
    psudobary=psudobary+Side%Node(iNode)%np%x
  END DO
  psudobary=psudobary/REAL(Side%nNodes)
  ! Treat periodic boundary sides
  VV_loc=0.
  IF(ASSOCIATED(Side%BC)) THEN
    IF(Side%BC%BCType .EQ. 1) THEN ! periodic
      VV_loc=VV(:,abs(Side%BC%BCAlphaInd))*sign(1,Side%BC%BCAlphaInd)  ! Displacement vector for periodic side
      SELECT CASE(VVpolicy)
      CASE(1) !ignore
        VV_loc=0.
      CASE(2)
        IF(Side%BC%BCalphaInd .GT. 0) VV_loc=0.
      CASE(3) !always
      END SELECT
    END IF
  END IF
  psudobary=psudobary+VV_loc  ! (New) pseudo-barycenter
  CALL getIdx(searchmesh,psudoBary,idx)
  IF(idxok(searchMesh,idx)) THEN
    ! Indices idx lie inside limits of search mesh: insert side
    CALL getNewToObject(newToObject)
    newToObject%Side=>Side
    CALL insertToObject(searchMesh%sm(idx(1),idx(2),idx(3))%ToObject,newToObject)
    searchMesh%objectCount=searchMesh%objectCount+1
  ELSEIF(mode .EQ. 3) THEN
    ! Indices idx do not lie inside limits: try nodes of side
    DO iNode=1,Side%nNodes
      CALL getIdx(searchMesh,Side%Node(iNode)%np%x+VV_loc,idx)  ! Calculate search mesh indices of node iNode
      IF(idxok(searchMesh,idx)) THEN
        ! If idx inside limits: insert side
        CALL getNewToObject(newToObject)
        newToObject%Side=>Side
        CALL insertToObject(searchMesh%sm(idx(1),idx(2),idx(3))%ToObject,newToObject)  
        searchMesh%objectCount=searchMesh%objectCount+1
        RETURN
      END IF
    END DO
    success=.FALSE.
    RETURN
  ELSE
    success=.FALSE.
  END IF
CASE(4) ! at first node ! fatal
    CALL getIdx(searchmesh,Side%Node(1)%np%x,idx)  ! Calculate search mesh indices of first node of side
    IF(idxok(searchMesh,idx)) THEN
      ! Insert side into search mesh
      CALL getNewToObject(newToObject)
      newToObject%Side=>Side
      CALL insertToObject(searchMesh%sm(idx(1),idx(2),idx(3))%ToObject,newToObject) ! fill Local side search mesh
      searchMesh%objectCount=searchMesh%objectCount+1
    ELSE
      CALL abort(__STAMP__,'search mesh error')
    END IF
CASE(5) ! if at least one node lies inside searchmesh, the side is inserted 
  DO iNode=1,Side%nNodes
    CALL getIdx(searchmesh,Side%Node(iNode)%np%x,idx) ! Calculate search mesh indices of node iNode
    IF(idxok(searchMesh,idx)) THEN  ! Check if indices lie inside the limits of the search mesh
      alreadyIn=.FALSE.
      ToObject=>searchMesh%sm(idx(1),idx(2),idx(3))%ToObject
      DO WHILE(ASSOCIATED(ToObject))
        IF(ASSOCIATED(ToObject%side)) THEN
          ! Compare sides
          IF(ASSOCIATED(ToObject%Side,Side)) THEN
            alReadyIn=.TRUE.
          END IF
        END IF
        ToObject=>ToObject%nextToObject
      END DO
      IF(.NOT. alreadyIn) THEN
        ! Insert side into search mesh
        CALL getNewToObject(newToObject)
        newToObject%Side=>Side
        CALL insertToObject(searchMesh%sm(idx(1),idx(2),idx(3))%ToObject,newToObject) ! fill Local side search mesh
        searchMesh%objectCount=searchMesh%objectCount+1
      END IF
    END IF
  END DO
END SELECT
END SUBROUTINE insertSide


SUBROUTINE insertAllMeshSides(searchMesh,mode)
!===================================================================================================================================
! Insert all sides in search mesh
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSearchMesh),POINTER,INTENT(INOUT) :: searchMesh   ! ?
INTEGER,INTENT(IN)        :: mode       ! See subroutine insert side
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER       :: Elem    ! Local element pointer
TYPE(tSide),POINTER       :: Side    ! Local side pointer
LOGICAL                   :: success   ! ?
!===================================================================================================================================
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    CALL insertSide(searchMesh,Side,mode,1,success)
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO
END SUBROUTINE insertAllMeshSides


FUNCTION idxok(searchMesh,idx)
!===================================================================================================================================
! Check if indices "idx" are locatet within the limits of "searchMesh": 1,searchMesh%nMesh(:)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSearchMesh),POINTER,INTENT(IN) :: searchMesh   ! ?
INTEGER,INTENT(IN)        :: idx(3)     ! Indices
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                   :: idxok      ! .TRUE. = indices within limits
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
IF ( (idx(1) .GE. 1) .AND. (idx(2) .GE. 1) .AND. (idx(1) .LE. searchMesh%nMesh(1)) .AND. (idx(2) .LE.  searchMesh%nMesh(2)) &
   .AND.(idx(3) .GE. 1) .AND. (idx(3) .LE. searchMesh%nMesh(3))) THEN
  idxok=.TRUE.
ELSE
  idxok=.FALSE.
END IF
END FUNCTION idxok


SUBROUTINE getIdx(searchmesh,coord,idx,testin)
!===================================================================================================================================
! Calculates Cartesian search mesh indices of a point with coordinates "coord"
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSearchMesh),POINTER,INTENT(IN) :: searchmesh  
REAL,INTENT(IN)              :: coord(3)    ! Point coordinates
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(OUT)          :: idx(3)      ! Indices i,j,k
LOGICAL,INTENT(OUT),OPTIONAL :: testin      ! .TRUE. = Indices inside search mesh limits
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
idx=1
idx(1:3)=INT((coord-searchMesh%xmin)/searchMesh%dx)+1
IF (PRESENT(testin)) testin=idxok(searchMesh,idx)
END SUBROUTINE getIdx


SUBROUTINE insertToObject(firstToObject,ToObject)
!===================================================================================================================================
! Insert object "toObject" in local (in search mesh) list that starts with "firstToObject". Object will be inserted before 
! "firstToObject".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tToObject),POINTER,INTENT(INOUT) :: firstToObject ! First object in local list of search objects
TYPE(tToObject),POINTER,INTENT(IN)    :: ToObject        ! ? 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
ToObject%nextToObject=>firstToObject
firstToObject=>ToObject
END SUBROUTINE insertToObject


SUBROUTINE insertNextToObject(prevToObject,ToObject)
!===================================================================================================================================
! Insert object "ToObject" in local (in search mesh) list. "ToObject" will be inserted behind "prevToObject".
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tToObject),POINTER,INTENT(IN) :: prevToObject ! First object in local list of search objects
TYPE(tToObject),POINTER,INTENT(IN) :: ToObject       ! ? 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tToObject),POINTER :: nextToObject   ! ?
!===================================================================================================================================
nextToObject=>prevToObject%nextToObject
prevToObject%nextToObject=>ToObject
ToObject%nextToObject=>nextToObject
END SUBROUTINE insertNextToObject


SUBROUTINE deleteToObjectList(firstToObject,withObjects)
!===================================================================================================================================
! Delete local (in search mesh) ToObject list.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tToObject),POINTER,INTENT(INOUT) :: firstToObject ! First object in local list of search objects
LOGICAL,INTENT(IN)                    :: withObjects   ! Delete mesh objects that are associated with "ToObject"
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tToObject),POINTER               :: sweepToObject,ToObject   ! ?
!===================================================================================================================================
sweepToObject=>firstToObject
DO WHILE(ASSOCIATED(sweepToObject))
  ToObject=>sweepToObject
  sweepToObject=>sweepToObject%nextToObject
  CALL deleteToObject(firstToObject,ToObject,withObjects)
END DO
NULLIFY(firstToObject)
END SUBROUTINE deleteToObjectList

END MODULE MOD_Search
