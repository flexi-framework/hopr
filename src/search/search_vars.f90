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
MODULE MOD_Search_Vars
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tNode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
TYPE tSearchMesh ! search mesh handle
  TYPE(tToObjectPtr),POINTER          ::       sm(:,:,:)              ! sm(nMesh(1),nMesh(2),nMesh(3)) Objects in the mesh
  TYPE(tToObject), POINTER            ::       actualToObject         ! pointer to current object
  REAL                                ::       Mid(3)
  REAL                                ::       Radius
  REAL                                ::       xmin(3)                ! Search mesh extend min
  REAL                                ::       xmax(3)                ! Search mesh extend max
  REAL                                ::       dx(3)                  ! step size for search
  INTEGER                             ::       nMesh(3)               ! number of point in each direction
  INTEGER                             ::       nSearch(3)             ! (-,+) Search area in each direction
  INTEGER                             ::       n                      ! number of points in search area
  INTEGER                             ::       nred                   ! reduced search area
  INTEGER                             ::       actualInd(3)           ! current search index
  INTEGER                             ::       idx_ur(3)              ! starting point index
  INTEGER                             ::       idxcounter             ! position in search area
  INTEGER,POINTER                     ::       mapping(:,:)           ! mapping(n,nDim) mMapping from search area point index 
  !                                                                   ! to search mesh point indices
  INTEGER,ALLOCATABLE                 ::       sortidx(:)             ! sortidx(n) search area index of current point 
  !                                                                   ! in search area
  INTEGER                             ::       objectCount            ! Counter for objects in search mesh
END TYPE tSearchMesh

TYPE tToObject ! list of objects for search mesh
  TYPE(tElem),POINTER                 ::       elem                   ! pointer to element
  TYPE(tSide),POINTER                 ::       side                   ! pointer to side
  TYPE(tNode),POINTER                 ::       node                   ! pointer to node
  TYPE(tToObject),POINTER             ::       nextToObject           ! pointer to next obj
END TYPE tToObject
      
TYPE tToObjectPtr
  TYPE(tToObject),POINTER             ::       ToObject               ! pointer to search object list
END TYPE tToObjectPtr

REAL                 ::       RefineSideSearch       ! factor(>=1) to increase the number of search elemens in 
                                                     ! side search mesh, used for rasterfahndung. increase
                                                     ! yields faster connect, but higher memory requirements
INTEGER              ::       nElemsNodeSearch       ! number of elems in node search mesh. is used for routine
                                                     ! get unique node

INTERFACE getNewToObject
   MODULE PROCEDURE getNewToObject
END INTERFACE

INTERFACE deleteToObject
   MODULE PROCEDURE deleteToObject
END INTERFACE
!===================================================================================================================================

CONTAINS
SUBROUTINE getNewToObject(ToObject)
!===================================================================================================================================
! Allocate and initialize new "ToObject"
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tToObject),POINTER,INTENT(OUT) :: ToObject ! Object in search mesh
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
ALLOCATE(ToObject)
NULLIFY(ToObject%nextToObject,ToObject%Elem,ToObject%Side,ToObject%Node)
END SUBROUTINE getNewToObject


SUBROUTINE deleteToObject(firstToObject,ToObject,withObjects)
!===================================================================================================================================
! Delete object "ToObject" from search mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:deleteElem,deleteSide,deleteNode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tToObject),POINTER,INTENT(OUT) :: firstToObject ! First object in local list of search objects
TYPE(tToObject),POINTER,INTENT(OUT) :: ToObject      ! object to be deleted
LOGICAL,INTENT(IN),OPTIONAL         :: withObjects   ! Delete mesh objects that are associated with "ToObject"
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tToObject),POINTER             :: prevToObject    ! ? 
!===================================================================================================================================
IF(.NOT.ASSOCIATED(ToObject))RETURN
! Disconnect ToObject from list
IF(ASSOCIATED(ToObject,firstToObject)) THEN
   firstToObject=>ToObject%nextToObject
ELSE
   prevToObject=>firstToObject
   DO WHILE(ASSOCIATED(prevToObject))
     IF(ASSOCIATED(prevToObject%nextToObject,ToObject)) EXIT
     prevToObject=>prevToObject%nextToObject
   END DO
   prevToObject%nextToObject=>ToObject%nextToObject
END IF
NULLIFY(ToObject%nextToObject)
IF(PRESENT(withObjects)) THEN
  IF(WithObjects) THEN
    IF(ASSOCIATED(ToObject%Elem)) CALL deleteElem(ToObject%Elem,ToObject%Elem)
    IF(ASSOCIATED(ToObject%Side)) CALL deleteSide(ToObject%Side,ToObject%Side)
    IF(ASSOCIATED(ToObject%Node)) CALL deleteNode(ToObject%Node)
  END IF
END IF
DEALLOCATE(ToObject)
NULLIFY(ToObject)
END SUBROUTINE deleteToObject

END MODULE MOD_Search_Vars
