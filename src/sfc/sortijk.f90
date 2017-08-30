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

MODULE MOD_SortIJK
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
TYPE tBox
  REAL(KIND=8)     :: mini(3)
  INTEGER(KIND=8)  :: intfact
  REAL(KIND=8)     :: spacing(3)
END TYPE tBox
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE SortElemsByCoords
  MODULE PROCEDURE SortElemsByCoords
END INTERFACE


PUBLIC::SortElemsByCoords
!===================================================================================================================================

CONTAINS
SUBROUTINE SortElemsByCoords(nElems,ElemBary,nElems_IJK,Elem_IJK)
!===================================================================================================================================
! The Element barycenters are sorted in x,y,z directions, to find out about structured cartesian parts of the mesh 
! (for example extrusion in z direction...), Element indices /i,j,k) for each element are set, 
!  where unstructured 3D meshes have indices (1:nElems,1,1)
! and 2.5D unstructured mesh for example (1:nElems/nz,1,1:nz) and fully structured domains have (1:nElems/ny/nz,1:ny,1:nz) 
!===================================================================================================================================
! MODULES
USE MOD_globals,ONLY:UNIT_stdOut,Timer
USE MOD_Output_Vars,ONLY:DebugVisu
USE MOD_Mesh_Vars,ONLY:SpaceQuandt
USE MOD_SortingTools,ONLY:Qsort1DoubleInt1Pint
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                         :: nElems  ! ?
REAL,INTENT(IN)                           :: ElemBary(nElems,3)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(OUT)                       :: nElems_IJK(3)  ! ?
INTEGER,INTENT(OUT)                       :: Elem_IJK(nElems,3)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL  :: lower(3)  ! ?
REAL  :: upper(3)  ! ?
type(tBox) :: Box  ! ?
INTEGER :: iElem  ! ?
INTEGER(KIND=8) :: IntCoords(nElems,3)  ! ?
INTEGER(KIND=8) :: IntList(nElems)  ! ?
INTEGER         :: IDList(nElems)  ! ?
INTEGER         :: dir,nStructDirs,counter  ! ?
INTEGER         :: nElemsMax,nElemsMin  ! ?
INTEGER         :: ii,jj,kk   ! ?
LOGICAL         :: structDir(3)  ! ?
INTEGER         :: tol  ! ?
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A,A,A)')'SORT ELEMENTS BY COORDINATES ...'
IF(nElems.EQ.1)THEN
  nElems_IJK=1
  Elem_IJK=1
  WRITE(UNIT_stdOut,'(A)')'... DONE'
  RETURN
END IF
! Determine extreme verticies for bounding box
lower(1) = minval(ElemBary(:,1))
lower(2) = minval(ElemBary(:,2))
lower(3) = minval(ElemBary(:,3))
upper(1) = maxval(ElemBary(:,1))
upper(2) = maxval(ElemBary(:,2))
upper(3) = maxval(ElemBary(:,3))
lower=lower-0.1*MAXVAL(upper-lower)
upper=upper+0.1*MAXVAL(upper-lower)
CALL setBoundingBox2(Box,lower,upper)

DO iElem=1,nElems
  IntCoords(iElem,:) = NINT((ElemBary(iElem,:)-box%mini)*box%spacing)
END DO
WRITE(*,'(a,E11.3)')'   REAL tolerance    : ',SpaceQuandt*PP_RealTolerance
tol=NINT(SpaceQuandt*PP_RealTolerance*MAXVAL(box%spacing)) !in the integer space
structDir=.FALSE.
! now find out, which directions are structured (only possible if they are oriented to cartesian coordinate system!)
DO dir=1,3
  ! Now sort the elements by coordinate directions
DO iElem=1,nElems
  IDList(iElem)      = iElem
END DO
IntList=IntCoords(:,dir)
CALL QSort1DoubleInt1PInt(IntList, IDList)
  nElemsMin=nElems
  nElemsMax=0
  counter=1
  DO iElem=2,nElems
    IF(ABS(IntCoords(IDlist(iElem),dir)-IntCoords(IDlist(iElem-1),dir)).GT.tol)THEN
      nElemsMin=MIN(nElemsMin,counter)
      nElemsMax=MAX(nElemsMax,counter)
      counter=1
    ELSE
      counter=counter+1
    END IF
  END DO
  WRITE(*,*)'  dir,min,max:',dir,nElemsMin,nElemsMax
  IF(nElemsMax.NE.nElemsMin) THEN
    nElems_IJK(dir)=0 !not structured
    StructDir(dir)=.FALSE.
  ELSE
    nElems_IJK(dir)=nElemsMax
    StructDir(dir)=.TRUE.
  END IF
END DO !dir
nStructDirs=0
DO dir=1,3
  IF(StructDir(dir)) nStructdirs=nStructDirs+1
END DO

SELECT CASE(nStructDirs)
CASE(0)
  nElems_IJK=(/nElems,1,1/)
CASE(1)
  ! in this case, the structured direction has all elements of a 2d plane
  ! for example 3 elements in z direction and 30 elements overall: i=1,j=1,k=10
  ! this must be transformed to i=10,j=1,k=3
  DO dir=1,3
    IF(structDir(Dir)) EXIT
  END DO
  SELECT CASE(dir)
  CASE(1)
    nElems_IJK(1)=nElems/nElems_IJK(1)
    nElems_IJK(2)=nElems/nElems_IJK(1)
    nElems_IJK(3)=1
  CASE(2)
    nElems_IJK(2)=nElems/nElems_IJK(2)
    nElems_IJK(1)=nElems/nElems_IJK(2)
    nElems_IJK(3)=1
  CASE(3)
    nElems_IJK(3)=nElems/nElems_IJK(3)
    nElems_IJK(1)=nElems/nElems_IJK(3)
    nElems_IJK(2)=1
  END SELECT !dir 
CASE(2) ! checked
  ! in this case, one direction has only one element
  DO dir=1,3
    IF(.NOT.structDir(Dir)) EXIT
  END DO
  SELECT CASE(dir)
  CASE(1)
    nElems_IJK(1)=1
    nElems_IJK(2:3)=(/nElems_IJK(3),nElems_IJK(2)/)
  CASE(2)
    nElems_IJK(2)=1
    nElems_IJK(1:3:2)=(/nElems_IJK(3),nElems_IJK(1)/)
  CASE(3)
    nElems_IJK(3)=1
    nElems_IJK(1:2)=(/nElems_IJK(2),nElems_IJK(1)/)
  END SELECT !dir 
CASE(3)  
  !nx'=ny*nz
  !ny'=nx*nz
  !nz'=nx*ny
  !nx=sqrt(ny'/nx'*nz') =sqrt(nx*nz*nx*ny/(ny*nz))
  !ny=nz'/nx
  !nz=ny'/nx
  nElems_IJK(1)=NINT(SQRT(REAL(nElems_IJK(2)*nElems_IJK(3)/nElems_IJK(1))))
  nElems_IJK(2:3)=(/nElems_IJK(3),nElems_IJK(2)/)/nElems_IJK(1)
END SELECT !nstructdirs


!check if everythings okay
IF(PRODUCT(nElems_IJK).NE.nElems) THEN
  WRITE(UNIT_stdOut,*)'WARNING!!!!'
  WRITE(UNIT_stdOut,*)'Problem during sort elements by coordinate nElems /= nElems_I*Elems_J*nELems_K...',nElems,nElems_IJK
END IF

WRITE(*,'(A,I8)')  '   number of structured directions :',nStructDirs
WRITE(*,'(A,3I8)') '   number of Elements in i-j-k     :',nElems_IJK

! Now sort the elements, where z coordinate is first criteria, then y then x coordinate
DO iElem=1,nElems
  IDList(iElem)      = iElem
  IntList(iElem)     = (IntCoords(iElem,3)*box%intfact+IntCoords(iElem,2))*box%intfact + IntCoords(iElem,1)
END DO
CALL Qsort1DoubleInt1Pint(IntList, IDList)

! and finally set the ijk indices

iElem=0
DO kk=1,nElems_IJK(3)
  DO jj=1,nElems_IJK(2)
    DO ii=1,nElems_IJK(1)
      iElem=iElem+1
      Elem_IJK(IDlist(iElem),:)=(/ii,jj,kk/) 
    END DO
  END DO
END DO
IF(DebugVisu)THEN
OPEN(UNIT=101,FILE='sorting.dat',STATUS='REPLACE')
WRITE(101,'(a)')'VARIABLES= "X" "Y" "Z"'
WRITE(101,'(a,I4,a,I4,a,I4,a)')'ZONE T="sortIJK", I=',nElems_IJK(1),', J=',nElems_IJK(2),', K=',nElems_IJK(3),', DATAPACKING=POINT'
iElem=0
DO kk=1,nElems_IJK(3)
  DO jj=1,nElems_IJK(2)
    DO ii=1,nElems_IJK(1)
      iElem=iElem+1
      WRITE(101,*) ElemBary(IDlist(iElem),:)
    END DO
  END DO
END DO
CLOSE(101)
END IF
CALL Timer(.FALSE.)
END SUBROUTINE SortElemsByCoords

SUBROUTINE setBoundingBox2(box, mini, maxi)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(KIND=8),INTENT(IN)      :: mini(3)  ! ?
REAL(KIND=8),INTENT(IN)      :: maxi(3)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tBox),INTENT(OUT)  :: box  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(KIND=8)      :: blen(3)  ! ?
INTEGER           :: nBits  ! ?
!===================================================================================================================================
box%mini = mini
blen = maxi - mini
nbits = (bit_size(box%intfact)-1)/3 
box%intfact = 2**nbits-1
box%spacing = REAL(box%intfact)/blen
END SUBROUTINE setBoundingBox2

END MODULE MOD_sortIJK
