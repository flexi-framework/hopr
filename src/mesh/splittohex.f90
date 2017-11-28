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
MODULE MOD_SplitToHex
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE SplitElementsToHex
  MODULE PROCEDURE SplitElementsToHex
END INTERFACE

INTERFACE SplitAllHexa
  MODULE PROCEDURE SplitAllHexa
END INTERFACE

INTERFACE SplitHexaByBoxes
  MODULE PROCEDURE SplitHexaByBoxes
END INTERFACE

PUBLIC:: SplitElementsToHex,SplitAllHexa,SplitHexaByBoxes
!===================================================================================================================================

CONTAINS
SUBROUTINE SplitAllHexa(nFine)
!===================================================================================================================================
! call routines to split all hexa mesh by a factor
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY:tElem,FirstElem
USE MOD_Mesh_Vars, ONLY:useCurveds
USE MOD_Mesh_Basis,ONLY:ElemGeometry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)            :: nFine  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER           :: Elem  ! ?
INTEGER                       :: maxInd,ii  ! ?
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A,I8,A)')' SPLIT EACH HEXA INTO nFineHexa^3=' ,nFine**3,' HEXA ...'
IF(useCurveds) CALL abort(__STAMP__, &
          'SplitAllHex cannot be used with curved elements up to now!')

!check first if all elements are hexa
Elem=>firstElem      !get maxInd of nodes
DO WHILE(ASSOCIATED(Elem))
  IF(Elem%nNodes .NE. 8) THEN
     WRITE(*,*)'non-hexaheral element found, nNodes= ',Elem%nNodes
     WRITE(*,*)' ===>> NO nFine APPLIED!!!'
     WRITE(UNIT_stdOut,*)'...DONE!'
     RETURN
  END IF
  Elem=>Elem%nextElem
ENDDO
maxInd=0
Elem=>firstElem      !get maxInd of nodes
DO WHILE(ASSOCIATED(Elem))
  DO ii=1,Elem%nNodes
    IF(Elem%node(ii)%np%ind.GT.maxInd) maxInd=Elem%node(ii)%np%ind 
  END DO
  Elem=>Elem%nextElem
ENDDO
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
    CALL SplitHexa8(Elem,nFine,maxInd)
  Elem=>Elem%nextElem
END DO
CALL Timer(.FALSE.)
END SUBROUTINE SplitAllHexa

SUBROUTINE SplitHexaByBoxes()
!===================================================================================================================================
! call routines to split only hexas with its barycenter within a splitbox. next splitbox is applied on top
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY:tElem,FirstElem
USE MOD_Mesh_Vars, ONLY:useCurveds
USE MOD_Mesh_Vars, ONLY:nSplitBoxes,SplitBoxes
USE MOD_Mesh_Basis,ONLY:ElemGeometry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER           :: Elem  
INTEGER                       :: maxInd,ii,iBox 
REAL                          :: xBary(3)
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A,I8,A)')' SPLIT HEXA BY SPLITBOXES ...'
IF(useCurveds) CALL abort(__STAMP__, &
          'SplitAllHex cannot be used with curved elements up to now!')

!check first if all elements are hexa
Elem=>firstElem      !get maxInd of nodes
DO WHILE(ASSOCIATED(Elem))
  IF(Elem%nNodes .NE. 8) THEN
     WRITE(*,*)'non-hexaheral element found, nNodes= ',Elem%nNodes
     WRITE(*,*)' ===>> NO SPLITBOX APPLIED!!!'
     WRITE(UNIT_stdOut,*)'...DONE!'
     RETURN
  END IF
  Elem=>Elem%nextElem
ENDDO
DO iBox=1,nSplitBoxes
  ASSOCIATE(xmin=>SplitBoxes(:,1,iBox),xmax=>SplitBoxes(:,2,iBox))
  maxInd=0
  Elem=>firstElem      !get maxInd of nodes
  DO WHILE(ASSOCIATED(Elem))
    DO ii=1,Elem%nNodes
      IF(Elem%node(ii)%np%ind.GT.maxInd) maxInd=Elem%node(ii)%np%ind 
    END DO
    Elem=>Elem%nextElem
  ENDDO
  Elem=>firstElem
  DO WHILE(ASSOCIATED(Elem))
    xBary=0.
    DO ii=1,Elem%nNodes
      xBary=xBary+Elem%Node(ii)%np%x
    END DO
    xBary=xBary/REAL(Elem%nNodes)
    IF(     ((xBary(1).GT.xmin(1)).AND.(xBary(1).LT.xmax(1))) &
       .AND.((xBary(2).GT.xmin(2)).AND.(xBary(2).LT.xmax(2))) &
       .AND.((xBary(3).GT.xmin(3)).AND.(xBary(3).LT.xmax(3))) )THEN
      CALL SplitHexa8(Elem,2,maxInd)
    END IF
    Elem=>Elem%nextElem
  END DO
  END ASSOCIATE !xmin,xmax
END DO !iBox
CALL Timer(.FALSE.)
END SUBROUTINE SplitHexaByBoxes

SUBROUTINE SplitElementsToHex()
!===================================================================================================================================
! call routines to split tetra, prism and hexa to hexas. no pyramids allowed!! 
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars, ONLY:tElem,FirstElem
USE MOD_Mesh_Vars, ONLY:useCurveds
USE MOD_Mesh_Basis,ONLY:ElemGeometry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER           :: Elem  ! ?
INTEGER                       :: maxInd,ii,counter(3)  ! ?
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')' SPLIT ELEMENTS...'
maxInd=0
IF(useCurveds) CALL abort(__STAMP__, &
          'SplitToHex cannot be used with curved elements up to now!')

Elem=>firstElem      !get maxInd of nodes
DO WHILE(ASSOCIATED(Elem))
  DO ii=1,Elem%nNodes
    IF(Elem%node(ii)%np%ind.GT.maxInd) maxInd=Elem%node(ii)%np%ind 
  END DO
  Elem=>Elem%nextElem
ENDDO
counter=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  SELECT CASE(Elem%nNodes)
  CASE(4) !Tetra
    counter(1)=counter(1)+1
    CALL SplitTetraToHexa(Elem,maxInd)
  CASE(5) !Pyra
    CALL abort(__STAMP__, &
          'pyramidal element found, cannot be splitted to hexa!')
  CASE(6) !linear bilinear prism
    counter(2)=counter(2)+1
    CALL SplitPentaToHexa(Elem,2,maxInd)
  CASE(8) !linear trilinear hexa
    counter(3)=counter(3)+1
    CALL SplitHexa8(Elem,2,maxInd)
  END SELECT
  Elem=>Elem%nextElem
END DO
WRITE(UNIT_stdOut,'(I8,A,I8,A)') counter(1), ' Tetrahedra splitted to ', 4* counter(1),' hexas'
WRITE(UNIT_stdOut,'(I8,A,I8,A)') counter(2), ' Pentas splitted to     ', 6* counter(2),' hexas'
WRITE(UNIT_stdOut,'(I8,A,I8,A)') counter(3), ' Hexas splitted to      ', 8* counter(3),' hexas'
CALL Timer(.FALSE.)
END SUBROUTINE SplitElementsToHex

SUBROUTINE SplitHexa8(Elem,M,maxInd)
!===================================================================================================================================
! split a trilinear hexa in 8
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,  ONLY:tElem,tElemPtr,tSide,tSidePtr,tNodePtr
USE MOD_Mesh_Vars,  ONLY:FirstElem
USE MOD_Mesh_Vars,  ONLY:deleteElem
USE MOD_Mesh_Vars,  ONLY:getNewNode,copyBC
USE MOD_Mesh_Basis, ONLY:getNewHexa
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(INOUT)            :: Elem  ! ?
INTEGER, INTENT(IN)            :: M    ! number of elements to refine
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)          :: maxInd  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER       :: aElem,nextElem,prevElem  ! ?
TYPE(tSide),POINTER       :: Side
TYPE(tElemPtr)            :: NewElems(0:M-1,0:M-1,0:M-1)  ! ?
TYPE(tNodePtr)            :: NewNodes(0:M,0:M,0:M)  ! ?
TYPE(tSidePtr)            :: NewSides(0:M-1,0:M-1,0:M-1,6)  ! ?
INTEGER                   :: i,j,k,l,locSideID  ! ?
REAL                      :: sM  ! ?
!===================================================================================================================================
sM=1./REAL(M)
!corner nodes by cgns standard
NewNodes(0,0,0)%np=>Elem%Node(1)%np
NewNodes(M,0,0)%np=>Elem%Node(2)%np
NewNodes(M,M,0)%np=>Elem%Node(3)%np
NewNodes(0,M,0)%np=>Elem%Node(4)%np
NewNodes(0,0,M)%np=>Elem%Node(5)%np
NewNodes(M,0,M)%np=>Elem%Node(6)%np
NewNodes(M,M,M)%np=>Elem%Node(7)%np
NewNodes(0,M,M)%np=>Elem%Node(8)%np

!midnodes at edges
DO i=0,M,M ; DO j=0,M,M; DO l=1,M-1
  CALL GetNewNode(NewNodes(l,i,j)%np,ind=maxInd)
  CALL GetNewNode(NewNodes(i,l,j)%np,ind=maxInd)
  CALL GetNewNode(NewNodes(i,j,l)%np,ind=maxInd)
  NewNodes(l,i,j)%np%x=sM*(REAL(M-l)*NewNodes(0,i,j)%np%x + REAL(l)*NewNodes(M,i,j)%np%x) 
  NewNodes(i,l,j)%np%x=sM*(REAL(M-l)*NewNodes(i,0,j)%np%x + REAL(l)*NewNodes(i,M,j)%np%x) 
  NewNodes(i,j,l)%np%x=sM*(REAL(M-l)*NewNodes(i,j,0)%np%x + REAL(l)*NewNodes(i,j,M)%np%x) 
END DO; END DO; END DO
! midnodes at faces
DO i=0,M,M; DO k=1,M-1; DO l=1,M-1 
  CALL GetNewNode(NewNodes(i,k,l)%np,ind=maxInd)
  CALL GetNewNode(NewNodes(k,i,l)%np,ind=maxInd)
  CALL GetNewNode(NewNodes(k,l,i)%np,ind=maxInd)
  NewNodes(i,k,l)%np%x=sM*(REAL(M-k)*NewNodes(i,0,l)%np%x + REAL(k)*NewNodes(i,M,l)%np%x) 
  NewNodes(k,i,l)%np%x=sM*(REAL(M-k)*NewNodes(0,i,l)%np%x + REAL(k)*NewNodes(M,i,l)%np%x) 
  NewNodes(k,l,i)%np%x=sM*(REAL(M-k)*NewNodes(0,l,i)%np%x + REAL(k)*NewNodes(M,l,i)%np%x) 
END DO; END DO; END DO
! inner element nodes
DO i=1,M-1; DO j=1,M-1; DO k=1,M-1 
  CALL GetNewNode(NewNodes(i,j,k)%np,ind=maxInd)
  NewNodes(i,j,k)%np%x=sM*(REAL(M-i)*NewNodes(0,j,k)%np%x + REAL(i)*NewNodes(M,j,k)%np%x) 
END DO; END DO; END DO;

!build new elements
DO i=0,M-1; DO j=0,M-1; DO k=0,M-1
  CALL GetNewHexa(NewElems(i,j,k)%EP,Elem%zone, &
                  NewNodes(i  ,j  ,k  )%np,NewNodes(i+1,j  ,k  )%np,NewNodes(i+1,j+1,k  )%np,NewNodes(i  ,j+1,k  )%np, &
                  NewNodes(i  ,j  ,k+1)%np,NewNodes(i+1,j  ,k+1)%np,NewNodes(i+1,j+1,k+1)%np,NewNodes(i  ,j+1,k+1)%np  )
  LocSideID=0
  Side=>NewElems(i,j,k)%EP%firstSide
  DO WHILE(ASSOCIATED(Side))
    LocSideID=LocSideID+1
    NewSides(i,j,k,LocSideID)%sp=>Side
    Side=>Side%nextElemSide
  END DO
END DO; END DO; END DO

!boundary conditions
locSideID=0
Side=>Elem%firstSide
DO WHILE(ASSOCIATED(Side))
  locSideID=locSideID+1
  IF(ASSOCIATED(Side%BC))THEN
    DO i=0,M-1; DO j=0,M-1
      SELECT CASE(locSideID)
      CASE(1) ! zeta_minus
        CALL copyBC(Side,NewSides(i,j,0,LocSideID)%sp)
      CASE(6) ! zeta plus
        CALL copyBC(Side,NewSides(i,j,M-1,LocSideID)%sp)
      CASE(2) ! eta_minus
        CALL copyBC(Side,NewSides(i,0,j,LocSideID)%sp)
      CASE(4) ! eta plus
        CALL copyBC(Side,NewSides(i,M-1,j,LocSideID)%sp)
      CASE(5) ! xi minus
        CALL copyBC(Side,NewSides(0,i,j,LocSideID)%sp)
      CASE(3) ! xi plus
        CALL copyBC(Side,NewSides(M-1,i,j,LocSideID)%sp)
      END SELECT
    END DO; END DO !i,j=0,1
  END IF
  Side=>Side%nextElemSide
END DO  ! WHILE(ASSOCIATED(Side))

!inline new elements in element list
aElem=>Elem
nextElem=>Elem%nextElem
prevElem=>Elem%prevElem
DO i=0,M-1; DO j=0,M-1; DO k=0,M-1
  aElem%nextElem          => NewElems(i,j,k)%EP
  aElem%nextElem%prevElem => aElem
  aElem                   => aElem%nextElem
END DO; END DO; END DO
!take out elem
IF(ASSOCIATED(nextElem)) THEN
  aElem%nextElem=>nextElem
  nextElem%prevElem=>aElem
END IF
CALL DeleteElem(firstElem,Elem)
Elem=>aElem

IF(ASSOCIATED(prevElem))THEN
  NewElems(0,0,0)%ep%prevElem=>prevElem
  prevElem%nextElem=>NewElems(0,0,0)%ep
ELSE
  firstElem=>NewElems(0,0,0)%ep
END IF
END SUBROUTINE SplitHexa8

SUBROUTINE SplitTetraToHexa(Elem,maxInd)
!===================================================================================================================================
! split a trilinear hexa in 8
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,  ONLY:tElem,tElemPtr,tSide,tNodePtr,tSidePtr
USE MOD_Mesh_Vars,  ONLY:FirstElem
USE MOD_Mesh_Vars,  ONLY:getNewNode,copyBC
USE MOD_Mesh_Vars,  ONLY:deleteElem
USE MOD_Mesh_Basis, ONLY:getNewHexa
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(INOUT)            :: Elem  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)          :: maxInd  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER       :: aElem,nextElem,prevElem  ! ?
TYPE(tSide),POINTER       :: Side  ! ?
TYPE(tSidePtr)            :: NewSides(4,6)  ! ?
TYPE(tElemPtr)            :: NewElems(4)  ! ?
TYPE(tNodePtr)            :: NewNodes(15)  ! ?
INTEGER                   :: i,locSideID  ! ?
!===================================================================================================================================
!corner nodes by cgns standard
NewNodes(1)%np=>Elem%Node(1)%np
NewNodes(2)%np=>Elem%Node(2)%np
NewNodes(3)%np=>Elem%Node(3)%np
NewNodes(4)%np=>Elem%Node(4)%np

!nodes on edges
CALL GetNewNode(NewNodes( 5)%np,ind=maxInd)
CALL GetNewNode(NewNodes( 6)%np,ind=maxInd)
CALL GetNewNode(NewNodes( 7)%np,ind=maxInd)
CALL GetNewNode(NewNodes( 8)%np,ind=maxInd)
CALL GetNewNode(NewNodes( 9)%np,ind=maxInd)
CALL GetNewNode(NewNodes(10)%np,ind=maxInd)
NewNodes( 5)%np%x=0.5*(NewNodes(1)%np%x+NewNodes(2)%np%x)
NewNodes( 6)%np%x=0.5*(NewNodes(2)%np%x+NewNodes(3)%np%x)
NewNodes( 7)%np%x=0.5*(NewNodes(1)%np%x+NewNodes(3)%np%x)
NewNodes( 8)%np%x=0.5*(NewNodes(1)%np%x+NewNodes(4)%np%x)
NewNodes( 9)%np%x=0.5*(NewNodes(2)%np%x+NewNodes(4)%np%x)
NewNodes(10)%np%x=0.5*(NewNodes(3)%np%x+NewNodes(4)%np%x)

!nodes on faces
CALL GetNewNode(NewNodes(11)%np,ind=maxInd)
CALL GetNewNode(NewNodes(12)%np,ind=maxInd)
CALL GetNewNode(NewNodes(13)%np,ind=maxInd)
CALL GetNewNode(NewNodes(14)%np,ind=maxInd)
NewNodes(11)%np%x=(NewNodes(1)%np%x+NewNodes(2)%np%x+NewNodes(3)%np%x)/3.
NewNodes(12)%np%x=(NewNodes(1)%np%x+NewNodes(2)%np%x+NewNodes(4)%np%x)/3.
NewNodes(13)%np%x=(NewNodes(2)%np%x+NewNodes(3)%np%x+NewNodes(4)%np%x)/3.
NewNodes(14)%np%x=(NewNodes(1)%np%x+NewNodes(3)%np%x+NewNodes(4)%np%x)/3.

!node inside
CALL GetNewNode(NewNodes(15)%np,ind=maxInd)
NewNodes(15)%np%x=(NewNodes(1)%np%x+NewNodes(2)%np%x+NewNodes(3)%np%x+NewNodes(4)%np%x)*0.25

!build new elements, corner node is always origin
  CALL GetNewHexa(NewElems(1)%EP,Elem%zone, &
                        NewNodes( 1)%np,NewNodes( 5)%np,NewNodes(11)%np,NewNodes( 7)%np, &
                        NewNodes( 8)%np,NewNodes(12)%np,NewNodes(15)%np,NewNodes(14)%np  )
  CALL GetNewHexa(NewElems(2)%EP,Elem%zone, &
                        NewNodes( 2)%np,NewNodes( 6)%np,NewNodes(11)%np,NewNodes( 5)%np,&
                        NewNodes( 9)%np,NewNodes(13)%np,NewNodes(15)%np,NewNodes(12)%np)
  CALL GetNewHexa(NewElems(3)%EP,Elem%zone, &
                        NewNodes( 3)%np,NewNodes( 7)%np,NewNodes(11)%np,NewNodes( 6)%np, &
                        NewNodes(10)%np,NewNodes(14)%np,NewNodes(15)%np,NewNodes(13)%np  )
  CALL GetNewHexa(NewElems(4)%EP,Elem%zone, &
                        NewNodes( 4)%np,NewNodes( 8)%np,NewNodes(14)%np,NewNodes(10)%np, &
                        NewNodes( 9)%np,NewNodes(12)%np,NewNodes(15)%np,NewNodes(13)%np  )

DO i=1,4
  LocSideID=0
  Side=>NewElems(i)%EP%firstSide
  DO WHILE(ASSOCIATED(Side))
    LocSideID=LocSideID+1
    NewSides(i,LocSideID)%sp=>Side
    Side=>Side%nextElemSide
  END DO
END DO
!boundary conditions
locSideID=0
Side=>Elem%firstSide
DO WHILE(ASSOCIATED(Side))
  locSideID=locSideID+1
  IF(ASSOCIATED(Side%BC))THEN
    SELECT CASE(locSideID)
    CASE(1) ! 1-2-3 side
      CALL copyBC(Side,NewSides(1,1)%sp) !Elem 1, locSideID 1
      CALL copyBC(Side,NewSides(2,1)%sp) !Elem 2, locSideID 1
      CALL copyBC(Side,NewSides(3,1)%sp) !Elem 3, locSideID 1
    CASE(2) ! 1-2-4 side
      CALL copyBC(Side,NewSides(1,2)%sp) !Elem 1, locSideID 2
      CALL copyBC(Side,NewSides(2,5)%sp) !Elem 2, locSideID 5
      CALL copyBC(Side,NewSides(4,2)%sp) !Elem 4, locSideID 2
    CASE(3) ! 2-3-4 side
      CALL copyBC(Side,NewSides(2,2)%sp) !Elem 2, locSideID 2
      CALL copyBC(Side,NewSides(3,5)%sp) !Elem 3, locSideID 5
      CALL copyBC(Side,NewSides(4,5)%sp) !Elem 4, locSideID 5
    CASE(4) ! 3-1-4 side
      CALL copyBC(Side,NewSides(3,2)%sp) !Elem 3, locSideID 2
      CALL copyBC(Side,NewSides(1,5)%sp) !Elem 1, locSideID 5
      CALL copyBC(Side,NewSides(4,1)%sp) !Elem 4, locSideID 1
    END SELECT
  END IF
  Side=>Side%nextElemSide
END DO  ! WHILE(ASSOCIATED(Side))

!inline new elements in element list
aElem=>Elem
nextElem=>Elem%nextElem
prevElem=>Elem%prevElem
DO i=1,4
  aElem%nextElem          => NewElems(i)%EP
  aElem%nextElem%prevElem => aElem
  aElem                   => aElem%nextElem
END DO
!take out elem
IF(ASSOCIATED(nextElem)) THEN
  aElem%nextElem=>nextElem
  nextElem%prevElem=>aElem
END IF
CALL DeleteElem(firstElem,Elem)
Elem=>aElem

IF(ASSOCIATED(prevElem))THEN
  NewElems(1)%ep%prevElem=>prevElem
  prevElem%nextElem=>NewElems(1)%ep
ELSE
  firstElem=>NewElems(1)%ep
END IF

END SUBROUTINE SplitTetraToHexa


SUBROUTINE SplitPentaToHexa(Elem,M,maxInd)
!===================================================================================================================================
! split a trilinear hexa in 8
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,  ONLY:tElem,tElemPtr,tSide,tSidePtr,tNodePtr
USE MOD_Mesh_Vars,  ONLY:FirstElem
USE MOD_Mesh_Vars,  ONLY:getNewNode,copyBC
USE MOD_Mesh_Vars,  ONLY:deleteElem
USE MOD_Mesh_Basis, ONLY:getNewHexa
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(INOUT)            :: Elem  ! ?
INTEGER, INTENT(IN)            :: M    ! number of elements to refine zeta direction
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER,INTENT(INOUT)          :: maxInd  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER       :: aElem,nextElem,prevElem  ! ?
TYPE(tSide),POINTER       :: Side  ! ?
TYPE(tElemPtr)            :: NewElems(3,0:M-1)  ! ?
TYPE(tSidePtr)            :: NewSides(3,0:M-1,6)  ! ?
TYPE(tNodePtr)            :: NewNodes(7,0:M)  ! ?
INTEGER                   :: i,l,locSideID  ! ?
REAL                      :: sM  ! ?
!===================================================================================================================================
sM=1./REAL(M)
!corner nodes by cgns standard
NewNodes(1,0)%np=>Elem%Node(1)%np
NewNodes(2,0)%np=>Elem%Node(2)%np
NewNodes(3,0)%np=>Elem%Node(3)%np
NewNodes(1,M)%np=>Elem%Node(4)%np
NewNodes(2,M)%np=>Elem%Node(5)%np
NewNodes(3,M)%np=>Elem%Node(6)%np

!midnodes at edges and triangle face
DO l=0,M,M
  CALL GetNewNode(NewNodes(4,l)%np,ind=maxInd)
  CALL GetNewNode(NewNodes(5,l)%np,ind=maxInd)
  CALL GetNewNode(NewNodes(6,l)%np,ind=maxInd)
  CALL GetNewNode(NewNodes(7,l)%np,ind=maxInd)
  NewNodes(4,l)%np%x=0.5*(NewNodes(1,l)%np%x+NewNodes(2,l)%np%x)
  NewNodes(5,l)%np%x=0.5*(NewNodes(2,l)%np%x+NewNodes(3,l)%np%x)
  NewNodes(6,l)%np%x=0.5*(NewNodes(1,l)%np%x+NewNodes(3,l)%np%x)
  NewNodes(7,l)%np%x=(NewNodes(1,l)%np%x+NewNodes(2,l)%np%x+NewNodes(3,l)%np%x)/3.
END DO
DO l=1,M-1
  DO i=1,7
    CALL GetNewNode(NewNodes(i,l)%np,ind=maxInd)
    NewNodes(i,l)%np%x=sM*(REAL(M-l)*NewNodes(i,0)%np%x + REAL(l)*NewNodes(i,M)%np%x) 
  END DO
END DO

!build new elements
DO l=0,M-1
  CALL GetNewHexa(NewElems(1,l)%EP,Elem%zone, &
                        NewNodes(1,l  )%np,NewNodes(4,l  )%np,NewNodes(7,l  )%np,NewNodes(6,l  )%np, &
                        NewNodes(1,l+1)%np,NewNodes(4,l+1)%np,NewNodes(7,l+1)%np,NewNodes(6,l+1)%np  )
  CALL GetNewHexa(NewElems(2,l)%EP,Elem%zone, &
                        NewNodes(2,l  )%np,NewNodes(5,l  )%np,NewNodes(7,l  )%np,NewNodes(4,l  )%np, &
                        NewNodes(2,l+1)%np,NewNodes(5,l+1)%np,NewNodes(7,l+1)%np,NewNodes(4,l+1)%np  )
  CALL GetNewHexa(NewElems(3,l)%EP,Elem%zone, &
                        NewNodes(3,l  )%np,NewNodes(6,l  )%np,NewNodes(7,l  )%np,NewNodes(5,l  )%np, &
                        NewNodes(3,l+1)%np,NewNodes(6,l+1)%np,NewNodes(7,l+1)%np,NewNodes(5,l+1)%np  )

  DO i=1,3
    LocSideID=0
    Side=>NewElems(i,l)%EP%firstSide
    DO WHILE(ASSOCIATED(Side))
      LocSideID=LocSideID+1
      NewSides(i,l,LocSideID)%sp=>Side
      Side=>Side%nextElemSide
    END DO
  END DO
END DO
!boundary conditions
locSideID=0
Side=>Elem%firstSide
DO WHILE(ASSOCIATED(Side))
  locSideID=locSideID+1
  IF(ASSOCIATED(Side%BC))THEN
    SELECT CASE(locSideID)
    CASE(1) ! 1-2-4 side
      DO l=0,M-1
        CALL copyBC(Side,NewSides(1,l,2)%sp) !Elem 1, locSideID 2
        CALL copyBC(Side,NewSides(2,l,5)%sp) !Elem 2, locSideID 5
      END DO
    CASE(2) ! 2-3-5 side
      DO l=0,M-1
        CALL copyBC(Side,NewSides(2,l,2)%sp) !Elem 2, locSideID 2
        CALL copyBC(Side,NewSides(3,l,5)%sp) !Elem 3, locSideID 5
      END DO
    CASE(3) ! 3-1-6 side
      DO l=0,M-1
        CALL copyBC(Side,NewSides(3,l,2)%sp) !Elem 3, locSideID 2
        CALL copyBC(Side,NewSides(1,l,5)%sp) !Elem 1, locSideID 5
      END DO
    CASE(4) ! 1-2-3 side
      CALL copyBC(Side,NewSides(1,0,1)%sp) !Elem 1, locSideID 1
      CALL copyBC(Side,NewSides(2,0,1)%sp) !Elem 2, locSideID 1
      CALL copyBC(Side,NewSides(3,0,1)%sp) !Elem 3, locSideID 1
    CASE(5) ! 4-5-6 side
      CALL copyBC(Side,NewSides(1,M-1,6)%sp) !Elem 1, locSideID 6
      CALL copyBC(Side,NewSides(2,M-1,6)%sp) !Elem 2, locSideID 6
      CALL copyBC(Side,NewSides(3,M-1,6)%sp) !Elem 3, locSideID 6
    END SELECT
  END IF
  Side=>Side%nextElemSide
END DO  ! WHILE(ASSOCIATED(Side))

!inline new elements in element list
aElem=>Elem
nextElem=>Elem%nextElem
prevElem=>Elem%prevElem
DO i=1,3; DO l=0,M-1
  aElem%nextElem          => NewElems(i,l)%EP
  aElem%nextElem%prevElem => aElem
  aElem                   => aElem%nextElem
END DO; END DO
!take out elem
IF(ASSOCIATED(nextElem)) THEN
  aElem%nextElem=>nextElem
  nextElem%prevElem=>aElem
END IF
CALL DeleteElem(firstElem,Elem)
Elem=>aElem

IF(ASSOCIATED(prevElem))THEN
  NewElems(1,0)%ep%prevElem=>prevElem
  prevElem%nextElem=>NewElems(1,0)%ep
ELSE
  firstElem=>NewElems(1,0)%ep
END IF
END SUBROUTINE SplitPentaToHexa

END MODULE MOD_SplitToHex
