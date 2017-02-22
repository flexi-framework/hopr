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
#include "hopr.h"
MODULE MOD_zcorrection
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
INTERFACE OrientElemsToZ 
  MODULE PROCEDURE OrientElemsToZ
END INTERFACE

INTERFACE zcorrection 
  MODULE PROCEDURE zcorrection
END INTERFACE

PUBLIC::OrientElemsToZ,zcorrection
!===================================================================================================================================

CONTAINS
SUBROUTINE OrientElemsToZ()
!===================================================================================================================================
! called when using zcorrection, orientates zeta  of element with z direction 
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:FirstElem,tElem,tNodePtr,tSidePtr,tSide
USE MOD_Mesh_Vars,ONLY:N
USE MOD_Mesh_Vars,ONLY:ElemSideMapping
USE MOD_Basis_Vars,ONLY:HexaMapInv
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER         :: Elem  ! ?
TYPE(tSide),POINTER         :: Side  ! ?
TYPE(tNodePtr),POINTER      :: Nodes(:),CurvedNodes(:,:,:)  ! ?
TYPE(tSidePtr),POINTER      :: Sides(:)  ! ?
REAL                        :: scalprod,dir(3,3)  ! ?
INTEGER                     :: whichdir,switch,Map(8,2,3),MapCurved(3,0:N,0:N,0:N,2,3),MapSides(6,2,3)  ! ?
INTEGER                     :: iNode,iSide,i,j,k,im,jm,km   ! ?
INTEGER                     :: counter(2,3)  ! ?
LOGICAL                     :: found
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)') ' ORIENT ELEMENTS IN Z FOR 3D->2D compatibility...'
counter=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  IF(Elem%nNodes.NE.8) THEN
    WRITE(UNIT_stdOut,*) 'WARNING! Rotation only implemented for hexas, skipping!'
    RETURN
  END IF
  Side=>Elem%firstSide
  DO iSide=1,6
    IF(.NOT.ASSOCIATED(Side)) STOP 'Problem Side Pointer not associated'
    IF(ASSOCIATED(Side%CurvedNode))THEN
      WRITE(UNIT_stdOut,*)  &
           'WARNING! Rotation only implemented for linear hexas and hexas only with volume curvednodes (no sides curved)!'
      RETURN
    END IF
  END DO
  Elem=>Elem%nextElem
END DO
ALLOCATE(Nodes(8),Sides(6),CurvedNodes(0:N,0:N,0:N))
!Mapping old node position to new one
! Map(iNode,switch,wichdir)
!whichdir=1,switch=1/2
Map(:,1,1)=(/1,5,6,2,4,8,7,3/)
Map(:,2,1)=(/6,2,1,5,7,3,4,8/)
!whichdir=2,switch=1/2
Map(:,1,2)=(/1,4,8,5,2,3,7,6/)
Map(:,2,2)=(/6,7,3,2,5,8,4,1/)
!whichdir=3,switch=1/2
Map(:,1,3)=(/1,2,3,4,5,6,7,8/)
Map(:,2,3)=(/8,7,6,5,4,3,2,1/)
DO i=0,N
  DO j=0,N
    DO k=0,N
      !whichdir=1,switch=1/2
      MapCurved(:,i,j,k,1,1)=(/j,k,i/)
      MapCurved(:,i,j,k,2,1)=(/N-j,k,N-i/)
      !whichdir=2,switch=1/2
      MapCurved(:,i,j,k,1,2)=(/k,i,j/)
      MapCurved(:,i,j,k,2,2)=(/N-k,i,N-j/)
      !whichdir=3,switch=1/2
      MapCurved(:,i,j,k,1,3)=(/i,j,k/)
      MapCurved(:,i,j,k,2,3)=(/i,N-j,N-k/)
    END DO !k
  END DO !j
END DO !i
!whichdir=1,switch=1/2
MapSides(:,1,1)=(/5,1,4,6,2,3/)
MapSides(:,2,1)=(/3,1,2,6,4,5/)
!whichdir=2,switch=1/2
MapSides(:,1,2)=(/2,5,6,3,1,4/)
MapSides(:,2,2)=(/4,5,1,3,6,2/)
!whichdir=3,switch=1/2
MapSides(:,1,3)=(/1,2,3,4,5,6/)
MapSides(:,2,3)=(/6,4,3,2,5,1/)

Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  dir(1,:)=NORMALIZE(Elem%node(2)%np%x-Elem%node(1)%np%x,3)
  dir(2,:)=NORMALIZE(Elem%node(4)%np%x-Elem%node(1)%np%x,3)
  dir(3,:)=NORMALIZE(Elem%node(5)%np%x-Elem%node(1)%np%x,3)
  found=.FALSE.
  DO whichdir=1,3
    scalprod=SUM(dir(whichdir,:)*(/0.,0.,1./))
    IF(ABS(scalprod).GT.0.95) THEN
      found=.TRUE.     
      EXIT
    END IF
  END DO
  IF(.NOT.found) THEN
    STOP "Element with no axis aligned with z direction found, cannot perform z orientation!"
  END IF
  switch=2
  IF(scalprod.GT.0) switch=1 !else is set to 2
  counter(switch,whichdir)=counter(switch,whichdir)+1
  DO iNode=1,8
    Nodes(iNode)%np=>Elem%Node(iNode)%np
  END DO
  DO iNode=1,8
    Elem%Node(Map(iNode,switch,whichdir))%np=>Nodes(iNode)%np
  END DO
  Side=>Elem%firstSide
  DO iSide=1,6
    Sides(iSide)%sp=>Side
    Side=>Side%nextElemSide
  END DO
  DO iSide=1,6
    NULLIFY(Sides(iSide)%sp%nextElemSide)
  END DO
  !Reassign Side pointers 
  im=MapSides(1,switch,whichdir)
  Elem%firstSide=>Sides(im)%sp
  Side=>Elem%firstSide
  Side%locSide=1
  DO iSide=2,6
    im=MapSides(iSide,switch,whichdir)
    Side%nextElemSide=>Sides(im)%sp
    Side%nextElemSide%locSide=iSide
    Side=>Side%nextElemSide
  END DO
  !Reassign Side Node pointers
  Side=>Elem%firstSide
  DO iSide=1,6
    DO iNode=1,4
      Side%Node(iNode)%np=>Elem%Node(ElemSideMapping(Elem%nNodes,iSide,iNode))%np
    END DO !iNode
    Side=>Side%nextElemSide
  END DO

  IF(ASSOCIATED(Elem%CurvedNode))THEN
    DO i=0,N
      DO j=0,N
        DO k=0,N
          CurvedNodes(i,j,k)%np=>Elem%curvedNode(HexaMapInv(i,j,k))%np
        END DO
      END DO
    END DO
    DO i=0,N
      DO j=0,N
        DO k=0,N
          im=MapCurved(1,i,j,k,switch,whichdir)
          jm=MapCurved(2,i,j,k,switch,whichdir)
          km=MapCurved(3,i,j,k,switch,whichdir)
          Elem%curvedNode(HexaMapInv(im,jm,km))%np=>CurvedNodes(i,j,k)%np
        END DO
      END DO
    END DO
  END IF!Curvednodes
  Elem=>Elem%nextElem
END DO
!WRITE(*,*)'direction,switch,number of elements'
!DO whichdir=1,3
!  DO switch=1,2
!    WRITE(*,*)whichdir,switch,counter(switch,whichdir)
!  END DO
!END DO
!check
!counter=0
!Elem=>firstElem
!DO WHILE(ASSOCIATED(Elem))
!  dir(1,:)=NORMALIZE(Elem%node(2)%np%x-Elem%node(1)%np%x,3)
!  dir(2,:)=NORMALIZE(Elem%node(4)%np%x-Elem%node(1)%np%x,3)
!  dir(3,:)=NORMALIZE(Elem%node(5)%np%x-Elem%node(1)%np%x,3)
!  DO whichdir=1,3
!    scalprod=SUM(dir(whichdir,:)*(/0.,0.,1./))
!    IF(ABS(scalprod).GT.0.5) EXIT
!  END DO
!  switch=2
!  IF(scalprod.GT.0) switch=1 !else is set to 2
!  counter(switch,whichdir)=counter(switch,whichdir)+1
!  Elem=>Elem%nextElem
!END DO
!WRITE(*,*)'direction,switch,number of elements after orientation'
!DO whichdir=1,3
!  DO switch=1,2
!    WRITE(*,*)whichdir,switch,counter(switch,whichdir)
!  END DO
!END DO

CALL Timer(.FALSE.)
END SUBROUTINE OrientElemsToZ
SUBROUTINE zcorrection()
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:FirstElem,tElem,tSide
USE MOD_Mesh_Vars,ONLY:zstart,dz,zLength,zPeriodic,nElemsZ
USE MOD_Mesh_Vars,ONLY:BoundaryType
USE MOD_Mesh_Tolerances,ONLY:SAMEPOINT
USE MOD_Mesh_Vars,ONLY:N,nMeshElems
USE MOD_Basis_Vars,ONLY:HexaMapInv
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER         :: Elem,lastElem  ! ?
TYPE(tSide),POINTER         :: Side,zminusSide,zplusSide  ! ?
INTEGER                     :: whichdirArr(nMeshElems),orientArr(nMeshElems)  ! ?
INTEGER                     :: Map(4,2),SideMap(2)  ! ?
INTEGER                     :: p,q,l,l1,l2,iSide  ! ?
INTEGER                     :: switch, switch2, zcounter, whichdir  ! ?
REAL                        :: scalprod,dir(3,3), displ1(3,0:N), displ2(3,0:N), displ(3,0:N)  ! ?
LOGICAL                     :: firstLayer  ! ?
LOGICAL                     :: dominant  ! ?
LOGICAL                     :: found     ! ?
INTEGER                     :: iNode,fNode,nPeriodicSides  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)') ' PERFORMING Z CORRECTION...'
CALL Timer(.TRUE.)
dz=zLength/REAL(nElemsZ)
DO l=0,N
 displ1(:,l) = REAL(l)  /REAL(N) * (/0.,0.,dz/)
 displ2(:,l) = REAL(N-l)/REAL(N) * (/0.,0.,dz/)
END DO

Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  IF(Elem%nNodes.NE.8) THEN
    WRITE(UNIT_stdOut,*) 'WARNING! zcorrection only implemented for hexas, skipping z correction!'
    RETURN
  END IF
  Elem=>Elem%nexTElem
END DO

IF(zPeriodic) nPeriodicSides=0
! set elem%tmp and all nodes%tmp =-1
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Elem%tmp=-1
  DO l=1,8
    Elem%node(l)%np%tmp=-1
  END DO
  ! set curved nodes %tmp =-1
  DO l=1,Elem%nCurvedNodes
    Elem%curvedNode(l)%np%tmp=-1
  END DO
  Elem=>Elem%nexTElem
END DO

whichdirArr=0
orientArr  =2
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  dir(1,:)=NORMALIZE(Elem%node(2)%np%x-Elem%node(1)%np%x,3)
  dir(2,:)=NORMALIZE(Elem%node(4)%np%x-Elem%node(1)%np%x,3)
  dir(3,:)=NORMALIZE(Elem%node(5)%np%x-Elem%node(1)%np%x,3)
  found=.FALSE.
  DO whichdir=1,3
    scalprod=SUM(dir(whichdir,:)*(/0.,0.,1./))
    IF(ABS(scalprod).GT.0.5) THEN
      found=.TRUE.
      EXIT
    END IF
  END DO
  IF(.NOT.found) THEN
    STOP "Element with no axis aligned with z direction found, cannot perform z correction!"
  END IF
  whichdirArr(Elem%ind)=whichdir
  IF(scalprod.GT.0) orientArr(Elem%ind)=1 !else is set to 2
  Elem=>Elem%nexTElem
END DO

! main correction loop
lastElem=>firstElem
DO WHILE(ASSOCIATED(lastElem))
  ! check if already marked
  IF(lastElem%tmp.EQ.0) THEN
    lastElem=>lastElem%nextElem
    CYCLE
  END IF
  Elem=>lastElem
  whichdir=whichdirArr(Elem%ind)
  SELECT CASE(whichdir)
    CASE(1)
      SideMap=(/5,3/)
    CASE(2)
      SideMap=(/2,4/)
    CASE(3)
      SideMap=(/1,6/)
  END SELECT

  ! switch that determines postive or negative orientation
  switch=orientArr(Elem%ind)
  ! find local side in -z direction
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%LocSide .EQ. SideMap(orientArr(Elem%ind))) EXIT
    Side=>Side%nextElemSide
  END DO

  ! go to cell with the bc associated in (-z)
  DO WHILE(.NOT. ASSOCIATED(Side%BC))
    Side=>Side%connection
    Elem=>Side%Elem
    Side => oppSide(Side)
  END DO !while not associated BC
  zminusSide=>Side
  ! now we are at a side and the corresponding cell at z=zmin, start the correction from here  
  ! correction  of nodes all elements in opposite direction (+z)
  Side => oppSide(Side)
  Elem=>Side%Elem
  zcounter=1
  firstLayer=.TRUE.
  DO  !correct element
    whichdir=whichdirArr(Elem%ind)
    switch=orientArr(Elem%ind)
   
    !Map for corner nodes (CGNS standard) 
    SELECT CASE(whichdir)
    CASE(1) !dir1 is periodic direction
      Map(:,1)=(/1,5,8,4/) !ximinus
      Map(:,2)=(/2,6,7,3/) !xiplus
    CASE(2) !dir2 is periodic direction
      Map(:,1)=(/1,2,6,5/) !etaminus
      Map(:,2)=(/4,3,7,8/) !etaplus
    CASE(3) !dir3 is periodic direction
      Map(:,1)=(/1,2,3,4/) !zetaminus
      Map(:,2)=(/5,6,7,8/) !zetaplus
    END SELECT
    ! set z coordinate of side at zmin to zstart
    IF(firstLayer)THEN
      DO l=1,4
        IF(Elem%node(Map(l,switch))%np%tmp .NE. 0) THEN
          Elem%node(Map(l,switch))%np%x(3)=zstart
          Elem%node(Map(l,switch))%np%tmp=0
        END IF
      END DO
      !curved
      IF(ASSOCIATED(Elem%CurvedNode))THEN
        iSide=0
        IF(switch.EQ.2) iSide=N
        DO p=0,N; DO q=0,N
          SELECT CASE(whichdir)
          CASE(1)
            Elem%curvedNode(HexaMapInv(iSide,p,q))%np%x(3)=zstart
            Elem%curvedNode(HexaMapInv(iSide,p,q))%np%tmp=0
          CASE(2)
            Elem%curvedNode(HexaMapInv(p,iSide,q))%np%x(3)=zstart
            Elem%curvedNode(HexaMapInv(p,iSide,q))%np%tmp=0
          CASE(3)
            Elem%curvedNode(HexaMapInv(p,q,iSide))%np%x(3)=zstart
            Elem%curvedNode(HexaMapInv(p,q,iSide))%np%tmp=0
          END SELECT
        END DO; END DO
      END IF !curved
    END IF ! first layer

    ! correction of corner points
    switch2=3-switch  !switch2=2 if switch=1, else switch2=1 if switch=2
    ! set x,y coordinates of local z+ face to those of local z- face
    DO l=1,4
      IF(Elem%node(Map(l,switch2))%np%tmp .NE. 0) THEN
        Elem%node(Map(l,switch2))%np%x=Elem%node(Map(l,switch))%np%x+(/0.,0.,dz/)
        Elem%node(Map(l,switch2))%np%tmp=0
      END IF
    END DO

    !curved
    IF(ASSOCIATED(Elem%CurvedNode))THEN
      IF(switch.EQ.1)THEN !scalprod is positive 
        l1=1; l2=N;   iSide=0; displ=displ1
      ELSE !scalprod<0
        l1=0; l2=N-1; iSide=N; displ=displ2
      END IF
      DO p=0,N; DO q=0,N
        DO l=l1,l2
          SELECT CASE(whichdir)
          CASE(1)
            IF(Elem%curvedNode(HexaMapInv(l,p,q))%np%tmp .NE. 0) THEN
               Elem%curvedNode(HexaMapInv(l,p,q))%np%x=Elem%curvedNode(HexaMapInv(iSide,p,q))%np%x + displ(:,l)
               Elem%curvedNode(HexaMapInv(l,p,q))%np%tmp=0
            END IF
          CASE(2)
            IF(Elem%curvedNode(HexaMapInv(p,l,q))%np%tmp .NE. 0) THEN
               Elem%curvedNode(HexaMapInv(p,l,q))%np%x=Elem%curvedNode(HexaMapInv(p,iSide,q))%np%x + displ(:,l) 
               Elem%curvedNode(HexaMapInv(p,l,q))%np%tmp=0
            END IF
          CASE(3)
            IF(Elem%curvedNode(HexaMapInv(p,q,l))%np%tmp .NE. 0) THEN
               Elem%curvedNode(HexaMapInv(p,q,l))%np%x=Elem%curvedNode(HexaMapInv(p,q,iSide))%np%x + displ(:,l)
               Elem%curvedNode(HexaMapInv(p,q,l))%np%tmp=0
            END IF
          END SELECT
        END DO
      END DO; END DO ! q,p
    END IF !curved

    !mark corrected element
    Elem%tmp=0
    firstLayer= .FALSE.
    IF(ASSOCIATED(Side%BC)) THEN
      IF(zcounter .NE. nElemsZ) STOP "Specified nElemsZ not correct."
      zplusSide=>Side
      IF(zPeriodic)THEN
        nPeriodicSides=nPeriodicSides+2
        !set BC to periodic between the two sides
        zminusSide%BC%BCType=1
        zminusSide%BC%BCalphaInd=1
        BoundaryType(zminusSide%BC%BCIndex,1)=1
        BoundaryType(zminusSide%BC%BCIndex,4)=1
        zplusSide%BC%BCType=1
        zplusSide%BC%BCalphaInd=-1
        BoundaryType(zplusSide%BC%BCIndex,1)=1
        BoundaryType(zplusSide%BC%BCIndex,4)=-1
        zminusSide%connection=>zplusSide
        zplusSide%connection=>zminusSide
        !adjustOrientednodes
        fNode=0
        DO iNode=1,zminusSide%nNodes
          IF(SAMEPOINT(zminusSide%Node(1)%np%x(:)+(/0.,0.,zLength/),zplusSide%Node(iNode)%np%x)) THEN
            fNode=iNode
          END IF
        END DO
        IF(fNode.EQ.0)STOP 'Problem with zPeriodic'
        dominant=.FALSE.
        IF(zminusSide%Elem%ind .LT. zplusSide%Elem%ind) dominant=.TRUE.
        IF(dominant) THEN
          DO iNode=1,zplusSide%nNodes
            zminusSide%orientedNode(iNode)%np=>zminusSide%Node(iNode)%np
            zplusSide%orientedNode(iNode)%np=>zplusSide%Node(fNode)%np
            fNode=fNode-1
            IF(fNode .LT.1) fNode=fNode+zplusSide%nNodes
          END DO
        ELSE
          DO iNode=1,zplusSide%nNodes
            zminusSide%orientedNode(iNode)%np=>zminusSide%Node(fNode)%np
            zplusSide%orientedNode(iNode)%np=>zplusSide%Node(iNode)%np
            fNode=fNode-1
            IF(fNode .LT.1) fNode=fNode+zplusSide%nNodes
          END DO
        END IF
      END IF !zPeriodic
      EXIT
    ELSE
      ! jump to next Element
      Side=>Side%connection
      Elem=>Side%Elem
      Side => oppSide(Side)
      zcounter=zcounter+1
    END IF
     
  END DO !correct element
  lastElem=>lastElem%nextElem
END DO
IF(zPeriodic) WRITE(UNIT_StdOut,*)'Number of Periodic sides built:' ,nPeriodicSides
CALL Timer(.FALSE.)
END SUBROUTINE zcorrection

FUNCTION oppSide(Side)
!===================================================================================================================================
! find oppsite side within the local cell
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tSide
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tSide),POINTER,INTENT(IN) :: Side  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tSide),POINTER            :: OppSide  ! ?
! normals on a Meshnode
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: oppSideMap(1:6)  ! ?
!===================================================================================================================================
oppSideMap(:) = (/6,4,5,2,3,1/)
OppSide=>Side%Elem%firstSide
DO WHILE(ASSOCIATED(OppSide))
  IF(OppSide%LocSide .EQ. oppSideMap(Side%LocSide)) EXIT
  OppSide=>OppSide%nextElemSide
END DO
END FUNCTION oppSide 


END MODULE MOD_zcorrection
