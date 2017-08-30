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
MODULE MOD_Mesh_Tools
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
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE CountSplines
  MODULE PROCEDURE CountSplines
END INTERFACE

INTERFACE NetVisu
  MODULE PROCEDURE NetVisu
END INTERFACE

INTERFACE BCVisu
  MODULE PROCEDURE BCVisu
END INTERFACE

INTERFACE chkSpl_surf
  MODULE PROCEDURE chkSpl_surf
END INTERFACE

INTERFACE chkSpl_vol
  MODULE PROCEDURE chkSpl_vol
END INTERFACE

INTERFACE SetTempMarker
  MODULE PROCEDURE SetTempMarker
END INTERFACE

INTERFACE CheckMortarWaterTight
  MODULE PROCEDURE CheckMortarWaterTight
END INTERFACE

PUBLIC::CountSplines,NetVisu,BCVisu
PUBLIC::chkspl_surf,chkspl_vol
PUBLIC::SetTempMarker
PUBLIC::CheckMortarWaterTight
!===================================================================================================================================

CONTAINS
SUBROUTINE CountSplines()
!===================================================================================================================================
! Count all Splines in MESH, only for DEBUGGING!!!
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:nBoundarySplines,nPeriodicSplines,nInnerSplines
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER            :: Elem  ! ?
TYPE(tSide),POINTER            :: Side  ! ?
INTEGER                        :: splineCounter(3) ! 1: Boundary splines, 2: periodic splines, 3: 
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, SAVE                  :: CallCount=0   ! Counter for calls of this subroutine, used for 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL                        :: hasChanged  ! ?
!===================================================================================================================================
IF(.NOT. Logging) RETURN
splineCounter(:)    = 0
hasChanged          = .FALSE.
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(ASSOCIATED(Side%CurvedNode)) THEN
      IF(ASSOCIATED(Side%BC)) THEN
        IF(Side%BC%BCType.EQ.1) THEN !Periodic, assuming that BCType=1=periodic
          splineCounter(2)=splineCounter(2)+1
        ELSE !Boundary
          splineCounter(1)=splineCounter(1)+1
        ENDIF
      ELSE !inner Spline
        splineCounter(3)=splineCounter(3)+1
      ENDIF
    ENDIF !spline
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO
IF(splineCounter(1) .NE. nBoundarySplines)THEN
  hasChanged=.TRUE.
  nBoundarySplines=splineCounter(1)
ENDIF
IF(splineCounter(2) .NE. nPeriodicSplines)THEN
  hasChanged=.TRUE.
  nPeriodicSplines=splineCounter(2)
ENDIF
IF(splineCounter(3) .NE. nInnerSplines)THEN
  hasChanged=.TRUE.
  nInnerSplines=splineCounter(3)
ENDIF
IF(hasChanged .OR. (CallCount .EQ. 0))THEN
  CallCount = CallCount+1
  WRITE(UNIT_logOut,'(132("~"))')
  WRITE(UNIT_logOut,'(A,I12,A)',ADVANCE='NO')'> countSplines(',CallCount,')'
  WRITE(UNIT_logOut,'(A,I12)',ADVANCE='NO')' | #Boundary Splines: ',splineCounter(1)
  WRITE(UNIT_logOut,'(A,I12)',ADVANCE='NO')' | #Periodic Splines: ',splineCounter(2)
  WRITE(UNIT_logOut,'(A,I12,A)')' | #Inner Splines:    ',splineCounter(3),' <'
  WRITE(UNIT_logOut,'(132("~"))')
END IF ! hasChanged
END SUBROUTINE CountSplines


SUBROUTINE NetVisu()
!===================================================================================================================================
! Debug visualization of the Grid. Needs only grid points and connectivity.
! writes a hybrid mesh in Tecplot format.
! the debug mesh is not the original mesh, as every element is shrinked to its barycenter. the shrinking can be controlled
! with the paremeter baryscale, default=1.0 = no shrinking (note: baryscale<=1.0)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Output   ,ONLY:Visualize
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER             :: Elem  ! ?
REAL                            :: BaryScale,BaryCoords(3)  ! ?
INTEGER                         :: nVal,iNode,iElem,nElems  ! ?
INTEGER                         :: NodeMap(1:8,4:8)  ! ?
CHARACTER(LEN=255)              :: FileString  ! ?
CHARACTER(LEN=255)              :: VarNames(3)  ! ?
REAL,ALLOCATABLE                :: Coord(:,:,:),Solution(:,:,:)  ! ?

PARAMETER (BaryScale   = 1.0)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'WRITING THE DEBUGMESH...'
CALL Timer(.TRUE.)
! Count elements and nodes and DoF
nElems      = 0
Elem=>firstElem
DO WHILE(ASSOCIATED(elem))
  nElems=nElems+1
  elem=>elem%nextElem
END DO
WRITE(UNIT_stdOut,*)'  #Elements ',nElems
filestring=TRIM(ProjectName)//'_'//'Debugmesh'
nVal=3
VarNames(1)='elemind'
VarNames(2)='zone'
VarNames(3)='Jacobian'

NodeMap=0
!mapping from iNode to i,j,k [0;1]
NodeMap(:,4)=(/1,2,3,3,4,4,4,4/) !tetra
NodeMap(:,5)=(/1,2,4,3,5,5,5,5/) !pyra
NodeMap(:,6)=(/1,2,3,3,4,5,6,6/) !prism
NodeMap(:,8)=(/1,2,4,3,5,6,8,7/) !prism

ALLOCATE(Coord(3,1:8,nElems))
ALLOCATE(Solution(nVal,1:8,nElems))

iElem=0
Elem=>firstElem
DO WHILE(ASSOCIATED(elem))
  iElem=iElem+1
  BaryCoords=0.
  DO iNode=1,Elem%nNodes
    BaryCoords=BaryCoords+Elem%Node(iNode)%np%x/REAL(Elem%nNodes)
  END DO
  DO iNode=1,8
    Coord(:,iNode,iElem)=(Elem%Node(NodeMap(iNode,Elem%nNodes))%np%x-BaryCoords)*BaryScale+BaryCoords
  END DO
  Solution(1,:,ielem)=Elem%ind
  Solution(2,:,ielem)=Elem%zone
  Solution(3,:,ielem)=Elem%detT
  Elem=>elem%nextElem
END DO
CALL Visualize(3,nVal,1,nElems,VarNames,Coord,Solution,FileString)

DEALLOCATE(Coord,Solution)

WRITE(UNIT_stdOut,'(3X,A,A)')' Mesh visualized for debug purposes in file : ',TRIM(filestring)
CALL Timer(.FALSE.)

END SUBROUTINE netVisu


SUBROUTINE BCVisu()
!===================================================================================================================================
! Debug visualization of the boundary conditions in Tecplot format.
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Output   ,ONLY:Visualize
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER             :: Elem  ! ?
TYPE(tSide),POINTER             :: Side   ! ?
REAL                            :: BaryScale,BaryCoords(3)  ! ?
INTEGER                         :: nVal,iNode,iBCside,nBCSides  ! ?
INTEGER                         :: NodeMap(1:4,3:4)  ! ?
CHARACTER(LEN=255)              :: FileString  ! ?
CHARACTER(LEN=255)              :: VarNames(6)  ! ?
REAL,ALLOCATABLE                :: Coord(:,:,:),Solution(:,:,:)  ! ?

PARAMETER (BaryScale   = 1.0)
!===================================================================================================================================
! ----------------------------------------- Visualize boundary conditions ----------------------------------------------------
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'WRITING THE BC MESH...'
CALL Timer(.TRUE.)
! Count boundary elements
nBCSides = 0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(ASSOCIATED(Side%BC))THEN
      nBCSides = nBCSides+1
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO
WRITE(UNIT_stdOut,*)'  #BCSides ',nBCSides
filestring=TRIM(ProjectName)//'_'//'Debugmesh_BC'
! Open/Create file and file header for 2D/3D Debugmesh.
nVal=6
Varnames(1)='BCIndex'
Varnames(2)='BCType'
Varnames(3)='BCCurveIndex'
Varnames(4)='BCState'
Varnames(5)='BCAlpha'
Varnames(6)='ElemID'

NodeMap=0
!mapping from iNode to i,j,k [0;1]
NodeMap(:,3)=(/1,2,3,3/) !tria
NodeMap(:,4)=(/1,2,4,3/) !quad

ALLOCATE(Coord(3,1:4,nBCSides))
ALLOCATE(Solution(nVal,1:4,nBCSides))
iBCSide=0
! Write nodes
Elem=>firstElem
DO WHILE(ASSOCIATED(elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(ASSOCIATED(Side%BC))THEN
      iBCSide=iBCside+1
      BaryCoords=0.
      DO iNode=1,Side%nNodes
        BaryCoords=BaryCoords+Side%Node(iNode)%np%x/REAL(Side%nNodes)
        Side%Node(iNode)%np%tmp=0
      END DO
      DO iNode=1,4
        Coord(:,iNode,iBCSide)=(Side%Node(NodeMap(iNode,Side%nNodes))%np%x-BaryCoords)*BaryScale+BaryCoords
      END DO
      Solution(1,:,iBCside)=Side%BC%BCindex
      Solution(2,:,iBCside)=Side%BC%BCType
      Solution(3,:,iBCside)=Side%CurveIndex
      Solution(4,:,iBCside)=Side%BC%BCState
      Solution(5,:,iBCside)=Side%BC%BCAlphaInd
      Solution(6,:,iBCside)=Elem%ind
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>elem%nextElem
END DO

CALL Visualize(2,nVal,1,nBCSides,VarNames,Coord,Solution,FileString)
DEALLOCATE(Coord,Solution)

WRITE(UNIT_stdOut,'(3X,A,A)')' Boundary mesh visualized for debug purposes in file : ',TRIM(FileString)
CALL Timer(.FALSE.)
END SUBROUTINE BCVisu


SUBROUTINE chkSpl_Surf()
!===================================================================================================================================
! Checks the Spline via Debug output of the Splines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars ,ONLY:N
USE MOD_Basis_Vars,ONLY:nVisu
USE MOD_Basis_Vars,ONLY:Vdm_Visu_tria,D_Visu_tria
USE MOD_Basis_Vars,ONLY:Vdm_Visu_quad,D_Visu_quad
USE MOD_Basis_Vars,ONLY:VisuTriaMapInv,VisuQuadMapInv
USE MOD_Mesh_Vars ,ONLY:tElem,tSide
USE MOD_Mesh_Vars ,ONLY:FirstElem
USE MOD_Output    ,ONLY:Visualize
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: Elem  ! ?
TYPE(tSide),POINTER             :: Side  ! ?
INTEGER                         :: i,j,l,nCurveds,iNode,nNodes,iSide  ! ?
INTEGER                         :: Nplot,Nplot_p1,Nplot_p1_2,Nplot_2  ! ?
INTEGER                         :: nVal   ! ?
CHARACTER(LEN=255)              :: FileString  ! ?
CHARACTER(LEN=255),ALLOCATABLE  :: VarNames(:)  ! ?
REAL,ALLOCATABLE                :: xNode(:,:),x(:,:),xt1(:,:),xt2(:,:)  ! ?
REAL,ALLOCATABLE                :: Coord(:,:,:),Values(:,:,:)   ! ?
REAL                            :: normal(3)  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'WRITE CURVED SURFACE VISUALIZATION...'
CALL Timer(.TRUE.)

nCurveds=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF((ASSOCIATED(Side%CurvedNode)).AND.(Side%CurveIndex.GT.0))THEN
      nCurveds=nCurveds+1
    END IF
    Side=>Side%nextElemSide
  END DO !Side
  Elem=>Elem%nextElem
END DO !Elem
IF(nCurveds.EQ.0)THEN
  WRITE(UNIT_stdOut,'(A)')'Nothing to visualize: No Boundaries with CurveIndex>0 found.'
  CALL Timer(.FALSE.)
  RETURN
ELSE
  WRITE(UNIT_stdOut,'(A,I8)')'   #CurvedSurfaces: ', nCurveds
END IF

FileString=TRIM(ProjectName)//'_'//'SplineSurf'

Nplot=nVisu !number of equidistant points for visualization mesh
Nplot_2=Nplot**2
Nplot_p1=Nplot+1
Nplot_p1_2=Nplot_p1*Nplot_p1

ALLOCATE(xNode((N+1)**2,3))
ALLOCATE(x(  Nplot_p1_2,3))
ALLOCATE(xt1(Nplot_p1_2,3))
ALLOCATE(xt2(Nplot_p1_2,3))

nVal=5
ALLOCATE(VarNames(nVal))
VarNames(1)='CurveIndex'
VarNames(2)='elemind'
VarNames(3)='nvecX'
VarNames(4)='nvecY'
VarNames(5)='nvecZ'

ALLOCATE(Coord(    3,Nplot_p1_2,nCurveds))
ALLOCATE(Values(nVal,Nplot_p1_2,nCurveds))
iSide=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF((ASSOCIATED(Side%CurvedNode)).AND.(Side%CurveIndex.GT.0))THEN
      iSide=iSide+1
      nNodes=Side%nCurvedNodes
      x=0.
      xt1=0.
      xt2=0.
      xNode=0.
      DO iNode=1,nNodes
        IF(ASSOCIATED(Side%curvedNode(iNode)%np))THEN
          xNode(iNode,:)=Side%curvedNode(iNode)%np%x
        ELSE
          CALL abort(__STAMP__, &
           'Curved node array has not the right size.',999,999.) ! check required due to intel compiler error (02)
        END IF
      END DO
      SELECT CASE(Side%nNodes)
      CASE(3)
        x  = MATMUL(Vdm_Visu_tria,xNode(1:nNodes,:))
        xt1= MATMUL(D_Visu_tria(:,:,1),xNode(1:nNodes,:))
        xt2= MATMUL(D_Visu_tria(:,:,2),xNode(1:nNodes,:))
        DO j=0,Nplot
          DO i=0,Nplot
            l=VisuTriaMapInv(i,j)
            Coord(:,1+i+(Nplot+1)*j,iSide)=x(l,:)
            normal=NORMALIZE(CROSS(xt1(l,:),xt2(l,:)),3)
            Values(3:5,1+i+(Nplot+1)*j,iSide)=-normal
          END DO
        END DO
      CASE(4)
        x  = MATMUL(Vdm_Visu_quad,xNode(1:nNodes,:))
        xt1= MATMUL(D_Visu_quad(:,:,1),xNode(1:nNodes,:))
        xt2= MATMUL(D_Visu_quad(:,:,2),xNode(1:nNodes,:))
        DO j=0,Nplot
          DO i=0,Nplot
            l=VisuQuadMapInv(i,j)
            Coord(:,1+i+(Nplot+1)*j,iSide)=x(l,:)
            normal=NORMALIZE(CROSS(xt1(l,:),xt2(l,:)),3)
            Values(3:5,1+i+(Nplot+1)*j,iSide)=-normal
          END DO
        END DO
      END SELECT
      Values(1,:,iSide)=Side%CurveIndex
      Values(2,:,iSide)=Elem%ind
    END IF !CurveIndex >0
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO

CALL Visualize(2,nVal,Nplot,nCurveds,VarNames,Coord,Values,FileString)

DEALLOCATE(VarNames,Coord,Values,xNode,x,xt1,xt2)
CALL Timer(.FALSE.)
END SUBROUTINE chkSpl_surf


SUBROUTINE chkSpl_Vol()
!===================================================================================================================================
! Checks the Spline via Debug output of the Splines
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars ,ONLY:tElem
USE MOD_Mesh_Vars ,ONLY:FirstElem
USE MOD_Mesh_Vars ,ONLY:N
USE MOD_Basis_Vars,ONLY:Vdm_Visu_Hexa,D_Visu_Hexa
USE MOD_Basis_Vars,ONLY:VisuHexaMapInv
USE MOD_Basis_Vars,ONLY:nVisu
USE MOD_Output    ,ONLY:Visualize
USE MOD_Output_vars,ONLY:Visu_sJ_limit
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: Elem  ! ?
INTEGER                         :: i,j,k,l,nCurveds,iNode,nNodes,iElem  ! ?
INTEGER                         :: Nplot,Nplot_p1,Nplot_p1_3  ! ?
INTEGER                         :: nVal   ! ?
CHARACTER(LEN=255)              :: FileString  ! ?
CHARACTER(LEN=255),ALLOCATABLE  :: VarNames(:)  ! ?
REAL,ALLOCATABLE                :: xNode(:,:),x(:,:),xt1(:,:),xt2(:,:),xt3(:,:),Jac(:)  ! ?
REAL,ALLOCATABLE                :: Coord(:,:,:),Values(:,:,:)   ! ?
REAL                            :: smaxJac  ! ?
INTEGER,ALLOCATABLE             :: ElemMap(:),ElemMapInv(:)  ! for visu_sJ_limit
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'WRITE CURVED VOLUME VISUALIZATION...'
CALL Timer(.TRUE.)


! count curved Elems
nCurveds=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  IF(ASSOCIATED(Elem%CurvedNode)) nCurveds=nCurveds+1
  Elem=>Elem%nextElem
END DO
IF(nCurveds.EQ.0)THEN
  WRITE(UNIT_stdOut,'(A)')'Nothing to visualize: No Curved volumes found.'
  CALL Timer(.FALSE.)
  RETURN
ELSE
  WRITE(UNIT_stdOut,'(A,I8)')'   #CurvedElements: ', nCurveds
END IF


FileString=TRIM(ProjectName)//'_'//'SplineVol'

Nplot=nVisu !number of equidistant points for visualization mesh
Nplot_p1=Nplot+1
Nplot_p1_3=Nplot_p1*Nplot_p1*Nplot_p1

ALLOCATE(xNode((N+1)**3,3))
ALLOCATE(x(  Nplot_p1_3,3))
ALLOCATE(xt1(Nplot_p1_3,3))
ALLOCATE(xt2(Nplot_p1_3,3))
ALLOCATE(xt3(Nplot_p1_3,3))
ALLOCATE(Jac(Nplot_p1_3))

nVal=4
ALLOCATE(VarNames(nVal))
VarNames(1)='elemind'
VarNames(2)='Jacobian'
VarNames(3)='scaledJacobian'
VarNames(4)='scaledJacElem'

ALLOCATE(Coord(    3,Nplot_p1_3,nCurveds))
ALLOCATE(Values(nVal,Nplot_p1_3,nCurveds))

iElem=0
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  IF(ASSOCIATED(Elem%CurvedNode)) THEN
    SELECT CASE(Elem%nNodes)
    CASE(8)
      iElem=iElem+1
      nNodes=Elem%nCurvedNodes
      x=0.
      xt1=0.
      xt2=0.
      xNode=0.
      DO iNode=1,nNodes
        IF(ASSOCIATED(Elem%curvedNode(iNode)%np))THEN
          xNode(iNode,:)=Elem%curvedNode(iNode)%np%x
        ELSE
          CALL abort(__STAMP__, &
           'Curved node array has not the right size.',999,999.) ! check required due to intel compiler error (02)
        END IF
      END DO
      x  = MATMUL(Vdm_Visu_Hexa,xNode(1:nNodes,:))
      xt1= MATMUL(D_Visu_Hexa(:,:,1),xNode(1:nNodes,:))
      xt2= MATMUL(D_Visu_Hexa(:,:,2),xNode(1:nNodes,:))
      xt3= MATMUL(D_Visu_Hexa(:,:,3),xNode(1:nNodes,:))
      DO iNode=1,Nplot_p1_3
        Jac(iNode)=SUM(xt1(iNode,:)*CROSS(xt2(iNode,:),xt3(iNode,:)))
      END DO
      sMaxJac=1./MAXVAL(Jac)
      DO k=0,Nplot
        DO j=0,Nplot
          DO i=0,Nplot
            l=VisuHexaMapInv(i,j,k)
            Coord(:,l,iElem)=x(l,:)
            Values(2,l,iElem)=Jac(l)
            Values(3,l,iElem)=Jac(l)*sMaxJac
          END DO
        END DO
      END DO
      Values(1,:,iElem)=Elem%ind
      Values(4,:,iElem)=MINVAL(Values(3,:,iElem))
    END SELECT
  END IF
  Elem=>Elem%nextElem
END DO
IF(Visu_sJ_limit.LT.1.)THEN
  ALLOCATE(ElemMap(nCurveds))
  ElemMap=0
  j=0
  DO iElem=1,nCurveds
    IF(MINVAL(Values(3,:,iElem)).LE.Visu_sJ_limit)THEN
      j=j+1
      ElemMap(iElem)=j
    END IF
  END DO
  IF(j.GT.0)THEN
    WRITE(UNIT_stdOut,'(A,I12,A,F6.3,A)')'Found',j,' elements with sJ<= ',Visu_sJ_limit,' to visualize...'
    ALLOCATE(ElemMapInv(j))
    j=0
    DO iElem=1,nCurveds
      IF(ElemMap(iElem).GT.0) THEN
        j=j+1
        ElemMapInv(j)=iElem
      END IF
    END DO
    CALL Visualize(3,nVal,Nplot,j,VarNames,Coord(:,:,ElemMapInv),Values(:,:,ElemMapInv),FileString)
    DEALLOCATE(ElemMapInv)
  END IF !j>0
  DEALLOCATE(ElemMap)
ELSE
  CALL Visualize(3,nVal,Nplot,nCurveds,VarNames,Coord,Values,FileString)
END IF !Visu_sJ_limit < 1

DEALLOCATE(VarNames,Coord,Values,xNode,x,xt1,xt2,xt3,Jac)
CALL Timer(.FALSE.)
END SUBROUTINE chkSpl_Vol


SUBROUTINE SetTempMarker(Elem,value,whichMarker)
!===================================================================================================================================
! Set temp markers equal to zero for a given element (corner nodes, sides, edges, curveds)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:tElem,tSide,tEdge,N
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER,INTENT(IN) :: Elem  ! ?
INTEGER,INTENT(IN)             :: value ! value to set tmp to
LOGICAL,INTENT(IN),OPTIONAL    :: whichMarker(8) ! 1/2: elem corner/curved, 3/4: side corner/curved, 5/6: edge corner/curved
                                                 ! 7: oriented nodes, 8: connection
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tSide),POINTER         :: Side  ! ?
TYPE(tEdge),POINTER         :: Edge  ! ?
INTEGER                     :: iNode,i ! ?
LOGICAL                     :: whichMarkerLoc(8)
!===================================================================================================================================
whichMarkerLoc=.TRUE.
whichMarkerLoc(7:8)=.FALSE.
IF(PRESENT(whichMarker)) whichMarkerLoc=whichMarker

IF(whichMarkerLoc(1))THEN
  DO iNode=1,Elem%nNodes
    Elem%Node(iNode)%np%tmp = value
  END DO !iNodes
END IF
IF(whichMarkerLoc(2).AND.ASSOCIATED(Elem%CurvedNode))THEN
  DO iNode=1,Elem%nCurvedNodes
    Elem%curvedNode(iNode)%np%tmp = value
  END DO
END IF

Side=>Elem%firstSide
DO WHILE(ASSOCIATED(Side))
  DO iNode=1,Side%nNodes
    IF(WhichMarkerLoc(3)) Side%Node(iNode)%np%tmp=value
    IF(WhichMarkerLoc(7)) Side%OrientedNode(iNode)%np%tmp=value
    IF(ASSOCIATED(Side%edge(iNode)%edp))THEN
      Edge=>Side%edge(iNode)%edp
      IF(WhichMarkerLoc(5))THEN
        Edge%Node(1)%np%tmp=value
        Edge%Node(2)%np%tmp=value
      END IF
      IF(WhichMarkerLoc(6).AND.ASSOCIATED(Edge%CurvedNode))THEN
        DO i=1,N+1
          Edge%curvedNode(i)%np%tmp=value
        END DO
      END IF
    END IF
  END DO

  IF(whichMarkerLoc(4).AND.ASSOCIATED(Side%CurvedNode))THEN
    DO iNode=1,Side%nCurvedNodes
      Side%curvedNode(iNode)%np%tmp=value
    END DO
  END IF
  !periodic side, only /= for connect!
  IF(Side%tmp2.GT.0.AND.ASSOCIATED(Side%Connection).AND.WhichMarkerLoc(8))THEN
    !dummy side found
    DO iNode=1,Side%nNodes
      Side%connection%Node(iNode)%np%tmp=value
    END DO
  END IF
  Side=>Side%nextElemSide
END DO !associated(side)
END SUBROUTINE SetTempMarker


SUBROUTINE checkMortarWatertight()
!===================================================================================================================================
! Checks if surface normals of mortars are defined such that the mesh is freestream preserving =^ watertight
! builds a 1D basis to change equidistant -> gauss points (0:N_GP) and then use tensor-product gauss 
! for differentiation and integration. n_GP=N should be exact, since normal vector is of degree (2*N-1 ,2*N-1)
! since its a dot product of two polynomials of degree (N-1,N) * (N,N-1)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars ,ONLY: tElem,tSide,VV
USE MOD_Mesh_Vars ,ONLY: FirstElem
USE MOD_Mesh_Vars ,ONLY: N
USE MOD_Mesh_Basis,ONLY: PackGeo,isOriented,getFlip
USE MOD_Mesh_Tolerances,ONLY: SAMEPOINT
USE MOD_Basis1D   ,ONLY: LegendreGaussNodesAndWeights
USE MOD_Basis1D   ,ONLY: PolynomialDerivativeMatrix
USE MOD_Basis1D   ,ONLY: BarycentricWeights
USE MOD_Basis1D   ,ONLY: InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: Elem  ! ?
TYPE(tSide),POINTER             :: Side,master,slave  ! ?
INTEGER                         :: p,q,N_GP
INTEGER                         :: flip,pf,qf
REAL                            :: xEq(0:N),wBaryEq(0:N),dir
REAL,ALLOCATABLE                :: xGP(:),wGP(:),DGP(:,:),VdmEqToGP(:,:)
REAL                            :: xGeoSide(3,0:N,0:N),xGeoSide2(3,0:N,0:N),XCheck(3),vec(3)
REAL                            :: NsurfBig(3),NsurfSmall(3,4)
REAL                            :: NsurfErr,MaxNsurfErr
INTEGER                         :: waterTight,nMortars
LOGICAL                         :: fail
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'CHECK IF MORTARS ARE WATERTIGHT...'
CALL Timer(.TRUE.)
DO p=0,N
  xEq(p)=2.*REAL(p)/REAL(N) -1.
END DO !i
CALL BarycentricWeights(N,xEq,wBaryEq)
N_GP=N
ALLOCATE(xGP(0:N_GP),wGP(0:N_GP),DGP(0:N_GP,0:N_GP),VdmEqToGP(0:N_GP,0:N))
CALL LegendreGaussNodesAndWeights(N_GP,xGP,wGP)
CALL PolynomialDerivativeMatrix(N_GP,xGP,DGP)
CALL InitializeVandermonde(N,N_GP,wBaryEq,xEq,xGP,VdmEqtoGP)

!IF(SUM(ABS(xGP-MATMUL(VdmEqToGP,xEq))).GT.1.0E-12) STOP 'PROBLEM WITH VANDERMONDE'
!IF(SUM(ABS(MATMUL(DGP,(1+0*xGP)))).GT.1.0E-12) STOP 'PROBLEM WITH DMATRIX (constant)'
!IF(ABS(SUM(ABS(MATMUL(DGP,xGP)))-REAL(N+1)).GT.1.0E-12) STOP 'PROBLEM WITH DMATRIX( linear)'
!IF(ABS(REAL(2*N+1)*SUM((xGP-0.1)**REAL(2*N)*wGP)-(0.9)**REAL(2*N+1)+(-1.1)**REAL(2*N+1)).GT.1.0E-12) &
!   STOP 'PROBLEM WITH INTEGRATION OF N'

WaterTight=0
nMortars=0
maxNsurfErr=0.

! check mortar sides only
Elem=>firstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%MortarType.GT.0) THEN
      ! check mortar sides

      nMortars=nMortars+1
      CALL PackGeo(N,Side,XgeoSide)

      ! In case of periodic BCs always add displacement vector. Pay attention to direction!
      IF(ASSOCIATED(Side%BC))THEN
        IF(Side%BC%BCType.EQ.1.AND.Side%BC%BCalphaInd.NE.0)THEN
          dir=1.*SIGN(1,Side%BC%BCalphaInd)
          DO q=0,N; DO p=0,N
            XGeoSide(:,p,q)=XGeoSide(:,p,q) + dir*VV(:,ABS(Side%BC%BCalphaInd))
          END DO; END DO
        END IF
      END IF

      NsurfBig=EvalNsurf(XgeoSide)
      NsurfSmall=0.
      DO p=1,Side%nMortars 
        CALL PackGeo(N,Side%MortarSide(p)%sp,XgeoSide)
        NsurfSmall(:,p)=EvalNsurf(XgeoSide)
      END DO 
      NsurfErr= ABS(NsurfBig(1)-SUM(NsurfSmall(1,:))) &
               +ABS(NsurfBig(2)-SUM(NsurfSmall(2,:))) &
               +ABS(NsurfBig(3)-SUM(NsurfSmall(3,:)))
      IF(NsurfErr.GT.1.0E-12) THEN
        ERRWRITE(*,*) &
                 '================> Mortar is not watertight, ERROR=',NsurfErr,' >1.0E-12, big side corners:'
        ERRWRITE(*,*)'Nsurf: ', NsurfBig 
        ERRWRITE(*,*)'   P1: ', Side%OrientedNode(1)%np%x
        ERRWRITE(*,*)'   P2: ', Side%OrientedNode(2)%np%x
        ERRWRITE(*,*)'   P3: ', Side%OrientedNode(3)%np%x
        ERRWRITE(*,*)'   P4: ', Side%OrientedNode(4)%np%x

        DO p=1,Side%nMortars
          ERRWRITE(*,*) &
                   '   ===> small side:',p
          ERRWRITE(*,*)'Nsurf : ', Nsurfsmall(:,p)
          ERRWRITE(*,*)'    P1: ', Side%MortarSide(p)%sp%OrientedNode(1)%np%x
          ERRWRITE(*,*)'    P2: ', Side%MortarSide(p)%sp%OrientedNode(2)%np%x
          ERRWRITE(*,*)'    P3: ', Side%MortarSide(p)%sp%OrientedNode(3)%np%x
          ERRWRITE(*,*)'    P4: ', Side%MortarSide(p)%sp%OrientedNode(4)%np%x
        END DO 
        WaterTight=WaterTight+1
        maxNsurfErr=max(maxNsurfErr,NsurfErr)
      END IF ! diffsurf>0

    ELSEIF(Side%MortarType.EQ.0)THEN

      ! check ALL periodic non-mortar (!) sides as they could potentially be affected by mortars
      IF(ASSOCIATED(Side%BC))THEN
        IF(Side%BC%BCType.EQ.1.AND.Side%BC%BCalphaInd.GT.0)THEN
          IF(isOriented(Side))THEN
            master=>Side; slave=> Side%connection
          ELSE
            master=>Side%connection; slave=>Side
          END IF
          flip=getFlip(master)

          CALL PackGeo(N,master,XgeoSide)
          CALL PackGeo(N,slave ,XgeoSide2)

          ! Compute periodic displacement of barycenters (cannot rely on VV, in case postdeform is used)
          DO p=1,3
            vec(p) = SUM(XgeoSide2(p,:,:)-XGeoSide(p,:,:))/(N+1)**2
          END DO

          !! In case of periodic BCs always add displacement vector. Pay attention to direction!
          !dir=1.*SIGN(1,master%BC%BCalphaInd)
          !vec=dir*VV(:,ABS(Side%BC%BCalphaInd))

          fail=.FALSE.
          DO q=0,N; DO p=0,N
            SELECT CASE(flip)
            CASE(0)
              pf=p; qf=q;
            CASE(1)
              pf=q; qf=p;
            CASE(2)
              pf=N-p; qf=q;
            CASE(3)
              pf=N-q; qf=N-p;
            CASE(4)
              pf=p; qf=N-q;
            END SELECT  

            XCheck=XGeoSide(:,p,q) + vec


            IF ( .NOT.SAMEPOINT(XCheck,XGeoSide2(:,pf,qf)) )THEN
              print*,'---------------'
              print*,Xcheck
              print*,XGeoSide(:,pf,qf)
              print*,XGeoSide2(:,pf,qf)
              fail=.TRUE.
              EXIT
            END IF
          END DO; END DO
          IF(fail)THEN
            ERRWRITE(*,*) &
                     '================> Periodic connection is not watertight! Error: ',ABS(XCheck-XGeoSide2(:,pf,qf))
            ERRWRITE(*,*)'   P1: ', Side%OrientedNode(1)%np%x
            ERRWRITE(*,*)'   P2: ', Side%OrientedNode(2)%np%x
            ERRWRITE(*,*)'   P3: ', Side%OrientedNode(3)%np%x
            ERRWRITE(*,*)'   P4: ', Side%OrientedNode(4)%np%x
            WaterTight=WaterTight+1
          END IF
        END IF
      END IF

    END IF !MortarType>0
    Side=>Side%nextElemSide
  END DO !while Side associated
  Elem=>Elem%nextElem
END DO !while Elem associated


IF(waterTight.GT.0) THEN
  WRITE(UNIT_stdOut,*)'  ERROR: ', waterTight,' Mortar sides of ',nMortars, &
                    ' are not watertight, MaxError= ',MaxNsurfErr
  CALL abort(__STAMP__, &
         'ERROR: mortars not watertight !!!',waterTight,maxNsurfErr)
ELSE
  WRITE(UNIT_stdOut,*)' ==> all mortars are watertight ;)'
END IF

DEALLOCATE(xGP,wGP,DGP,VdmEqToGP)

CALL Timer(.FALSE.)

CONTAINS

  FUNCTION EvalNsurf(Xgeo_in) RESULT(Nsurf)
    USE MOD_ChangeBasis,ONLY:ChangeBasis2D
    !uses N,wGP, VdmEqToGP and DGP from main routine
    IMPLICIT NONE
    !-------------------------------------------------------------
    !INPUT/ OUTPUT VARIABLES
    REAL,INTENT(IN)    :: Xgeo_in(3,0:N,0:N)
    REAL               :: NSurf(3) !non-normalized normal integrated over Side
    !-------------------------------------------------------------
    !LOCAL VARIABLES
    REAL    :: XgeoGP  (3,0:N_GP,0:N_GP)
    REAL    :: dXdxiGP (3)
    REAL    :: dXdetaGP(3)
    REAL    :: nvec(3)
    INTEGER :: i,j,l
    !-------------------------------------------------------------
    CALL ChangeBasis2D(3,N,N_GP,VdmEqToGP,XGeo_in,XGeoGP)
    Nsurf=0.
    DO j=0,N_GP
      DO i=0,N_GP
        !evaluate derivative at i,j point
        dXdxiGP =0.
        dXdetaGP=0.
        DO l=0,N_GP
          dXdxiGP (:) = dXdxiGP (:) + DGP(i,l)*XgeoGP(:,l,j)
          dXdetaGP(:) = dXdetaGP(:) + DGP(j,l)*XgeoGP(:,i,l)
        END DO
        nVec=CROSS(dXdxiGP(:),dXdetaGP(:))
        Nsurf(:)=Nsurf(:)+wGP(i)*wGP(j)*nVec(:)
      END DO !i 
    END DO !j 
  END FUNCTION EvalNsurf

END SUBROUTINE checkMortarWaterTight


END MODULE MOD_Mesh_Tools
