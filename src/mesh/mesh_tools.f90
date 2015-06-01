#include "defines.f90"
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

INTERFACE PostDeform
  MODULE PROCEDURE PostDeform
END INTERFACE

PUBLIC::CountSplines,NetVisu,BCVisu
PUBLIC::chkspl_surf,chkspl_vol
PUBLIC::PostDeform
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
WRITE(UNIT_stdOut,'(132("~"))')

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
WRITE(UNIT_stdOut,'(132("~"))')
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
VarNames(3)='nx'
VarNames(4)='ny'
VarNames(5)='nz'

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

nVal=3
ALLOCATE(VarNames(nVal))
VarNames(1)='elemind'
VarNames(2)='Jacobian'
VarNames(3)='scaledJacobian'

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
    END SELECT
  END IF
  Elem=>Elem%nextElem
END DO

CALL Visualize(3,nVal,Nplot,nCurveds,VarNames,Coord,Values,FileString)

DEALLOCATE(VarNames,Coord,Values,xNode,x,xt1,xt2,xt3,Jac)
CALL Timer(.FALSE.)
END SUBROUTINE chkSpl_Vol


SUBROUTINE PostDeform()
!===================================================================================================================================
! input x,y,z node coordinates are transformed by a smooth (!) mapping to new x,y,z coordinates 
!===================================================================================================================================
!MODULE INPUT VARIABLES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:firstElem,tElem,MeshPostDeform
!MODULE OUTPUT VARIABLES
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER          :: aElem   ! ?
REAL                         :: x_loc(3)! ?
INTEGER                      :: iNode   ! ?
!===================================================================================================================================
IF(MeshPostDeform.EQ.0) RETURN
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,'(A)')'POST DEFORM THE MESH...'
CALL Timer(.TRUE.)
aElem=>FirstElem
DO WHILE(ASSOCIATED(aElem))
  DO iNode=1,aElem%nNodes
    aElem%Node(iNode)%np%tmp=-1
  END DO
  DO iNode=1,aElem%nCurvedNodes
    aElem%CurvedNode(iNode)%np%tmp=-1
  END DO
  aElem=>aElem%nextElem
END DO ! WHILE(ASSOCIATED(aElem))
aElem=>FirstElem
DO WHILE(ASSOCIATED(aElem))
  DO iNode=1,aElem%nNodes
    IF(aElem%Node(iNode)%np%tmp.EQ.-1)THEN
      x_loc(:)=aElem%Node(iNode)%np%x(:)
      aElem%Node(iNode)%np%x(:)=PostDeformFunc(x_loc)
      aElem%Node(iNode)%np%tmp=0
    END IF
  END DO
  DO iNode=1,aElem%nCurvedNodes
    IF(aElem%CurvedNode(iNode)%np%tmp.EQ.-1)THEN
      x_loc(:)=aElem%CurvedNode(iNode)%np%x(:)
      aElem%CurvedNode(iNode)%np%x(:)=PostDeformFunc(x_loc)
      aElem%CurvedNode(iNode)%np%tmp=0
    END IF
  END DO
  aElem=>aElem%nextElem
END DO ! WHILE(ASSOCIATED(aElem))
CALL Timer(.FALSE.)
END SUBROUTINE PostDeform



FUNCTION PostDeformFunc(X) RESULT(xout)
!===================================================================================================================================
! input x,y,z node coordinates are transformed by a smooth (!) mapping to new x,y,z coordinates 
!===================================================================================================================================
!MODULE INPUT VARIABLES
USE MOD_Globals
USE MOD_Mesh_Vars,ONLY:MeshPostDeform,PostDeform_R0  
!MODULE OUTPUT VARIABLES
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: x(3) ! contains original xyz coords
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL            :: Xout(3)    ! contains new XYZ position 
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL            :: dx(3),xi(3),alpha,dx1(3),dx2(3)
!===================================================================================================================================
dx=0.
SELECT CASE(MeshPostDeform)
CASE(1) ! 2D box, x,y in [-1,1]^2, to cylinder with radius PostDeform_R0
        ! all points outside [-1,1]^2 will be mapped directly to a circle (p.e. 2,2 => sqrt(0.5)*PostDeform_R0*(2,2) )
        ! inside [-1,1]^2 and outside [-0.5,0.5]^2 there will be a blending from a circle to a square
        ! the inner square [-0.5,0.5]^2 will be a linear blending of the bounding curves
  IF((ABS(x(1)).LT.0.5).AND.(ABS(x(2)).LT.0.5))THEN !inside [-0.5,0.5]^2
    !right side at x=0.5
    dx1(1)=0.5*SQRT(2.)*COS(0.25*Pi*x(2)/0.5)-0.5
    dx1(2)=(0.5*SQRT(2.)*SIN(0.25*Pi*x(2)/0.5)-x(2))
    !upper side at y=0.5
    dx2(1)=0.5*SQRT(2.)*SIN(0.25*Pi*x(1)/0.5)-x(1)
    dx2(2)=0.5*SQRT(2.)*COS(0.25*Pi*x(1)/0.5)-0.5
    alpha=0.5
    dx(1:2)=alpha*2*(dx1(1:2)*(/x(1),ABS(x(1))/)+dx2(1:2)*(/ABS(x(2)),x(2)/))
  ELSE !outside [-0.5,0.5]^2
    IF(ABS(x(2)).LT.ABS(x(1)))THEN !left and right quater
      dx(1)=x(1)*SQRT(2.)*COS(0.25*Pi*x(2)/x(1))-x(1)
      dx(2)=x(1)*SQRT(2.)*SIN(0.25*Pi*x(2)/x(1))-x(2)
    ELSEIF(ABS(x(2)).GE.ABS(x(1)))THEN !upper and lower quater
      dx(1)=x(2)*SQRT(2.)*SIN(0.25*Pi*x(1)/x(2))-x(1)
      dx(2)=x(2)*SQRT(2.)*COS(0.25*Pi*x(1)/x(2))-x(2)
    END IF
    IF((ABS(x(1)).LE.1).AND.(ABS(x(2)).LE.1))THEN
      !alpha=(1.-PRODUCT(1.-x(1:2)**2)) !only <1 inside [-1,1]^2
      alpha=MAX(ABS(x(1)),ABS(x(2)))
    ELSE !outside [-1,1]^2
      alpha=1.
    END IF
    dx(1:2)=alpha*dx(1:2)
  END IF
  xout(1:2)=PostDeform_R0*SQRT(0.5)*(x(1:2)+dx(1:2))
  xout(3)=x(3)
END SELECT

END FUNCTION PostDeformFunc

END MODULE MOD_Mesh_Tools
