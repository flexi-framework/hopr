#include "defines.f90"
MODULE MOD_Mesh_Jacobians
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
INTERFACE CheckJacobians
  MODULE PROCEDURE CheckJacobians
END INTERFACE

PUBLIC:: CheckJacobians
!===================================================================================================================================

CONTAINS
SUBROUTINE CheckJacobians()
!===================================================================================================================================
! Solve the Jacobian determinant and sharp boundaries for J_min and J_max
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars ,ONLY:FirstElem,nMeshElems
USE MOD_Mesh_Vars ,ONLY:tElem
USE MOD_Basis_Vars,ONLY:Vdm_analyze_Hexa,D_analyze_Hexa
USE MOD_Mesh_Vars,ONLY:N
USE MOD_Basis_Vars,ONLY:nAnalyze
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: nNodes,iNode,iElem,i  ! ?
INTEGER              :: nAnalyze_3  ! ?
REAL                 :: maxJac,minJac  ! ?
REAL,ALLOCATABLE     :: xNode(:,:),x(:,:),xt1(:,:),xt2(:,:),xt3(:,:),Jac(:)  ! ?
REAL                 :: scaledJac(nMeshElems)  ! ?
INTEGER              :: scaledJacStat(0:10)  ! ?
TYPE(tElem),POINTER  :: aElem  ! ?
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)') 'Check Elements Jacobian...'
nAnalyze_3=nAnalyze**3
nNodes=(N+1)**3
ALLOCATE(xNode(nNodes,3))
ALLOCATE(x(nAnalyze_3,3))
ALLOCATE(xt1(nAnalyze_3,3))
ALLOCATE(xt2(nAnalyze_3,3))
ALLOCATE(xt3(nAnalyze_3,3))
ALLOCATE(Jac(nAnalyze_3))

scaledJac=1.
iElem=0
! Allocate Curved Node Coordinates, Array for Jacobians_Determinants J, Bezier Coefficients B
! Loop over all Elements
aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  IF (aElem%nNodes.NE. 8) THEN
    WRITE(*,*)' WARNING: checkJacobian only for Hexas...',aElem%nnodes
    aElem=>aElem%nextElem
    CYCLE
  END IF ! Only Elements with nNodes
  iElem=iElem+1 ! Element Counter
  ! check if Element is curved
  IF(aElem%nCurvedNodes.GT.0) THEN ! curved Element
    ! check if there are the right number of curved Nodes
    IF (aElem%nCurvedNodes.NE.nNodes) THEN
      WRITE(*,*) 'wrong Number of curved Nodes in Elem:',iElem
      WRITE(*,*) 'should be:',nNodes,',but is:',aElem%nCurvedNodes
      STOP
    END IF

    ! fill Node coordinates array X with curved nodes
    DO iNode=1,nNodes
      xNode(iNode,1)=aElem%CurvedNode(iNode)%np%x(1)
      xNode(iNode,2)=aElem%CurvedNode(iNode)%np%x(2)
      xNode(iNode,3)=aElem%CurvedNode(iNode)%np%x(3)
    END DO

    ! compute Jacobian at Analyzenodes
    x  = MATMUL(Vdm_Analyze_Hexa,xNode(1:nNodes,:))
    xt1= MATMUL(D_Analyze_Hexa(:,:,1),xNode(1:nNodes,:))
    xt2= MATMUL(D_Analyze_Hexa(:,:,2),xNode(1:nNodes,:))
    xt3= MATMUL(D_Analyze_Hexa(:,:,3),xNode(1:nNodes,:))
    DO iNode=1,nAnalyze_3
      Jac(iNode)=SUM(xt1(iNode,:)*CROSS(xt2(iNode,:),xt3(iNode,:)))
    END DO
    maxJac=MAXVAL(ABS(Jac))
    minJac=MINVAL(Jac)
   
    scaledJac(iElem)=minJac/maxJac
    ! Check that Nodal Jacobians are positive
    !DO iNode=1,nAnalyze_3
    !  IF (Jac(iNode)/maxJac.LE.0.01) THEN
    !     WRITE(*,*)'Negative Jacobian, iELem=',iElem,'Jac=',Jac(iNode),'maxJac',maxJac
    !  END IF
    !END DO
  END IF
  aElem=>aElem%nextElem
END DO
IF (iElem.NE.nMeshElems) THEN
  WRITE(*,*) 'iElem /= nMeshElems',iElem,nMeshElems
END IF
DEALLOCATE(x,xNode,xt1,xt2,xt3,Jac)

! Error Section
IF(ANY(scaledJac.LE.0.01))THEN
  ALLOCATE(X(1,3))
  OPEN(UNIT=100,FILE='Jacobian_Error.out',STATUS='UNKNOWN',ACTION='WRITE')
  WRITE(100,*) 'Corrupt Elems Found, i.e. the Jacobian is negative, '
  WRITE(100,*)
  WRITE(100,'(4(A12,4X))') 'iElem','xBary','yBary','zBary'
  WRITE(100,*)
  iElem=0
  aElem=>firstElem
  DO WHILE(ASSOCIATED(aElem))
    iElem=iElem+1
    IF(scaledJac(iElem).LT.0.01) THEN
      X(:,:)=0.
      DO iNode=1,aElem%nNodes
        X(1,:)=X(1,:)+aElem%node(inode)%np%x
      END DO
      X(1,:)=X(1,:)/REAL(aElem%nNodes)
      WRITE(100,'(I12,3(4X,F12.8))') iElem, X(1,:)
    END IF
    aElem=>aElem%nextElem
  END DO
  CLOSE(UNIT=100)
  WRITE(*,*) 'WARNING!!! Elements with negative Jacobian found, see Jacobian_Error.out for details.'
  DEALLOCATE(X)
END IF
scaledJacStat(:)=0

DO iElem=1,nMeshElems
  i=CEILING(MAX(0.,scaledJac(iElem)*10))
  scaledJacStat(i)=scaledJacStat(i)+1 
END DO
WRITE(Unit_StdOut,'(A)') ' Number of element with scaled Jacobians ranging between:'
WRITE(Unit_StdOut,'(A)') '   <  0.0  <  0.1  <  0.2  <  0.3  <  0.4  <  0.5  <  0.6  <  0.7  <  0.8  <  0.9  <  1.0 '
DO i=0,10
  WRITE(Unit_StdOut,'(I6,X,A1)',ADVANCE='NO')scaledJacStat(i),'|'
END DO
WRITE(Unit_StdOut,'(A1)')' '
 

CALL Timer(.FALSE.)
WRITE(UNIT_StdOut,'(132("="))')
WRITE(UNIT_StdOut,*)

END SUBROUTINE CheckJacobians

END MODULE MOD_Mesh_Jacobians
