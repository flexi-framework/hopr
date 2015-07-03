#include "defines.f90"
MODULE MOD_CartMesh
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
TYPE tCartesianMesh ! provides data structure for Cartesian mesh
  REAL                                ::       Corner(3,8)            ! Corner nodes of zone
  REAL                                ::       l0(3)                  ! first length (+/-) = direction
  REAL                                ::       factor(3)              ! strech factor (+/-) = direction
  INTEGER                             ::       BCIndex(6)             ! index in UserdefinedBoundaries
  INTEGER                             ::       elemType               ! Element type in zone
  INTEGER                             ::       nElems(3)              ! number of elements in zone
END TYPE tCartesianMesh

TYPE tCartesianMeshPtr
  TYPE(tCartesianMesh),POINTER        ::       CM                     ! cartesian mesh pointer
END TYPE tCartesianMeshPtr

! Public Part ----------------------------------------------------------------------------------------------------------------------
TYPE(tCartesianMeshPtr),POINTER       ::       CartMeshes(:)          ! cartMeshes(nZones) pointer to Cartesian mesh information

PUBLIC::tCartesianMesh,CartMeshes,tCartesianMeshPtr

INTERFACE CartesianMesh
  MODULE PROCEDURE CartesianMesh
END INTERFACE

INTERFACE GetNewHexahedron
  MODULE PROCEDURE GetNewHexahedron
END INTERFACE

PUBLIC::CartesianMesh, GetNewHexahedron
!===================================================================================================================================

CONTAINS
SUBROUTINE GetNewTetrahedron(CornerNode)
!===================================================================================================================================
! Build new tetrahedron for cartesian mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:getNewElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNodePtr),INTENT(IN)                  :: CornerNode(8)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: aElem   ! ?
!===================================================================================================================================
! 1. Tetrahedron
CALL GetNewElem(aElem)
aElem%nNodes=4
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(1)%NP
aElem%Node(2)%NP=>CornerNode(3)%NP
aElem%Node(3)%NP=>CornerNode(4)%NP
aElem%Node(4)%NP=>CornerNode(5)%NP
! 2. Tetrahedron
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=4
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(1)%NP
aElem%Node(2)%NP=>CornerNode(2)%NP
aElem%Node(3)%NP=>CornerNode(3)%NP
aElem%Node(4)%NP=>CornerNode(5)%NP
! 3. Tetrahedron
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=4
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(3)%NP
aElem%Node(2)%NP=>CornerNode(5)%NP
aElem%Node(3)%NP=>CornerNode(7)%NP
aElem%Node(4)%NP=>CornerNode(8)%NP
! 4. Tetrahedron
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=4
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(3)%NP
aElem%Node(2)%NP=>CornerNode(5)%NP
aElem%Node(3)%NP=>CornerNode(6)%NP
aElem%Node(4)%NP=>CornerNode(7)%NP
! 5. Tetrahedron
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=4
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(2)%NP
aElem%Node(2)%NP=>CornerNode(3)%NP
aElem%Node(3)%NP=>CornerNode(5)%NP
aElem%Node(4)%NP=>CornerNode(6)%NP
! 6. Tetrahedron
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=4
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(3)%NP
aElem%Node(2)%NP=>CornerNode(4)%NP
aElem%Node(3)%NP=>CornerNode(5)%NP
aElem%Node(4)%NP=>CornerNode(8)%NP

! Add elements to list
IF(.NOT.ASSOCIATED(FirstElem))THEN
  FirstElem=>aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
ELSE
  aElem%nextElem          => FirstElem
  aElem%nextElem%prevElem => aElem
  FirstElem          => aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
END IF
NULLIFY(aElem)

END SUBROUTINE GetNewTetrahedron


SUBROUTINE GetNewPyramid(CornerNode)
!===================================================================================================================================
! Build new pyramid for cartesian mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem,tNode
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:getNewElem,getNewNode
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNodePtr),INTENT(IN)                  :: CornerNode(8)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: aElem  ! ?
TYPE(tNode),POINTER             :: CenterNode  ! ?
REAL                            :: xm(3)   ! ?
INTEGER                         :: i  ! ?
!===================================================================================================================================
! Get center node
xm=0.
DO i=1,8
  xm=xm+CornerNode(i)%np%x
END DO
xm=0.125*xm
CALL GetNewNode(CenterNode)
CenterNode%x=xm

! 1. Pyramid
CALL GetNewElem(aElem)
aElem%nNodes=5
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(1)%NP
aElem%Node(2)%NP=>CornerNode(2)%NP
aElem%Node(3)%NP=>CornerNode(3)%NP
aElem%Node(4)%NP=>CornerNode(4)%NP
aElem%Node(5)%NP=>CenterNode
! 2. Pyramid
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=5
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(1)%NP
aElem%Node(2)%NP=>CornerNode(5)%NP
aElem%Node(3)%NP=>CornerNode(6)%NP
aElem%Node(4)%NP=>CornerNode(2)%NP
aElem%Node(5)%NP=>CenterNode
! 3. Pyramid
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=5
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(2)%NP
aElem%Node(2)%NP=>CornerNode(6)%NP
aElem%Node(3)%NP=>CornerNode(7)%NP
aElem%Node(4)%NP=>CornerNode(3)%NP
aElem%Node(5)%NP=>CenterNode
! 4. Pyramid
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=5
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(1)%NP
aElem%Node(2)%NP=>CornerNode(4)%NP
aElem%Node(3)%NP=>CornerNode(8)%NP
aElem%Node(4)%NP=>CornerNode(5)%NP
aElem%Node(5)%NP=>CenterNode
! 5. Pyramid
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=5
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(5)%NP
aElem%Node(2)%NP=>CornerNode(8)%NP
aElem%Node(3)%NP=>CornerNode(7)%NP
aElem%Node(4)%NP=>CornerNode(6)%NP
aElem%Node(5)%NP=>CenterNode
! 6. Pyramid
CALL GetNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=5
ALLOCATE(aElem%node(aElem%nnodes))
aElem%Node(1)%NP=>CornerNode(7)%NP
aElem%Node(2)%NP=>CornerNode(8)%NP
aElem%Node(3)%NP=>CornerNode(4)%NP
aElem%Node(4)%NP=>CornerNode(3)%NP
aElem%Node(5)%NP=>CenterNode

! Add elements to list
IF(.NOT.ASSOCIATED(FirstElem))THEN
  FirstElem=>aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
ELSE
  aElem%nextElem          => FirstElem
  aElem%nextElem%prevElem => aElem
  FirstElem          => aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
END IF
NULLIFY(aElem)

END SUBROUTINE GetNewPyramid


SUBROUTINE GetNewPrism(CornerNode)
!===================================================================================================================================
! Build new prism for cartesian mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:getNewElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNodePtr),INTENT(IN)                  :: CornerNode(8)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: aElem   ! ?
!===================================================================================================================================
CALL getNewElem(aElem)
aElem%nNodes=6
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(1)%NP
aElem%Node(2)%NP=>CornerNode(4)%NP
aElem%Node(3)%NP=>CornerNode(5)%NP
aElem%Node(4)%NP=>CornerNode(2)%NP
aElem%Node(5)%NP=>CornerNode(3)%NP
aElem%Node(6)%NP=>CornerNode(6)%NP
CALL getNewElem(aElem%nextElem)
aElem%nextElem%prevElem => aElem
aElem                   => aElem%nextElem
aElem%nNodes=6
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CornerNode(4)%NP
aElem%Node(2)%NP=>CornerNode(8)%NP
aElem%Node(3)%NP=>CornerNode(5)%NP
aElem%Node(4)%NP=>CornerNode(3)%NP
aElem%Node(5)%NP=>CornerNode(7)%NP
aElem%Node(6)%NP=>CornerNode(6)%NP

! Add elements to list
IF(.NOT.ASSOCIATED(FirstElem))THEN
  FirstElem=>aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
ELSE
  aElem%nextElem          => FirstElem
  aElem%nextElem%prevElem => aElem
  FirstElem          => aElem
  DO WHILE(ASSOCIATED(FirstElem%prevElem))
    FirstElem=>FirstElem%prevElem
  END DO
END IF
NULLIFY(aElem)

END SUBROUTINE GetNewPrism


SUBROUTINE GetNewHexahedron(CornerNode)
!===================================================================================================================================
! Build new hexahedron for cartesian mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:getNewElem
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tNodePtr),INTENT(IN)                  :: CornerNode(8)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: aElem   ! ?
INTEGER                         :: i  ! ?
!===================================================================================================================================
CALL getNewElem(aElem)
aElem%nNodes=8
ALLOCATE(aElem%Node(aElem%nNodes))
DO i=1,8
  aElem%Node(i)%NP=>CornerNode(i)%NP
END DO

! Add elements to list
IF(.NOT.ASSOCIATED(FirstElem))THEN
  FirstElem=>aElem
ELSE
  aElem%nextElem          => FirstElem
  aElem%nextElem%prevElem => aElem
  FirstElem          => aElem
END IF
NULLIFY(aElem)
END SUBROUTINE GetNewHexahedron


SUBROUTINE CartesianMesh()
!===================================================================================================================================
! Builds cartesian mesh. Called by fillMesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:nZones,BoundaryType
USE MOD_Mesh_Vars,ONLY:getNewSide,getNewNode,getNewBC
USE MOD_Mesh_Vars,ONLY:deleteSide,deleteNode
USE MOD_Mesh_Basis,ONLY:CreateSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! Some parameters set in INI-File:
! nZones                    : Number of mesh zones
! CartMeshes(iCm)%cm%nElems : Number of elements in zone icm
! CartMesh%Corner                : Corner points of zone
! CartMesh%elemType              : Element type in zone
! CartMesh%BCtype                : Type of boundary condition for zone boundaries
! CartMesh%BCstate               : Boundary state corresponding to boundary type
! CartMesh%BCalphaInd            : Boundary vector index for periodic zone boundaries
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!   
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tCartesianMesh),POINTER :: CartMesh  ! ?
TYPE(tNodePtr)               :: CornerNode(8)   ! ?
TYPE(tNodePtr),POINTER       :: Mnodes(:,:,:)   ! ?
TYPE(tElem),POINTER          :: aElem   ! ?
TYPE(tSide),POINTER          :: aSide  ! ?
REAL                         :: x1(3),x2(3),e1(3)  ! ?
REAL                         :: Co(3,8) ,dx(3),fac(3)   ! ?
REAL                         :: F, dF  ! ?
INTEGER                      :: iZone,i_Dim,iSide  ! ?
INTEGER                      :: i,l,m,n  ! ?
INTEGER                      :: nElems(3),NodeInd  ! ?
INTEGER                      :: ne(3)  ! ?
LOGICAL                      :: onBnd  ! ?
CHARACTER(LEN=32)            :: formatstr
!===================================================================================================================================
WRITE(UNIT_stdOut,*)'Building cartesian mesh'
NodeInd=0
DO iZone=1,nZones
  WRITE(UNIT_stdOut,*)'  build zone no.:' ,iZone
  CartMesh=>CartMeshes(iZone)%CM
  nElems=CartMesh%nElems
  ! Stretching elements with cartmesh%factor and/or carmesh%l0
  DO i_Dim=1,3
    SELECT CASE(i_Dim)
    CASE(1)
      e1(i_Dim)=CartMesh%Corner(i_Dim,2)-CartMesh%Corner(i_Dim,1) ! e1 container for cartmesh-dimensions
    CASE(2)
      e1(i_Dim)=CartMesh%Corner(i_Dim,4)-CartMesh%Corner(i_Dim,1)
    CASE(3)
      e1(i_Dim)=CartMesh%Corner(i_Dim,5)-CartMesh%Corner(i_Dim,1)
    END SELECT
    IF(ABS(CartMesh%l0(i_Dim)) .LT. PP_RealTolerance ) THEN !l0 = 0 = inactive
      IF(ABS(CartMesh%factor(i_Dim)) .LT. PP_RealTolerance ) THEN !fac= 0 = equidistant
        fac(i_Dim)=1.       
      ELSE ! stretched, (nElem,f) given
        fac(i_Dim)=(ABS(CartMesh%factor(i_Dim)))**(SIGN(1.,CartMesh%factor(i_Dim))) ! sign for direction 
      END IF
    ELSE !l0 active
      dx(i_Dim)=e1(i_Dim)/ABS(cartMesh%l0(i_Dim)) ! l/l0
      IF(dx(i_Dim) .LT. (1.-PP_RealTolerance)) THEN ! l0 > l
        !CALL abort(__STAMP__,  & 
         stop  'stretching error, length l0 longer than grid region, in direction '  !,i_Dim,999.)
      END IF
      IF(ABS(CartMesh%factor(i_Dim)) .LT. PP_RealTolerance ) THEN ! fac=0 , (nElem,l0) given, fac calculated
        ! 
      ELSE ! (factor, l0) given, change nElems
        fac(i_Dim)=(ABS(CartMesh%factor(i_Dim)))**(SIGN(1.,CartMesh%factor(i_Dim)*CartMesh%l0(i_Dim))) ! sign for direction 
        IF(fac(i_Dim) .NE. 1.) THEN
          CartMesh%nElems(i_Dim)=NINT(LOG(1.-dx(i_Dim)*(1.-fac(i_Dim))) / LOG(fac(i_Dim))) !nearest Integer 
        ELSE
          CartMesh%nElems(i_Dim)=NINT(dx(i_Dim))
        END IF
        IF (CartMesh%nElems(i_Dim) .LT. 1) CartMesh%nElems(i_Dim)=1
        WRITE(UNIT_stdOut,*)'   -Element number in dir',i_Dim,'changed from', &
                              nElems(i_Dim),'to',CartMesh%nElems(i_Dim) 
        nElems(i_Dim)=CartMesh%nElems(i_Dim) 
      END IF
      IF(nElems(i_Dim).EQ. 1) THEN 
        fac(i_Dim)=1.0          
      ELSEIF(nElems(i_Dim).EQ. 2) THEN 
        fac(i_Dim)=dx(i_Dim)-1.
      ELSE !nElems > 2
        fac(i_Dim)=dx(i_Dim)/nElems(i_Dim) !start value for Newton iteration
        IF (abs(fac(i_Dim)-1.) .GT. PP_RealTolerance) THEN ! NEWTON iteration, only if not equidistant case
          F=1.
          DO WHILE (ABS(F) .GT. PP_RealTolerance)
            F = fac(i_Dim)**nElems(i_Dim) + dx(i_Dim)*(1.-fac(i_Dim)) -1. ! non-linear function
            dF= nElems(i_Dim)*fac(i_Dim)**(nElems(i_Dim)-1) -dx(i_Dim)  !dF/dfac
            fac(i_Dim)= fac(i_Dim) - F/dF
          END DO
          fac(i_Dim)=fac(i_Dim)**SIGN(1.,CartMesh%l0(i_Dim)) ! sign for direction
        END IF
      END IF
      WRITE(UNIT_stdOut,*)'   -stretching factor in dir',i_Dim,'is now', fac(i_Dim)  
    END IF

    IF( ABS((nElems(i_Dim)-1.)*LOG(fac(i_Dim))/LOG(10.)) .GE. 4. )                      &
      !CALL abort(__STAMP__, &
       stop  'stretching error, length ratio > 1.0E4 in direction '  !,i_Dim,999.)
    IF(ABS(fac(i_Dim)-1.) .GT. PP_RealTolerance ) THEN
      dx(i_Dim)=(1.-fac(i_Dim))/(1.-fac(i_Dim)**nElems(i_Dim)) ! first length to start
    ELSE !equidistant case
      dx(i_Dim)=1./ nElems(i_Dim)
    END IF
  END DO  !i_Dim=1,3
  IF(ABS(SUM(fac(:))-3.).GT. PP_RealTolerance) THEN
    WRITE(UNIT_stdOut,*) '   -STRETCHING INFORMATION'
    WRITE(formatstr,'(A5,I1,A8)')'(A19,',3,'(I9,5X))'
    WRITE(UNIT_stdOut,formatstr) '      -nElems(:) : ', CartMesh%nElems(:)
    WRITE(formatstr,'(A5,I1,A11)')'(A19,',3,'(G13.7,1X))'
    WRITE(UNIT_stdOut,formatstr) '      -l0(:)     : ', e1(:)*dx(:)
    WRITE(UNIT_stdOut,formatstr) '      -factor(:) : ', fac(:)
  END IF

  ! The low / up stuff is a remainder of the mpi version and could be replaced...
  ne(:) =CartMesh%nElems
  ALLOCATE(Mnodes(0:ne(1),0:ne(2),0:ne(3)))
  Co=CartMesh%Corner
  ! Build nodes
  e1(:)=dx(:)/fac(:) !container for dx
  !start position x2 in unit square grid
  x2(:)=0.
  ! calculate node postions, x1(:) coordinates in unit square grid
  x1(1)=x2(1)
  dx(1)=e1(1)
  DO l=0,ne(1)
!    x1(1)=REAL(l)/REAL(nSplit(1))
    x1(2)=x2(2)
    dx(2)=e1(2) 
    DO m=0,ne(2)
!      x1(2)=REAL(m)/REAL(nSplit(2))
      x1(3)=x2(3)
      dx(3)=e1(3)
      DO n=0,ne(3)
!        x1(3)=REAL(n)/REAL(nSplit(3))
        CALL GetNewNode(Mnodes(l,m,n)%np)
        NodeInd=NodeInd+1
        Mnodes(l,m,n)%np%ind=NodeInd 
        Mnodes(l,m,n)%np%x=Co(:,1)+x1(1)*(Co(:,2)-Co(:,1))+x1(2)*(Co(:,4)-Co(:,1))+x1(3)*(Co(:,5)-Co(:,1))+ &
                         (Co(:,3)-Co(:,4)-Co(:,2)+Co(:,1))*x1(1)*x1(2)+(Co(:,6)-Co(:,5)-Co(:,2)+Co(:,1))*x1(1)*x1(3)+&
                            (Co(:,8)-Co(:,4)-Co(:,5)+Co(:,1))*x1(2)*x1(3)+(Co(:,7)-Co(:,1)+Co(:,4)+Co(:,5)+Co(:,2)-  &
                                  Co(:,8)-Co(:,3)-Co(:,6))*x1(1)*x1(2)*x1(3)
        dx(3)=dx(3)*fac(3)
        x1(3)=x1(3)+dx(3)
      END DO
      dx(2)=dx(2)*fac(2)
      x1(2)=x1(2)+dx(2)
    END DO
    dx(1)=dx(1)*fac(1)
    x1(1)=x1(1)+dx(1)
  END DO
  DO l=1,ne(1)
    DO m=1,ne(2)
      DO n=1,ne(3)
        CornerNode(1)%np=>Mnodes(l-1,m-1,n-1)%np
        CornerNode(2)%np=>Mnodes(l,m-1,n-1)%np
        CornerNode(3)%np=>Mnodes(l,m,n-1)%np
        CornerNode(4)%np=>Mnodes(l-1,m,n-1)%np
        CornerNode(5)%np=>Mnodes(l-1,m-1,n)%np
        CornerNode(6)%np=>Mnodes(l,m-1,n)%np
        CornerNode(7)%np=>Mnodes(l,m,n)%np
        CornerNode(8)%np=>Mnodes(l-1,m,n)%np
        SELECT CASE(CartMesh%ElemType)
        CASE(104) ! Tetrahedron
          CALL GetNewTetrahedron(CornerNode)
        CASE(105) ! Pyramid
          CALL GetNewPyramid(CornerNode)
        CASE(106) ! Dreiecks-prismon
          CALL GetNewPrism(CornerNode)
        CASE(108) ! Hexaeder
          CALL GetNewHexahedron(CornerNode)
        CASE DEFAULT
          CALL abort(__STAMP__,&
            'The specified element type is not known. Valid types: 104,105,106,108',CartMesh%ElemType)
        END SELECT
      END DO !n
    END DO !m
  END DO !l
  !for Boundary Conditions
  DO l=0,ne(1)
    DO m=0,ne(2)
      DO n=0,ne(3)
        Mnodes(l,m,n)%np%tmp=0 
        IF(n.EQ.0) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+1 !zeta minus
        IF(m.EQ.0) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+20 !eta minus
        IF(l.EQ.ne(1) ) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+300 !xi plus
        IF(m.EQ.ne(2) ) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+4000  !eta plus
        IF(l.EQ.0) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+50000 !xi minus
        IF(n.EQ.ne(3) ) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+600000 !zeta plus
      END DO !n
    END DO !m
  END DO !l

  ! Set element zone
  aElem=>FirstElem
  DO WHILE(ASSOCIATED(aElem))
    IF(aElem%Zone .EQ. 0)THEN
      aElem%Zone=iZone
      CALL CreateSides(aElem,.TRUE.)
    END IF
    aElem=>aElem%nextElem
  END DO

  aElem=>FirstElem
  DO WHILE(ASSOCIATED(aElem))
    IF(aElem%Zone .EQ. iZone)THEN
      aSide=>aElem%firstSide
      DO WHILE(ASSOCIATED(aSide))
        DO iSide=1,6
          onBnd=.TRUE.
          DO i=1, aSide%nNodes
            IF((MOD(aSide%Node(i)%np%tmp,10**iSide)/(10**(iSide-1))).NE.iSide)THEN
              onBnd=.FALSE.
              EXIT ! Loop
            END IF
          END DO
          IF(CartMesh%BCIndex(iSide).EQ.0) onBnd=.FALSE.
          IF(onBnd)THEN
            CALL getNewBC(aSide%BC)
            aSide%BC%BCIndex    = CartMesh%BCIndex(iSide)
            aSide%BC%BCType     = BoundaryType(aSide%BC%BCIndex,1)
            aSide%curveIndex    = BoundaryType(aSide%BC%BCIndex,2)
            aSide%BC%BCstate    = BoundaryType(aSide%BC%BCIndex,3)
            aSide%BC%BCalphaInd = BoundaryType(aSide%BC%BCIndex,4)
          END IF
        END DO  ! iSide=1,6
        aSide=>aSide%nextElemSide
      END DO ! WHILE(ASSOCIATED(aSide))
    END IF  ! aElem%Zone .EQ. iZone
    aElem=>aElem%nextElem
  END DO ! WHILE(ASSOCIATED(aElem))
  DEALLOCATE(Mnodes)
END DO ! iZone
END SUBROUTINE CartesianMesh

END MODULE MOD_CartMesh
