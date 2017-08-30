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
MODULE MOD_Mesh
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
INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

INTERFACE FillMesh
  MODULE PROCEDURE FillMesh
END INTERFACE

PUBLIC::InitMesh,FillMesh
!===================================================================================================================================

CONTAINS
SUBROUTINE InitMesh()
!===================================================================================================================================
! Generates mesh data.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars
USE MOD_Output_Vars,ONLY:OutputInitDone
USE MOD_ReadInTools
USE MOD_CartMesh,ONLY:CartMeshes
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                           :: buf(24)  ! ?
INTEGER                        :: i,j,tmpInt  ! ?
INTEGER,ALLOCATABLE            :: BC_TypeDummy(:,:)  ! ?
CHARACTER(LEN=255),ALLOCATABLE :: BC_NameDummy(:)  ! ?
CHARACTER(LEN=255)             :: DefStr  ! ?
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT MESH...'
IF(.NOT.OutputInitDone)THEN
  CALL abort(__STAMP__,  &
             'ERROR: InitOutput has to be called before InitMesh!')
END IF

! Curved
useCurveds      = GETLOGICAL('useCurveds','.FALSE.')   ! Use / reconstruct spline boundaries
IF(useCurveds)THEN
  N = GETINT('BoundaryOrder','4')-1     ! Get spline order / type
ELSE
  N=1
END IF

! Mesh mode: 0=HDF5 input, 1=Internal cartesian, 2=Gambit, 3=CGNS, 4=ANSA, 5=GMSH
MeshMode=GETINT('Mode')
! Get number of mesh zones in domain
nZones  =GETINT('nZones')
meshIsAlreadyCurved=.FALSE.
IF (MeshMode .EQ. 1) THEN
  ! ---------- INTERNAL CARTESIAN MESH ---------------------------------------------------------------------------------------------
  IF(useCurveds) THEN
    InnerElemStretch = GETLOGICAL('InnerElemStretch','.TRUE.')
  ELSE
    InnerElemStretch = .FALSE.
  END IF
  IF(InnerElemStretch) MeshIsAlreadyCurved=.TRUE.
  ALLOCATE(CartMeshes(nZones))
  DO i=1,nZones
    ALLOCATE(CartMeshes(i)%CM)
    buf=GETREALARRAY('Corner',24)  ! Corner nodes of zone
    DO j=1,8
      CartMeshes(i)%CM%Corner(:,j)=buf((j-1)*3+1:j*3)
    END DO
    CartMeshes(i)%CM%BCIndex =GETINTARRAY('BCIndex',6)  ! Boundary type 
    CartMeshes(i)%CM%nElems  =GETINTARRAY('nElems',3)   ! No of elems in each direction
    CartMeshes(i)%CM%l0      =GETREALARRAY('l0',3,'0.,0.,0.') ! first length (+/-) = direction
    CartMeshes(i)%CM%factor  =GETREALARRAY('factor',3,'0.,0.,0.') ! stretch factor (+/-) = direction
    CartMeshes(i)%CM%ElemType=GETINT('elemtype') ! Element type
    CartMeshes(i)%CM%meshTemplate=GETINT('meshTemplate','1') ! Element type
  END DO
ELSEIF (MeshMode .EQ. 11) THEN
  ! ---------- INTERNAL 1 BLOCK CURVED CARTESIAN MESH -----------------------------------------------------------------
  stretchType = GETINTARRAY('stretchType',3,'0,0,0')
  fac = GETREALARRAY('fac',3,'1.,1.,1.')
  WRITE(DefStr,'(E21.11,A1,E21.11,A1,E21.11)')fac(1),',',fac(2),',',fac(3)
  fac2 = GETREALARRAY('fac2',3,TRIM(DefStr))
  DxmaxToDxmin = GETREALARRAY('DxmaxToDxmin',3,'1.,1.,1.')
  !read in Mesh Parameters
  CurvedMeshType = GETINT('Meshtype','1')
  BCIndex =GETINTARRAY('BCIndex',6)  ! Boundary type 
  nElems  =GETINTARRAY('nElems',3)   ! No of elems in each direction
  ! rescale the factor if innercell stretching is used, so the same factor can be supplied 
  ! in parameter file for cell and innercell stretching
  DO i=1,3
    IF(stretchType(i).EQ.6) THEN
      fac(i)=fac(i)**(nElems(i))
      fac2(i)=fac2(i)**(nElems(i))
    END IF
  END DO
  SELECT CASE(CurvedmeshType)
  CASE(1,11) ! cartesian domain X0,Y0,Z0,DX,DY,DZ
    X0 = GETREALARRAY('X0',3,'0.,0.,0.')
    DX = GETREALARRAY('DX',3,'1.,1.,1.')
  CASE(2) ! trilinear with 8 points
    XP(:,1) = GETREALARRAY('GeometricXP1',3,'0.,0.,0.')
    XP(:,2) = GETREALARRAY('GeometricXP2',3,'1.,0.,0.')
    XP(:,3) = GETREALARRAY('GeometricXP3',3,'1.,1.,0.')
    XP(:,4) = GETREALARRAY('GeometricXP4',3,'0.,1.,0.')
    XP(:,5) = GETREALARRAY('GeometricXP5',3,'0.,0.,1.')
    XP(:,6) = GETREALARRAY('GeometricXP6',3,'1.,0.,1.')
    XP(:,7) = GETREALARRAY('GeometricXP7',3,'1.,1.,1.')
    XP(:,8) = GETREALARRAY('GeometricXP8',3,'0.,1.,1.')
  CASE(3) ! general curved
    ! Used for curved meshes (i.e. meshmode = 3) and determines the curved mapping
    WhichMapping = GETINT('WhichMapping','1')
    SELECT CASE(WhichMapping)
    CASE(1,2) ! special case trilinear, equal to meshmode 2, needs 8 points!! if 2...adds bubble function on the side
      DX = GETREALARRAY('DX',3,'1.,1.,1.')
      XP(:,1) = GETREALARRAY('GeometricXP1',3,'0.,0.,0.')
      XP(:,2) = GETREALARRAY('GeometricXP2',3,'1.,0.,0.')
      XP(:,3) = GETREALARRAY('GeometricXP3',3,'1.,1.,0.')
      XP(:,4) = GETREALARRAY('GeometricXP4',3,'0.,1.,0.')
      XP(:,5) = GETREALARRAY('GeometricXP5',3,'0.,0.,1.')
      XP(:,6) = GETREALARRAY('GeometricXP6',3,'1.,0.,1.')
      XP(:,7) = GETREALARRAY('GeometricXP7',3,'1.,1.,1.')
      XP(:,8) = GETREALARRAY('GeometricXP8',3,'0.,1.,1.')
    CASE(3,4) ! half/full cylinder ...needs r_0 = radius of the cylinder and r_inf=radius of the domain, dy=extension in y dir
      R_0   = GETREAL('R_0','0.5')
      R_INF = GETREAL('R_INF','10.')
      DY    = GETREAL('DZ','2.')
      IF (WhichMapping.EQ.3) PHI   = GETREAL('PHI','180.') !angle of segment (def.: half cylinder, i.e. 180deg)
    CASE(5) ! SINE BUMP  
      DX    = GETREALARRAY('DX',3,'4.,1.,1.') ! half length, width and height of domain
      R_0   = GETREAL('R_0','0.1') ! height of bump
      R_INF = GETREAL('R_INF','1.') ! length of bump
    END SELECT ! WhichMapping
  END SELECT !MEshType
  meshisAlreadyCurved=.TRUE. 
ELSE
  ! ---------- EXTERNAL MESH -------------------------------------------------------------------------------------------------------
  nMeshFiles=1
  SELECT CASE(MeshMode)
  CASE(-1,0)   ! HDF5 input file (-1: old for conversion)
    nMeshFiles=GETINT('nMeshFiles','1') ! Number of mesh files: each mesh file = one zone
    meshIsAlreadyCurved=.TRUE.
  CASE(2)   ! Gambit 
    nMeshFiles=GETINT('nMeshFiles','1') ! Number of mesh files: each mesh file = one zone
    useBinary =GETLOGICAL('useBinary','.FALSE.')  ! Read binary mesh file (.halo files) instead of ASCII
  CASE(3)   ! CGNS mesh
    meshIsAlreadyCurved = GETLOGICAL('meshIsAlreadyCurved','.FALSE.')  ! build curveds by agglomeration (only block structured)
    nskipZ=GETINT('nskipZ','1') !skip every nskip point (=1, nothing is skipped)
    WRITE(DefStr,*) N
    NBlock=GETINT('NBlock',TRIM(DefStr)) !initial polynomial degree of block structured mesh
    nSkip=1 ! TODO: integrate into reader or remove
    nMeshFiles=GETINT('nMeshFiles','1') ! Number of mesh files: each mesh file = one zone
    BugFix_ANSA_CGNS=GETLOGICAL('BugFix_ANSA_CGNS','.FALSE.')
  CASE(5)   ! GMSH file
    meshIsAlreadyCurved=.TRUE.
  END SELECT
  ALLOCATE(MeshFileName(nMeshFiles))
  DO i=1,nMeshFiles
    MeshFileName(i)=GETSTR('FileName') ! Read file name
  END DO
END IF

! Geometry
postScale=GETLOGICAL('postScaleMesh','.FALSE.') ! apply scaling either after readin or before output
MeshScale=GETREAL('meshScale','1.0')           ! scaling factor applied to node coordinates during read in
doScale  = (ABS(MeshScale-1.).GT.PP_RealTolerance)
SpaceQuandt=GETREAL('SpaceQuandt','0.1')

! Curved
IF(useCurveds) THEN
  ! If curved mesh is read in (e.g. CGNS/GMSH/HDF5) curved boundaries can be rebuilt using the linear mesh
  rebuildCurveds=.FALSE.
  IF(meshIsAlreadyCurved) rebuildCurveds=GETLOGICAL('rebuildCurveds','.FALSE.')

  IF((.NOT.meshIsAlreadyCurved).OR.rebuildCurveds)THEN
    curvingMethod=GETINT('curvingMethod','-1') !Type of normals to be used
    SELECT CASE(curvingMethod)
    CASE(-1) !do nothing
    CASE(1) !normals from CAD
      N=3 
      NormalsType=GETINT('NormalsType','1')
      SELECT CASE(normalsType)
      CASE(1) ! normals by reconstruction
      CASE(2) ! CAD normals
        NormalVectFile=GETSTR('NormalVectFile')
        minNormalAngle=GETREAL('minNormalAngle','2.')
        minNormalAngle=COS(minNormalAngle*PI/180.)
      CASE(3) ! exact normals
        N=3
        tmpInt=GETINT('nExactNormals')
        ALLOCATE(ExactNormals(tmpInt*2))
        ExactNormals=GETINTARRAY('exactNormals',tmpInt*2) !(Curvedindex,exactnormaltype,...)
        ! >0 if an analytical normal can be used (see curved.f90,exactNormals)
      CASE DEFAULT
        CALL abort(__STAMP__,&
          'No normal type specified.')
      END SELECT
    CASE(3) ! refined surface elements (only ANSA readin)
      SplitElemFile=GETSTR('SplitElemFile')
      IF(INDEX(TRIM(SplitElemFile),'.cgns').NE.0)THEN
        SplitMeshMode=4
      ELSE
        SplitMeshMode=3
      END IF
    CASE(4) ! ICEM spectral elements
      SpecElemFile=GETSTR('specelemfile')
    CASE DEFAULT
      CALL abort(__STAMP__,&
        'Specified curving method does not exist. =1: NormalVectors, =3: SplitElemFile, =4: SpecElemFile')
    END SELECT
  END IF
  doExactSurfProjection=GETLOGICAL('doExactSurfProjection','.FALSE.')
  IF(doExactSurfProjection)THEN
    nExactSurfFuncs=CNTSTR('exactSurfFunc','0')
    ALLOCATE(ExactSurfFunc(2*nExactSurfFuncs))
    ExactSurfFunc(:)=GETINTARRAY('exactSurfFunc',2*nExactSurfFuncs) !(Curvedindex,exactsurffunc,...)
  END IF

  ! If domain is curved, try to uncurve it and leave only the sides with BCs speciefied curved
  ! -1: deactivated, 0: only boundary is curved, 1: only first element is curved,
  ! 2-n: first n layers from the boundary are curved
  nCurvedBoundaryLayers=GETINT('nCurvedBoundaryLayers','-1')

END IF !usecurveds

BoundaryOrder=N+1

! Boundaries
nUserDefinedBoundaries=CNTSTR('BoundaryName','0')
IF(nUserDefinedBoundaries.NE.CNTSTR('BoundaryType','0'))&
  CALL abort(__STAMP__,&
             'The number of boundary names and boundary types has to be identical.')
IF(nUserDefinedBoundaries .GT. 0)THEN
  ALLOCATE(BoundaryName(nUserDefinedBoundaries))
  ALLOCATE(BoundaryType(nUserDefinedBoundaries,4))
  DO i=1,nUserDefinedBoundaries
    BoundaryName(i)  =GETSTR('BoundaryName')
    BoundaryType(i,:)=GETINTARRAY('BoundaryType',4)
  END DO
END IF

! Get displacement vector for periodic boundaries
nVV=CNTSTR('VV','0')
IF(nVV .GT. 0)THEN
  ALLOCATE(VV(3,nVV))
  DO i=1,nVV
    VV(:,i)=GETREALARRAY('vv',3)
  END DO
END IF  !(nVV .GT. 0)

! 2.5D mesh
MeshDim=3
IF((MeshMode .EQ. 2) .OR. (MeshMode .EQ. 3))THEN
 ! 2.5D mesh: convert 2D mesh to 3D mesh (gambit and cgns mesh only)
  MeshDim=GETINT('MeshDim','3') 
END IF
!for SpecMesh only 2D available
IF((MeshMode .EQ. 6)) MeshDim=2

IF(MeshDim .EQ. 2)THEN
  zLength=GETREAL('zLength')
  nElemsZ=GETINT('nElemsZ')
  dz      = zLength/DBLE(nElemsZ)
  zLength = dz*DBLE(nElemsZ)
  lowerZ_BC=GETINTARRAY('lowerZ_BC',4)
  upperZ_BC=GETINTARRAY('upperZ_BC',4)
  ALLOCATE(BC_NameDummy(nUserDefinedBoundaries))
  ALLOCATE(BC_TypeDummy(nUserDefinedBoundaries,4))
  BC_NameDummy=BoundaryName
  BC_TypeDummy=BoundaryType
  DEALLOCATE(BoundaryName,BoundaryType)
  ! Add lower and upper z boundary
  nUserDefinedBoundaries=nUserDefinedBoundaries+2
  ALLOCATE(BoundaryName(nUserDefinedBoundaries))
  ALLOCATE(BoundaryType(nUserDefinedBoundaries,4))
  BoundaryName(1:nUserDefinedBoundaries-2    )=BC_NameDummy(:)
  BoundaryType(1:nUserDefinedBoundaries-2,1:4)=BC_TypeDummy(:,1:4)
  BoundaryName(nUserDefinedBoundaries-1      )='LowerZ_BC'
  BoundaryType(nUserDefinedBoundaries-1,1:4  )=lowerZ_BC
  BoundaryName(nUserDefinedBoundaries        )='UpperZ_BC'
  BoundaryType(nUserDefinedBoundaries  ,1:4  )=upperZ_BC
  lowerZ_BC_Ind=nUserDefinedBoundaries-1
  upperZ_BC_Ind=nUserDefinedBoundaries
  DEALLOCATE(BC_NameDummy,BC_TypeDummy)
END IF  ! MeshDim=2

! zcorrection
doZcorrection=GETLOGICAL('doZcorrection','.FALSE.')
OrientZ=GETLOGICAL('OrientZ','.FALSE.')
IF(doZcorrection)THEN
  zLength=GETREAL('zLength')
  nElemsZ=GETINT('nElemsZ')
  zstart=GETREAL('zstart')
  zPeriodic=GETLOGICAL('zPeriodic','.FALSE.')
END IF
!splitting of elements
SplitToHex=GETLOGICAL('SplitToHex','.FALSE.')   ! split all elements to hexa
nFineHexa=GETINT('nFineHexa','1')               ! split all hexa by a factor 

nSplitBoxes=CNTSTR('SplitBox','0')
ALLOCATE(SplitBoxes(3,2,nSplitBoxes))
DO i=1,nSplitBoxes
  SplitBoxes(:,:,i)=RESHAPE(GETREALARRAY('SplitBox',6),(/3,2/))
END DO !nSplitBoxes


! for mortarmeshes ensure that small mortar geometry is identical to big mortar geometry
! does not work for periodic mortars, will be set true by default for postdeform!
doRebuildMortarGeometry=GETLOGICAL('doRebuildMortarGeometry','.TRUE.')

meshPostDeform=GETINT('MeshPostDeform','0')
IF(meshPostDeform.GT.0) THEN
  PostDeform_useGL=GETLOGICAL('PostDeform_useGL','.TRUE.')
  PostDeform_R0=GETREAL('PostDeform_R0','1.')
  PostDeform_Lz=GETREAL('PostDeform_Lz','1.')
  PostDeform_sq=GETREAL('PostDeform_sq','0.')
  PostDeform_Rtorus=GETREAL('PostDeform_Rtorus','-1.') !from cyl-> torus
END IF !PostDeform
postConnect=GETINT('postConnect','0')

! Connect
ConformConnect=GETLOGICAL('ConformConnect','.TRUE.') ! Fast connect for conform mesh

! Elem Check
checkElemJacobians=GETLOGICAL('checkElemJacobians','.TRUE.')
jacobianTolerance=GETREAL('jacobianTolerance','1.E-16')


! Initialize variables
NULLIFY(FirstElem,FirstSplitElem)
! Element side mappings following the CGNS standard
nNodesElemSideMapping = 0
ElemSideMapping       = 0
NodeCount      = 0
SideCount      = 0
ElemCount      = 0
! Number of nodes per side
nNodesElemSideMapping(4,:) = (/3,3,3,3,0,0/)
nNodesElemSideMapping(5,:) = (/4,3,3,3,3,0/)
nNodesElemSideMapping(6,:) = (/4,4,4,3,3,0/)
nNodesElemSideMapping(8,:) = (/4,4,4,4,4,4/)
! Map nodes to sides
! Sides of tetrahedron
ElemSideMapping(4,1,:) = (/1,3,2,0/)
ElemSideMapping(4,2,:) = (/1,2,4,0/)
ElemSideMapping(4,3,:) = (/2,3,4,0/)
ElemSideMapping(4,4,:) = (/3,1,4,0/)
! Sides of pyramid
ElemSideMapping(5,1,:) = (/1,4,3,2/)
ElemSideMapping(5,2,:) = (/1,2,5,0/)
ElemSideMapping(5,3,:) = (/2,3,5,0/)
ElemSideMapping(5,4,:) = (/3,4,5,0/)
ElemSideMapping(5,5,:) = (/4,1,5,0/)
! Sides of prism
ElemSideMapping(6,1,:) = (/1,2,5,4/)
ElemSideMapping(6,2,:) = (/2,3,6,5/)
ElemSideMapping(6,3,:) = (/3,1,4,6/)
ElemSideMapping(6,4,:) = (/1,3,2,0/)
ElemSideMapping(6,5,:) = (/4,5,6,0/)
! Sides of hexahedron
ElemSideMapping(8,1,:) = (/1,4,3,2/)
ElemSideMapping(8,2,:) = (/1,2,6,5/)
ElemSideMapping(8,3,:) = (/2,3,7,6/)
ElemSideMapping(8,4,:) = (/3,4,8,7/)
ElemSideMapping(8,5,:) = (/1,5,8,4/)
ElemSideMapping(8,6,:) = (/5,6,7,8/)

MeshInitDone=.TRUE.
WRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
WRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh


SUBROUTINE fillMesh()
!===================================================================================================================================
! Generates mesh data.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars
USE MOD_zcorrection,      ONLY: zcorrection
USE MOD_zcorrection,      ONLY: OrientElemsToZ
USE MOD_SplitToHex,       ONLY: SplitElementsToHex,SplitAllHexa,SplitHexaByBoxes
USE MOD_Output_Vars,      ONLY: DebugVisu,DebugVisuLevel
USE MOD_Curved,           ONLY: SplitToSpline,ReconstructNormals,getExactNormals,deleteDuplicateNormals,readNormals
USE MOD_Curved,           ONLY: create3dSplines,curvedEdgesToSurf,curvedSurfacesToElem
USE MOD_Curved,           ONLY: buildCurvedElementsFromVolume,buildCurvedElementsFromBoundarySides
USE MOD_Curved,           ONLY: readNormals
USE MOD_Curved,           ONLY: ProjectToExactSurfaces
USE MOD_Curved,           ONLY: RebuildMortarGeometry
USE MOD_Mesh_Basis,       ONLY: BuildEdges,ElemGeometry,FindElemTypes
USE MOD_Mesh_Connect,     ONLY: Connect
USE MOD_Mesh_Connect,     ONLY: Connect2DMesh
USE MOD_GlobalUniqueNodes,ONLY: GlobalUniqueNodes
USE MOD_CartMesh,         ONLY: CartesianMesh
USE MOD_CurvedCartMesh,   ONLY: CurvedCartesianMesh
USE MOD_Mesh_Tools,       ONLY: CountSplines,Netvisu,BCvisu,chkspl_surf,chkspl_vol
USE MOD_Mesh_Tools,       ONLY: CheckMortarWaterTight
USE MOD_Mesh_PostDeform,  ONLY: PostDeform
USE MOD_Output_HDF5,      ONLY: WriteMeshToHDF5
USE MOD_Mesh_Jacobians,   ONLY: CheckJacobians
USE MOD_Readin_ANSA
USE MOD_Readin_CGNS
USE MOD_Readin_Gambit
USE MOD_Readin_GMSH
USE MOD_Readin_HDF5
USE MOD_Readin_ICEM
USE MOD_Readin_SpecMesh2D
USE MOD_Output_Vars,ONLY:useSpaceFillingCurve
USE MOD_Output_HDF5,      ONLY: SpaceFillingCurve
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
! Some parameters set in INI-File:
! MeshMode      : Mesh mode: 1-internal, cartesian  2-mesh from Gambit file  3-mesh from CGNS file  4-mesh from ANSA file
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER :: Elem  ! ?
TYPE(tSide),POINTER :: Side    ! ?
LOGICAL             :: curvedFound  ! ?
INTEGER             :: iElem  ! ?
!===================================================================================================================================
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(132("="))')

IF(.NOT.useCurveds) curvingMethod = -1
nMeshElems=0
NULLIFY(FirstElem)

! read mesh
SELECT CASE (MeshMode)
  CASE(0)
    CALL ReadMeshFromHDF5(MeshFileName(1),.TRUE.) ! meshfile
    meshIsAlreadyCurved=.TRUE.
  CASE(100)
    CALL ReadMeshFromHDF5(MeshFileName(1),.FALSE.) ! meshfile
    meshIsAlreadyCurved=.TRUE.
  CASE(1)
    CALL CartesianMesh()  ! Build cartesian mesh
  CASE(11)
    CALL CurvedCartesianMesh()  ! Build cartesian mesh
    meshIsAlreadyCurved=.TRUE.
  CASE(2)
    CALL readGambit()     ! Read .neu file (Gambit)
    IF(MeshDim .EQ. 2) CALL fill25DMesh() ! Build 3D mesh
  CASE(3)
    CALL readCGNSmesh()   ! Read CGNS mesh
    IF(MeshDim .EQ. 2) CALL fill25DMesh() ! Build 3D mesh
  CASE(4)
    CALL readStar()       ! Read Star file (ANSA)
  CASE(5)
    CALL readGMSH()       ! Read .MSH file (GMSH)
  CASE(6)
    MeshDim=3 !overwrite, build first 3D element layer in readin
    CALL readSpecMesh2D()   
    meshIsAlreadyCurved=.TRUE.
    CALL fill25DMesh() 
  CASE DEFAULT
    CALL abort(__STAMP__, &
      'Not known how to construct mesh')
END SELECT

! apply meshscale after readin
IF(doScale.AND..NOT.postScale) CALL ApplyMeshScale(FirstElem)

! get element trafo and ensure right-handed coordinate system
Elem=>FirstElem
DO WHILE (ASSOCIATED(Elem))
  CALL ElemGeometry(Elem)
  Elem=>Elem%nextElem
END DO

IF(SplitToHex)     THEN
  CALL SplitElementsToHex()
  AdaptedMesh=.TRUE.
END IF
IF(nFineHexa.GT.1) THEN
  CALL SplitAllHexa(nFineHexa)
  AdaptedMesh=.TRUE.
END IF
IF(nSplitBoxes.GT.0) THEN
  CALL SplitHexaByBoxes()
  AdaptedMesh=.TRUE.
END IF

! Count elements 
nMeshElems=0
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  nMeshElems=nMeshElems+1
  Elem=>Elem%nextElem
END DO
WRITE(UNIT_stdOut,*)'Number of Elements: ',nMeshElems
ALLOCATE(Elems(nMeshElems))
iElem=0
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  iElem=iElem+1
  Elems(iElem)%ep=>Elem
  Elem%ind=iElem
  Elem=>Elem%nextElem
END DO

! Set local and global mesh maxDX (will be used in ALL following searchmeshes unless redefined)
maxDX=0.
DO iElem=1,nMeshElems
  Side=>Elems(iElem)%ep%firstSide
  DO WHILE(ASSOCIATED(Side))
    !determine local max dx
    maxDX=MAX(maxDX,ABS(Side%Node(2)%np%x-Side%Node(1)%np%x))
    maxDX=MAX(maxDX,ABS(Side%Node(3)%np%x-Side%Node(2)%np%x))
    maxDX=MAX(maxDX,ABS(Side%Node(1)%np%x-Side%Node(3)%np%x))
    IF(Side%nNodes.EQ.4) maxDX=MAX(maxDX,ABS(Side%Node(4)%np%x-Side%Node(3)%np%x))
    Side=>Side%nextElemSide
  END DO
END DO !iElem

! WRITE some visualization, build connectivity and build edge data structure
IF(DebugVisu)THEN
  CALL netVisu() ! visualize mesh
  CALL BCVisu()  ! visualize BC
END IF
IF(useCurveds.AND.Logging)   CALL CountSplines()  ! In case of restart there can be splines

IF(OrientZ) CALL OrientElemsToZ() 

IF(MeshMode .GT. 0)THEN
  CALL Connect(reconnect=.FALSE.,deletePeriodic=.FALSE.)                           ! Create connection between elements
  IF(useCurveds.AND.Logging) CALL CountSplines()  ! In case of restart there can be splines
END IF
CALL buildEdges()

! check if sides to be curved exist
curvedFound=.FALSE.
DO iElem=1,nMeshElems
  Side=>Elems(iElem)%ep%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%curveIndex.GT.0) curvedFound=.TRUE.
    Side=>Side%nextElemSide
  END DO
END DO !iElem
IF(.NOT.curvedFound) curvingMethod=-1

IF(useCurveds)THEN
  IF(meshIsAlreadyCurved.AND..NOT.rebuildCurveds)THEN
    ! if curved nodes have been read in and mesh should not be modified, just distribute nodes
    IF(nCurvedBoundaryLayers.GE.0)THEN
      CALL buildCurvedElementsFromBoundarySides(nCurvedBoundaryLayers)
    ELSE
      CALL buildCurvedElementsFromVolume()
    END IF
  ELSE
    ! Build curved mesh
    SELECT CASE(curvingMethod)
    CASE(1)
      SELECT CASE(normalsType)
      CASE(1)
        CALL ReconstructNormals()
      CASE(2)
        CALL readNormals()                ! Read normal vector file
      CASE(3)
        CALL getExactNormals()
      END SELECT
      CALL deleteDuplicateNormals()
      CALL create3DSplines()              ! Reconstruct curved boundaries
      CALL curvedEdgesToSurf(keepExistingCurveds=.FALSE.)
    CASE(3) ! STAR/ANSA: generate curved mesh from subgrid
      SELECT CASE(SplitMeshMode)
      CASE(3)
        CALL readStar_split(FirstSplitElem,SplitElemFile)
      CASE(4)
        CALL ReadCGNSSurfaceMesh(FirstSplitElem,SplitElemFile)
      CASE DEFAULT
        CALL abort(__STAMP__, &
             'Splitted meshes can only be read in star cd format or cgns!')
      END SELECT
      CALL Connect2DMesh(FirstSplitElem)
      IF(doScale.AND..NOT.postScale) CALL ApplyMeshScale(FirstSplitElem)
      CALL splitToSpline()
      CALL curvedEdgesToSurf(keepExistingCurveds=.FALSE.)
    CASE(4) !!!only in combination with icem cgns meshes!!!! 
      CALL readSpecEdges()
      CALL curvedEdgesToSurf(keepExistingCurveds=.FALSE.)
      ! jetzt muessen bei allen randbedingungen die curved sind, die Edge%curvedNode(:) besetzt sein
    CASE DEFAULT
      ! don't curve
      IF(.NOT.curvedFound) WRITE(UNIT_stdOut,*) 'Curved sides have not been specified, no curving can be performed.'
      WRITE(UNIT_stdOut,*)'No curving has been performed.'
    END SELECT
    CALL curvedSurfacesToElem()
  END IF

  CALL CountSplines()

END IF ! useCurveds

! make all nodes unique
CALL GlobalUniqueNodes(.TRUE.)

! check if sides with mortars exist
mortarFound=.FALSE.
DO iElem=1,nMeshElems
  Side=>Elems(iElem)%ep%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%MortarType.GT.0) THEN
      mortarFound=.TRUE.
      EXIT !do loop
    END IF
    Side=>Side%nextElemSide
  END DO
  IF(mortarFound) EXIT !do loop
END DO !iElem

IF(doExactSurfProjection) CALL ProjectToExactSurfaces()
! get element types
CALL FindElemTypes()
! correct displacement in z (e.g. for periodic/2.5D)
IF(doZcorrection) CALL zCorrection()

CALL CheckNodeConnectivity()

! prepare sorting by space filling curve
! NOTE: SpaceFillingcurve is not used, if existing hdf5 mesh is read in and the sorting should stay identical
IF(useSpaceFillingCurve)THEN
  CALL SpaceFillingCurve(nMeshElems)
END IF

CALL PostDeform()

SELECT CASE(postConnect)
CASE(0) !do nothing
CASE(1) !reconnect all sides
  CALL Connect(reconnect=.TRUE.,deletePeriodic=.FALSE.)                           ! Create connection between elements
CASE(2) !reconnect all sides, delete periodic connections (sides on top by postdeform)
  CALL Connect(reconnect=.TRUE.,deletePeriodic=.TRUE.)                           ! Create connection between elements
END SELECT !postConnect

IF(mortarFound) THEN
  IF(doRebuildMortarGeometry) CALL RebuildMortarGeometry()
  !after rebuild , mortars should be fine, but checking is better:
  CALL CheckMortarWaterTight()
END IF !mortarFound

! apply meshscale before output (default)
IF(doScale.AND.postScale) CALL ApplyMeshScale(FirstElem)
IF(DebugVisu) THEN
  CALL netVisu()  ! visualize mesh
  CALL BCVisu()   ! visualize BC
  IF(useCurveds) THEN
    IF(DebugVisuLevel.GE.1) CALL chkSpl_Surf()
    IF(DebugVisuLevel.GE.2) CALL chkSpl_Vol()
  END IF
END IF
IF(checkElemJacobians) CALL CheckJacobians()

IF(useCurveds .AND. Logging) CALL CountSplines()   ! In case of restart there can be splines
CALL WriteMeshToHDF5(TRIM(ProjectName)//'_mesh.h5')
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.FALSE.)
END SUBROUTINE fillMesh


SUBROUTINE fill25DMesh()
!===================================================================================================================================
! Convert 2D mesh into 3D mesh. Input: linear mesh -> first+second layer (one element), indices of second not set
!                                      curved mesh -> first+second layer (one element)+curved element nodes, all indices set
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tSide,FirstElem
USE MOD_Mesh_Vars,ONLY:MeshDim,SpaceQuandt
USE MOD_Mesh_Vars,ONLY:DZ,zLength,lowerZ_BC,upperZ_BC,lowerZ_BC_Ind,upperZ_BC_Ind,n2dNodes
USE MOD_Mesh_Vars,ONLY:getNewBC,getNewElem,getNewNode
USE MOD_Mesh_Basis,ONLY:createSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER  :: Elem,newElem,FirstElem_loc,lastElem,firstNewElem,lastNewElem  ! ?
TYPE(tSide),POINTER  :: Side,TempSide  ! ?
REAL                 :: zPos  ! ?
INTEGER              :: iNode,BCData(8,5),nNodes  ! ?
LOGICAL              :: LowerBCSide,UpperBCSide,ConnectionSide,copyBC  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')'Fill 2.5D mesh...'

! First we have to set the indices of the new nodes created during mesh readin (only for linear meshes, has to be done already for
! curved ones
Elem => FirstElem
DO WHILE(ASSOCIATED(Elem))
  IF(Elem%nCurvedNodes.EQ.0)THEN
    nNodes=Elem%nNodes/2  ! Only Hexahedrons or prisms in 2.5D case -> 8/6 nodes
    DO iNode=nNodes+1,Elem%nNodes ! nodes 1 to Elem%nNodes/2 on lower z-layer, other nodes on upper z-layer
      Elem%Node(iNode)%NP%Ind=Elem%Node(iNode-nNodes)%NP%Ind+n2dNodes
      IF(Elem%Node(iNode)%NP%x(3) .EQ. 0.)&
        WRITE(UNIT_stdOut,*)'Corrupt 2.5D elem found in fill25DMesh!' ! Security check
    END DO
  END IF
  Elem=>Elem%nextElem
END DO

MeshDim = 3  ! only 3D elements
zPos         = DZ
NULLIFY(newElem,lastElem,firstNewElem,lastNewElem)

FirstElem_loc => FirstElem
Elem          => FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    LowerBCSide    = .TRUE.
    ConnectionSide = .TRUE.
    UpperBCSide    = .TRUE.
    DO iNode=1,Side%nNodes
      IF(ABS( Side%Node(iNode)%np%x(3))          .GT. (PP_MeshTolerance*SpaceQuandt)) LowerBCSide    = .FALSE.
      IF(ABS((Side%Node(iNode)%np%x(3)-zLength)) .GT. (PP_MeshTolerance*SpaceQuandt)) UpperBCSide    = .FALSE.
      IF(ABS((Side%Node(iNode)%np%x(3)-zPos))    .GT. (PP_MeshTolerance*SpaceQuandt)) ConnectionSide = .FALSE.
    END DO
    IF(.NOT.(LowerBCSide.OR.UpperBCSide.OR.ConnectionSide))THEN
      Side=>Side%nextElemSide
      CYCLE ! nothing to do, continue
    END IF

    IF(LowerBCSide)THEN
      CALL getNewBC(Side%BC)
      Side%BC%BCtype     = lowerZ_BC(1) !assign BC to side
      Side%curveIndex    = lowerZ_BC(2) 
      Side%BC%BCstate    = lowerZ_BC(3)
      Side%BC%BCalphaInd = lowerZ_BC(4)
      Side%BC%BCIndex    = lowerZ_BC_Ind
    END IF
    IF(UpperBCSide)THEN
      CALL getNewBC(Side%BC)
      Side%BC%BCtype     = upperZ_BC(1) !assign BC to side
      Side%curveIndex    = upperZ_BC(2) 
      Side%BC%BCstate    = upperZ_BC(3)
      Side%BC%BCalphaInd = upperZ_BC(4)
      Side%BC%BCIndex    = upperZ_BC_Ind
    END IF

    IF(ConnectionSide .AND. .NOT.UpperBCSide)THEN
      ! Copy elements and add z-displacement, always build new nodes, multiple nodes will be removed in rasterfahndung
      CALL getNewElem(newElem)
      newElem%Zone  = Elem%Zone
      newElem%nNodes= Elem%nNodes
      ALLOCATE(newElem%Node(newElem%nNodes))
      newElem%nCurvedNodes=Elem%nCurvedNodes

      IF(Elem%nCurvedNodes.GT.0)THEN
        ALLOCATE(newElem%curvedNode(newElem%nCurvedNodes))
        DO iNode=1,newElem%nCurvedNodes
          CALL getNewNode(newElem%curvedNode(iNode)%np,1)
          newElem%curvedNode(iNode)%np%x  =Side%Elem%curvedNode(iNode)%np%x+(/0.,0.,DZ/)   ! Open side
          newElem%curvedNode(iNode)%np%ind=Side%Elem%curvedNode(iNode)%np%ind+n2dNodes     ! Open side
        END DO
        DO iNode=1,newElem%nNodes
          CALL getNewNode(newElem%Node(iNode)%np,1)
          newElem%Node(iNode)%np%x  =Side%Elem%Node(iNode)%np%x+(/0.,0.,DZ/)  ! Open side
          newElem%Node(iNode)%np%ind=Side%Elem%Node(iNode)%np%ind+n2dNodes    ! Open side
        END DO
      ELSE
        DO iNode=1,newElem%nNodes
          CALL getNewNode(newElem%Node(iNode)%np,1)
          IF(iNode .LE. Side%nNodes)THEN 
            newElem%Node(iNode)%np%x  =Side%Elem%Node(iNode+Side%nNodes)%np%x   ! Connected side
            ! Use elem nodes, side node order could be wrong
            newElem%Node(iNode)%np%ind=Side%Elem%Node(iNode+Side%nNodes)%np%ind ! Connected side
          ELSE
            newElem%Node(iNode)%np%x  =Side%Elem%Node(iNode)%np%x+(/0.,0.,DZ/)  ! Open side
            newElem%Node(iNode)%np%ind=Side%Elem%Node(iNode)%np%ind+n2dNodes    ! Open side
          END IF
        END DO
      END IF
      CALL createSides(newElem,.TRUE.)

      ! Copy BC
      BCData=0
      TempSide=>Elem%firstSide
      DO WHILE(ASSOCIATED(TempSide))
        IF(ASSOCIATED(TempSide%BC))THEN
          copyBC=.FALSE.
          IF(Elem%nNodes .EQ. 6)THEN
            IF(TempSide%LocSide .LT. 4) copyBC=.TRUE.
          ELSE
            IF((TempSide%LocSide .GT. 1) .AND. (TempSide%LocSide .LT. 6)) copyBC=.TRUE.
          END IF
          IF(copyBC)THEN
            BCData(TempSide%LocSide,1)=TempSide%BC%BCType
            BCData(TempSide%LocSide,2)=TempSide%curveIndex
            BCData(TempSide%LocSide,3)=TempSide%BC%BCState
            BCData(TempSide%LocSide,4)=TempSide%BC%BCalphaInd
            BCData(TempSide%LocSide,5)=TempSide%BC%BCIndex
          END IF
        END IF
        TempSide=>TempSide%nextElemSide
      END DO
      TempSide=>newElem%firstSide
      DO WHILE(ASSOCIATED(TempSide))
        IF(BCData(TempSide%LocSide,1) .NE. 0)THEN
          CALL getNewBC(TempSide%BC)
          TempSide%BC%BCType     = BCData(TempSide%LocSide,1)
          TempSide%curveIndex    = BCData(TempSide%LocSide,2)
          TempSide%BC%BCState    = BCData(TempSide%LocSide,3)
          TempSide%BC%BCalphaInd = BCData(TempSide%LocSide,4)
          TempSide%BC%BCIndex    = BCData(TempSide%LocSide,5)
        END IF
        TempSide=>TempSide%nextElemSide
      END DO 
      IF(.NOT.ASSOCIATED(firstNewElem))THEN
        firstNewElem => newElem
        lastNewElem  => newElem
      ELSE
        newElem%nextElem => firstNewElem
        firstNewElem     => newElem
      END IF
    END IF
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
  IF(.NOT.ASSOCIATED(Elem))THEN
    IF(ASSOCIATED(lastElem))THEN
      lastElem%nextElem  => FirstElem
      FirstElem%prevElem => lastElem
      FirstElem          => FirstElem_loc
    END IF
    IF(ASSOCIATED(firstNewElem))THEN
      FirstElem_loc => firstNewElem
      lastElem      => lastNewElem
      NULLIFY(firstNewElem,lastNewElem)
      Elem => FirstElem_loc  ! Next layer
      zPos =zPos+DZ
    END IF
  END IF
END DO

CALL Timer(.FALSE.)
END SUBROUTINE fill25DMesh

SUBROUTINE CheckNodeConnectivity()
!===================================================================================================================================
! Eliminates multiple nodes, checks periodic boundary conditions and connects elements to their neighbours.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:tElem,tSide,tEdge,FirstElem
USE MOD_Mesh_Vars, ONLY:N
USE MOD_Basis_Vars,ONLY:QuadMap
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER       :: Elem                                                       ! Local element pointers
TYPE(tSide),POINTER       :: Side                                                      ! Local side pointers
TYPE(tEdge),POINTER       :: Edge                                                      ! Local side pointers
INTEGER                   :: iNode,i
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')' =======================CHECK NODE CONNECTIVITY ========================================='
WRITE(UNIT_stdOut,'(A)')'###### CHECK CORNER NODES'
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  DO iNode=1,Elem%nNodes
    Elem%node(iNode)%np%tmp=-888
  END DO !iNodes
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    DO iNode=1,Side%nNodes
      IF(Side%node(iNode)%np%tmp.NE.-888)THEN
        WRITE(UNIT_stdOut,'(A,I8)') 'Corner nodes der Seiten /= Corner nodes vom Element',Elem%ind
      END IF
    END DO 
    DO i=1,Side%nNodes
      Edge=>Side%edge(i)%edp
      IF(Edge%Node(1)%np%tmp.NE.-888) WRITE(UNIT_stdOut,'(A,I8)')'Edge node 1 /= Corner nodes vom Element',Elem%ind
      IF(Edge%Node(2)%np%tmp.NE.-888) WRITE(UNIT_stdOut,'(A,I8)')'Edge node 2 /= Corner nodes vom Element',Elem%ind
    END DO
    Side=>Side%nextElemSide
  END DO !associated(side)
  Elem=>Elem%nextElem
END DO

WRITE(UNIT_stdOut,'(A)')'###### CHECK EDGE CURVED NODES'
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    IF(Side%nCurvedNodes.NE.0) THEN
      DO iNode=1,Side%ncurvedNodes
        Side%curvednode(iNode)%np%tmp=-777
      END DO 
      DO iNode=1,Side%nNodes
        IF(Side%node(iNode)%np%tmp.NE.-777)&
          WRITE(UNIT_stdOut,'(A,I8,A,I8)') 'Side Corner node /= Curved surface corner nodes, iNode=',iNode,'of element',Elem%ind
      END DO 
      DO i=1,Side%nNodes
        Edge=>Side%edge(i)%edp
        IF(ASSOCIATED(Edge%CurvedNode))THEN
          DO iNode=1,N+1
            IF(Edge%CurvedNode(iNode)%np%tmp.NE.-777)&
              WRITE(UNIT_stdOut,'(A,I8,A,I8)')    'Curved edge node /= Curved surface nodes, iNode=',iNode,'of element',Elem%ind
          END DO
        END IF
      END DO
    END IF
    Side=>Side%nextElemSide
  END DO !associated(side)
  Elem=>Elem%nextElem
END DO

WRITE(UNIT_stdOut,'(A)')'###### CHECK SURFACE CURVED NODES'
Elem=>FirstElem
DO WHILE(ASSOCIATED(Elem))
  IF(Elem%nCurvedNodes.NE.0) THEN
    DO iNode=1,Elem%nCurvedNodes
      Elem%curvednode(iNode)%np%tmp=-999
    END DO !iNodes
    DO iNode=1,Elem%nNodes
      IF(Elem%node(iNode)%np%tmp.NE.-999)&
       WRITE(UNIT_stdOut,'(A,I8,A,I8)') 'Corner nodes of element /= Curved volume corner nodes, iNode=',iNode,'of element',Elem%ind
    END DO !iNodes
    Side=>Elem%firstSide
    DO WHILE(ASSOCIATED(Side))
      DO iNode=1,Side%ncurvedNodes
        IF(Side%curvednode(iNode)%np%tmp.NE.-999)&
         WRITE(UNIT_stdOut,'(A,2I8,A,I8)')'Curved surface node /= Curved volume nodes, p,q=',QuadMap(iNode,:),'of element',Elem%ind
      END DO 
      DO i=1,Side%nNodes
        Edge=>Side%edge(i)%edp
        IF(ASSOCIATED(Edge%CurvedNode))THEN
          DO iNode=1,N+1
            IF(Edge%CurvedNode(iNode)%np%tmp.NE.-999)&
              WRITE(UNIT_stdOut,'(A,I8,A,I8)') 'Curved edge node /= Curved volume nodes, iNode=',iNode,'of element',Elem%ind
          END DO
        END IF
      END DO
      Side=>Side%nextElemSide
    END DO !associated(side)
  END IF
  Elem=>Elem%nextElem
END DO
END SUBROUTINE CheckNodeConnectivity

SUBROUTINE ApplyMeshScale(StartElem)
!===================================================================================================================================
! Apply meshscale to existing mesh (called either directly after readin or before output)
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars, ONLY:tElem,tSide,tEdge
USE MOD_Mesh_Vars, ONLY:N,MeshScale
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tElem),POINTER,INTENT(IN) :: StartElem  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tElem),POINTER       :: Elem  ! ?
TYPE(tSide),POINTER       :: Side  ! ?
TYPE(tEdge),POINTER       :: Edge  ! ?
INTEGER                   :: iNode,iEdge  ! ?
!===================================================================================================================================
! make all nodes unique
Elem=>StartElem
DO WHILE(ASSOCIATED(Elem))
  ! elem nodes
  DO iNode=1,Elem%nNodes
    IF(Elem%node(iNode)%np%tmp.NE.-11)THEN
      Elem%node(iNode)%np%x=Elem%node(iNode)%np%x*MeshScale
      Elem%node(iNode)%np%tmp=-11
    END IF
  END DO
  ! curved elem nodes
  DO iNode=1,Elem%nCurvedNodes
    IF(Elem%curvedNode(iNode)%np%tmp.NE.-11)THEN
      Elem%curvedNode(iNode)%np%x=Elem%curvedNode(iNode)%np%x*MeshScale
      Elem%curvedNode(iNode)%np%tmp=-11
    END IF
  END DO
  ! sides
  Side=>Elem%firstSide
  DO WHILE(ASSOCIATED(Side))
    ! side nodes
    DO iNode=1,Side%nNodes
      IF(side%node(iNode)%np%tmp.NE.-11)THEN
        side%node(iNode)%np%x=side%node(iNode)%np%x*MeshScale
        side%node(iNode)%np%tmp=-11
      END IF
    END DO
    ! curved side nodes
    DO iNode=1,Side%nCurvedNodes
      IF(side%curvedNode(iNode)%np%tmp.NE.-11)THEN
        side%curvedNode(iNode)%np%x=side%curvedNode(iNode)%np%x*MeshScale
        side%curvedNode(iNode)%np%tmp=-11
      END IF
    END DO
    ! edges
    DO iEdge=1,Side%nNodes
      IF(.NOT.ASSOCIATED(Side%edge(iEdge)%edp)) CYCLE
      Edge=>Side%edge(iEdge)%edp
      ! edge nodes
      DO iNode=1,2
        IF(.NOT.ASSOCIATED(edge%node(iNode)%np)) CYCLE
        IF(edge%Node(iNode)%np%tmp.NE.-11)THEN
          edge%node(iNode)%np%x=edge%node(iNode)%np%x*MeshScale
          edge%node(iNode)%np%tmp=-11
        END IF
      END DO
      ! curved edge nodes
      IF(ASSOCIATED(edge%curvedNode))THEN
        DO iNode=1,N+1
          IF(edge%curvedNode(iNode)%np%tmp.NE.-11)THEN
            edge%curvedNode(iNode)%np%x=edge%curvedNode(iNode)%np%x*MeshScale
            edge%curvedNode(iNode)%np%tmp=-11
          END IF
        END DO
      END IF
    END DO
    Side=>Side%nextElemSide
  END DO
  Elem=>Elem%nextElem
END DO
END SUBROUTINE ApplyMeshScale

END MODULE MOD_Mesh
