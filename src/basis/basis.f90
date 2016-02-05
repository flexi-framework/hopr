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
#include "defines.f90"
MODULE MOD_Basis
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
INTERFACE fillBasisMapping
   MODULE PROCEDURE fillBasisMapping
END INTERFACE

INTERFACE InitBasis
   MODULE PROCEDURE InitBasis
END INTERFACE

INTERFACE ChangeBasisHexa 
  MODULE PROCEDURE ChangeBasisHexa
END INTERFACE

PUBLIC:: InitBasis
PUBLIC:: ChangeBasisHexa
!===================================================================================================================================

CONTAINS
SUBROUTINE InitBasis()
!===================================================================================================================================
! Initialize...
!===================================================================================================================================
! MODULES
USE MOD_Basis_Vars
USE MOD_Mesh_Vars,ONLY:MeshInitDone,N
USE MOD_ReadInTools
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
CHARACTER(LEN=4) :: tmpstr  ! ?
!===================================================================================================================================
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A)') ' INIT BASIS...'
IF(.NOT.MeshInitDone)THEN
  CALL abort(__STAMP__,  &
             'ERROR: InitMesh has to be called before InitBasis!',999,999.)
END IF
WRITE(tmpstr,'(I4)')N
nVisu=GETINT('nVisu',tmpstr)
WRITE(tmpstr,'(I4)')N+2
nAnalyze=GETINT('nAnalyze',tmpstr)
IF(nVisu.LT.1)THEN
  CALL abort(__STAMP__,'nVisu has to be >= 1')
END IF
CALL fillBasisMapping()
WRITE(UNIT_stdOut,'(A)')' INIT BASIS DONE!'
WRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitBasis


SUBROUTINE fillBasisMapping()
!===================================================================================================================================
! determines the hierarchical relation from iAns to the position of the testfunction in the basis triangle, tetrahedra 
! example: linear monomial basis function in 2D: nAns=3: 1,x,y
! iAns=1: 1, iAns=2: x, iAns=3: y
!===================================================================================================================================
! MODULES
USE MOD_Basis_Vars,ONLY:nVisu,nAnalyze
USE MOD_Basis_Vars,ONLY:TriaMap,TriaMapInv,VisuTriaMap,VisuTriaMapInv,Vdm_visu_Tria,D_visu_Tria
USE MOD_Basis_Vars,ONLY:QuadMap,QuadMapInv,VisuQuadMap,VisuQuadMapInv,Vdm_visu_Quad,D_visu_Quad
USE MOD_Basis_Vars,ONLY:TetraMap,TetraMapInv !,Vdm_visu_Tetra,D_visu_Tetra
USE MOD_Basis_Vars,ONLY:PyraMap,PyraMapInv   !,Vdm_visu_Pyra,D_visu_Pyra
USE MOD_Basis_Vars,ONLY:PrismMap,PrismMapInv !,Vdm_visu_Prism,D_visu_Prism
USE MOD_Basis_Vars,ONLY:HexaMap,HexaMapInv,Vdm_visu_Hexa,D_visu_Hexa,VisuHexaMap,VisuHexaMapInv
USE MOD_Basis_Vars,ONLY:MapSideToVol
USE MOD_Basis_Vars,ONLY:Vdm_Analyze_Hexa,D_Analyze_Hexa
USE MOD_Basis_Vars,ONLY:EdgeToTria,EdgeToQuad
USE MOD_Mesh_Vars,ONLY:N
USE MOD_QuadBasis
USE MOD_TriaBasis
USE MOD_TetraBasis
USE MOD_PyraBasis
USE MOD_PrismBasis
USE MOD_HexaBasis
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: i,j,p,q,nNodes   ! ?
!===================================================================================================================================

! Build mappings for each element type
CALL getBasisMappingTria(N,nNodes,TriaMap,TriaMapInv)
CALL getBasisMappingQuad(N,nNodes,QuadMap,QuadMapInv)
CALL getBasisMappingTetra(N,nNodes,TetraMap,TetraMapInv)
CALL getBasisMappingPyra(N,nNodes,PyraMap,PyraMapInv)
CALL getBasisMappingPrism(N,nNodes,PrismMap,PrismMapInv)

! basis mappings used for visualization (often number of visu points != order) + mapping tria->quad(tensorproduct) required
CALL getBasisMappingQuad(nVisu,nNodes,VisuTriaMap,VisuTriaMapInv)
CALL getBasisMappingQuad(nVisu,nNodes,VisuQuadMap,VisuQuadMapInv)
CALL getBasisMappingHexa(N,nNodes,HexaMap,HexaMapInv)
CALL getBasisMappingHexa(nVisu,nNodes,VisuHexaMap,VisuHexaMapInv)

! get Vandermonde and D matrices for visualization
CALL getTriaBasis(N,nVisu+1,Vdm_visu_Tria,D_visu_Tria)
CALL getQuadBasis(N,nVisu+1,Vdm_visu_Quad,D_visu_Quad)
!CALL getTetraBasis(N,nVisu+1,Vdm_visu_Tetra,D_visu_Tetra)
CALL getHexaBasis(N,nVisu+1,Vdm_visu_Hexa,D_visu_Hexa)
CALL getHexaBasis(N,nAnalyze,Vdm_Analyze_Hexa,D_Analyze_Hexa)

! map edges of triangle surface counterclockwise points to BoundaBasisMappingInv(i,j)
ALLOCATE(EdgeToTria(3,N+1))
!first Side
j=0
DO i=0,N
  edgeToTria(1,i+1)=TriaMapInv(i,j)
END DO
!Second Side
DO i=0,N
  edgeToTria(2,i+1)=TriaMapInv(N-i,i)
END DO
!last Side
i=0
DO j=0,N
  edgeToTria(3,j+1)=TriaMapInv(i,N-j)
END DO
! map edges of quadrangle counterclockwise points to BoundaBasisMappingInv(i,j)
ALLOCATE(EdgeToQuad(4,N+1))
!first Side
j=0
DO i=0,N
  EdgeToQuad(1,i+1)=QuadMapInv(i,j)
END DO
!Second Side
i=N
DO j=0,N
  EdgeToQuad(2,j+1)=QuadMapInv(i,j)
END DO
!Third Side
j=N
DO i=0,N
  EdgeToQuad(3,i+1)=QuadMapInv(N-i,j)
END DO
!last Side
i=0
DO j=0,N
  EdgeToQuad(4,j+1)=QuadMapInv(i,N-j) 
END DO

!side mappings
ALLOCATE(MapSideToVol((N+1)**2,8,4:8)) !sidenodes,locSide,Elem%nNodes
MapSideToVol=-1 !not associated = -1!!
DO q=0,N
  DO p=0,N-q
      MapSideToVol(TriaMapInv(p,q),1,4)=TetraMapInv(q,p,0)
      MapSideToVol(TriaMapInv(p,q),2,4)=TetraMapInv(p,0,q)
      MapSideToVol(TriaMapInv(p,q),3,4)=TetraMapInv(N-q-p,p,q)
      MapSideToVol(TriaMapInv(p,q),4,4)=TetraMapInv(0,q,p)
  END DO !p
END DO !q
!pyra  
DO q=0,N
  DO p=0,N
    MapSideToVol(QuadMapInv(p,q),1,5)=PyraMapInv(q,p,0)
  END DO !p
END DO !q
DO q=0,N
  DO p=0,N-q
    MapSideToVol(TriaMapInv(p,q),2,5)=PyraMapInv(p,0,q)
    MapSideToVol(TriaMapInv(p,q),3,5)=PyraMapInv(N-q,p,q)
    MapSideToVol(TriaMapInv(p,q),4,5)=PyraMapInv(N-q-p,N-q,q)
    MapSideToVol(TriaMapInv(p,q),5,5)=PyraMapInv(0,q,p)
  END DO !p
END DO !q
!prism
DO q=0,N
  DO p=0,N
      MapSideToVol(QuadMapInv(p,q),1,6)=PrismMapInv(p,0,q)
      MapSideToVol(QuadMapInv(p,q),2,6)=PrismMapInv(N-p,p,q)
      MapSideToVol(QuadMapInv(p,q),3,6)=PrismMapInv(0,q,p)
  END DO !p
END DO !q
DO q=0,N
  DO p=0,N-q
      MapSideToVol(triaMapInv(p,q),4,6)=PyraMapInv(p,q,N)
      MapSideToVol(triaMapInv(p,q),5,6)=PrismMapInv(q,p,0)
  END DO !p
END DO !q
!hexa
DO q=0,N
  DO p=0,N
    MapSideToVol(QuadMapInv(p,q),1,8)=HexaMapInv(q,p,0)
    MapSideToVol(QuadMapInv(p,q),2,8)=HexaMapInv(p,0,q)
    MapSideToVol(QuadMapInv(p,q),3,8)=HexaMapInv(N,p,q)
    MapSideToVol(QuadMapInv(p,q),4,8)=HexaMapInv(N-p,N,q)
    MapSideToVol(QuadMapInv(p,q),5,8)=HexaMapInv(0,q,p)
    MapSideToVol(QuadMapInv(p,q),6,8)=HexaMapInv(p,q,N)
  END DO !p
END DO !q
END SUBROUTINE fillBasisMapping


SUBROUTINE ChangeBasisHexa(Dim1,N_In,N_Out,Vdm,X3D_In,X3D_Out)
!===================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
! to another 3D tensor product node positions (number of nodes N_out+1) 
! defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!  xi is defined in the 1DrefElem xi=[-1,1]
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: Dim1,N_In,N_Out
REAL,INTENT(IN)     :: X3D_In(1:Dim1,0:N_In,0:N_In,0:N_In)
REAL,INTENT(IN)     :: Vdm(0:N_Out,0:N_In)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)    :: X3D_Out(1:Dim1,0:N_Out,0:N_Out,0:N_Out)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
INTEGER             :: iN_In,jN_In,kN_In,iN_Out,jN_Out,kN_Out
REAL                :: X3D_Buf1(1:Dim1,0:N_Out,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:Dim1,0:N_Out,0:N_Out,0:N_In) ! second intermediate results from 1D interpolations
!===================================================================================================================================
X3D_buf1=0.
! first direction iN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO iN_In=0,N_In
      DO iN_Out=0,N_Out
        X3D_Buf1(:,iN_Out,jN_In,kN_In)=X3D_Buf1(:,iN_Out,jN_In,kN_In)+Vdm(iN_Out,iN_In)*X3D_In(:,iN_In,jN_In,kN_In)
      END DO
    END DO
  END DO
END DO
X3D_buf2=0.
! second direction jN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO jN_Out=0,N_Out
      DO iN_Out=0,N_Out
        X3D_Buf2(:,iN_Out,jN_Out,kN_In)=X3D_Buf2(:,iN_Out,jN_Out,kN_In)+Vdm(jN_Out,jN_In)*X3D_Buf1(:,iN_Out,jN_In,kN_In)
      END DO
    END DO
  END DO
END DO
X3D_Out=0.
! last direction kN_In
DO kN_In=0,N_In
  DO kN_Out=0,N_Out
    DO jN_Out=0,N_Out
      DO iN_Out=0,N_Out
        X3D_Out(:,iN_Out,jN_Out,kN_Out)=X3D_Out(:,iN_Out,jN_Out,kN_Out)+Vdm(kN_Out,kN_In)*X3D_Buf2(:,iN_Out,jN_Out,kN_In)
      END DO
    END DO
  END DO
END DO
END SUBROUTINE ChangeBasisHexa

END MODULE MOD_Basis
