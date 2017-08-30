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

INTERFACE GetNodesAndWeights
   MODULE PROCEDURE GetNodesAndWeights
END INTERFACE

INTERFACE GetVandermonde
   MODULE PROCEDURE GetVandermonde
END INTERFACE

PUBLIC:: InitBasis
PUBLIC:: GetNodesAndWeights
PUBLIC:: GetVandermonde
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
WRITE(tmpstr,'(I4)')3*(N+1)
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
USE MOD_Basis_Vars,ONLY:nVisu
USE MOD_Basis_Vars,ONLY:TriaMap,TriaMapInv,VisuTriaMap,VisuTriaMapInv,Vdm_visu_Tria,D_visu_Tria
USE MOD_Basis_Vars,ONLY:QuadMap,QuadMapInv,VisuQuadMap,VisuQuadMapInv,Vdm_visu_Quad,D_visu_Quad
USE MOD_Basis_Vars,ONLY:TetraMap,TetraMapInv !,Vdm_visu_Tetra,D_visu_Tetra
USE MOD_Basis_Vars,ONLY:PyraMap,PyraMapInv   !,Vdm_visu_Pyra,D_visu_Pyra
USE MOD_Basis_Vars,ONLY:PrismMap,PrismMapInv !,Vdm_visu_Prism,D_visu_Prism
USE MOD_Basis_Vars,ONLY:HexaMap,HexaMapInv,Vdm_visu_Hexa,D_visu_Hexa,VisuHexaMap,VisuHexaMapInv
USE MOD_Basis_Vars,ONLY:MapSideToVol
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


SUBROUTINE GetNodesAndWeights(N_in,NodeType_in,xIP,wIP,wIPBary)
!==================================================================================================================================
! Compute 1D nodes and weights and build the Vandermonde-Matrix
!==================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_Basis1D, ONLY: LegendreGaussNodesAndWeights,LegGaussLobNodesAndWeights,ChebyshevGaussNodesAndWeights
USE MOD_Basis1D, ONLY: ChebyGaussLobNodesAndWeights,BarycentricWeights
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: N_in         ! Number of 1D points
CHARACTER(LEN=*),INTENT(IN)        :: NodeType_in  ! Type of 1D points
REAL,INTENT(OUT)                   :: xIP(0:N_in)
REAL,INTENT(OUT),OPTIONAL          :: wIP(0:N_in)
REAL,INTENT(OUT),OPTIONAL          :: wIPBary(0:N_in)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i
!==================================================================================================================================
IF(PRESENT(wIP))THEN
  SELECT CASE(TRIM(NodeType_in))
  CASE('GAUSS')
    CALL LegendreGaussNodesAndWeights( N_in,xIP,wIP)
  CASE('GAUSS-LOBATTO')
    CALL LegGaussLobNodesAndWeights(   N_in,xIP,wIP)
  CASE('CHEBYSHEV-GAUSS')
    CALL ChebyshevGaussNodesAndWeights(N_in,xIP,wIP)
  CASE('CHEBYSHEV-GAUSS-LOBATTO')
    CALL ChebyGaussLobNodesAndWeights( N_in,xIP,wIP)
  CASE('VISU')
    DO i=0,N_in
      xIP(i) = 2.*REAL(i)/REAL(N_in) - 1.
    END DO
    ! Trapez rule for integration !!!
    wIP(:) = 2./REAL(N_in)
    wIP(0) = 0.5*wIP(0)
    wIP(N_in) = 0.5*wIP(N_in)
  CASE('VISU_INNER')
    DO i=0,N_in
      xIP(i) = 1./REAL(N_in+1)+2.*REAL(i)/REAL(N_in+1) - 1.
    END DO
    ! first order intergration !!!
    wIP=2./REAL(N_in+1)
  CASE DEFAULT
    CALL Abort(__STAMP__,&
      'NodeType "'//TRIM(NodeType_in)//'" in GetNodesAndWeights not found!')
  END SELECT
ELSE
  SELECT CASE(TRIM(NodeType_in))
  CASE('GAUSS')
    CALL LegendreGaussNodesAndWeights(N_in,xIP)
  CASE('GAUSS-LOBATTO')
    CALL LegGaussLobNodesAndWeights(N_in,xIP)
  CASE('CHEBYSHEV-GAUSS')
    CALL ChebyshevGaussNodesAndWeights(N_in,xIP)
  CASE('CHEBYSHEV-GAUSS-LOBATTO')
    CALL ChebyGaussLobNodesAndWeights(N_in,xIP)
  CASE('VISU')
    DO i=0,N_in
      xIP(i) = 2.*REAL(i)/REAL(N_in) - 1.
    END DO
  CASE('VISU_INNER')
    DO i=0,N_in
      xIP(i) = 1./REAL(N_in+1)+2.*REAL(i)/REAL(N_in+1) - 1.
    END DO
  CASE DEFAULT
    CALL Abort(__STAMP__,&
      'NodeType "'//TRIM(NodeType_in)//'" in GetNodesAndWeights not found!')
  END SELECT
END IF !present wIP
IF(PRESENT(wIPBary)) CALL BarycentricWeights(N_in,xIP,wIPBary)
END SUBROUTINE GetNodesAndWeights


SUBROUTINE GetVandermonde(N_in,NodeType_in,N_out,NodeType_out,Vdm_In_Out,Vdm_Out_In,modal)
!==================================================================================================================================
! Compute 1D nodes and weights and build the Vandermonde-Matrix
!==================================================================================================================================
! MODULES
USE MOD_Basis1D, ONLY:BarycentricWeights,InitializeVandermonde,BuildLegendreVdm
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                 :: N_in,N_out                 ! Number of 1D input points / output points
CHARACTER(LEN=*),INTENT(IN)        :: NodeType_in,NodeType_out   ! Type of 1D input points / output points
LOGICAL,INTENT(IN),OPTIONAL        :: modal                      !
REAL,INTENT(OUT)                   :: Vdm_In_out(0:N_out,0:N_in)
REAL,INTENT(OUT),OPTIONAL          :: Vdm_Out_In(0:N_in,0:N_out)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i
REAL                               :: xIP_in(0:N_in)
REAL                               :: xIP_out(0:N_out)
REAL                               :: wBary_in(0:N_in)
REAL                               ::  Vdm_Leg_in( 0:N_in,0:N_in)
REAL                               :: sVdm_Leg_in( 0:N_in,0:N_in)
REAL                               ::  Vdm_Leg_out(0:N_out,0:N_out)
REAL                               :: sVdm_Leg_out(0:N_out,0:N_out)
LOGICAL                            :: modalLoc
!==================================================================================================================================
modalLoc=.FALSE.
IF(PRESENT(modal)) modalLoc=modal

! Check if change Basis is needed
IF((TRIM(NodeType_out).EQ.TRIM(NodeType_in)).AND.(N_in.EQ.N_out))THEN
  Vdm_In_Out=0.
  DO i=0,N_in
    Vdm_In_out(i,i)=1.
  END DO
  IF(PRESENT(Vdm_Out_In))THEN
    Vdm_Out_In=0.
    DO i=0,N_Out
      Vdm_Out_In(i,i)=1.
    END DO
  END IF
ELSE
  ! Input points
  CALL GetNodesAndWeights(N_in,NodeType_in,xIP_in)
  CALL BarycentricWeights(N_in,xIP_in,wBary_in)
  ! Output points
  CALL GetNodesAndWeights(N_out,NodeType_out,xIP_out)

  IF(modalLoc)THEN
    CALL BuildLegendreVdm(N_In, xIP_in, Vdm_Leg_in, sVdm_Leg_in)
    CALL BuildLegendreVdm(N_Out,xIP_out,Vdm_Leg_out,sVdm_Leg_out)
  END IF

  IF((N_Out.LT.N_In).AND.modalLoc)THEN
    Vdm_In_Out=MATMUL(Vdm_Leg_Out(0:N_Out,0:N_Out),sVdm_Leg_In(0:N_Out,0:N_In))
  ELSE
    CALL InitializeVandermonde(N_in,N_out,wBary_in,xIP_in,xIP_out,Vdm_In_Out)
  END IF
  IF(PRESENT(Vdm_Out_In))THEN
    IF((N_In.LT.N_Out).AND.modalLoc)THEN
      Vdm_Out_In=MATMUL(Vdm_Leg_In(0:N_In,0:N_In),sVdm_Leg_Out(0:N_In,0:N_Out))
    ELSE
      CALL BarycentricWeights(N_out,xIP_out,wBary_in)
      CALL InitializeVandermonde(N_out,N_in,wBary_in,xIP_out,xIP_in,Vdm_Out_In)
    END IF
  END IF
END IF
END SUBROUTINE GetVandermonde


END MODULE MOD_Basis
