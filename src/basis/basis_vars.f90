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
MODULE MOD_Basis_Vars
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

! Vandermande and D matrices
REAL,ALLOCATABLE,TARGET        :: VdM_visu_Tria(:,:)          
REAL,ALLOCATABLE,TARGET        :: D_visu_Tria(:,:,:)          
REAL,ALLOCATABLE,TARGET        :: VdM_visu_Quad(:,:)          
REAL,ALLOCATABLE,TARGET        :: D_visu_Quad(:,:,:)          
REAL,ALLOCATABLE,TARGET        :: VdM_visu_Tetra(:,:)          
REAL,ALLOCATABLE,TARGET        :: D_visu_Tetra(:,:,:)          
REAL,ALLOCATABLE,TARGET        :: VdM_visu_Pyra(:,:)          
REAL,ALLOCATABLE,TARGET        :: D_visu_Pyra(:,:,:)          
REAL,ALLOCATABLE,TARGET        :: VdM_visu_Prism(:,:)          
REAL,ALLOCATABLE,TARGET        :: D_visu_Prism(:,:,:)          
REAL,ALLOCATABLE,TARGET        :: VdM_visu_Hexa(:,:)          
REAL,ALLOCATABLE,TARGET        :: D_visu_Hexa(:,:,:)          

! Tensorproduct mappings + inverse mappings for all elements
INTEGER,ALLOCATABLE,TARGET     :: TriaMap(:,:) 
INTEGER,ALLOCATABLE,TARGET     :: TriaMapInv(:,:)
INTEGER,ALLOCATABLE,TARGET     :: QuadMap(:,:) 
INTEGER,ALLOCATABLE,TARGET     :: QuadMapInv(:,:)
INTEGER,ALLOCATABLE,TARGET     :: TetraMap(:,:) 
INTEGER,ALLOCATABLE,TARGET     :: TetraMapInv(:,:,:)
INTEGER,ALLOCATABLE,TARGET     :: PyraMap(:,:) 
INTEGER,ALLOCATABLE,TARGET     :: PyraMapInv(:,:,:)
INTEGER,ALLOCATABLE,TARGET     :: PrismMap(:,:) 
INTEGER,ALLOCATABLE,TARGET     :: PrismMapInv(:,:,:)
INTEGER,ALLOCATABLE,TARGET     :: HexaMap(:,:) 
INTEGER,ALLOCATABLE,TARGET     :: HexaMapInv(:,:,:)

INTEGER                        :: nVisu                  ! number of 1D Points -1  for visualization
INTEGER                        :: nAnalyze               ! number of 1D Points -1  for analysis (Jacobian...)
INTEGER,ALLOCATABLE            :: VisuTriaMap(:,:) 
INTEGER,ALLOCATABLE            :: VisuTriaMapInv(:,:)
INTEGER,ALLOCATABLE            :: VisuQuadMap(:,:) 
INTEGER,ALLOCATABLE            :: VisuQuadMapInv(:,:)
INTEGER,ALLOCATABLE            :: VisuHexaMap(:,:) 
INTEGER,ALLOCATABLE            :: VisuHexaMapInv(:,:,:)

INTEGER,ALLOCATABLE            :: edgeToTria(:,:)        ! mapping from edges of a triangle to surface
INTEGER,ALLOCATABLE            :: edgeToQuad(:,:)        ! mapping from edges of a quadrangle to surface
INTEGER,ALLOCATABLE            :: MapSideToVol(:,:,:)    ! iVolNode=MapSideToVol(iSideNode,Side%locSide,Elem%nNodes)
!===================================================================================================================================
END MODULE MOD_Basis_Vars
