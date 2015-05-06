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
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
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
REAL,ALLOCATABLE,TARGET        :: VdM_analyze_Hexa(:,:)          
REAL,ALLOCATABLE,TARGET        :: D_analyze_Hexa(:,:,:)          

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
!===================================================================================================================================
END MODULE MOD_Basis_Vars
