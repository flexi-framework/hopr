MODULE MOD_Output_Vars
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
LOGICAL                     :: DebugVisu                  ! set .TRUE. for debug visualization (INPUT)
INTEGER                     :: DebugVisuLevel             !=0, only linear mesh, =1 + surfspline (default), =2  +volspline
REAL                        :: Visu_sJ_limit              ! limit to visualize only curved elements with sJ<=Visu_sJ_limit
INTEGER                     :: outputFormat               !=0: VTK, =1 tecplot ascii, =2 CGNS 
CHARACTER(LEN=100)          :: sfc_type                   ! morton or hilbert
LOGICAL                     :: doSortIJK
LOGICAL                     :: useSpaceFillingCurve
LOGICAL                     :: OutputInitDone
!===================================================================================================================================
END MODULE MOD_Output_Vars
