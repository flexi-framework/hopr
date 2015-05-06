MODULE MOD_SpaceFillingCurve
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Output_Vars,ONLY:sfc_type
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
TYPE tBox
  REAL(KIND=8) :: mini(3)
  INTEGER      :: nbits
  REAL(KIND=8) :: spacing(3)
END TYPE tBox
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE SortElemsBySpaceFillingCurve
  MODULE PROCEDURE SortElemsBySpaceFillingCurve
END INTERFACE

INTERFACE
   FUNCTION evalhilbert(disc,nbits,ndims) Result(indx)
       INTEGER(KIND=8),INTENT(IN)  :: disc(3) ! ?
       INTEGER(KIND=8)  :: indx ! ?
       INTEGER(KIND=4),INTENT(IN)  :: nbits ! ?
       INTEGER(KIND=4),INTENT(IN)  :: nDims ! ?
   END FUNCTION evalhilbert
END INTERFACE

PUBLIC::SortElemsBySpaceFillingCurve
PUBLIC::EVAL_MORTON,EVAL_MORTON_ARR
!===================================================================================================================================

CONTAINS
SUBROUTINE SortElemsBySpaceFillingCurve(nElems,ElemBary,IDList,whichBoundBox)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_globals,ONLY:UNIT_stdOut
USE MOD_SortingTools,ONLY:Qsort1DoubleInt1Pint
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)                         :: nElems  ! ?
REAL,INTENT(IN)                            :: ElemBary(nElems,3)  ! ?
INTEGER,INTENT(IN)                         :: whichBoundBox  !=1: each direction independant, =2: cube
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(INOUT)                       :: IDlist(nElems)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tBox) :: SCBox  ! ?
REAL       :: lower(3)  ! ?
REAL       :: upper(3)  ! ?
INTEGER    :: iElem  ! ?
INTEGER(KIND=8) :: IntList(nElems)  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,A,A)')'SORT ELEMENTS ON SPACE FILLING CURVE, TYPE ',TRIM(sfc_type),' ...'
IF(nElems.GT.1)THEN
  ! Determine extreme verticies for bounding box
  SELECT CASE(whichBoundBox)
  CASE(1)
    lower(1) = minval(ElemBary(:,1))
    lower(2) = minval(ElemBary(:,2))
    lower(3) = minval(ElemBary(:,3))
    upper(1) = maxval(ElemBary(:,1))
    upper(2) = maxval(ElemBary(:,2))
    upper(3) = maxval(ElemBary(:,3))
  CASE(2)
  lower=MINVAL(ElemBary)
  upper=MAXVAL(ElemBary)
  END SELECT
  CALL setBoundingBox(SCBox,lower,upper)
  
  DO iElem=1,nElems
    IntList(iElem) = COORD2INT(SCBox, ElemBary(iElem,:))
  END DO
  
  ! Now sort the elements according to their index
  ! on the space filling curve.
  CALL Qsort1DoubleInt1Pint(IntList, IDList)
ELSE
  IDList=1
END IF  
WRITE(UNIT_stdOut,'(A)')'... DONE'
END SUBROUTINE SortElemsBySpaceFillingCurve

FUNCTION COORD2INT(Box, Coord) RESULT(ind)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
TYPE(tBox),INTENT(IN)       :: Box  ! ?
REAL(KIND=8),INTENT(IN)     :: Coord(3)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=8)             :: ind  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8)             :: disc(3)  ! ?
!===================================================================================================================================
! compute the integer discretization in each
! direction
disc = NINT((coord-box%mini)*box%spacing)

! Map the three coordinates on a single integer
! value.
SELECT CASE(sfc_type)
CASE('morton')
  ind = EVAL_MORTON(disc,box%nBits)
CASE('hilbert')
  ind = evalhilbert(disc(1:3),box%nbits,3)
END SELECT
END FUNCTION COORD2INT 

FUNCTION EVAL_MORTON(intcoords,nBits) RESULT(ind)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)  :: intcoords(3)  ! ?
INTEGER,INTENT(IN)          :: nBits  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=8)             :: ind  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: dir,i
!===================================================================================================================================
  ind = 0
  DO dir=1,3
    DO i=0,nbits-1
      ! Interleave the three directions, start
      ! from position 0 with counting the bits.
      IF(BTEST(intcoords(dir),i))  ind = IBSET(ind,3*i-dir+3)
    END DO
  END DO
END FUNCTION EVAL_MORTON

SUBROUTINE EVAL_MORTON_ARR(ind,intcoords,nP,nBits)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER(KIND=8),INTENT(IN)  :: intcoords(3,nP)  ! ?
INTEGER,INTENT(IN)          :: nBits,nP  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(OUT) :: ind(nP)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER         :: dir,i
!===================================================================================================================================
ind = 0
DO dir=1,3
  DO i=0,nbits-1
    ! Interleave the three directions, start
    ! from position 0 with counting the bits.
    ind = MERGE(IBSET(ind,3*i-dir+3),ind,BTEST(intcoords(dir,:),i))
  END DO
END DO
END SUBROUTINE EVAL_MORTON_ARR

SUBROUTINE setBoundingBox(box,mini,maxi)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL(KIND=8),INTENT(IN)      :: mini(3)  ! ?
REAL(KIND=8),INTENT(IN)      :: maxi(3)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
TYPE(tBox),INTENT(OUT)        :: box  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(KIND=8)      :: blen(3)  ! ?
INTEGER(KIND=8)   :: intfact  ! ?
!===================================================================================================================================
box%mini = mini
blen = maxi - mini
box%nbits = (bit_size(intfact)-1) / 3
intfact = 2**box%nbits-1
box%spacing = REAL(intfact)/blen
END SUBROUTINE setBoundingBox

END MODULE MOD_SpaceFillingCurve
