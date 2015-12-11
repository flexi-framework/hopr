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
!*****************************************************************************!
!*
!*      $RCSfile: tools.f90,v $
!*
!*      $Revision: 1.1 $
!*
!*      $Date: 2007/04/10 15:37:06 $
!*
!*      $Author: iaghinde $
!*
!*****************************************************************************!

MODULE MOD_SortingTools
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------

  INTERFACE QSort1int
     MODULE PROCEDURE QSort1int
  END INTERFACE

  INTERFACE QSort1int1Pint
     MODULE PROCEDURE QSort1int1Pint
  END INTERFACE

  INTERFACE QSort1Doubleint1Pint
     MODULE PROCEDURE QSort1DoubleInt1Pint
  END INTERFACE

  INTERFACE QSort1int2double
     MODULE PROCEDURE QSort1int2double
  END INTERFACE

  INTERFACE QSort2int
     MODULE PROCEDURE QSort2int
  END INTERFACE

  INTERFACE QSort3int
     MODULE PROCEDURE QSort3int
  END INTERFACE

  INTERFACE QSort4int
     MODULE PROCEDURE QSort4int
  END INTERFACE

  INTERFACE QSort2Real
     MODULE PROCEDURE QSort2Real
  END INTERFACE

  INTERFACE QSort3Real
     MODULE PROCEDURE QSort3Real
  END INTERFACE

  INTERFACE MSortNInt
     MODULE PROCEDURE MSortNInt
  END INTERFACE

  INTERFACE MSortNLong
     MODULE PROCEDURE MSortNLong
  END INTERFACE

  INTERFACE MergeSortInt
     MODULE PROCEDURE MergeSortInt
  END INTERFACE

  INTERFACE MergeSortLong
     MODULE PROCEDURE MergeSortLong
  END INTERFACE

  !----------------------------------------------------------------------------
  PUBLIC :: QSort1int
  PUBLIC :: QSort1int1Pint
  PUBLIC :: QSort1Doubleint1Pint
  PUBLIC :: QSort1int2double
  PUBLIC :: QSort2int
  PUBLIC :: QSort3int
  PUBLIC :: QSort4int
  PUBLIC :: QSort2Real
  PUBLIC :: QSort3Real
  PUBLIC :: MSortNInt
  PUBLIC :: MSortNLong
  PUBLIC :: MergeSortInt
  PUBLIC :: MergeSortLong
  !----------------------------------------------------------------------------
  PRIVATE :: Partition1int
  PRIVATE :: Partition1int1Pint
  PRIVATE :: Partition1DoubleInt1Pint
  PRIVATE :: Partition1int2double
  PRIVATE :: Partition2int
  PRIVATE :: Partition3int
  PRIVATE :: Partition4int
  PRIVATE :: Partition2Real
  PRIVATE :: Partition3Real
  !----------------------------------------------------------------------------

CONTAINS
RECURSIVE SUBROUTINE Qsort1Int(A)
!===================================================================================================================================
! QSort1int: (c) by Mark Haas
!  Uses the Quicksort algorithm to sort an INTEGER table with one relevant columns.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
   INTEGER,INTENT(INOUT)            :: A(:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: marker  ! ?
!===================================================================================================================================

   IF(SIZE(A) .GT. 1) THEN
      CALL Partition1Int(A,marker)
      CALL Qsort1Int(A(:marker-1))
      CALL Qsort1Int(A(marker:))
   END IF
 RETURN
END SUBROUTINE Qsort1Int
SUBROUTINE Partition1Int(A,marker)
!===================================================================================================================================
!  Sorting routine used by QSort1int above. This routine is PRIVATE 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
   INTEGER,INTENT(INOUT)            :: A(:)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
   INTEGER,INTENT(OUT)              :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: i,j  ! ?
   INTEGER                          :: temp,x  ! ?
!===================================================================================================================================
   x = A(1)
   i = 0
   j  = SIZE(A)+1

   DO
      j = j-1
      DO
         IF (A(j) .LE. x)THEN
            EXIT
         END IF
         j = j-1
      END DO
      i = i+1
      DO
         IF(A(i) .GE. x)THEN
            EXIT
         END IF
         i = i+1
      END DO
      IF (i .LT. j) THEN
         ! exchange A(i) and A(j)
         temp = A(i)
         A(i)  = A(j)
         A(j)  = temp
      ELSE IF (i .EQ. j) THEN
         marker = i+1
         RETURN
      ELSE
         marker = i
         RETURN
      END IF
   END DO
 RETURN
END SUBROUTINE Partition1Int


RECURSIVE SUBROUTINE Qsort1Int1PInt(A,P)
!===================================================================================================================================
! QSort1int: (c) by Mark Haas
!  Uses the Quicksort algorithm to sort an INTEGER table with one relevant columns, along with a passive onedimensional array.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
   INTEGER,INTENT(INOUT)            :: A(:)  ! active array, will be sorted
   INTEGER,INTENT(INOUT)            :: P(:)  ! passive array, is sorted like A
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: marker  ! ?
!===================================================================================================================================

   IF(SIZE(A) .GT. 1) THEN
      CALL Partition1Int1Pint(A,P,marker)
      CALL Qsort1Int1Pint(A(:marker-1),P(:marker-1))
      CALL Qsort1Int1Pint(A(marker:),P(marker:))
   END IF
 RETURN
END SUBROUTINE Qsort1Int1PInt


SUBROUTINE Partition1Int1PInt(A,P,marker)
!===================================================================================================================================
!  Sorting routine used by QSort1int above. This routine is PRIVATE 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
   INTEGER,INTENT(INOUT)            :: A(:) ! active array, will be sorted
   INTEGER,INTENT(INOUT)            :: P(:) ! passive array, is sorted like A
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
   INTEGER,INTENT(OUT)              :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: i,j  ! ?
   INTEGER                          :: temp,x  ! ?
!===================================================================================================================================
   x = A(1)
   i = 0
   j  = SIZE(A)+1

   DO
      j = j-1
      DO
         IF (A(j) .LE. x)THEN
            EXIT
         END IF
         j = j-1
      END DO
      i = i+1
      DO
         IF(A(i) .GE. x)THEN
            EXIT
         END IF
         i = i+1
      END DO
      IF (i .LT. j) THEN
         ! exchange A(i) and A(j), P(i) and P(j)
         temp = A(i)
         A(i)  = A(j)
         A(j)  = temp
         temp = P(i)
         P(i)  = P(j)
         P(j)  = temp
      ELSE IF (i .EQ. j) THEN
         marker = i+1
         RETURN
      ELSE
         marker = i
         RETURN
      END IF
   END DO
 RETURN
END SUBROUTINE Partition1Int1PInt


RECURSIVE SUBROUTINE Qsort1DoubleInt1PInt(A,P)
!===================================================================================================================================
! QSort1int: (c) by Mark Haas
!  Uses the Quicksort algorithm to sort an INTEGER table with one relevant columns, along with a passive onedimensional array.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
   INTEGER(KIND=8),INTENT(INOUT)    :: A(:) ! active array, will be sorted
   INTEGER,INTENT(INOUT)            :: P(:) ! passive array, is sorted like A
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: marker  ! ?
!===================================================================================================================================

   IF(SIZE(A) .GT. 1) THEN
      CALL Partition1DoubleInt1Pint(A,P,marker)
      CALL Qsort1DoubleInt1Pint(A(:marker-1),P(:marker-1))
      CALL Qsort1DoubleInt1Pint(A(marker:),P(marker:))
   END IF
 RETURN
END SUBROUTINE Qsort1DoubleInt1PInt


SUBROUTINE Partition1DoubleInt1PInt(A,P,marker)
!===================================================================================================================================
!  Sorting routine used by QSort1int above. This routine is PRIVATE 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
   INTEGER(KIND=8),INTENT(INOUT)    :: A(:) ! active array, will be sorted
   INTEGER,INTENT(INOUT)            :: P(:) ! passive array, is sorted like A
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
   INTEGER,INTENT(OUT)              :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: i,j  ! ?
   INTEGER(KIND=8)                  :: temp,x  ! ?
   INTEGER                          :: ptemp  ! ?
!===================================================================================================================================
   x = A(1)
   i = 0
   j  = SIZE(A)+1

   DO
      j = j-1
      DO
         IF (A(j) .LE. x)THEN
            EXIT
         END IF
         j = j-1
      END DO
      i = i+1
      DO
         IF(A(i) .GE. x)THEN
            EXIT
         END IF
         i = i+1
      END DO
      IF (i .LT. j) THEN
         ! exchange A(i) and A(j), P(i) and P(j)
         temp = A(i)
         A(i)  = A(j)
         A(j)  = temp
         ptemp = P(i)
         P(i)  = P(j)
         P(j)  = ptemp
      ELSE IF (i .EQ. j) THEN
         marker = i+1
         RETURN
      ELSE
         marker = i
         RETURN
      END IF
   END DO
 RETURN
END SUBROUTINE Partition1DoubleInt1PInt


RECURSIVE SUBROUTINE Qsort1Int2double(A)
!===================================================================================================================================
! QSort1int: (c) by Mark Haas
!  Uses the Quicksort algorithm to sort an INTEGER(KIND=8) table with one relevant columns.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
   INTEGER(KIND=8),INTENT(INOUT),DIMENSION(:,:)         :: A  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: marker  ! ?
!===================================================================================================================================

   IF(SIZE(A,1) .GT. 1) THEN
      CALL Partition1Int2double(A,marker)
      CALL Qsort1Int2double(A(:marker-1,:))
      CALL Qsort1Int2double(A(marker:,:))
   END IF
 RETURN
END SUBROUTINE Qsort1Int2double


SUBROUTINE Partition1Int2double(A,marker)
!===================================================================================================================================
!  Sorting routine used by QSort1int above. This routine is PRIVATE
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
   INTEGER(KIND=8),INTENT(INOUT),DIMENSION(:,:)         :: A  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
   INTEGER,INTENT(OUT)              :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
   INTEGER                          :: i,j  ! ?
   INTEGER(KIND=8)                  :: temp(size(A,2))  ! ?
   INTEGER(KIND=8)                  :: x      ! pivot point
!===================================================================================================================================
   x = A(1,1)
   i = 0
   j  = SIZE(A,1)+1

   DO
      j = j-1
      DO
         IF (A(j,1) .LE. x)THEN
            EXIT
         END IF
         j = j-1
      END DO
      i = i+1
      DO
         IF(A(i,1) .GE. x)THEN
            EXIT
         END IF
         i = i+1
      END DO
      IF (i .LT. j) THEN
         ! exchange A(i) and A(j)
         temp(:) = A(i,:)
         A(i,:)  = A(j,:)
         A(j,:)  = temp(:)
      ELSE IF (i .EQ. j) THEN
         marker = i+1
         RETURN
      ELSE
         marker = i
         RETURN
      END IF
   END DO
 RETURN
END SUBROUTINE Partition1Int2double
RECURSIVE SUBROUTINE Qsort2int(A)
!===================================================================================================================================
!  Uses the Quicksort algorithm to sort an INTEGER table with two relevant columns. There may be an arbitrary number of additional
!  columns that will be moved around at the same time but do not influence the sorting process (except for slowing things down of
!  course!)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
INTEGER, INTENT(INOUT), DIMENSION(:,:) :: A ! first dimension is sorted, second dimension are the columns
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iq  ! ?
!===================================================================================================================================
IF(size(A,1) .GT. 1) THEN
   CALL Partition2int(A, iq)
   CALL Qsort2int(A(:iq-1,:))
   CALL Qsort2int(A(iq:,:))
END IF
END SUBROUTINE Qsort2int


SUBROUTINE Partition2int(A, marker)
!===================================================================================================================================
!  Sorting routine used by QSort2int above. This routine is PRIVATE
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
INTEGER, INTENT(INOUT), DIMENSION(:,:) :: A  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)                   :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER   :: i, j  ! ?
INTEGER   :: temp(size(A,2))  ! ?
INTEGER   :: x(2)              ! pivot point
!===================================================================================================================================
x(:) = A(1,1:2)
i= 0
j= size(A,1) + 1
DO
   j = j-1
   DO
      IF ((A(j,1) .LT. x(1)).OR.(((A(j,1).EQ.x(1)).AND.(A(j,2) .LE. x(2))))) EXIT
      j = j-1
   END DO
   i = i+1
   DO
      if ((A(i,1) .GT. x(1)).OR.(((A(i,1).EQ.x(1)).AND.(A(i,2) .GE. x(2))))) EXIT
      i = i+1
   END DO
   IF (i .LT. j) THEN
      ! exchange A(i) and A(j)
      temp(:) = A(i,:)
      A(i,:) = A(j,:)
      A(j,:) = temp(:)
   ELSE IF (i.EQ.j) THEN
      marker = i+1
      RETURN
   ELSE
      marker = i
      RETURN
   END IF
END DO
END SUBROUTINE Partition2int

RECURSIVE SUBROUTINE Qsort3int(A)
!===================================================================================================================================
! QSort3int:
!  Uses the Quicksort algorithm to sort an INTEGER table with three relevant columns.
!  There may be an arbitrary number of additional columns that will be moved around at the same
!  time but do not influence the sorting process (except for slowing things down of course!)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
  INTEGER, INTENT(INOUT), DIMENSION(:,:) :: A  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iq  ! ?
!===================================================================================================================================

  IF(size(A,1) .GT. 1) THEN
     CALL Partition3int(A, iq)
     CALL Qsort3int(A(:iq-1,:))
     CALL Qsort3int(A(iq:,:))
  END IF
END SUBROUTINE Qsort3int
SUBROUTINE Partition3int(A, marker)
!===================================================================================================================================
!  Sorting routine used by QSort3int above. This routine is PRIVATE
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
  INTEGER, INTENT(INOUT), DIMENSION(:,:) :: A  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER, INTENT(OUT)                   :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: i, j  ! ?
  INTEGER   :: temp(size(A,2))  ! ?
  INTEGER   :: x(3)              ! pivot point
!===================================================================================================================================
  x(:) = A(1,1:3)
  i= 0
  j= size(A,1) + 1

  DO
     j = j-1
     DO
        IF ((A(j,1) .LT. x(1)).OR.(((A(j,1).EQ.x(1)).AND.(A(j,2) .LT. x(2)))).OR. &
           (((A(j,1).EQ.x(1)).AND.(A(j,2).EQ.x(2)).AND.(A(j,3) .LE. x(3))))) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        if ((A(i,1) .GT. x(1)).OR.(((A(i,1).EQ.x(1)).AND.(A(i,2) .GT. x(2)))).OR. &
           (((A(i,1).EQ.x(1)).AND.(A(i,2).EQ.x(2)).AND.(A(i,3) .GE. x(3))))) EXIT
        i = i+1
     END DO
     IF (i .LT. j) THEN
        ! exchange A(i) and A(j)
        temp(:) = A(i,:)
        A(i,:) = A(j,:)
        A(j,:) = temp(:)
     ELSE IF (i.EQ.j) THEN
        marker = i+1
        RETURN
     ELSE
        marker = i
        RETURN
     END IF
  END DO
 RETURN
END SUBROUTINE Partition3int

RECURSIVE SUBROUTINE Qsort4int(A)
!===================================================================================================================================
! QSort3int:
!  Uses the Quicksort algorithm to sort an INTEGER table with four relevant columns.
!  There may be an arbitrary number of additional columns that will be moved around at the same
!  time but do not influence the sorting process (except for slowing things down of course!)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
INTEGER, INTENT(INOUT), DIMENSION(:,:) :: A  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iq  ! ?
!===================================================================================================================================
IF(size(A,1) .GT. 1) THEN
   CALL Partition4int(A, iq)
   CALL Qsort4int(A(:iq-1,:))
   CALL Qsort4int(A(iq:,:))
END IF
END SUBROUTINE Qsort4int
SUBROUTINE Partition4int(A, marker)
!===================================================================================================================================
!  Sorting routine used by QSort4int above. This routine is PRIVATE
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
INTEGER, INTENT(INOUT), DIMENSION(:,:) :: A   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)                   :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER   :: i,j                ! ?
INTEGER   :: temp(size(A,2))        ! ?
INTEGER   :: x(4)              ! pivot point
!===================================================================================================================================
x(:) = A(1,1:4)
i= 0
j= size(A,1) + 1

DO
  j = j-1
  DO
    IF(A(j,1) < x(1))THEN
      EXIT
    ELSE IF(A(j,1) == x(1))THEN
      IF(A(j,2) < x(2))THEN
        EXIT
      ELSE IF(A(j,2) == x(2))THEN
        IF(A(j,3) < x(3))THEN
          EXIT
        ELSE IF(A(j,3) == x(3))THEN
          IF(A(j,4) <= x(4))THEN
            EXIT
          END IF
        END IF
      END IF
    END IF
    j = j-1
  END DO
  i = i+1
  DO
    IF(A(i,1) > x(1))THEN
      EXIT
    ELSE IF(A(i,1) == x(1))THEN
      IF(A(i,2) > x(2))THEN
        EXIT
      ELSE IF(A(i,2) == x(2))THEN
        IF(A(i,3) > x(3))THEN
          EXIT
        ELSE IF(A(i,3) == x(3))THEN
          IF(A(i,4) >= x(4))THEN
            EXIT
          END IF
        END IF
      END IF
    END IF
    i = i+1
  END DO
  IF (i < j) THEN
     ! exchange A(i) and A(j)
     temp(:) = A(i,:)
     A(i,:) = A(j,:)
     A(j,:) = temp(:)
  ELSE IF (i == j) THEN
     marker = i+1
     RETURN
  ELSE
     marker = i
     RETURN
  END IF
END DO
END SUBROUTINE Partition4int
RECURSIVE SUBROUTINE Qsort2Real(A)
!===================================================================================================================================
! QSort2Real:
!  Uses the Quicksort algorithm to sort a REAL table with two relevant columns.
!  There may be an arbitrary number of additional columns that will be moved around at the same
!  time but do not influence the sorting process (except for slowing things down of course!)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
  REAL, INTENT(INOUT), DIMENSION(:,:) :: A    ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iq   ! ?
!===================================================================================================================================
  IF(size(A,1) > 1) THEN
     CALL Partition2Real(A, iq)
     CALL Qsort2Real(A(:iq-1,:))
     CALL Qsort2Real(A(iq:,:))
  END IF
END SUBROUTINE Qsort2Real

SUBROUTINE Partition2Real(A, marker)
!===================================================================================================================================
! Partition2Real:
!  Sorting routine used by QSort2Real above. This routine is PRIVATE.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
  REAL, INTENT(INOUT), DIMENSION(:,:) :: A   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER, INTENT(OUT)                :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: i, j     ! ?
  REAL      :: temp(size(A,2))   ! ?
  REAL      :: x(2)              ! pivot point
!===================================================================================================================================
  x(:) = A(1,1:2)
  i= 0
  j= size(A,1) + 1

  DO
     j = j-1
     DO
        IF ((A(j,1) < x(1)).OR.(((A(j,1) == x(1)).AND.(A(j,2) <= x(2))))) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        if ((A(i,1) > x(1)).OR.(((A(i,1) == x(1)).AND.(A(i,2) >= x(2))))) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        temp(:) = A(i,:)
        A(i,:) = A(j,:)
        A(j,:) = temp(:)
     ELSE IF (i == j) THEN
        marker = i+1
        RETURN
     ELSE
        marker = i
        RETURN
     END IF
  END DO

END SUBROUTINE Partition2Real

RECURSIVE SUBROUTINE Qsort3Real(A)
!===================================================================================================================================
! QSort3Real:
!  Uses the Quicksort algorithm to sort a REAL table with three relevant columns.
!  There may be an arbitrary number of additional columns that will be moved around at the same
!  time but do not influence the sorting process (except for slowing things down of course!)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
  REAL, INTENT(INOUT), DIMENSION(:,:) :: A    ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER :: iq      ! ?
!===================================================================================================================================
  IF(size(A,1) > 1) THEN
     CALL Partition3Real(A, iq)
     CALL Qsort3Real(A(:iq-1,:))
     CALL Qsort3Real(A(iq:,:))
  END IF
END SUBROUTINE Qsort3Real

SUBROUTINE Partition3Real(A, marker)
!===================================================================================================================================
! Partition3Real:
!  Sorting routine used by QSort3Real above. This routine is PRIVATE.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
  REAL, INTENT(INOUT), DIMENSION(:,:) :: A      ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
  INTEGER, INTENT(OUT)                :: marker   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
  INTEGER   :: i, j   ! ?
  REAL      :: temp(size(A,2))    ! ?
  REAL      :: x(3)              ! pivot point
!===================================================================================================================================
  x(:) = A(1,1:3)
  i= 0
  j= size(A,1) + 1

  DO
     j = j-1
     DO
        IF ((A(j,1) < x(1)).OR.(((A(j,1) == x(1)).AND.(A(j,2) < x(2)))).OR. &
           (((A(j,1) == x(1)).AND.(A(j,2) == x(2)).AND.(A(j,3) <= x(3))))) EXIT
        j = j-1
     END DO
     i = i+1
     DO
        if ((A(i,1) > x(1)).OR.(((A(i,1) == x(1)).AND.(A(i,2) > x(2)))).OR. &
           (((A(i,1) == x(1)).AND.(A(i,2) == x(2)).AND.(A(i,3) >= x(3))))) EXIT
        i = i+1
     END DO
     IF (i < j) THEN
        ! exchange A(i) and A(j)
        temp(:) = A(i,:)
        A(i,:) = A(j,:)
        A(j,:) = temp(:)
     ELSE IF (i == j) THEN
        marker = i+1
        RETURN
     ELSE
        marker = i
        RETURN
     END IF
  END DO

END SUBROUTINE Partition3Real


RECURSIVE SUBROUTINE MSortNInt(A,N,ind)
!===================================================================================================================================
! Sorts array of integers
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: N,ind
INTEGER,INTENT(INOUT) :: A(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i,first,sA
!===================================================================================================================================
sA=SIZE(A,2)
IF((N.LE.0).OR.(sA.EQ.1)) RETURN
IF(N > SIZE(A,1)-ind+1) STOP 'Number of entries to sort greater then array dimension!'

CALL MergeSortInt(A,ind)
IF(N.EQ.1) RETURN

first=1
DO i=1,sA-1
  IF(A(ind,i).NE.A(ind,i+1))THEN
    CALL MSortNInt(A(:,first:i),N-1,ind+1)
    first=i+1
  END IF
END DO
IF(first.NE.sA) CALL MSortNInt(A(:,first:sA),N-1,ind+1)

END SUBROUTINE MSortNInt

RECURSIVE SUBROUTINE MergeSortInt(A,ind)
!===================================================================================================================================
! Sorts array of integers
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: ind
INTEGER,INTENT(INOUT) :: A(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: nA,nB,tmp(SIZE(A,1)),nTotal
!===================================================================================================================================
nTotal=SIZE(A,2)
IF(nTotal.LT.2) RETURN
IF(nTotal.EQ.2)THEN
  IF(A(ind,1).GT.A(ind,2))THEN
    tmp    = A(:,1)
    A(:,1) = A(:,2)
    A(:,2) = tmp
  ENDIF
  RETURN
ENDIF
nA=(nTotal+1)/2
CALL MergeSortInt(A(:,1:nA),ind)
nB=nTotal-nA
CALL MergeSortInt(A(:,nA+1:nTotal),ind)
! Performed first on lowest level
IF(A(ind,nA).GT.A(ind,nA+1)) CALL DoMergeInt(A,nA,nB,SIZE(A,1),ind)
END SUBROUTINE MergeSortInt

SUBROUTINE DoMergeInt(A,nA,nB,nInd,ind)
!===================================================================================================================================
! Merge subarrays
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: nA ! ?
INTEGER,INTENT(IN)    :: nB ! ?
INTEGER,INTENT(IN)    :: nInd ! ?
INTEGER,INTENT(IN)    :: ind ! ?
INTEGER,INTENT(INOUT) :: A(nInd,nA+nB) ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k
INTEGER :: part1(nInd,nA),part2(nInd,nB)
!===================================================================================================================================
part1(:,1:nA)=A(:,1:nA)
part2(:,1:nB)=A(:,nA+1:nA+nB)
i=1; j=1; k=1;
DO WHILE((i.LE.nA).AND.(j.LE.nB))
  IF(part1(ind,i).LE.part2(ind,j))THEN
    A(:,k)=part1(:,i)
    i=i+1
  ELSE
    A(:,k)=part2(:,j)
    j=j+1
  ENDIF
  k=k+1
END DO
j=nA-i
A(:,k:k+nA-i)=part1(:,i:nA)
END SUBROUTINE DoMergeInt



RECURSIVE SUBROUTINE MSortNLong(A,N,ind)
!===================================================================================================================================
! Sorts array of integers
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: N,ind
INTEGER(KIND=8),INTENT(INOUT) :: A(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: i,first,sA
!===================================================================================================================================
sA=SIZE(A,2)
IF((N.LE.0).OR.(sA.EQ.1)) RETURN
IF(N > SIZE(A,1)-ind+1) STOP 'Number of entries to sort greater then array dimension!'

CALL MergeSortLong(A,ind)
IF(N.EQ.1) RETURN

first=1
DO i=1,sA-1
  IF(A(ind,i).NE.A(ind,i+1))THEN
    CALL MSortNLong(A(:,first:i),N-1,ind+1)
    first=i+1
  END IF
END DO
IF(first.NE.sA) CALL MSortNLong(A(:,first:sA),N-1,ind+1)

END SUBROUTINE MSortNLong
RECURSIVE SUBROUTINE MergeSortLong(A,ind)
!===================================================================================================================================
! Sorts array of integers
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: ind
INTEGER(KIND=8),INTENT(INOUT) :: A(:,:)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER(KIND=8)       :: tmp(SIZE(A,1))
INTEGER               :: nA,nB,nTotal
!===================================================================================================================================
nTotal=SIZE(A,2)
IF(nTotal.LT.2) RETURN
IF(nTotal.EQ.2)THEN
  IF(A(ind,1).GT.A(ind,2))THEN
    tmp    = A(:,1)
    A(:,1) = A(:,2)
    A(:,2) = tmp
  ENDIF
  RETURN
ENDIF
nA=(nTotal+1)/2
CALL MergeSortLong(A(:,1:nA),ind)
nB=nTotal-nA
CALL MergeSortLong(A(:,nA+1:nTotal),ind)
! Performed first on lowest level
IF(A(ind,nA).GT.A(ind,nA+1)) CALL DoMergeLong(A,nA,nB,SIZE(A,1),ind)
END SUBROUTINE MergeSortLong
SUBROUTINE DoMergeLong(A,nA,nB,nInd,ind)
!===================================================================================================================================
! Merge subarrays
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)    :: nA ! ?
INTEGER,INTENT(IN)    :: nB ! ?
INTEGER,INTENT(IN)    :: nInd ! ?
INTEGER,INTENT(IN)    :: ind ! ?
INTEGER(KIND=8),INTENT(INOUT) :: A(nInd,nA+nB)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,j,k
INTEGER(KIND=8) :: part1(nInd,nA),part2(nInd,nB)
!===================================================================================================================================
part1(:,1:nA)=A(:,1:nA)
part2(:,1:nB)=A(:,nA+1:nA+nB)
i=1; j=1; k=1;
DO WHILE((i.LE.nA).AND.(j.LE.nB))
  IF(part1(ind,i).LE.part2(ind,j))THEN
    A(:,k)=part1(:,i)
    i=i+1
  ELSE
    A(:,k)=part2(:,j)
    j=j+1
  ENDIF
  k=k+1
END DO
j=nA-i
A(:,k:k+nA-i)=part1(:,i:nA)
END SUBROUTINE DoMergeLong


END MODULE MOD_SortingTools
