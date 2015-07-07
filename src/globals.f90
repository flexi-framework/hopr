#include "defines.f90"
MODULE MOD_Globals
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
#ifndef F03
USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
REAL                        :: PI                         ! container for Pi=3.141592...
REAL                        :: rPI                        ! rPI = SQRT(Pi)
! Unit numbers for file io
#ifndef F03
INTEGER, PARAMETER          :: UNIT_stdIn  = input_unit   ! Terminal input
INTEGER, PARAMETER          :: UNIT_stdOut = output_unit  ! Terminal output
INTEGER, PARAMETER          :: UNIT_errOut = error_unit   ! For error output
#else
INTEGER, PARAMETER          :: UNIT_stdIn  = 5            ! Terminal input
INTEGER, PARAMETER          :: UNIT_stdOut = 6            ! Terminal output
INTEGER, PARAMETER          :: UNIT_errOut = 0            ! For error output
#endif
INTEGER, PARAMETER          :: UNIT_logOut = 100          ! For logging
CHARACTER(LEN=255)          :: ProjectName                ! necessary data for in/output and name used to generate filenames from
LOGICAL                     :: Logging                    ! Set .TRUE. to activate logging function for each processor



INTERFACE Abort
   MODULE PROCEDURE Abort
END INTERFACE

INTERFACE Timer
   MODULE PROCEDURE Timer
END INTERFACE

INTERFACE IS_NAN
   MODULE PROCEDURE IS_NAN
END INTERFACE

INTERFACE CROSS
   MODULE PROCEDURE CROSS
END INTERFACE

INTERFACE NORMALIZE
   MODULE PROCEDURE NORMALIZE
END INTERFACE

INTERFACE getDet3
   MODULE PROCEDURE getDet3
END INTERFACE

INTERFACE getInv3
   MODULE PROCEDURE getInv3
END INTERFACE

INTERFACE getInverse
   MODULE PROCEDURE getInverse
END INTERFACE

!===================================================================================================================================

CONTAINS
SUBROUTINE Abort(SourceFile,SourceLine,CompDate,CompTime,ErrorMessage,IntInfoOpt,RealInfoOpt)
!===================================================================================================================================
! Terminate program correctly if an error has occurred (important in MPI mode!).
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)                  :: SourceFile      ! Source file where error has occurred
INTEGER,INTENT(IN)                           :: SourceLine      ! Line in source file
CHARACTER(LEN=*),INTENT(IN)                  :: CompDate        ! Compilation date
CHARACTER(LEN=*),INTENT(IN)                  :: CompTime        ! Compilation time
CHARACTER(LEN=*),INTENT(IN)                  :: ErrorMessage    ! Error message
INTEGER,OPTIONAL,INTENT(IN)                  :: IntInfoOpt      ! Error info (integer)
REAL,OPTIONAL,INTENT(IN)                     :: RealInfoOpt     ! Error info (real)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!   There is no way back!
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                           :: IntInfo         ! Error info (integer)
REAL                              :: RealInfo        ! Error info (real)
!===================================================================================================================================
IntInfo  = MERGE(IntInfoOpt ,999 ,PRESENT(IntInfoOpt) )
RealInfo = MERGE(RealInfoOpt,999.,PRESENT(RealInfoOpt))
WRITE(UNIT_stdOut,*)
WRITE(UNIT_stdOut,*)'_____________________________________________________________________________'
WRITE(UNIT_stdOut,*)'Program abort caused in File : ',TRIM(SourceFile),' Line ',SourceLine
WRITE(UNIT_stdOut,*)'This file was compiled at ',TRIM(CompDate),'  ',TRIM(CompTime)
WRITE(UNIT_stdOut,'(A10,A)',ADVANCE='NO')'Message: ',TRIM(ErrorMessage)
IF(IntInfo  .NE. 999 ) WRITE(UNIT_stdOut,'(I8)',ADVANCE='NO')IntInfo
IF(RealInfo .NE. 999.) WRITE(UNIT_stdOut,'(E16.8)')RealInfo
WRITE(UNIT_stdOut,*)
WRITE(UNIT_stdOut,'(A,A,A,I6.6,A)')'See ',TRIM(ProjectName),'_ERRORS.out for more details'
WRITE(UNIT_stdOut,*)
STOP 0001
END SUBROUTINE Abort

SUBROUTINE Timer(start,unit_in)
!===================================================================================================================================
! Get average and max time needed for a task. You always have to call this routine twice. The first time you call it
! ("start"=.TRUE.) the start time a new TimerObject will be created
! will be stored in "startTime", "start" is set to .FALSE. (no barrier in MPI case). The second time, the end time will be measured
! and the difference between start and end time will be calculated. In the MPI case, the maximum and the average of this time is
! calculated over all processors (blocking / synchronizing). "start" is reset to .TRUE.
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, OPTIONAL,INTENT(IN)              :: unit_in            ! UNIT for output. If not set, UNIT_stdOut will be used and
                                                                ! only proc 0 writes messages to terminal. Otherwise all procs
                                                                ! write to UNIT given by unit_in.
LOGICAL,INTENT(IN)                        :: start              ! Start / stop switch
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!   Returns average and max time and corresponding processors
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE tTimerObject
  TYPE(tTimerObject),POINTER   :: nextTimer                     ! Pointer to next timer object
  REAL                         :: StartTime                     ! Time of first call
END TYPE tTimerObject
TYPE(tTimerObject),POINTER,SAVE:: lastTimer    ! ?
TYPE(tTimerObject),POINTER     :: TimerPtr  ! ?
REAL                           :: sendBuffer(2)                 ! Just a buffer
INTEGER, SAVE                  :: nTimerObjects=0    ! ?
INTEGER                        :: output_unit                   ! UNIT for output
!===================================================================================================================================
IF(start)THEN
  ALLOCATE(TimerPtr)
  NULLIFY(TimerPtr%nextTimer)
  ! Get start time
  CALL CPU_TIME(TimerPtr%StartTime)
  IF(nTimerObjects .EQ. 0)THEN
    lastTimer=>TimerPtr
  ELSE
    TimerPtr%nextTimer=>lastTimer
    lastTimer=>TimerPtr
  END IF
  nTimerObjects=nTimerObjects+1
  NULLIFY(TimerPtr)
ELSE
  output_unit=UNIT_stdOut
  IF(PRESENT(unit_in)) output_unit=unit_in
  ! Get end time / local difference between end time and start time
  TimerPtr=>lastTimer
  CALL CPU_TIME(sendBuffer(1))
  sendBuffer(2)=sendBuffer(1)-TimerPtr%startTime
  WRITE(output_unit,'(A,25("   ."),A,F0.3,A)')'DONE!','   [ time: ',sendBuffer(2),'s ]'
  lastTimer=>TimerPtr%nextTimer
  NULLIFY(TimerPtr%nextTimer)
  DEALLOCATE(TimerPtr)
  NULLIFY(TimerPtr)
  nTimerObjects=nTimerObjects-1
END IF
END SUBROUTINE Timer


FUNCTION IS_NAN(x)
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
#ifdef IEEE_ISNAN 
USE IEEE_ARITHMETIC
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: x              ! Floating point value to be tested
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
LOGICAL :: is_nan         ! Number or not a number?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
#ifndef IEEE_ISNAN
REAL    :: testvalue  ! ?
#endif
!===================================================================================================================================
#ifdef IEEE_ISNAN
is_nan=IEEE_IS_NAN(x)
#else
is_nan=.FALSE.
testvalue=0.
! Too large
IF(ABS(x) .GT. HUGE(testvalue))THEN
  ERRWRITE(*,*) 'Floating overflow detected!'
  ERRWRITE(*,*) 'x=',x
  is_nan=.TRUE.
ENDIF
! Not defined
IF(x .NE. x)THEN 
  ERRWRITE(*,*) 'Floating invalid detected!'
  ERRWRITE(*,*) 'x=',x
  is_nan=.TRUE.
ENDIF
#endif
END FUNCTION IS_NAN


PURE FUNCTION NORMALIZE(v1,nVal)
!===================================================================================================================================
! normalizes a nDim vector with respect to the eucledian norm
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN)  :: nVal ! ?
REAL,INTENT(IN)     :: v1(nVal) ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: normalize(nVal) ! ? 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
normalize=v1/SQRT(SUM(v1*v1))
END FUNCTION NORMALIZE


PURE FUNCTION CROSS(v1,v2)
!===================================================================================================================================
! computes the cross product of to 3 dimensional vectpors: cross=v1 x v2
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN) :: v1(3)    ! ? 
REAL,INTENT(IN) :: v2(3)    ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL            :: cross(3) ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
cross=(/v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3),v1(1)*v2(2)-v1(2)*v2(1)/)
END FUNCTION CROSS
FUNCTION getDet3(Mat)
!===================================================================================================================================
! compute determinant of 3x3 matrix
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Mat(3,3) ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: getDet3 ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
getDet3=   ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * Mat(3,3) &
         + ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * Mat(3,1) &
         + ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * Mat(3,2)
END FUNCTION getDet3

FUNCTION getInv3(Mat,sDet_in)
!===================================================================================================================================
! compute inverse of 3x3 matrix, needs sDet=1/det(Mat)
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)  :: Mat(3,3) ! ?
REAL,INTENT(IN),OPTIONAL    :: sDet_in ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL             :: getInv3(3,3) ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL             :: sDet ! ?
!===================================================================================================================================
IF(PRESENT(sDet_in))THEN
  sDet=sDet_in
ELSE
  sDet=1./getDet3(Mat)
END IF
getInv3(1,1) = ( Mat(2,2) * Mat(3,3) - Mat(2,3) * Mat(3,2) ) * sDet
getInv3(1,2) = ( Mat(1,3) * Mat(3,2) - Mat(1,2) * Mat(3,3) ) * sDet
getInv3(1,3) = ( Mat(1,2) * Mat(2,3) - Mat(1,3) * Mat(2,2) ) * sDet
getInv3(2,1) = ( Mat(2,3) * Mat(3,1) - Mat(2,1) * Mat(3,3) ) * sDet
getInv3(2,2) = ( Mat(1,1) * Mat(3,3) - Mat(1,3) * Mat(3,1) ) * sDet
getInv3(2,3) = ( Mat(1,3) * Mat(2,1) - Mat(1,1) * Mat(2,3) ) * sDet
getInv3(3,1) = ( Mat(2,1) * Mat(3,2) - Mat(2,2) * Mat(3,1) ) * sDet
getInv3(3,2) = ( Mat(1,2) * Mat(3,1) - Mat(1,1) * Mat(3,2) ) * sDet
getInv3(3,3) = ( Mat(1,1) * Mat(2,2) - Mat(1,2) * Mat(2,1) ) * sDet
END FUNCTION getInv3



FUNCTION GetInverse(dim1,A) RESULT(Ainv)
!===================================================================================================================================
! invert a matrix (dependant in LAPACK Routines) 
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER, INTENT(IN) :: dim1   ! size of matrix a
REAL,INTENT(IN)     :: A(dim1,dim1) ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL                :: Ainv(dim1,dim1) ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: IPIV(dim1),INFO,lwork ! ?
REAL               :: WORK(dim1*dim1)  ! ?
!===================================================================================================================================
! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  CALL DGETRF(dim1, dim1, Ainv, dim1, IPIV, INFO)

  IF (INFO /= 0) THEN
     STOP 'MATRIX IS NUMERICALLY SINGULAR!'
  END IF

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  lwork=dim1*dim1
  CALL DGETRI(dim1, Ainv, dim1, IPIV, WORK, lwork , INFO)

  IF (INFO /= 0) THEN
     STOP 'MATRIX INVERSION FAILED!'
  END IF
END FUNCTION GetInverse


END MODULE MOD_Globals
