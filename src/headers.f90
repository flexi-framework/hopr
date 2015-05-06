! cgnslib_f.h defines a parameter NULL which confilcts with the Fortran 95
! function NULL(). Use CGNS_header only for IO routines.
MODULE CGNS_header
  IMPLICIT NONE
  PUBLIC
  PRIVATE  :: Null
  SAVE
  
  ! CGNS header file
  INCLUDE 'cgnslib_f.h'

#  ifdef BLUEGENE
  ! Explicit interface for JUGENE
  interface cg_goto_f
     module procedure cg_goto_4
     module procedure cg_goto_6
     module procedure cg_goto_8
     module procedure cg_goto_10
     module procedure cg_goto_12
  end interface

contains

  subroutine cg_goto_4(FileName, Base, iError, Label)
    implicit none
    ! Arguments:
    integer :: FileName
    integer :: Base
    integer :: iError
    character(LEN=*) :: Label

    external :: cg_goto_f

    call cg_goto_f(FileName, Base, iError, Label)
  end subroutine cg_goto_4

  subroutine cg_goto_6(FileName, Base, iError, TypLabel, TypIndex, Label)
    implicit none
    ! Arguments:
    integer :: FileName
    integer :: Base
    integer :: iError
    integer :: TypIndex
    character(LEN=*) :: Label
    character(LEN=*) :: TypLabel

    external :: cg_goto_f

    call cg_goto_f(FileName, Base, iError, TypLabel, TypIndex, Label)
  end subroutine cg_goto_6

  subroutine cg_goto_8(FileName, Base, iError, &
                       TypLabel1, TypIndex1, TypLabel2, TypIndex2, Label)
    implicit none
    ! Arguments:
    integer :: FileName
    integer :: Base
    integer :: iError
    integer :: TypIndex1
    integer :: TypIndex2
    character(LEN=*) :: Label
    character(LEN=*) :: TypLabel1
    character(LEN=*) :: TypLabel2

    external :: cg_goto_f

    call cg_goto_f(FileName, Base, iError, &
                   TypLabel1, TypIndex1, TypLabel2, TypIndex2, Label)
  end subroutine cg_goto_8

  subroutine cg_goto_10(FileName, Base, iError, &
                       TypLabel1, TypIndex1, TypLabel2, TypIndex2, TypLabel3, TypIndex3, Label)
    implicit none
    ! Arguments:
    integer :: FileName
    integer :: Base
    integer :: iError
    integer :: TypIndex1
    integer :: TypIndex2
    integer :: TypIndex3
    character(LEN=*) :: Label
    character(LEN=*) :: TypLabel1
    character(LEN=*) :: TypLabel2
    character(LEN=*) :: TypLabel3

    external :: cg_goto_f

    call cg_goto_f(FileName, Base, iError, &
                   TypLabel1, TypIndex1, TypLabel2, TypIndex2, TypLabel3, TypIndex3, Label)
  end subroutine cg_goto_10

  subroutine cg_goto_12(FileName, Base, iError, &
                       TypLabel1, TypIndex1, TypLabel2, TypIndex2, TypLabel3, TypIndex3, TypLabel4, TypIndex4, Label)
    implicit none
    ! Arguments:
    integer :: FileName
    integer :: Base
    integer :: iError
    integer :: TypIndex1
    integer :: TypIndex2
    integer :: TypIndex3
    integer :: TypIndex4
    character(LEN=*) :: Label
    character(LEN=*) :: TypLabel1
    character(LEN=*) :: TypLabel2
    character(LEN=*) :: TypLabel3
    character(LEN=*) :: TypLabel4

    external :: cg_goto_f

    call cg_goto_f(FileName, Base, iError, &
                   TypLabel1, TypIndex1, TypLabel2, TypIndex2, TypLabel3, TypIndex3, TypLabel4, TypIndex4, Label)
  end subroutine cg_goto_12

#  endif
END MODULE CGNS_header
