PROGRAM test

!Implicit Treatment
IMPLICIT NONE

TYPE tMortonBox
  REAL(KIND=8) :: mini(3)
  INTEGER(KIND=4) :: nbits
  REAL(KIND=8) :: spacing(3)
END TYPE tMortonBox

INTERFACE
   FUNCTION evalhilbert(disc,nbits,ndims) Result(indx)
       INTEGER(KIND=8)  :: disc(3),indx
       INTEGER(KIND=4)  :: nbits,ndims
   END FUNCTION evalhilbert
END INTERFACE

INTEGER(KIND=4) :: nDims,nBits
INTEGER(KIND=8) :: disc(3),res
REAL            :: coord(3,8)
TYPE(tMortonBox) :: box 
INTEGER            :: ii

!CALL printer(3)

nDims = 3
nBits = 21
 coord(1:3,1)=(/0.,0.,0./)
 coord(1:3,2)=(/1.,0.,0./)
 coord(1:3,3)=(/1.,1.,0./)
 coord(1:3,4)=(/0.,1.,0./)
 coord(1:3,5)=(/0.,0.,1./)
 coord(1:3,6)=(/1.,0.,1./)
 coord(1:3,7)=(/1.,1.,1./)
 coord(1:3,8)=(/0.,1.,1./)

box%mini(1:3) = 0
box%nbits = 21
box%spacing(1:3) = (2**21)/4

box%spacing(1:3) = REAL(2**box%nbits-1)/(4-box%mini(1:3))

Do ii=1,8
    disc(1:3) =  NINT((coord(1:3,ii)-box%mini(1:3))*box%spacing(1:3))
    res = evalhilbert(disc(1:3),box%nbits,nDims)
    WRITE(*,*) res
End Do 

END PROGRAM test
