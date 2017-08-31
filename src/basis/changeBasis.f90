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
MODULE MOD_ChangeBasis
!==================================================================================================================================
! Changes a 2D or 3D Tensor Product Lagrange Points of Lagrange Basis of degree N_In to
! Lagrange points of a Lagrange Basis N_Out, using two
! arbitrary point disributions xi_In(0:N_In) and xi_Out(0:N_Out)
!==================================================================================================================================
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ChangeBasis3D
  MODULE PROCEDURE ChangeBasis3D_Single
  MODULE PROCEDURE ChangeBasis3D_Mult
END INTERFACE
!
INTERFACE ChangeBasis2D
  MODULE PROCEDURE ChangeBasis2D_Single
  MODULE PROCEDURE ChangeBasis2D_Mult
END INTERFACE

INTERFACE ChangeBasis3D_XYZ
  MODULE PROCEDURE ChangeBasis3D_XYZ
END INTERFACE

PUBLIC :: ChangeBasis3D
PUBLIC :: ChangeBasis2D
PUBLIC :: ChangeBasis3D_XYZ
!==================================================================================================================================

CONTAINS



SUBROUTINE ChangeBasis3D_Mult(nVar,nElems,NIn,NOut,Vdm,UIn,UOut,addToOutput)
!==================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
! to another 3D tensor product node positions (number of nodes N_out+1)
! defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!  xi is defined in the 1DrefElem xi=[-1,1]
!==================================================================================================================================
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: nVar,nElems,NIn,NOut
REAL,INTENT(IN)     :: UIn(nVar,0:NIn,0:NIn,0:NIn,nElems)
REAL,INTENT(IN)     :: Vdm(0:NOut,0:NIn)
LOGICAL,INTENT(IN)  :: addToOutput
REAL,INTENT(INOUT)  :: UOut(nVar,0:NOut,0:NOut,0:NOut,nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iI,jI,kI,iO,jO,kO,iElem,a,b,nVar2
REAL,ALLOCATABLE    :: UBuf1(:,:,:,:),UBuf2(:,:,:,:)
!==================================================================================================================================
nVar2=nVar*nElems
IF(nVar2.GT.2*nVar)THEN
  ALLOCATE(UBuf2(nVar2,0:NIn,0:NIn,0:NIn))
  ALLOCATE(UBuf1(nVar2,0:NOut,0:NIn,0:NIn))

  ! pack solution
  DO iElem=1,nElems
    a=nVar*(iElem-1)+1
    b=nVar*iElem
    DO kI=0,NIn; DO jI=0,NIn; DO iI=0,NIn
      Ubuf2(a:b,iI,jI,kI)=UIn(:,iI,jI,kI,iElem)
    END DO; END DO; END DO
  END DO

  ! first direction iI
  DO kI=0,NIn; DO jI=0,NIn
    DO iO=0,NOut
      UBuf1(:,iO,jI,kI)=Vdm(iO,0)*Ubuf2(:,0,jI,kI)
    END DO
    DO iI=1,NIn
      DO iO=0,NOut
        UBuf1(:,iO,jI,kI)=UBuf1(:,iO,jI,kI)+Vdm(iO,iI)*Ubuf2(:,iI,jI,kI)
      END DO
    END DO
  END DO; END DO

  DEALLOCATE(Ubuf2)
  ALLOCATE(UBuf2(nVar2,0:NOut,0:NOut,0:NIn))

  ! second direction jI
  DO kI=0,NIn
    DO jO=0,NOut; DO iO=0,NOut
      UBuf2(:,iO,jO,kI)=Vdm(jO,0)*UBuf1(:,iO,0,kI)
    END DO; END DO
    DO jI=1,NIn
      DO jO=0,NOut; DO iO=0,NOut
        UBuf2(:,iO,jO,kI)=UBuf2(:,iO,jO,kI)+Vdm(jO,jI)*UBuf1(:,iO,jI,kI)
      END DO; END DO
    END DO
  END DO

  DEALLOCATE(Ubuf1)
  ALLOCATE(UBuf1(nVar2,0:NOut,0:NOut,0:NOut))

  ! last direction kI
  DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
    Ubuf1(:,iO,jO,kO)=Vdm(kO,0)*UBuf2(:,iO,jO,0)
  END DO; END DO; END DO
  DO kI=1,NIn
    DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
      Ubuf1(:,iO,jO,kO)=Ubuf1(:,iO,jO,kO)+Vdm(kO,kI)*UBuf2(:,iO,jO,kI)
    END DO; END DO; END DO
  END DO

  ! unpack solution
  IF(addToOutput)THEN
    DO iElem=1,nElems
      a=nVar*(iElem-1)+1
      b=nVar*iElem
      DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
        UOut(:,iO,jO,kO,iElem)=UOut(:,iO,jO,kO,iElem)+Ubuf1(a:b,iO,jO,kO)
      END DO; END DO; END DO
    END DO
  ELSE
    DO iElem=1,nElems
      a=nVar*(iElem-1)+1
      b=nVar*iElem
      DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
        UOut(:,iO,jO,kO,iElem)=Ubuf1(a:b,iO,jO,kO)
      END DO; END DO; END DO
    END DO
  END IF
  DEALLOCATE(UBuf1,Ubuf2)

ELSE

  ALLOCATE(UBuf1(nVar,0:NOut,0:NIn,0:NIn))
  ALLOCATE(UBuf2(nVar,0:NOut,0:NOut,0:NIn))
  DO iElem=1,nElems
    ! first direction iI
    DO kI=0,NIn; DO jI=0,NIn
      DO iO=0,NOut
        UBuf1(:,iO,jI,kI)=Vdm(iO,0)*UIn(:,0,jI,kI,iElem)
      END DO
      DO iI=1,NIn
        DO iO=0,NOut
          UBuf1(:,iO,jI,kI)=UBuf1(:,iO,jI,kI)+Vdm(iO,iI)*UIn(:,iI,jI,kI,iElem)
        END DO
      END DO
    END DO; END DO
    ! second direction jI
    DO kI=0,NIn
      DO jO=0,NOut; DO iO=0,NOut
        UBuf2(:,iO,jO,kI)=Vdm(jO,0)*UBuf1(:,iO,0,kI)
      END DO; END DO
      DO jI=1,NIn
        DO jO=0,NOut; DO iO=0,NOut
          UBuf2(:,iO,jO,kI)=UBuf2(:,iO,jO,kI)+Vdm(jO,jI)*UBuf1(:,iO,jI,kI)
        END DO; END DO
      END DO
    END DO
    ! last direction kI
    IF(addToOutput)THEN
      DO kI=0,NIn
        DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
          UOut(:,iO,jO,kO,iElem)=UOut(:,iO,jO,kO,iElem)+Vdm(kO,kI)*UBuf2(:,iO,jO,kI)
        END DO; END DO; END DO
      END DO
    ELSE
      DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
        UOut(:,iO,jO,kO,iElem)=Vdm(kO,0)*UBuf2(:,iO,jO,0)
      END DO; END DO; END DO
      DO kI=1,NIn
        DO kO=0,NOut; DO jO=0,NOut; DO iO=0,NOut
          UOut(:,iO,jO,kO,iElem)=UOut(:,iO,jO,kO,iElem)+Vdm(kO,kI)*UBuf2(:,iO,jO,kI)
        END DO; END DO; END DO
      END DO
    END IF
  END DO
  DEALLOCATE(UBuf1,Ubuf2)

END IF
END SUBROUTINE ChangeBasis3D_Mult


SUBROUTINE ChangeBasis3D_Single(Dim1,N_In,N_Out,Vdm,X3D_In,X3D_Out)
!==================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
! to another 3D tensor product node positions (number of nodes N_out+1)
! defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!  xi is defined in the 1DrefElem xi=[-1,1]
!==================================================================================================================================
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: Dim1,N_In,N_Out
REAL,INTENT(IN)     :: X3D_In(1:Dim1,0:N_In,0:N_In,0:N_In)
REAL,INTENT(IN)     :: Vdm(0:N_Out,0:N_In)
REAL,INTENT(OUT)    :: X3D_Out(1:Dim1,0:N_Out,0:N_Out,0:N_Out)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iN_In,jN_In,kN_In,iN_Out,jN_Out,kN_Out
REAL                :: X3D_Buf1(1:Dim1,0:N_Out,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:Dim1,0:N_Out,0:N_Out,0:N_In) ! second intermediate results from 1D interpolations
!==================================================================================================================================
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
END SUBROUTINE ChangeBasis3D_Single



SUBROUTINE ChangeBasis2D_Mult(nVar,firstSide,iStart,iEnd,NIn,NOut,Vdm,UIn,UOut,addToOutput)
!==================================================================================================================================
! interpolate a 2D tensor product Lagrange basis defined by (Nin+1) 1D interpolation point positions xi_In(0:NIn)
! to another 2D tensor product node positions (number of nodes Nout+1)
! defined by (Nout+1) interpolation point  positions xi_Out(0:NOut)
!  xi is defined in the 1DrefElem xi=[-1,1]
!==================================================================================================================================
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: nVar,firstSide,iStart,iEnd,NIn,NOut
REAL,INTENT(IN)     :: UIn(nVar,0:NIn,0:NIn,firstSide:iEnd)
REAL,INTENT(IN)     :: Vdm(0:NOut,0:NIn)
LOGICAL,INTENT(IN)  :: addToOutput
REAL,INTENT(INOUT)  :: UOut(nVar,0:NOut,0:NOut,firstSide:iEnd)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iI,jI,iO,jO,iSide,a,b,nVar2,nLocSides
REAL,ALLOCATABLE    :: UBuf1(:,:,:),UBuf2(:,:,:)
!==================================================================================================================================
IF(iEnd.LT.iStart) RETURN

nLocSides=iEnd-iStart+1
nVar2=nVar*nLocSides
IF(nVar2.GT.2*nVar)THEN
  ALLOCATE(UBuf2(nVar2,0:NIn,0:NIn))
  ALLOCATE(UBuf1(nVar2,0:NOut,0:NIn))
  ! pack solution
  DO iSide=iStart,iEnd
    a=nVar*(iSide-iStart)+1
    b=nVar*(iSide-iStart+1)
    DO jI=0,NIn; DO iI=0,NIn
      Ubuf2(a:b,iI,jI)=UIn(:,iI,jI,iSide)
    END DO; END DO
  END DO

  ! first direction iI
  DO jI=0,NIn
    DO iO=0,NOut
      UBuf1(:,iO,jI)=Vdm(iO,0)*Ubuf2(:,0,jI)
    END DO
    DO iI=1,NIn
      DO iO=0,NOut
        UBuf1(:,iO,jI)=UBuf1(:,iO,jI)+Vdm(iO,iI)*Ubuf2(:,iI,jI)
      END DO
    END DO
  END DO

  DEALLOCATE(UBuf2)
  ALLOCATE(UBuf2(nVar2,0:NOut,0:NOut))

  ! second direction jI
  DO jO=0,NOut; DO iO=0,NOut
    UBuf2(:,iO,jO)=Vdm(jO,0)*UBuf1(:,iO,0)
  END DO; END DO
  DO jI=1,NIn
    DO jO=0,NOut; DO iO=0,NOut
      UBuf2(:,iO,jO)=UBuf2(:,iO,jO)+Vdm(jO,jI)*UBuf1(:,iO,jI)
    END DO; END DO
  END DO

  ! unpack solution
  IF(addToOutput)THEN
    DO iSide=iStart,iEnd
      a=nVar*(iSide-iStart)+1
      b=nVar*(iSide-iStart+1)
      DO jO=0,NOut; DO iO=0,NOut
        UOut(:,iO,jO,iSide)=UOut(:,iO,jO,iSide)+Ubuf2(a:b,iO,jO)
      END DO; END DO
    END DO
  ELSE
    DO iSide=iStart,iEnd
      a=nVar*(iSide-iStart)+1
      b=nVar*(iSide-iStart+1)
      DO jO=0,NOut; DO iO=0,NOut
        UOut(:,iO,jO,iSide)=Ubuf2(a:b,iO,jO)
      END DO; END DO
    END DO
  END IF
  DEALLOCATE(UBuf1,UBuf2)

ELSE

  ALLOCATE(UBuf1(nVar,0:NOut,0:NIn))
  DO iSide=iStart,iEnd
    ! first direction iI
    DO jI=0,NIn
      DO iO=0,NOut
        UBuf1(:,iO,jI)=Vdm(iO,0)*UIn(:,0,jI,iSide)
      END DO
      DO iI=1,NIn
        DO iO=0,NOut
          UBuf1(:,iO,jI)=UBuf1(:,iO,jI)+Vdm(iO,iI)*UIn(:,iI,jI,iSide)
        END DO
      END DO
    END DO

    ! second direction jI
    IF(addToOutput)THEN
      DO jI=0,NIn
        DO jO=0,NOut; DO iO=0,NOut
          UOut(:,iO,jO,iSide)=UOut(:,iO,jO,iSide)+Vdm(jO,jI)*UBuf1(:,iO,jI)
        END DO; END DO
      END DO
    ELSE
      DO jO=0,NOut; DO iO=0,NOut
        UOut(:,iO,jO,iSide)=Vdm(jO,0)*UBuf1(:,iO,0)
      END DO; END DO
      DO jI=1,NIn
        DO jO=0,NOut; DO iO=0,NOut
          UOut(:,iO,jO,iSide)=UOut(:,iO,jO,iSide)+Vdm(jO,jI)*UBuf1(:,iO,jI)
        END DO; END DO
      END DO
    END IF

  END DO
  DEALLOCATE(UBuf1)
END IF
END SUBROUTINE ChangeBasis2D_Mult


SUBROUTINE ChangeBasis2D_Single(Dim1,N_In,N_Out,Vdm,X2D_In,X2D_Out)
!==================================================================================================================================
! interpolate a 2D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
! to another 2D tensor product node positions (number of nodes N_out+1)
! defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!  xi is defined in the 1DrefElem xi=[-1,1]
!==================================================================================================================================
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: Dim1,N_In,N_Out
REAL,INTENT(IN)     :: X2D_In(1:Dim1,0:N_In,0:N_In)
REAL,INTENT(IN)     :: Vdm(0:N_Out,0:N_In)
REAL,INTENT(OUT)    :: X2D_Out(1:Dim1,0:N_Out,0:N_Out)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iN_In,jN_In,iN_Out,jN_Out
REAL                :: X2D_Buf1(1:Dim1,0:N_Out,0:N_In)  ! first intermediate results from 1D interpolations
!==================================================================================================================================
X2D_buf1=0.
! first direction iN_In
DO jN_In=0,N_In
  DO iN_In=0,N_In
    DO iN_Out=0,N_Out
      X2D_Buf1(:,iN_Out,jN_In)=X2D_Buf1(:,iN_Out,jN_In)+Vdm(iN_Out,iN_In)*X2D_In(:,iN_In,jN_In)
    END DO
  END DO
END DO
X2D_Out=0.
! second direction jN_In
DO jN_In=0,N_In
  DO jN_Out=0,N_Out
    DO iN_Out=0,N_Out
      X2D_Out(:,iN_Out,jN_Out)=X2D_Out(:,iN_Out,jN_Out)+Vdm(jN_Out,jN_In)*X2D_Buf1(:,iN_Out,jN_In)
    END DO
  END DO
END DO
END SUBROUTINE ChangeBasis2D_Single



SUBROUTINE ChangeBasis3D_XYZ(Dim1,N_In,N_Out,Vdm_xi,Vdm_eta,Vdm_zeta,X3D_In,X3D_Out)
!==================================================================================================================================
! interpolate a 3D tensor product Lagrange basis defined by (N_in+1) 1D interpolation point positions xi_In(0:N_In)
! to another 3D tensor product node positions (number of nodes N_out+1)
! defined by (N_out+1) interpolation point  positions xi_Out(0:N_Out)
!  xi is defined in the 1DrefElem xi=[-1,1]
!==================================================================================================================================
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: Dim1,N_In,N_Out
REAL,INTENT(IN)     :: X3D_In(1:Dim1,0:N_In,0:N_In,0:N_In)
REAL,INTENT(IN)     :: Vdm_xi(0:N_Out,0:N_In)
REAL,INTENT(IN)     :: Vdm_eta(0:N_Out,0:N_In)
REAL,INTENT(IN)     :: Vdm_zeta(0:N_Out,0:N_In)
REAL,INTENT(OUT)    :: X3D_Out(1:Dim1,0:N_Out,0:N_Out,0:N_Out)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iN_In,jN_In,kN_In,iN_Out,jN_Out,kN_Out
REAL                :: X3D_Buf1(1:Dim1,0:N_Out,0:N_In,0:N_In)  ! first intermediate results from 1D interpolations
REAL                :: X3D_Buf2(1:Dim1,0:N_Out,0:N_Out,0:N_In) ! second intermediate results from 1D interpolations
!==================================================================================================================================
X3D_buf1=0.
! first direction iN_In
DO kN_In=0,N_In
  DO jN_In=0,N_In
    DO iN_In=0,N_In
      DO iN_Out=0,N_Out
        X3D_Buf1(:,iN_Out,jN_In,kN_In)=X3D_Buf1(:,iN_Out,jN_In,kN_In)+Vdm_xi(iN_Out,iN_In)*X3D_In(:,iN_In,jN_In,kN_In)
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
        X3D_Buf2(:,iN_Out,jN_Out,kN_In)=X3D_Buf2(:,iN_Out,jN_Out,kN_In)+Vdm_eta(jN_Out,jN_In)*X3D_Buf1(:,iN_Out,jN_In,kN_In)
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
        X3D_Out(:,iN_Out,jN_Out,kN_Out)=X3D_Out(:,iN_Out,jN_Out,kN_Out)+Vdm_zeta(kN_Out,kN_In)*X3D_Buf2(:,iN_Out,jN_Out,kN_In)
      END DO
    END DO
  END DO
END DO
END SUBROUTINE ChangeBasis3D_XYZ

END MODULE MOD_ChangeBasis
