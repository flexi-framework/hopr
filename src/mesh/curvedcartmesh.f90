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
! Copyright (C) 2017 Claus-Dieter Munz <munz@iag.uni-stuttgart.de>
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
MODULE MOD_CurvedCartMesh
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE CurvedCartesianMesh
  MODULE PROCEDURE CurvedCartesianMesh
END INTERFACE

INTERFACE GetNewCurvedHexahedron
  MODULE PROCEDURE GetNewCurvedHexahedron
END INTERFACE

PUBLIC::CurvedCartesianMesh
PUBLIC::GetNewCurvedHexahedron
!===================================================================================================================================

CONTAINS


SUBROUTINE GetNewCurvedHexahedron(CurvedNode,Ngeo,Zone)
!===================================================================================================================================
! Build new hexahedron for cartesian mesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:getNewElem
USE MOD_Mesh_Basis,ONLY:CreateSides
USE MOD_Basis_Vars,ONLY:HexaMap
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)              :: NGeo  ! ?
INTEGER,INTENT(IN)              :: Zone  ! ?
TYPE(tNodePtr),INTENT(IN)       :: CurvedNode(0:Ngeo,0:NGeo,0:Ngeo)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
TYPE(tElem),POINTER             :: aElem   ! ?
INTEGER                         :: i  ! ?
!===================================================================================================================================
CALL getNewElem(aElem)
aElem%zone=Zone
aElem%nNodes=8
ALLOCATE(aElem%Node(aElem%nNodes))
aElem%Node(1)%NP=>CurvedNode(   0,   0,   0)%NP
aElem%Node(2)%NP=>CurvedNode(NGeo,   0,   0)%NP
aElem%Node(3)%NP=>CurvedNode(NGeo,NGeo,   0)%NP
aElem%Node(4)%NP=>CurvedNode(   0,NGeo,   0)%NP
aElem%Node(5)%NP=>CurvedNode(   0,   0,NGeo)%NP
aElem%Node(6)%NP=>CurvedNode(NGeo,   0,NGeo)%NP
aElem%Node(7)%NP=>CurvedNode(NGeo,NGeo,NGeo)%NP
aElem%Node(8)%NP=>CurvedNode(   0,NGeo,NGeo)%NP
CALL CreateSides(aElem,.TRUE.)
IF(Ngeo.GT.1)THEN !curved
  aElem%nCurvedNodes=(Ngeo+1)**3
  ALLOCATE(aElem%CurvedNode(aElem%nCurvedNodes))
  DO i=1,aElem%nCurvedNodes
    aElem%CurvedNode(i)%NP=>CurvedNode(HexaMap(i,1),HexaMap(i,2),HexaMap(i,3))%NP
  END DO
END IF !curved

! Add elements to list
IF(.NOT.ASSOCIATED(FirstElem))THEN
  FirstElem=>aElem
ELSE
  aElem%nextElem          => FirstElem
  aElem%nextElem%prevElem => aElem
  FirstElem          => aElem
END IF
NULLIFY(aElem)
END SUBROUTINE GetNewCurvedHexahedron


SUBROUTINE CurvedCartesianMesh()
!===================================================================================================================================
! Builds cartesian mesh. Called by fillMesh.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tNodePtr,tElem,tSide
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:BoundaryType
USE MOD_Mesh_Vars,ONLY:getNewSide,getNewNode,getNewBC
USE MOD_Mesh_Vars,ONLY:deleteSide,deleteNode
USE MOD_Mesh_Vars,ONLY:BoundaryOrder,nElems,BCIndex
USE MOD_Mesh_Vars,ONLY:stretchType,fac,DXMaxToDXMin
USE MOD_IO_HDF5  ,ONLY:Elem_IJK,nElems_IJK
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
TYPE(tNodePtr)               :: CurvedNode(0:BoundaryOrder-1,0:BoundaryOrder-1,0:BoundaryOrder-1)   ! ?
TYPE(tNodePtr),POINTER       :: Mnodes(:,:,:)   ! ?
TYPE(tElem),POINTER          :: aElem   ! ?
TYPE(tSide),POINTER          :: aSide  ! ?
REAL,ALLOCATABLE             :: X_Ngeo(:,:,:,:,:,:,:)  ! ?
REAL                         :: snp(3),s,ds  ! ?
INTEGER                      :: iZone,iSide  ! ?
INTEGER                      :: i,j,k,ii,jj,kk,l,m,n,iDir  ! ?
INTEGER                      :: NodeInd  ! ?
INTEGER                      :: np(3)  ! ?
INTEGER                      :: Ngeo  ! ?
INTEGER                      :: elemID  ! ?
LOGICAL                      :: onBnd  ! ?
TYPE tDR
  REAL,ALLOCATABLE :: d(:)   ! ?
END TYPE
TYPE(tDR)  :: DR(3) ! grid spacing in r,s,t reference domain, used for stretching
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
WRITE(UNIT_stdOut,*)'BUILDING CURVED CARTESIAN MESH...'
CALL Timer(.TRUE.)
Ngeo=BoundaryOrder-1
NodeInd=0
iZone=1
! The low / up stuff is a remainder of the mpi version and could be replaced...
np(:) =nElems*Ngeo
snp=1./REAL(np)
nElems_IJK=nElems
ALLOCATE(Elem_IJK(PRODUCT(nElems),3))

ALLOCATE(Mnodes(0:np(1),0:np(2),0:np(3)))
DO l=0,np(1)
  DO m=0,np(2)
    DO n=0,np(3)
            CALL GetNewNode(Mnodes(l,m,n)%np)
            NodeInd=NodeInd+1
            Mnodes(l,m,n)%np%ind=NodeInd 
    END DO !n
  END DO !m
END DO !l
ALLOCATE(X_Ngeo(3,0:Ngeo,0:Ngeo,0:Ngeo,nElems(1),nElems(2),nElems(3)))


DO iDir=1,3
  ALLOCATE(DR(iDir)%d(nElems(iDir)))
  DR(iDir)%d=0. ! Initialize
    SELECT CASE(stretchType(iDir))
    CASE(0,6,7,30) !equidistant or for inner cell stretchings
      DR(iDir)%d=1.
    CASE(1) ! factor in one direction
      DR(iDir)%d(1)=1.
      DO ii=2,nElems(iDir)
        DR(iDir)%d(ii)=DR(iDir)%d(ii-1)*fac(iDir)
      END DO
    CASE(2) ! length ratio
      fac(iDir)=DxMaxToDxMin(iDir)**(1./(nElems(iDir)-1.))
      DR(iDir)%d(1)=1.
      DO ii=2,nElems(iDir)
        DR(iDir)%d(ii)=DR(iDir)%d(ii-1)*fac(iDir)
      END DO
    CASE(3) ! bell shaped exp(-s^2), s=[-1,1]
            ! DR= 1+ (dxmax/dxmin-2)*(exp(-(s*fac)^2)-exp(-(1*fac)^2))/(exp(-(0*fac)^2)-exp(-(1*fac)^2)) 
            ! fac>1, more stretching outwards, fac<1 less stretching  
      s=-1
      ds=2./(nElems(iDir)-1)
      DO ii=1,nElems(iDir)
        DR(iDir)%d(ii)=1.+(DxMaxToDxmin(iDir)-1.)*(exp(-(s*fac(iDir))**2) &
                       - exp(-(fac(iDir))**2))/(exp(0.) - exp(-(fac(iDir))**2))
        s=s+ds
      END DO
    CASE(4) ! half bell shaped exp(-s^2), s=[-1,0]
            ! DR= 1+ (dxmax/dxmin-2)*(exp(-(s*fac)^2)-exp(-(1*fac)^2))/(exp(-(0*fac)^2)-exp(-(1*fac)^2)) 
            ! fac>1, more stretching outwards, fac<1 less stretching  
      ds=(DxMaxToDxmin(iDir)-1.)/(1.-DxMaxToDxmin(iDir)*exp((-1+1./REAL(nElems(iDir)))*fac(iDir)))
      DO ii=1,nElems(iDir)
        DR(iDir)%d(ii)=1.+ds*exp((-1.+REAL(ii)/REAL(nElems(iDir)))*fac(iDir))
      END DO
    CASE(5) ! get custom node distribution from files
     ! CALL getNodeDistrib(DR(iDir)%d,iDir)
    END SELECT
    !normalize to reference space [-1,1]
    DR(iDir)%d=DR(iDir)%d*(2./SUM(DR(iDir)%d))
END DO  ! iDir

CALL BuildRefDomain(NGeo,nElems,X_NGeo,DR(1)%d,DR(2)%d,DR(3)%d)

DEALLOCATE(DR(1)%d,DR(2)%d,DR(3)%d)


! Build nodes
CALL BuildPhysicalDomain(NGeo,nElems,X_Ngeo)
!start position x2 in unit square grid
! calculate node postions, x1(:) coordinates in unit square grid
ElemID=PRODUCT(nElems)
DO kk=1,nElems(3)
  DO jj=1,nElems(2)
    DO ii=1,nElems(1)
      DO k=0,Ngeo
        DO j=0,Ngeo
          DO i=0,Ngeo
            l=i+(ii-1)*Ngeo
            m=j+(jj-1)*Ngeo
            n=k+(kk-1)*Ngeo
            Mnodes(l,m,n)%np%x(:)=X_Ngeo(:,i,j,k,ii,jj,kk) 
            CurvedNode(i,j,k)%np=>Mnodes(l,m,n)%np
          END DO !i
        END DO !j
      END DO !k
      !ONLY HEXA
      CALL GetNewCurvedHexahedron(CurvedNode,BoundaryOrder-1,iZone)
      !set ijk here
      Elem_IJK(elemID,:)=(/ii,jj,kk/)
      ElemID=ElemID-1
    END DO !ii
  END DO !jj
END DO !kk
!for Boundary Conditions
DO l=0,np(1)
  DO m=0,np(2)
    DO n=0,np(3)
      Mnodes(l,m,n)%np%tmp=0 
      IF(n.EQ.0) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+1 !zeta minus
      IF(m.EQ.0) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+20 !eta minus
      IF(l.EQ.np(1) ) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+300 !xi plus
      IF(m.EQ.np(2) ) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+4000  !eta plus
      IF(l.EQ.0) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+50000 !xi minus
      IF(n.EQ.np(3) ) Mnodes(l,m,n)%np%tmp=Mnodes(l,m,n)%np%tmp+600000 !zeta plus
    END DO !n
  END DO !m
END DO !l


aElem=>FirstElem
DO WHILE(ASSOCIATED(aElem))
  aSide=>aElem%firstSide
  DO WHILE(ASSOCIATED(aSide))
    DO iSide=1,6
      onBnd=.TRUE.
      DO i=1, aSide%nNodes
        IF((MOD(aSide%Node(i)%np%tmp,10**iSide)/(10**(iSide-1))).NE.iSide)THEN
          onBnd=.FALSE.
          EXIT ! Loop
        END IF
      END DO
      IF(BCIndex(iSide).EQ.0) onBnd=.FALSE.
      IF(onBnd)THEN
        CALL getNewBC(aSide%BC)
        aSide%BC%BCIndex    = BCIndex(iSide)
        aSide%BC%BCType     = BoundaryType(aSide%BC%BCIndex,1)
        aSide%curveIndex    = BoundaryType(aSide%BC%BCIndex,2)
        aSide%BC%BCstate    = BoundaryType(aSide%BC%BCIndex,3)
        aSide%BC%BCalphaInd = BoundaryType(aSide%BC%BCIndex,4)
      END IF
    END DO  ! iSide=1,6
    aSide=>aSide%nextElemSide
  END DO ! WHILE(ASSOCIATED(aSide))
  aElem=>aElem%nextElem
END DO ! WHILE(ASSOCIATED(aElem))
DEALLOCATE(Mnodes)
CALL Timer(.FALSE.)
END SUBROUTINE CurvedCartesianMesh


SUBROUTINE BuildRefDomain(Nin,nElems,R_glob,DR1,DR2,DR3)
!===================================================================================================================================
! compute the coordinates (r,s,t) in the reference domain [-1,1]x[-1,1]x[-1,1] at the points Xi
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:fac,fac2,stretchType,DxMaxToDxMin
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: Nin,nElems(3)  ! ?
REAL,INTENT(IN)    :: DR1(nElems(1)),DR2(nElems(2)),DR3(nElems(3)) ! Ref Domain grid spacings
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: R_glob(3,0:Nin,0:Nin,0:Nin,nElems(1),nElems(2),nElems(3)) ! Ref Domain Coordinates r,s,t for all elems
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: ii,jj,kk  ! ?
INTEGER          :: i,j,k  ! ?
INTEGER          :: dd,ee  ! ?
REAL             :: r1,r2,r3      ! coordinates of the barycenter of the actual element in reference domain
REAL             :: tmpfac,Pi  ! ?
REAL             :: ef2,c1,c2  ! ?
REAL             :: xi(0:Nin)  ! ?
!===================================================================================================================================

!        x-------------x
!        |      |      |
!        |      |      |
!        -----(rst)-----  A
!        |      |      |  |ds
!  s     |      |      |  |
! A      x-------------x  v
! |      <-dr->
! O--> r
!
!regular
DO i=0,Nin
 xi(i)=2*REAL(i)/REAL(Nin)-1.
END DO

r3=-1.
DO kk=1,nElems(3)
  r2=-1.
  DO jj=1,nElems(2)
    r1=-1.
    DO ii=1,nElems(1)
      DO k=0,Nin
        DO j=0,Nin
          DO i=0,Nin
            R_glob(1,i,j,k,ii,jj,kk)=r1+0.5*(1.+Xi(i))*DR1(ii)
            R_glob(2,i,j,k,ii,jj,kk)=r2+0.5*(1.+Xi(j))*DR2(jj)
            R_glob(3,i,j,k,ii,jj,kk)=r3+0.5*(1.+Xi(k))*DR3(kk)
          END DO !i=0,N
        END DO !j=0,N
      END DO !k=0,N
      r1=r1+DR1(ii)
      END DO !ii=1,nElems(1)
    r2=r2+DR2(jj)
  END DO !jj=1,nElems(2)
  r3=r3+DR3(kk)
END DO !kk=1,nElems(3)

DO dd=1,2
  ee=1
  IF(dd.EQ.1)ee=2
  SELECT CASE(stretchType(dd))
  CASE(6) !stretching with factor
    IF(.NOT.((fac(dd).EQ.1.).AND.(fac2(dd).EQ.1))) THEN
      DO kk=1,nElems(3)
        DO jj=1,nElems(2)
          DO ii=1,nElems(1)
            ! loop over Chebyshev points
            DO k=0,Nin
              DO j=0,Nin
                DO i=0,Nin
                  r1 = R_glob(dd,i,j,k,ii,jj,kk) !coordinate in stretching direction
                  r2 = R_glob(ee,i,j,k,ii,jj,kk) !perpendicular coordinate
                  tmpfac=0.5*(fac2(dd)*(1+r2)+fac(dd)*(1-r2))
                  R_glob(dd,i,j,k,ii,jj,kk)=(2./(tmpfac-1.))*(tmpfac**(0.5*(r1+1.))-1.) -1.
                END DO !i=0:Nin
              END DO !j=0:Nin
            END DO !k=0:Nin
          END DO !ii=1,nElems(1)
        END DO !jj=1,nElems(2)
      END DO !kk=1,nElems(3)
    END IF
 
  CASE(7) ! sin  symmetric mesh
    Pi = 4.*ATAN(1.)
    IF(.NOT.((fac(dd).EQ.1.).AND.(fac2(dd).EQ.1))) THEN
      DO kk=1,nElems(3)
        DO jj=1,nElems(2)
          DO ii=1,nElems(1)
            ! loop over Chebyshev points
            DO k=0,Nin
              DO j=0,Nin
                DO i=0,Nin
                  r1 = R_glob(dd,i,j,k,ii,jj,kk) !coordinate in stretching direction
                  r2 = R_glob(ee,i,j,k,ii,jj,kk) !perpendicular coordinate
                  tmpfac=0.5*(fac2(dd)*(1+r2)+fac(dd)*(1-r2))
                  R_glob(dd,i,j,k,ii,jj,kk)= r1-2*(tmpfac-1.)/(tmpfac+1.)*SIN(Pi*r1)/(2.*Pi) 
                END DO !i=0:Nin
              END DO !j=0:Nin
            END DO !k=0:Nin
          END DO !ii=1,nElems(1)
        END DO !jj=1,nElems(2)
      END DO !kk=1,nElems(3)
    END IF
  CASE(30) ! bell shaped (like case(3), but with inner cellstretching), uses fac (strecngth of the bell shape) i
           ! and dxmaxtodxmin (size between cells in the middle to cells at the ends).
           ! its the normalized integral of the function in Case(3):
           ! r = \int_{-1}^s 1+(dxmaxtodxmin-1)/(1-e^(-f^2))*( e^(-s^2 f^2)-e^(-f^2) ) ds
           ! => r= x*[1-e^(-f^2)*(1+(dxmaxtodxmin-1)/(1-e^(-f^2)))] +(dxmaxtodxmin-1)/(1-e^(-f^2))*(\int_{-1}^s e^(-s^2 f^2)ds)
           ! the integral of the exponential function is the error function erf,  
           ! and the integral is the normalized, so that  r is in [-1,1]
    IF(.NOT.(DxMaxToDxMin(dd).EQ.1.)) THEN
      Pi = 4.*ATAN(1.)
      ef2= EXP(-fac(dd)**2.)
      c1=1.-ef2*(DxMaxToDxmin(dd)-1.)/(1.-ef2)
      c2=(DxMaxToDxMin(dd)-1.)/((1.-ef2)*(2.*fac(dd)))*SQRT(pi)
      DO kk=1,nElems(3)
        DO jj=1,nElems(2)
          DO ii=1,nElems(1)
            ! loop over Chebyshev points
            DO k=0,Nin
              DO j=0,Nin
                DO i=0,Nin
                  r1 = R_glob(dd,i,j,k,ii,jj,kk) !coordinate in stretching direction
                  R_glob(dd,i,j,k,ii,jj,kk)=((r1+1.)*c1+c2*(ERF(r1*fac(dd))-ERF(-fac(dd))))/(c1+c2*ERF(fac(dd)))-1. 
                END DO !i=0:Nin
              END DO !j=0:Nin
            END DO !k=0:Nin
          END DO !ii=1,nElems(1)
        END DO !jj=1,nElems(2)
      END DO !kk=1,nElems(3)
    END IF
  END SELECT !CASE(stretchType)
END DO !dd=1,2

END SUBROUTINE BuildRefDomain


SUBROUTINE BuildPhysicalDomain(Nin,nElems,X)
!===================================================================================================================================
! compute the coordinates (X,Y,Z) in the physical space by evaluating at r,s,t positions of reference points
!===================================================================================================================================
USE MOD_Mesh_Vars,ONLY:CurvedMeshType,WhichMapping,X0,DX,XP,fac
! MODULES
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
INTEGER,INTENT(IN) :: Nin,nElems(3) 
! We use the CGNS notation of the Hexeader: http://www.grc.nasa.gov/WWW/cgns/sids/conv.html#unst_hexa
                        !      P8------P7
                        !     /|      /|
                        !    P5------P6|
                        !    | P4----|-P3
                        !    |/      |/
                        !    P1------P2  
! Faces
                        !F1              P1,P4,P3,P2      
                        !F2              P1,P2,P6,P5     
                        !F3              P2,P3,P7,P6    
                        !F4              P3,P4,P8,P7   
                        !F5              P1,P5,P8,P4  
                        !F6              P5,P6,P7,P8
! Edges
                        !E1              P1,P2           
                        !E2              P2,P3          
                        !E3              P3,P4         
                        !E4              P4,P1        
                        !E5              P1,P5       
                        !E6              P2,P6      
                        !E7              P3,P7     
                        !E8              P4,P8    
                        !E9              P5,P6    
                        !E10             P6,P7    
                        !E11             P7,P8  
                        !E12             P8,P5  
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: X(3,0:Nin,0:Nin,0:Nin,nElems(1),nElems(2),nElems(3)) ! IN: contains the Ref Domain coords r,s,t
                                                                           ! OUT: contains XYZ position of the ref points
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER          :: ii,jj,kk  ! ?
INTEGER          :: i,j,k  ! ?
REAL             :: DXh(3),R_loc(3)  ! ?
REAL             :: EP(3,12),FP(3,6)  ! ?
!===================================================================================================================================
SELECT CASE(CurvedmeshType)
CASE(1) ! cartesian domain X0,Y0,Z0,DX,DY,DZ
  DXh(:)=0.5*DX(:)
  DO kk=1,nElems(3)
    DO jj=1,nElems(2)
      DO ii=1,nElems(1)
        ! loop over Chebyshev points
        DO k=0,Nin
          DO j=0,Nin
            DO i=0,Nin
              ! overwrite r,s,t with x,y,z
              X(:,i,j,k,ii,jj,kk)=X0(:)+(1.+X(:,i,j,K,ii,jj,kk))*DXh(:) 
            END DO !i=0:Nin
          END DO !j=0:Nin
        END DO !k=0:Nin
      END DO !ii=1,nElems(1)
    END DO !jj=1,nElems(2)
  END DO !kk=1,nElems(3)
CASE(2) !trilinear domain
  DO kk=1,nElems(3)
    DO jj=1,nElems(2)
      DO ii=1,nElems(1)
        ! loop over Chebyshev points
        DO k=0,Nin
          DO j=0,Nin
            DO i=0,Nin
              ! store the local coordinates r,s,t
              R_loc = X(:,i,j,k,ii,jj,kk)
              ! overwrite r,s,t with x,y,z
              X(:,i,j,k,ii,jj,kk)=0.125*( XP(:,1)*(1.-R_loc(1)) &
                                                                 *(1.-R_loc(2)) &
                                                                 *(1.-R_loc(3)) &
                                                         +XP(:,2)*(1.+R_loc(1)) &
                                                                 *(1.-R_loc(2)) &
                                                                 *(1.-R_loc(3)) &
                                                         +XP(:,3)*(1.+R_loc(1)) &
                                                                 *(1.+R_loc(2)) &
                                                                 *(1.-R_loc(3)) &
                                                         +XP(:,4)*(1.-R_loc(1)) &
                                                                 *(1.+R_loc(2)) &
                                                                 *(1.-R_loc(3)) &
                                                         +XP(:,5)*(1.-R_loc(1)) &
                                                                 *(1.-R_loc(2)) &
                                                                 *(1.+R_loc(3)) &
                                                         +XP(:,6)*(1.+R_loc(1)) &
                                                                 *(1.-R_loc(2)) &
                                                                 *(1.+R_loc(3)) &
                                                         +XP(:,7)*(1.+R_loc(1)) &
                                                                 *(1.+R_loc(2)) &
                                                                 *(1.+R_loc(3)) &
                                                         +XP(:,8)*(1.-R_loc(1)) &
                                                                 *(1.+R_loc(2)) &
                                                                 *(1.+R_loc(3)) )
            END DO !i=0:Nin                    
          END DO !j=0:Nin
        END DO !k=0:Nin
      END DO !ii=1,nElems(1)
    END DO !jj=1,nElems(2)
  END DO !kk=1,nElems(3)
CASE(3) ! general curved domain
  ! we have to compute the corner points using the mapping
  CALL EvaluateFace(XP(:,1),(/-1.,-1./),DX,XP,WhichMapping,1)
  CALL EvaluateFace(XP(:,2),(/+1.,-1./),DX,XP,WhichMapping,1)
  CALL EvaluateFace(XP(:,3),(/+1.,+1./),DX,XP,WhichMapping,1)
  CALL EvaluateFace(XP(:,4),(/-1.,+1./),DX,XP,WhichMapping,1)
  CALL EvaluateFace(XP(:,5),(/-1.,-1./),DX,XP,WhichMapping,6)
  CALL EvaluateFace(XP(:,6),(/+1.,-1./),DX,XP,WhichMapping,6)
  CALL EvaluateFace(XP(:,7),(/+1.,+1./),DX,XP,WhichMapping,6)
  CALL EvaluateFace(XP(:,8),(/-1.,+1./),DX,XP,WhichMapping,6)
  DO kk=1,nElems(3)
    DO jj=1,nElems(2)
      DO ii=1,nElems(1)
        ! loop over Chebyshev points
        DO k=0,Nin
          DO j=0,Nin
            DO i=0,Nin
              ! store the local coordinates r,s,t
              R_loc = X(:,i,j,k,ii,jj,kk)
              ! overwrite r,s,t with x,y,z
              ! As we have to subtract the contribution of the edges, the contribution of the 
              ! corner nodes is subtracted two times.... we just add another contribution of 
              ! corner points to compensate
              X(:,i,j,k,ii,jj,kk)=0.125*( XP(:,1)*(1.-R_loc(1)) &
                                                                 *(1.-R_loc(2)) &
                                                                 *(1.-R_loc(3)) &
                                                         +XP(:,2)*(1.+R_loc(1)) &
                                                                 *(1.-R_loc(2)) &
                                                                 *(1.-R_loc(3)) &
                                                         +XP(:,3)*(1.+R_loc(1)) &
                                                                 *(1.+R_loc(2)) &
                                                                 *(1.-R_loc(3)) &
                                                         +XP(:,4)*(1.-R_loc(1)) &
                                                                 *(1.+R_loc(2)) &
                                                                 *(1.-R_loc(3)) &
                                                         +XP(:,5)*(1.-R_loc(1)) &
                                                                 *(1.-R_loc(2)) &
                                                                 *(1.+R_loc(3)) &
                                                         +XP(:,6)*(1.+R_loc(1)) &
                                                                 *(1.-R_loc(2)) &
                                                                 *(1.+R_loc(3)) &
                                                         +XP(:,7)*(1.+R_loc(1)) &
                                                                 *(1.+R_loc(2)) &
                                                                 *(1.+R_loc(3)) &
                                                         +XP(:,8)*(1.-R_loc(1)) &
                                                                 *(1.+R_loc(2)) &
                                                                 *(1.+R_loc(3)) )
            ! Determine the coordinates of the projection of the ref Point on the 6 faces FP
            CALL EvaluateFace(FP(:,1),(/R_loc(1),R_loc(2)/),DX,XP,WhichMapping,1)
            CALL EvaluateFace(FP(:,2),(/R_loc(1),R_loc(3)/),DX,XP,WhichMapping,2)
            CALL EvaluateFace(FP(:,3),(/R_loc(2),R_loc(3)/),DX,XP,WhichMapping,3)
            CALL EvaluateFace(FP(:,4),(/R_loc(1),R_loc(3)/),DX,XP,WhichMapping,4)
            CALL EvaluateFace(FP(:,5),(/R_loc(2),R_loc(3)/),DX,XP,WhichMapping,5)
            CALL EvaluateFace(FP(:,6),(/R_loc(1),R_loc(2)/),DX,XP,WhichMapping,6)

            ! Add the contribution of the Face Points
            X(:,i,j,k,ii,jj,kk)=X(:,i,j,k,ii,jj,kk)&
                                                   +0.5*(   FP(:,1)*(1.-R_loc(3)) &
                                                           +FP(:,2)*(1.-R_loc(2)) &
                                                           +FP(:,3)*(1.+R_loc(1)) &
                                                           +FP(:,4)*(1.+R_loc(2)) &
                                                           +FP(:,5)*(1.-R_loc(1)) &
                                                           +FP(:,6)*(1.+R_loc(3))) 
            ! Determine the coordinates of the projection of the ref Point on the 12 Edges EP
            CALL EvaluateFace(EP(:, 1),(/R_loc(1),-1./),DX,XP,WhichMapping,1)
            CALL EvaluateFace(EP(:, 2),(/+1.,R_loc(2)/),DX,XP,WhichMapping,1)
            CALL EvaluateFace(EP(:, 3),(/R_loc(1),+1./),DX,XP,WhichMapping,1)
            CALL EvaluateFace(EP(:, 4),(/-1.,R_loc(2)/),DX,XP,WhichMapping,1)
            CALL EvaluateFace(EP(:, 5),(/-1.,R_loc(3)/),DX,XP,WhichMapping,2)
            CALL EvaluateFace(EP(:, 6),(/+1.,R_loc(3)/),DX,XP,WhichMapping,2)
            CALL EvaluateFace(EP(:, 7),(/+1.,R_loc(3)/),DX,XP,WhichMapping,4)
            CALL EvaluateFace(EP(:, 8),(/-1.,R_loc(3)/),DX,XP,WhichMapping,4)
            CALL EvaluateFace(EP(:, 9),(/R_loc(1),-1./),DX,XP,WhichMapping,6)
            CALL EvaluateFace(EP(:,10),(/+1.,R_loc(2)/),DX,XP,WhichMapping,6)
            CALL EvaluateFace(EP(:,11),(/R_loc(1),+1./),DX,XP,WhichMapping,6)
            CALL EvaluateFace(EP(:,12),(/-1.,R_loc(2)/),DX,XP,WhichMapping,6)

            ! Subtract the constribution of the Edge Points
            X(:,i,j,k,ii,jj,kk)=X(:,i,j,k,ii,jj,kk) &
                                                  -0.25*(   EP(:,1) *(1.-R_loc(2)) &
                                                                    *(1.-R_loc(3)) &
                                                           +EP(:,2) *(1.+R_loc(1)) &
                                                                    *(1.-R_loc(3)) &
                                                           +EP(:,3) *(1.+R_loc(2)) &
                                                                    *(1.-R_loc(3)) &
                                                           +EP(:,4) *(1.-R_loc(1)) &
                                                                    *(1.-R_loc(3)) &
                                                           +EP(:,5) *(1.-R_loc(1)) &
                                                                    *(1.-R_loc(2)) &
                                                           +EP(:,6) *(1.+R_loc(1)) &
                                                                    *(1.-R_loc(2)) &
                                                           +EP(:,7) *(1.+R_loc(1)) &
                                                                    *(1.+R_loc(2)) &
                                                           +EP(:,8) *(1.-R_loc(1)) &
                                                                    *(1.+R_loc(2)) &
                                                           +EP(:,9) *(1.-R_loc(2)) &
                                                                    *(1.+R_loc(3)) &
                                                           +EP(:,10)*(1.+R_loc(1)) &
                                                                    *(1.+R_loc(3)) &
                                                           +EP(:,11)*(1.+R_loc(2)) &
                                                                    *(1.+R_loc(3)) &
                                                           +EP(:,12)*(1.-R_loc(1)) &
                                                                    *(1.+R_loc(3))) 
            END DO !i=0:Nin                    
          END DO !j=0:Nin
        END DO !k=0:Nin
      END DO !ii=1,nElems(1)
    END DO !jj=1,nElems(2)
  END DO !kk=1,nElems(3)
CASE(11) ! use a global stretching function in all directions, input parameter fac=(1.,1.,1.) is no stretching
  IF(fac(1).GT.1.) THEN
    DXh(1)=DX(1)/(fac(1)-1.)
    DO kk=1,nElems(3); DO jj=1,nElems(2); DO ii=1,nElems(1)
      ! loop over Chebyshev points
      DO k=0,Nin;DO j=0,Nin;DO i=0,Nin
        ! overwrite r,s,t with x,y,z
        X(1,i,j,k,ii,jj,kk)=X0(1)+(fac(1)**(0.5*(1.+X(1,i,j,K,ii,jj,kk)))-1.)*DXh(1) 
      END DO; END DO; END DO 
    END DO; END DO; END DO 
  ELSE
    DXh(1)=0.5*DX(1)
    DO kk=1,nElems(3); DO jj=1,nElems(2); DO ii=1,nElems(1)
      ! loop over Chebyshev points
      DO k=0,Nin;DO j=0,Nin;DO i=0,Nin
        ! overwrite r,s,t with x,y,z
        X(1,i,j,k,ii,jj,kk)=X0(1)+(1.+X(1,i,j,K,ii,jj,kk))*DXh(1) 
      END DO; END DO; END DO 
    END DO; END DO; END DO 
  END IF
  IF(fac(2).GT.1.) THEN
    DXh(2)=DX(2)/(fac(2)-1.)
    DO kk=1,nElems(3); DO jj=1,nElems(2); DO ii=1,nElems(1)
      ! loop over Chebyshev points
      DO k=0,Nin;DO j=0,Nin;DO i=0,Nin
        ! overwrite r,s,t with x,y,z
        X(2,i,j,k,ii,jj,kk)=X0(2)+(fac(2)**(0.5*(1.+X(2,i,j,K,ii,jj,kk)))-1.)*DXh(2) 
      END DO; END DO; END DO 
    END DO; END DO; END DO 
  ELSE
    DXh(2)=0.5*DX(2)
    DO kk=1,nElems(3); DO jj=1,nElems(2); DO ii=1,nElems(1)
      ! loop over Chebyshev points
      DO k=0,Nin;DO j=0,Nin;DO i=0,Nin
        ! overwrite r,s,t with x,y,z
        X(2,i,j,k,ii,jj,kk)=X0(2)+(1.+X(2,i,j,K,ii,jj,kk))*DXh(2) 
      END DO; END DO; END DO 
    END DO; END DO; END DO 
  END IF
  IF(fac(3).GT.1.) THEN
    DXh(3)=DX(3)/(fac(3)-1.)
    DO kk=1,nElems(3); DO jj=1,nElems(2); DO ii=1,nElems(1)
      ! loop over Chebyshev points
      DO k=0,Nin;DO j=0,Nin;DO i=0,Nin
        ! overwrite r,s,t with x,y,z
        X(3,i,j,k,ii,jj,kk)=X0(3)+(fac(3)**(0.5*(1.+X(3,i,j,K,ii,jj,kk)))-1.)*DXh(3) 
      END DO; END DO; END DO 
    END DO; END DO; END DO 
  ELSE
    DXh(3)=0.5*DX(3)
    DO kk=1,nElems(3); DO jj=1,nElems(2); DO ii=1,nElems(1)
      ! loop over Chebyshev points
      DO k=0,Nin;DO j=0,Nin;DO i=0,Nin
        ! overwrite r,s,t with x,y,z
        X(3,i,j,k,ii,jj,kk)=X0(3)+(1.+X(3,i,j,K,ii,jj,kk))*DXh(3) 
      END DO; END DO; END DO 
    END DO; END DO; END DO 
  END IF
END SELECT
END SUBROUTINE BuildPhysicalDomain

SUBROUTINE EvaluateFace(X,Para,DX,XP,WhichMapping,WhichFace)
!===================================================================================================================================
! Determines the value of the Face Function with number whichFaceFunction of the volume mapping with number whichMapping
! gets 2D reference parameter Para (-1...1) and gives 3D physical value X
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:R_0,R_INF,DY,PHI
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: Para(2),XP(3,8),DX(3)  ! ?
INTEGER,INTENT(IN) :: WhichMapping,WhichFace  ! ?
REAL,INTENT(OUT)   :: X(3)  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3)   :: e1,e2,e3        ! kartesian unit vectors
REAL,DIMENSION(3)   :: Gamma1,Gamma2,Gamma3,Gamma4     ! values of the 4 1D curves
REAL,DIMENSION(3)   :: Corner1,Corner2,Corner3,Corner4 ! values of the 4 corner points of the 2D face
REAL                :: PI,Z_Value,Fac,Para2(2)  ! ?
!===================================================================================================================================
SELECT CASE(WhichMapping)
CASE(1) ! TriLinearMapping 
  SELECT CASE(WhichFace)
  CASE(1) ! Face 1, Zeta=-1, Para = Xi,Eta
    X = 0.25*( XP(:,1)*(1.-Para(1))*(1.-Para(2))& 
              +XP(:,2)*(1.+Para(1))*(1.-Para(2))& 
              +XP(:,3)*(1.+Para(1))*(1.+Para(2))& 
              +XP(:,4)*(1.-Para(1))*(1.+Para(2))) 
  CASE(2) ! Face 2, Eta=-1, Para = Xi,Zeta
    X = 0.25*( XP(:,1)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,2)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,5)*(1.-Para(1))*(1.+Para(2))&
              +XP(:,6)*(1.+Para(1))*(1.+Para(2)))
  CASE(3) ! Face 3, Xi=+1, Para = Eta, Zeta
    X = 0.25*( XP(:,2)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,3)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,6)*(1.-Para(1))*(1.+Para(2))&
              +XP(:,7)*(1.+Para(1))*(1.+Para(2)))

  CASE(4) ! Face 4, Eta=+1, Para = Xi, Zeta
    X = 0.25*( XP(:,3)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,4)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,7)*(1.+Para(1))*(1.+Para(2))&
              +XP(:,8)*(1.-Para(1))*(1.+Para(2)))
  CASE(5) ! Face 5, Xi=-1, Para = Eta, Zeta
    X = 0.25*( XP(:,1)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,4)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,5)*(1.-Para(1))*(1.+Para(2))&
              +XP(:,8)*(1.+Para(1))*(1.+Para(2)))
  CASE(6) ! Face 6, Zeta=+1, Xi, Eta
    X = 0.25*( XP(:,5)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,6)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,7)*(1.+Para(1))*(1.+Para(2))&
              +XP(:,8)*(1.-Para(1))*(1.+Para(2)))
  END SELECT
CASE(2) ! Vulgo No Tron
  SELECT CASE(WhichFace)
  CASE(1) ! Face 1, Zeta=-1, Para = Xi,Eta
    X = 0.25*( XP(:,1)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,2)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,3)*(1.+Para(1))*(1.+Para(2))&
              +XP(:,4)*(1.-Para(1))*(1.+Para(2)))&
              +(/0.,0.,DX(3)/)*(1.-Para(1)**2)*(1.-Para(2)**2) 
  CASE(2) ! Face 2, Eta=-1, Para = Xi,Zeta
    X = 0.25*( XP(:,1)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,2)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,5)*(1.-Para(1))*(1.+Para(2))&
              +XP(:,6)*(1.+Para(1))*(1.+Para(2)))&
              +(/0.,DX(2),0./)*(1.-Para(1)**2)*(1.-Para(2)**2) 

  CASE(3) ! Face 3, Xi=+1, Para = Eta, Zeta
    X = 0.25*( XP(:,2)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,3)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,6)*(1.-Para(1))*(1.+Para(2))&
              +XP(:,7)*(1.+Para(1))*(1.+Para(2)))&
              +(/DX(1),0.,0./)*(1.-Para(1)**2)*(1.-Para(2)**2) 

  CASE(4) ! Face 4, Eta=+1, Para = Xi, Zeta
    X = 0.25*( XP(:,3)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,4)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,7)*(1.+Para(1))*(1.+Para(2))&
              +XP(:,8)*(1.-Para(1))*(1.+Para(2)))&
              +(/0.,DX(2),0./)*(1.-Para(1)**2)*(1.-Para(2)**2) 
  CASE(5) ! Face 5, Xi=-1, Para = Eta, Zeta
    X = 0.25*( XP(:,1)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,4)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,5)*(1.-Para(1))*(1.+Para(2))&
              +XP(:,8)*(1.+Para(1))*(1.+Para(2)))&
              +(/DX(1),0.,0./)*(1.-Para(1)**2)*(1.-Para(2)**2) 
  CASE(6) ! Face 6, Zeta=+1, Xi, Eta
    X = 0.25*( XP(:,5)*(1.-Para(1))*(1.-Para(2))&
              +XP(:,6)*(1.+Para(1))*(1.-Para(2))&
              +XP(:,7)*(1.+Para(1))*(1.+Para(2))&
              +XP(:,8)*(1.-Para(1))*(1.+Para(2)))&
              +(/0.,0.,DX(3)/)*(1.-Para(1)**2)*(1.-Para(2)**2) 
  END SELECT
CASE(3) ! half cylinder, Kopriva Book pg. 263 (we change y and z direction!! Changed to segment: PHI=180deg for half cylinder)
  PI=ACOS(-1.)
  e1=(/1.,0.,0./)
  e2=(/0.,1.,0./)
  e3=(/0.,0.,1./)
  Fac=PHI/180. !scale azimuth of half cylinder
  Para2=Para !respective component might be scaled with Fac
  SELECT CASE(WhichFace)
  CASE(1) ! Face 1, Zeta=-1, Para = Xi,Eta
    ! Scale azimuth of half cylinder
    Para2(2)= Fac*(Para(2)+1.)-1.
    ! Determine the values of the face and move it to y=-dy
    ! GammaI 
    Gamma1 = -r_0 * COS(PI*(Para2(2) + 1.)*0.5)*e1 + r_0*sin(PI*(Para2(2) + 1.)*0.5)*e2
    Gamma2 = (r_0 + (r_inf-r_0)*(Para2(1)+1.)*0.5)*e1
    Gamma3 = -r_inf * COS(PI*(Para2(2) + 1.)*0.5)*e1 + r_inf*sin(PI*(Para2(2) + 1.)*0.5)*e2
    Gamma4 = (-r_0 - (r_inf-r_0)*(Para2(1)+1.)*0.5)*e1
    ! Determine the 4 corner points of the face
    Corner1 = -r_0 * COS(PI*(-1. + 1.)*0.5)*e1 + r_0*sin(PI*(-1. + 1.)*0.5)*e2 != Gamma1 mit Para = -1.: Para2=-1
    Corner2 = (r_0 + (r_inf-r_0)*(-1.+1.)*0.5)*e1 != Gamma2 mit Para = -1.: Para2=-1
    Corner3 = -r_inf * COS(PI*((2.*Fac-1.) + 1.)*0.5)*e1 + r_inf*sin(PI*((2.*Fac-1.) + 1.)*0.5)*e2 != Gamma3 mit Para = +1.: Para2=2*Fac-1
    Corner4 = (-r_0 - (r_inf-r_0)*((2.*Fac-1.)+1.)*0.5)*e1 != Gamma4 mit Para = +1.: Para2=2*Fac-1
    ! Now use the classic "2D" coons mapping to determine the face, Kopriva Book pg. 230 
    X = 0.5 * ( Gamma4*(1.-Para2(1)) + Gamma2*(1.+Para2(1)) + Gamma1*(1.-Para2(2)) + Gamma3*(1.+Para2(2)))&
      - 0.25* ( (1.-Para2(1))*(Corner1*(1.-Para2(2))+Corner4*(1.+Para2(2)))+(1.+Para2(1))*(Corner2*(1.-Para2(2))+Corner3*(1.+Para2(2))))
    X(3) = -dy
  CASE(2) ! Face 2, Eta=-1, Para = Xi,Zeta
    ! Scale azimuth of half cylinder
    Para2(1)= Fac*(Para(1)+1.)-1.
    ! X = Gamma1 + Zeta*dy*e2
    X = -r_0 * COS(PI*(Para2(1) + 1.)*0.5)*e1 + r_0*sin(PI*(Para2(1) + 1.)*0.5)*e2 + e3*dy*Para2(2)
    !X = (r_0 + (r_inf-r_0)*(Para2(1)+1.)*0.5)*e1 + e3*dy*Para2(2)
  CASE(3) ! Face 3, Xi=+1, Para = Eta, Zeta
    ! Do not scale Para, only rotate entire face
    ! X = Gamma2 + Zeta*dy*e2
    !X = r_inf * COS(PI*(Para2(1) + 1.)*0.5)*e1 + r_inf*sin(PI*(Para2(1) + 1.)*0.5)*e2 + e3*dy*Para2(2)
    !X = (r_0 + (r_inf-r_0)*(Para2(1)+1.)*0.5)*e1 + e3*dy*Para2(2)
    X = -(r_0 + (r_inf-r_0)*(Para2(1)+1.)*0.5)*e1*COS(PHI*PI/180.) + (r_0 + (r_inf-r_0)*(Para2(1)+1.)*0.5)*e2*SIN(PHI*PI/180.) + e3*dy*Para2(2)
  CASE(4) !  Face 4, Eta=+1, Para = Xi, Zeta
    ! Scale azimuth of half cylinder
    Para2(1)= Fac*(Para(1)+1.)-1.
    ! Similar to face 1, expect we have now r_inf instead of r_inf
    ! X = Gamma3 + Zeta*dy*e2
    X = - r_inf * COS(PI*(Para2(1) + 1.)*0.5)*e1 + r_inf*sin(PI*(Para2(1) + 1.)*0.5)*e2 + e3*dy*Para2(2)
    !X = (r_0 + (r_inf-r_0)*(Para2(1)+1.)*0.5)*e1 + e3*dy*Para2(2)
  CASE(5) ! Face 5, Xi=-1, Para = Eta, Zeta
    ! No scaling of azimuth needed, since this face is always at same position
    ! For the full cylinder Face5 = Face3 !!
    ! X = Gamma4 + Zeta*dy*e2
    X = (-r_0 - (r_inf-r_0)*(Para2(1)+1.)*0.5)*e1 + e3*dy*Para2(2)
    !X = r_0 * COS(PI*(Para2(1) + 1.)*0.5)*e1 + r_0*sin(PI*(Para2(1) + 1.)*0.5)*e2 + e3*dy*Para2(2)
  CASE(6) !  Face 6, Zeta=+1, Xi, Eta
    ! The same as Face 2, except that we have now +dy instead of -dy in y direction!
    ! Scale azimuth of half cylinder
    Para2(2)= Fac*(Para(2)+1.)-1.
    ! Determine the values of the face and move it to y=+dy
    ! GammaI 
    Gamma1 = -r_0 * COS(PI*(Para2(2) + 1.)*0.5)*e1 + r_0*sin(PI*(Para2(2) + 1.)*0.5)*e2
    Gamma2 = (r_0 + (r_inf-r_0)*(Para2(1)+1.)*0.5)*e1
    Gamma3 = -r_inf * COS(PI*(Para2(2) + 1.)*0.5)*e1 + r_inf*sin(PI*(Para2(2) + 1.)*0.5)*e2
    Gamma4 = (-r_0 - (r_inf-r_0)*(Para2(1)+1.)*0.5)*e1
    ! Determine the 4 corner points of the face
    Corner1 = -r_0 * COS(PI*(-1. + 1.)*0.5)*e1 + r_0*sin(PI*(-1. + 1.)*0.5)*e2 != Gamma1 mit Para = -1.: Para2=-1
    Corner2 = (r_0 + (r_inf-r_0)*(-1.+1.)*0.5)*e1 != Gamma2 mit Para = -1.: Para2=-1
    Corner3 = -r_inf * COS(PI*((2.*Fac-1.) + 1.)*0.5)*e1 + r_inf*sin(PI*((2.*Fac-1.) + 1.)*0.5)*e2 != Gamma3 mit Para = +1.: Para2=2*Fac-1
    Corner4 = (-r_0 - (r_inf-r_0)*((2.*Fac-1.)+1.)*0.5)*e1 != Gamma4 mit Para = +1.: Para2=2*Fac-1
    ! Now use the classic "2D" coons mapping to determine the face, Kopriva Book pg. 230 
    X = 0.5 * ( Gamma4*(1.-Para2(1)) + Gamma2*(1.+Para2(1)) + Gamma1*(1.-Para2(2)) + Gamma3*(1.+Para2(2)))&
      - 0.25* ( (1.-Para2(1))*(Corner1*(1.-Para2(2))+Corner4*(1.+Para2(2)))+(1.+Para2(1))*(Corner2*(1.-Para2(2))+Corner3*(1.+Para2(2))))
    X(3)=dy
  END SELECT !whichface

CASE(4) ! full cylinder, Kopriva Book pg. 263 (we change y and z direction!!)
  PI=2.*ACOS(-1.)
  e1=(/1.,0.,0./)
  e2=(/0.,1.,0./)
  e3=(/0.,0.,1./)
  SELECT CASE(WhichFace)
  CASE(1) ! Face 1, Zeta=-1, Para = Xi,Eta
    ! Determine the values of the face and move it to y=-dy
    ! GammaI 
    Gamma1 = +r_0 * COS(PI*(Para(2) + 1.)*0.5)*e1 - r_0*sin(PI*(Para(2) + 1.)*0.5)*e2
    Gamma2 = (r_0 + (r_inf-r_0)*(Para(1)+1.)*0.5)*e1
    Gamma3 = +r_inf * COS(PI*(Para(2) + 1.)*0.5)*e1 - r_inf*sin(PI*(Para(2) + 1.)*0.5)*e2
    Gamma4 = (r_0 + (r_inf-r_0)*(Para(1)+1.)*0.5)*e1
    ! Determine the 4 corner points of the face
    Corner1 = (/+r_0,0.,0./) != Gamma1 mit Para(1) = -1. =XP1 without y component
    Corner2 = (/+r_0,0.,0./) != Gamma2 mit Para(2) = -1. = XP2 without y component
    Corner3 = (/+r_inf,0.,0./) != Gamma3 mit Para(1) = +1. = XP6 without y component
    Corner4 = (/+r_inf,0.,0./) != Gamma4 mit Para(2) = +1. = XP7 without y component
    ! Now use the classic "2D" coons mapping to determine the face, Kopriva Book pg. 230 
    X = 0.5 * ( Gamma4*(1.-Para(1)) + Gamma2*(1.+Para(1)) + Gamma1*(1.-Para(2)) + Gamma3*(1.+Para(2)))&
      - 0.25* ( (1.-Para(1))*(Corner1*(1.-Para(2))+Corner4*(1.+Para(2)))+(1.+Para(1))*(Corner2*(1.-Para(2))+Corner3*(1.+Para(2))))
    X(3) = -dy
  CASE(2) ! Face 2, Eta=-1, Para = Xi,Zeta
    ! X = Gamma1 + Zeta*dy*e2
    X = r_0 * COS(PI*(Para(1) + 1.)*0.5)*e1 - r_0*sin(PI*(Para(1) + 1.)*0.5)*e2 + e3*dy*Para(2)
    !X = (r_0 + (r_inf-r_0)*(Para(1)+1.)*0.5)*e1 + e3*dy*Para(2)
  CASE(3) ! Face 3, Xi=+1, Para = Eta, Zeta
    ! X = Gamma2 + Zeta*dy*e2
    !X = r_inf * COS(PI*(Para(1) + 1.)*0.5)*e1 + r_inf*sin(PI*(Para(1) + 1.)*0.5)*e2 + e3*dy*Para(2)
    X = (r_0 + (r_inf-r_0)*(Para(1)+1.)*0.5)*e1 + e3*dy*Para(2)
  CASE(4) !  Face 4, Eta=+1, Para = Xi, Zeta
    ! Similar to face 1, expect we have now r_inf instead of r_inf
    ! X = Gamma3 + Zeta*dy*e2
    X = r_inf * COS(PI*(Para(1) + 1.)*0.5)*e1 - r_inf*sin(PI*(Para(1) + 1.)*0.5)*e2 + e3*dy*Para(2)
    !X = (r_0 + (r_inf-r_0)*(Para(1)+1.)*0.5)*e1 + e3*dy*Para(2)
  CASE(5) ! Face 5, Xi=-1, Para = Eta, Zeta
    ! For the full cylinder Face5 = Face3 !!
    ! X = Gamma4 + Zeta*dy*e2
    X = (r_0 + (r_inf-r_0)*(Para(1)+1.)*0.5)*e1 + e3*dy*Para(2)
    !X = r_0 * COS(PI*(Para(1) + 1.)*0.5)*e1 + r_0*sin(PI*(Para(1) + 1.)*0.5)*e2 + e3*dy*Para(2)
  CASE(6) !  Face 6, Zeta=+1, Xi, Eta
    ! The same as Face 2, except that we have now +dy instead of -dy in y direction!
    ! Determine the values of the face and move it to y=+dy
    ! GammaI 
    Gamma1 = r_0 * COS(PI*(Para(2) + 1.)*0.5)*e1 - r_0*sin(PI*(Para(2) + 1.)*0.5)*e2
    Gamma2 = (r_0 + (r_inf-r_0)*(Para(1)+1.)*0.5)*e1
    Gamma3 = r_inf * COS(PI*(Para(2) + 1.)*0.5)*e1 - r_inf*sin(PI*(Para(2) + 1.)*0.5)*e2
    Gamma4 = (r_0 + (r_inf-r_0)*(Para(1)+1.)*0.5)*e1
    ! Determine the 4 corner points of the face
    Corner1 = (/r_0,0.,0./) != Gamma1 mit Para(1) = -1.
    Corner2 = (/r_0,0.,0./) != Gamma2 mit Para(2) = -1.
    Corner3 = (/r_inf,0.,0./) != Gamma3 mit Para(1) = +1.
    Corner4 = (/r_inf,0.,0./) != Gamma4 mit Para(2) = +1.
    ! Now use the classic "2D" coons mapping to determine the face, Kopriva Book pg. 230 
    X = 0.5 * ( Gamma4*(1.-Para(1)) + Gamma2*(1.+Para(1)) + Gamma1*(1.-Para(2)) + Gamma3*(1.+Para(2)))&
      - 0.25* ( (1.-Para(1))*(Corner1*(1.-Para(2))+Corner4*(1.+Para(2)))+(1.+Para(1))*(Corner2*(1.-Para(2))+Corner3*(1.+Para(2))))
    X(3)=dy
 
  END SELECT !whichface

CASE(5) !SINE BUMP
  PI=ACOS(-1.)
  e1=(/1.,0.,0./)
  e2=(/0.,1.,0./)
  e3=(/0.,0.,1./)
  SELECT CASE(WhichFace)
  CASE(1) ! Face 1, Zeta=-1, Para = Xi,Eta
    X(1:2)=DX(1:2)*Para(1:2)
    !X(3)= R_0*0.5*(1.+cos(MIN(ABS(X(1))/R_inf,1.)*PI))
    X(3)= R_0*(1. - MIN(ABS(DX(1)*Para(1))/R_inf,1.)**2)**6
  CASE(2) ! Face 2, Eta=-1, Para = Xi,Zeta
    !gamma1= DX(1)*Para(1)*e1 + R_0*0.5*(1.+cos(MIN(ABS(DX(1)*Para(1))/R_inf,1.)*PI))*e3
    Z_Value=R_0*(1. - MIN(DX(1)/R_inf,1.)**2)**6
    gamma1= DX(1)*Para(1)*e1 + R_0*(1. - MIN(ABS(DX(1)*Para(1))/R_inf,1.)**2)**6*e3
    gamma2= DX(1)*e1+ (Z_Value*(1.-Para(2))*0.5 + DX(3)*(1+Para(2)))*e3
    gamma3= DX(1)*Para(1)*e1 + 2.*DX(3)*e3 
    gamma4=-DX(1)*e1+ (Z_Value*(1.-Para(2))*0.5 + DX(3)*(1+Para(2)))*e3
    ! Determine the 4 corner points of the face
    Corner1 = (/-DX(1),0.,Z_Value/) != Gamma1 mit Para(1) = -1.
    Corner2 = (/+DX(1),0.,Z_Value/) != Gamma2 mit Para(2) = -1.
    Corner3 = (/+DX(1),0.,2.*DX(3)/) != Gamma3 mit Para(1) = +1.
    Corner4 = (/-DX(1),0.,2.*DX(3)/) != Gamma4 mit Para(2) = +1.
    ! Now use the classic "2D" coons mapping to determine the face, Kopriva Book pg. 230 
    X = 0.5 * ( Gamma4*(1.-Para(1)) + Gamma2*(1.+Para(1)) + Gamma1*(1.-Para(2)) + Gamma3*(1.+Para(2)))&
      - 0.25* ( (1.-Para(1))*(Corner1*(1.-Para(2))+Corner4*(1.+Para(2)))+(1.+Para(1))*(Corner2*(1.-Para(2))+Corner3*(1.+Para(2))))
    X(2)=DX(2)
  CASE(3) ! Face 3, Xi=+1, Para = Eta, Zeta
    Z_Value=R_0*(1. - MIN(DX(1)/R_inf,1.)**2)**6
    X(1)=DX(1)
    X(2)=DX(2)*Para(1) 
    X(3)=Z_Value*(1.-Para(2))*0.5+DX(3)*(1.+Para(2))
  CASE(4) ! Face 4, Eta=+1, Para = Xi, Zeta
    !gamma1= DX(1)*Para(1)*e1 + R_0*0.5*(1.+cos(MIN(ABS(DX(1)*Para(1))/R_inf,1.)*PI))*e3
    Z_Value=R_0*(1. - MIN(DX(1)/R_inf,1.)**2)**6
    gamma1= DX(1)*Para(1)*e1 + R_0*(1. - MIN(ABS(DX(1)*Para(1))/R_inf,1.)**2)**6*e3
    gamma2= DX(1)*e1+ (Z_Value*(1.-Para(2))*0.5 + DX(3)*(1+Para(2)))*e3
    gamma3= DX(1)*Para(1)*e1 + 2.*DX(3)*e3 
    gamma4=-DX(1)*e1+ (Z_Value*(1.-Para(2))*0.5 + DX(3)*(1+Para(2)))*e3
    ! Determine the 4 corner points of the face
    Corner1 = (/-DX(1),0.,Z_Value/) != Gamma1 mit Para(1) = -1.
    Corner2 = (/+DX(1),0.,Z_Value/) != Gamma2 mit Para(2) = -1.
    Corner3 = (/+DX(1),0.,2.*DX(3)/) != Gamma3 mit Para(1) = +1.
    Corner4 = (/-DX(1),0.,2.*DX(3)/) != Gamma4 mit Para(2) = +1.
    ! Now use the classic "2D" coons mapping to determine the face, Kopriva Book pg. 230 
    X = 0.5 * ( Gamma4*(1.-Para(1)) + Gamma2*(1.+Para(1)) + Gamma1*(1.-Para(2)) + Gamma3*(1.+Para(2)))&
      - 0.25* ( (1.-Para(1))*(Corner1*(1.-Para(2))+Corner4*(1.+Para(2)))+(1.+Para(1))*(Corner2*(1.-Para(2))+Corner3*(1.+Para(2))))
    X(2)=-DX(2)
  CASE(5) ! Face 5, Xi=-1, Para = Eta, Zeta
    Z_Value=R_0*(1. - MIN(DX(1)/R_inf,1.)**2)**6
    X(1)=-DX(1)
    X(2)=DX(2)*Para(1) 
    X(3)=Z_Value*(1.-Para(2))*0.5+DX(3)*(1.+Para(2))
  CASE(6) ! Face 6, Zeta=+1, Xi, Eta
    X(1:2)=DX(1:2)*Para(1:2)
    X(3)  =2.*DX(3)
  END SELECT !whichFace

CASE(6) !quadratic duct (a la lounge)
  PI=ACOS(-1.)
  e1=(/1.,0.,0./)
  e2=(/0.,1.,0./)
  e3=(/0.,0.,1./)
  SELECT CASE(WhichFace)
  CASE(1) ! Face 1, Zeta=-1, Para = Xi,Eta
    X(1) = (1.+Para(1))*0.5
    X(2)=0.25*Para(2)
    X(3)= X(1)**2
  CASE(2) ! Face 2, Eta=-1, Para = Xi,Zeta
    X(1)=(1.+Para(1))*0.5
    gamma1= X(1)*e1 + X(1)**2*e3
    gamma2=(/1.,0.,1./)*(1.-Para(2))*0.5 + (/2./3.,0.,7./6./)*(1.+Para(2))*0.5
    X(1)=(1.+Para(1))*0.5*2./3.
    gamma3=X(1)*e1+(0.5+1.5*X(1)**2)*e3
    gamma4=0.25*(1+Para(2))*e3
    ! Determine the 4 corner points of the face
    Corner1 = (/0.,0.,0./) != Gamma1 mit Para(1) = -1.
    Corner2 = (/1.,0.,1./) != Gamma2 mit Para(2) = -1.
    Corner3 = (/2./3.,0.,7./6./) != Gamma3 mit Para(1) = +1.
    Corner4 = (/0.,0.,0.5/) != Gamma4 mit Para(2) = +1.
    ! Now use the classic "2D" coons mapping to determine the face, Kopriva Book pg. 230 
    X = 0.5 * ( Gamma4*(1.-Para(1)) + Gamma2*(1.+Para(1)) + Gamma1*(1.-Para(2)) + Gamma3*(1.+Para(2)))&
      - 0.25* ( (1.-Para(1))*(Corner1*(1.-Para(2))+Corner4*(1.+Para(2)))+(1.+Para(1))*(Corner2*(1.-Para(2))+Corner3*(1.+Para(2))))
    X(2)=0.25
  CASE(3) ! Face 3, Xi=+1, Para = Eta, Zeta
    X = (/1.,0.,1./)*(1.-Para(2))*0.5 + (/2./3.,0.,7./6./)*(1.+Para(2))*0.5 +0.25*Para(1)*e2
  CASE(4) ! Face 4, Eta=+1, Para = Xi, Zeta
    X(1)=(1.+Para(1))*0.5
    gamma1= X(1)*e1 + X(1)**2*e3
    gamma2=(/1.,0.,1./)*(1.-Para(2))*0.5 + (/2./3.,0.,7./6./)*(1.+Para(2))*0.5
    X(1)=(1.+Para(1))*0.5*2./3.
    gamma3=X(1)*e1+(0.5+1.5*X(1)**2)*e3
    gamma4=0.25*(1+Para(2))*e3
    ! Determine the 4 corner points of the face
    Corner1 = (/0.,0.,0./) != Gamma1 mit Para(1) = -1.
    Corner2 = (/1.,0.,1./) != Gamma2 mit Para(2) = -1.
    Corner3 = (/2./3.,0.,7./6./) != Gamma3 mit Para(1) = +1.
    Corner4 = (/0.,0.,0.5/) != Gamma4 mit Para(2) = +1.
    ! Now use the classic "2D" coons mapping to determine the face, Kopriva Book pg. 230 
    X = 0.5 * ( Gamma4*(1.-Para(1)) + Gamma2*(1.+Para(1)) + Gamma1*(1.-Para(2)) + Gamma3*(1.+Para(2)))&
      - 0.25* ( (1.-Para(1))*(Corner1*(1.-Para(2))+Corner4*(1.+Para(2)))+(1.+Para(1))*(Corner2*(1.-Para(2))+Corner3*(1.+Para(2))))
    X(2)=-0.25
  CASE(5) ! Face 5, Xi=-1, Para = Eta, Zeta
    X = 0.25*(1+Para(2))*e3+0.25*Para(1)*e2
  CASE(6) ! Face 6, Zeta=+1, Xi, Eta
    X(1) = (1.+Para(1))*0.5*2./3.
    X(2)=0.25*Para(2)
    X(3)= 1.5*X(1)**2+0.5
  END SELECT !whichFace

END SELECT
END SUBROUTINE EvaluateFace


END MODULE MOD_CurvedCartMesh
