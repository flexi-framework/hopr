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
MODULE MOD_Readin_ICEM
!===================================================================================================================================
! ?
!===================================================================================================================================
! MODULES
USE MOD_Globals
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE readSpecEdges
  MODULE PROCEDURE readSpecEdges
END INTERFACE

PUBLIC::readSpecEdges
!===================================================================================================================================

CONTAINS
SUBROUTINE readSpecEdges()
!===================================================================================================================================
! This subroutine reads the ASCII File exported from ICEM containing the Chebychev-Lobatto points. These are only curved edges.
! The corner nodes of the corresponding edge are found by their node index, then the additional curved nodes are associated to the
! edge%curvednode array and are transformed to equidistant nodes.
!===================================================================================================================================
! MODULES
USE MOD_Mesh_Vars,ONLY:tElem,tEdge,tNode,tNodePtr
USE MOD_Mesh_Vars,ONLY:FirstElem
USE MOD_Mesh_Vars,ONLY:N
USE MOD_Mesh_Vars,ONLY:SpecElemFile
USE MOD_Mesh_Vars,ONLY:getNewNode
USE MOD_Basis1D,ONLY:ChebyGaussLobNodesAndWeights,BarycentricWeights,InitializeVandermonde
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,ALLOCATABLE          :: EdgeIndex(:,:)  ! ?
INTEGER                      :: iEdge,iCounter,maxInd,helpVarI  ! ?
INTEGER                      :: read_stat,stat_old  ! ?
INTEGER                      :: ii,nCurves  ! ?
REAL,ALLOCATABLE             :: EdgeNodes(:,:,:)  ! ?
REAL,ALLOCATABLE             :: helpVar(:,:),xi_Out(:),xCL(:),wBary(:),VdM(:,:)  ! ?
CHARACTER (LEN=60)           :: filename  ! ?
CHARACTER (LEN=200)          :: dummyline,helpstring  ! ?
LOGICAL                      :: exist,edgefound  ! ?
TYPE(tElem), POINTER         :: aElem  ! ?
TYPE(tNode), POINTER         :: aNode,bNode  ! ?
TYPE(tNodePtr),ALLOCATABLE   :: Nodes(:)  ! ?
TYPE(tEdge), POINTER         :: aEdge  ! ?
!===================================================================================================================================
WRITE(UNIT_stdOut,'(132("~"))')
CALL Timer(.TRUE.)
WRITE(UNIT_stdOut,'(A)')'Reading Chebychev-Lobatto points from specElemFile ',TRIM(specElemFile),'...'
filename=TRIM(specElemFile)
ALLOCATE(xi_Out(0:N))
ALLOCATE(xCL(0:N))
ALLOCATE(wBary(0:N))
ALLOCATE(VdM(0:N,0:N))
! equidistant nodes
DO ii=0,N
  xi_Out(ii)=-1.+ii*2./N
END DO
CALL ChebyGaussLobNodesAndWeights(N,xCL)  
CALL BarycentricWeights(N,xCL,wBary)
CALL InitializeVandermonde(N,N,wBary,xCL,xi_out,Vdm)
DEALLOCATE(xi_Out,xCL,wBary)

aElem=>firstElem      !get maxInd of nodes
maxInd=0
DO WHILE(ASSOCIATED(aElem))
  DO ii=1,aElem%nNodes
    maxInd=MAX(maxInd,aElem%node(ii)%np%ind)
  END DO
  aElem=>aElem%nextElem
ENDDO

ALLOCATE(Nodes(maxInd))    !create nodes array according to index
DO ii=1,maxInd
  NULLIFY(Nodes(ii)%np)
END DO

aElem=>firstElem
DO WHILE(ASSOCIATED(aElem))
  DO ii=1,aElem%nNodes
    Nodes(aElem%node(ii)%np%ind)%np=>aElem%node(ii)%np
  END DO
  aElem=>aElem%nextElem
END DO

INQUIRE (FILE = Filename, EXIST = exist)
IF (exist) THEN     
      OPEN(   UNIT    = 105           , &
              FILE    = Filename      , &
              FORM    = 'FORMATTED'   , &
              ACTION  = 'READ'        , &
              STATUS  = 'OLD'         , &
              IOSTAT  = stat_old)
ELSE !  Not yet Specified
  WRITE (*,*) '####### Looking for Filename:', Filename
    CALL abort(__STAMP__,&
              'Spetral element file with Chebychev-Lobatto Nodes not found.')
END IF

! Find the number of curved edges (lines in file)
read_stat=0
iCounter=0
DO WHILE ( read_stat == 0 )
  READ (105,*, IOSTAT=read_stat) dummyline
  IF (read_stat == 0 .OR. TRIM(dummyline)=='') THEN     !!!!!
    iCounter = iCounter + 1
  END IF
  IF(MOD(iCounter-1,N).EQ.0.AND.dummyline.NE.'Edge')THEN
    WRITE(UNIT_stdOut,*)TRIM(filename),' in line',iCounter,' is not conform to boundaryOrder or corrupted'
    STOP
  ENDIF
END DO

! WRITE (*,*) 'Rows in file (iCounter):',iCounter !Debugging
nCurves=iCounter/N                !Now nCurves gives the number of curved edges 
ALLOCATE(EdgeIndex(nCurves,2)) 
ALLOCATE(helpVar(N-1,3)) 
ALLOCATE(EdgeNodes(nCurves,0:N,3)) 

REWIND (105)
DO iEdge=1,nCurves                     !read curves into EdgeIndex and EdgeNodes
  READ(105,'(A60)') dummyline
  helpstring=dummyline((index(dummyline,'('))+1:)
  helpstring=helpstring(1:(index(helpstring,')'))-1)
  read(helpstring,*,IOSTAT=read_stat)EdgeIndex(iEdge,1:2)
  IF (read_stat.NE.0.OR.dummyline(1:17).NE.'Edge definition (') THEN
    WRITE(UNIT_stdOut,*)TRIM(filename),' in line',N*(iEdge-1)+1,'does not have the right format'
    STOP
  END IF
  DO ii=1,N-1
    READ(105,*,IOSTAT=read_stat) EdgeNodes(iEdge,ii,1:3)
    IF (read_stat.ne.0) THEN
      WRITE(UNIT_stdOut,*)TRIM(filename),' does not have the right format'
      STOP
    END IF
  END DO
  IF (EdgeIndex(iEdge,1).GT.EdgeIndex(iEdge,2)) THEN         !reverse nodes and indices if not oriented
    helpvar(1:N-1,1:3)=EdgeNodes(iEdge,1:N-1,1:3)
    helpvari=EdgeIndex(iEdge,1)
    EdgeIndex(iEdge,1)=EdgeIndex(iEdge,2)
    EdgeIndex(iEdge,2)=helpvarI 
    DO ii=1,(N-1)
      EdgeNodes(iEdge,ii,1:3)=helpvar(N-ii,1:3)
    END DO
  END IF
END DO
DEALLOCATE(helpvar)

iCounter=0
DO iEdge=1,nCurves     ! fill edge%curvedNode(:) and transform them to equidistant
  aNode=>Nodes(EdgeIndex(iEdge,1))%np
  aEdge=>aNode%firstEdge
  edgefound=.FALSE.
  DO WHILE(ASSOCIATED(aEdge))
    bNode=>aEdge%node(2)%np
    IF (bNode%ind.EQ.EdgeIndex(iEdge,2)) THEN
      edgefound=.TRUE.
      ALLOCATE(aEdge%curvedNode(N+1))
      aEdge%curvedNode(1)%np=>aNode
      aEdge%curvedNode(N+1)%np=>bNode
      EdgeNodes(iEdge,0,1:3)=aNode%x(1:3)
      EdgeNodes(iEdge,N,1:3)=bNode%x(1:3)
      EdgeNodes(iEdge,:,:)=MATMUL(VdM,EdgeNodes(iEdge,:,:)) !transform to equidist
      iCounter=iCounter+1
      DO ii=2,N
        maxInd=maxInd+1
        CALL getNewNode(aEdge%curvedNode(ii)%np)
        aEdge%curvedNode(ii)%np%ind=maxInd
        aEdge%curvedNode(ii)%np%x(1:3)=EdgeNodes(iEdge,ii-1,1:3)              
      END DO
      EXIT
    ELSE
      aEdge=>aEdge%nextEdge
    END IF
  END DO
  IF(.NOT.edgefound) THEN
    WRITE(UNIT_stdOut,*)'Found no edge in mesh with indices',EdgeIndex(iEdge,:)
  END IF
END DO
maxInd=SIZE(Nodes(:))
DO ii=1,maxInd
  NULLIFY(Nodes(ii)%np)
END DO
DEALLOCATE(Nodes,EdgeNodes,EdgeIndex,VdM)

IF (iCounter.NE.nCurves) THEN
  WRITE(UNIT_stdOut,*)'#curvedNodes read not equal to #curvedNodes assigned!!!'
  WRITE(UNIT_stdOut,*)'Number of curved edges read     =',nCurves
  WRITE(UNIT_stdOut,*)'Number of curved edges assigned =',iCounter
END IF
CALL Timer(.FALSE.)
END SUBROUTINE readSpecEdges

END MODULE MOD_Readin_ICEM
