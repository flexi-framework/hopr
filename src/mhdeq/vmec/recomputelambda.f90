
SUBROUTINE RecomputeLambda()
!===================================================================================================================================
! recompute lambda, from interpolated R,Z coordinates. Algorithm taken from track code 
! uses the Jacobian: R*(dR/dtheta*dZ/drho-dR/drho*dZ/dtheta)
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_VMEC_Mappings, ONLY: nFluxVMEC,mn_mode_nyq,xm_nyq,xn_nyq,mn_mode,xm,xn
USE MOD_VMEC_Vars,     ONLY: rho,Rmnc_Spl,Zmns_Spl
USE SPLINE1_MOD,       ONLY: SPLINE1_EVAL 
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: max_m,max_n   !maximum  mode numbers
INTEGER           :: iNode,jNode,iMode,Mode00        
INTEGER           :: iGuess,iFlux
INTEGER           :: ntheta,nZeta  !number of points in theta/zeta
REAL,ALLOCATABLE  :: wtheta(:)     !integration weights
REAL,ALLOCATABLE  :: wzeta(:)      !integration weights
REAL              :: dtheta,dzeta  !equidistant distance between points
REAL              :: theta,zeta  !poloidal/toroidal angle
REAL              :: R,dRdU,dRdrho 
REAL              :: Z,dZdU,dZdrho 
REAL              :: rhom,drhom,splout(3) 
REAL              :: COSMN(mn_mode),SINMN(mn_mode) 
REAL,ALLOCATABLE  :: gsqrt(:,:,:)  !jacobian
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A)')'  RECOMPUTE LAMBDA ...'
Mode00=-1
DO iMode=1,mn_mode_nyq
  IF((xm_nyq(iMode).EQ.0).AND.(xn_nyq(iMode).EQ.0))THEN
    Mode00=iMode
    EXIT
  END IF
END DO
IF(Mode00.EQ.-1) STOP 'Problem with Mode00'

max_m=MAXVAL(ABS(xm_nyq))
ntheta=4*max_m+1
dtheta=2*Pi/REAL(ntheta-1)
ALLOCATE(wtheta(1:ntheta))
wtheta(1)=1.0
wtheta(ntheta)=1.0
DO iNode=2,ntheta-1 !Over poloidal nodes
  IF(MOD(iNode,2) == 0) THEN
    !Even node
    wtheta(iNode)=4.0
  ELSE
    !Odd mode
    wtheta(iNode)=2.0
  END IF
END DO !Over poloidal nodes
wtheta(:)=wtheta(:)*(2.0*pi/(3.0*REAL(ntheta-1)))

max_n=MAXVAL(ABS(xn_nyq))
nzeta=4*max_n+1
dzeta=2*Pi/REAL(nzeta-1)
ALLOCATE(wzeta(1:nzeta))
wzeta(1)=1.0
wzeta(nzeta)=1.0
DO jNode=2,nzeta-1 !Over toroidal nodes
  IF(MOD(jNode,2) == 0) THEN
    !Even node
    wzeta(jNode)=4.0
  ELSE
    !Odd mode
    wzeta(jNode)=2.0
  END IF
END DO !Over toroidal nodes
wzeta(:)=wzeta(:)*(2.0*pi/(3.0*REAL(nzeta-1)))

ALLOCATE(gsqrt(2:nFluxVMEC,1:ntheta,1:nzeta))
DO jNode=1,nzeta
  zeta=(jNode-1)*dzeta
  DO iNode=1,ntheta
    theta=(iNode-1)*dtheta
    CosMN(:)  = COS(xm(:) * theta - xn(:) * zeta)
    SinMN(:)  = SIN(xm(:) * theta - xn(:) * zeta) 
    DO iFlux=2,nFluxVMEC
      iGuess=iFlux
      R      =0.
      Z      =0.
      dRdu   =0.
      dRdrho =0.
      dZdu   =0.
      dZdrho =0.
      DO iMode=1,mn_mode
        IF(xmabs_nyq(iMode).EQ.0)THEN
          rhom=1.
          drhom=0.
        ELSEIF(xmabs_nyq(iMode).EQ.1)THEN
          rhom=rho(iFlux)
          drhom=1.
        ELSE
          rhom=rho(iFlux)**xmabs_nyq(iMode)
          drhom=xmabs_nyq(iMode)*rho(iFlux)**(xmabs_nyq(iMode)-1)
        END IF
        CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho(iFlux),rho,Rmnc_Spl(:,:,iMode),iGuess,splout) 
        R    = R    + rhom*splout(1)*CosMN(iMode)
        !dR/dtheta (cos(m*theta-n*zeta))=-m*sin(m*theta-n*zeta)
        dRdU = dRdU - rhom*splout(1)*SinMN(iMode)*xm(iMode) 
        !dR/drho = sum_mn [ rho**m * dR_mn/drho + R_mn *d(rho**m)/drho ] * cos(m*theta-n*zeta)
        dRdrho=dRdrho+(rhom*splout(2)+splout(1)*drhom)*CosMN(iMode)
        CALL SPLINE1_EVAL((/1,1,0/), nFluxVMEC,rho(iFlux),rho,Zmns_Spl(:,:,iMode),iGuess,splout) 
        Z    = Z    + rhom*splout(1)*SinMN(iMode)
        !dZ/dtheta (sin(m*theta-n*zeta))=m*cos(m*theta-n*zeta)
        dZdU = dZdU + rhom*splout(1)*CosMN(iMode)*xm(iMode) 
        !dZ/dzeta  (sin(m*theta-n*zeta))=-n*cos(m*theta-n*zeta)
        dZdrho=dZdrho+(rhom*splout(2)+splout(1)*drhom)*SinMN(iMode)
        !Jacobian: R*(dR/dtheta*dZ/drho-dR/drho*dZ/dtheta)
      END DO !iMode=1,mn_mode
      gsqrt(iFlux,iNode,jNode)=R*(dRdU*dZdrho - dRdrho*dZdU)
    END DO! iFlux=2,nFluxVMEC
  END DO! iNode=1,ntheta
END DO! jNode=1,nzeta
!CALCULATE LAMBDA

      DO iMode=1,mn_mode_nyq !Over modes

        IF(iMode /= Mode00) THEN

!Loop over radial grid
          DO iFlux=2,nFluxVMEC

!Loop over theta grid
            DO iNode=1,ntheta_3d !Over poloidal nodes

              g1(iNode)=gsqrt_3d(iFlux,iNode,1)/rcyl_3d(iFlux,iNode,1)**2
              theta=(iNode-1)*dtheta_3d
              g2(iNode)=g1(iNode)*COS(xm_nyq(iMode)*theta)

            ENDDO !Over poloidal nodes

!Theta averages to find lambda coefficients
            lam_3d(1,i,k)=2.0_rspec*SUM(wtheta_3d*g2)/SUM(wtheta_3d*g1)        &
     &                    /REAL(m_3d(k),rspec) !Over poloidal nodes

          ENDDO !Over radial nodes

          lam_3d(1,2:nrho_3d,k)=lam_3d(1,2:nrho_3d,k)                          &
     &                          /rho_3d(2:nrho_3d)**mabs_3d(k)

!Make a parabolic fit to axis
          lam_3d(1,1,k)=(lam_3d(1,2,k)*rho_3d(3)**2                            &
     &                  -lam_3d(1,3,k)*rho_3d(2)**2)                           &
     &                  /(rho_3d(3)**2-rho_3d(2)**2)
 
!-------------------------------------------------------------------------------
!Spline the Lmn coefficients for internal storage (not-a-knot edge BC)
!-------------------------------------------------------------------------------
          CALL SPLINE1_FIT(nrho_3d,rho_3d,lam_3d(:,:,k),                       &
     &                     K_BC1=k_bc1,                                        &
     &                     K_BCN=k_bcn)

        ENDIF

      ENDDO !Over modes



WRITE(UNIT_stdOut,'(A)')'  ... DONE'
END SUBROUTINE RecomputeLambda 
