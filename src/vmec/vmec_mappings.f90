! ----------------------------------------------------------------------
! Version 3.75
!
!> compute smoothing of VMEC r,z
!> recalculate magnetic field
!> compute mapping for euterpe
! ----------------------------------------------------------------------
! last changed <03.12.15 by mzb>                           M. Borchardt
! ----------------------------------------------------------------------

MODULE MOD_VMEC_Mappings

  IMPLICIT NONE
  PUBLIC

!> high precision real numbers
  INTEGER, PARAMETER :: hp = SELECTED_REAL_KIND(10, 20)
!> Pi
  REAL(KIND = hp), PARAMETER :: Pi    = 3.1415926535897932384626433832795_hp
!> 2 Pi
  REAL(KIND = hp), PARAMETER :: TwoPi = 6.2831853071795864769252867665590_hp
!> mu_0 (Permeability)
  REAL(KIND = hp), PARAMETER :: mu0   = 4 * Pi * 1E-7_hp         ! V s / (A m)

!> version number
  REAL(KIND = hp), PARAMETER :: version = 3.75_hp

!> local node number
  INTEGER :: me_rank

!> number of nodes
  INTEGER :: me_size

!> MPI error
  INTEGER :: ierr

!> write mapping for Euterpe (default: .TRUE.)
  LOGICAL :: mapEuterpe = .TRUE.

!> use sqrt(s) as radial coordinate (default: .FALSE.)
  LOGICAL :: useRho = .FALSE.

!> write values for interpolation using arcus tangens
  LOGICAL :: useArcTan = .FALSE.

!> write data for perpendicular gyro ring
  LOGICAL :: gyroperp = .FALSE.

!> maximal number of iteration for GMRES
  INTEGER, PARAMETER :: itSolMaxIter = 800

!> minimal number of iteration for GMRES
  INTEGER, PARAMETER :: minIter = 10

!> small number
  REAL(KIND = hp), PARAMETER :: tiny = 1E-20_hp

!> number of integration points in u direction
  INTEGER :: intPointsU = 128

!> number of integration points in v direction
  INTEGER :: intPointsV = 128

!> use secant method instead of newton method
  LOGICAL :: useSecant = .FALSE.

!> convergence criterion for the Newton method
  REAL(KIND = hp) :: newtonEps = 1E-8_hp

!> maximum number of newton iterations
  INTEGER :: newtonMax = 100

!> maximal number of recursions in FindInterSection
!> per recursion relaxation is changed by factor of relax
  INTEGER :: maxRec = 50

!> number of grid points in R direction
  INTEGER :: rGridPoints = 51

!> number of grid points in z direction
  INTEGER :: zGridPoints = 51

!> number of flux surfaces to store
  INTEGER :: nFlux = -1

!> last extrapolated "flux" surface
  REAL(KIND = hp) :: lastExtS = 1.1_hp

!> number of flux surfaces outside of s = 1
  INTEGER :: nExtraFlux = -1

!> number of theta values to store
  INTEGER :: nTheta = 128

!> number of points in u direction for searching ray starting point
  INTEGER :: nPointsU = 512

!> number of grid points in u direction
  INTEGER :: nu = 100

!> number of grid points in v direction
  INTEGER :: nv = 100

!> number of interpolation points for theta grid (default: 200)
  INTEGER :: nThetaInt = 200

!> number of interpolation points for phi grid (default: 200)
  INTEGER :: nPhiInt = 200

!> number of phi cuts
  INTEGER :: nPhiCuts = 64

!> starting cut
  INTEGER :: startCut = -1

!> ending cut
  INTEGER :: endCut = -1

!> number geometrical angles for mapping from geometrical angle to PEST angle
!> theta
  INTEGER :: nGeomAngle = 360

!> relaxation factor
  REAL(KIND = hp) :: relax = .97_hp

!> distance until values is extrapolated (with respect to s = 1 flux surface)
  REAL(KIND = hp) :: constFac = 2.5_hp

!> number of inner flux surfaces for mapping near the center
  INTEGER :: nInnerS

!> number of angles per flux surfaces for mapping near the center
  INTEGER :: nInTheta = 400

!> flux surface label until which "polar"-like coordinates will be calculated
  REAL(KIND = hp) :: inSVal = .1

!> grid boundary low r
  REAL(KIND = hp) :: lowBR

!> grid boundary high r
  REAL(KIND = hp) :: highBR

!> grid boundary low z
  REAL(KIND = hp) :: lowBZ

!> grid boundary high z
  REAL(KIND = hp) :: highBZ

!> major radius of the device
  REAL(KIND = hp) :: rMajor

!> starting index for polar grid
  INTEGER :: stInd
  
!> number of flux surfaces in VMEC output
  INTEGER :: nFluxVMEC

!> number of flux surfaces in Mapping
  INTEGER :: radius

!> number of modes (2D)
  INTEGER :: mn_mode

!> number of modes (Nyquist, 2D)
  INTEGER :: mn_mode_nyq

!> number of modes for magnetic axis (1D)
  INTEGER :: nAxis

!> number of field periods
  INTEGER :: nfp

!> dimension of s
  INTEGER :: nS

!> poloidal mode number
  INTEGER :: mPol

!> toroidal mode number
  INTEGER :: nTor

!> mnmax
  INTEGER :: mnmax

!> mnmax_nyq
  INTEGER :: mnmax_nyq

!> B_0
  REAL(KIND = hp) :: b0

!> volume
  REAL(KIND = hp) :: volume

!> signum of sqrtG
  INTEGER :: signgs

  !! vector parameter
!> poloidal mode numbers
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: xm

!> toroidal mode numbers
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: xn

!> poloidal mode numbers (Nyquist)
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: xm_nyq

!> toroidal mode numbers (Nyquist)
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: xn_nyq

!> iota on half mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: iotas

!> pressure on full mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: presf

!> toroidal flux on full mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: phi

!> d(phi)/ds: Toroidal flux derivative on full mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: phipf

!> i: Toroidal current on full mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: icur

!> j: poloidal current on full mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: jcur

!> R (cosine components on full mesh)
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: rmnc

!> z (sine components on full mesh)
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: zmns

!> lambda (sine components (read on half mesh, interpolated on full mesh))
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: lmns

!> kappa
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: kmns

!> jacobian (cosine components (read on half mesh, interpolated on full
!> mesh))
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gmnc

!> dsqrtG/ds
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dgdSmnc

!> |B} (cosine components (read on half mesh, interpolated on full mesh))
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: bmnc

  !! B^u, B^v (cosine components (read on half mesh, interpolated on full
  !! mesh))
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: bsupumnc, bsupvmnc

  !! dlambda/du, dlambda/dv (cosine components on full mesh)
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: lUmnc, lVmnc

  !! dkappa/du, dkappa/dv (cosine components on full mesh)
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: kUmnc, kVmnc

  !! dlambda/du2, dlambda/(du dv), dlambda/dv2 (sine components on full mesh)
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: lU2mns, lUVmns, lV2mns

  !! cosine components of cylindrical R in PEST coordinates (full mesh)
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: cylRmnc

  !! sine components of cylindrical z in PEST coordinates (full mesh)
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: cylZmns

  !! fourier components of the derivatives of R
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dRdSmnc, dRdUmns, dRdVmns

  !! fourier components of the derivatives of z
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dZdSmns, dZdUmnc, dZdVmnc

  !! fourier components of the second derivatives of r, z
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: d2RdS2mnc, d2ZdS2mns
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: d2RdUdSmns, d2ZdUdSmnc
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: d2RdVdSmns, d2ZdVdSmnc
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: d2RdU2mnc, d2ZdU2mns
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: d2RdUdVmnc, d2ZdUdVmns
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: d2RdV2mnc, d2ZdV2mns

  !! fourier components of grad s
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradSRmnc, gradSPmns
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradSZmns

  !! fourier components of grad theta
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradTRmns, gradTPmnc
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradTZmnc

  !! fourier components of |B|
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: absBmnc

  !! fourier components of |B| * |sqrtG|
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: absBsmnc

  !! fourier components of d|B|/ds
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dBdSmnc

  !! fourier components of d(|B| * |sqrtG|)/ds
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dBdSsmnc

  !! fourier components of b
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: bRmns, bZmnc, bPmnc

  !! fourier components of dB/du, dB/dv
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dBdUmns, dBdVmns

  !! fourier components of scaled dB/du, dB/dv
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dBdUsmns, dBdVsmns

  !! fourier components of grad B
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradBRmnc, gradBPmns
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradBZmns

  !! fourier components of b x grad B / |B|
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: bgradBRmns, bgradBPmnc
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: bgradBZmnc

  !! fourier components of div b
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: divBmns

  !! fourier components of rot B parallel
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: rotBparAbsmnc

  !! fourier components of rot B perpendicular
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: rotBperpRmns, rotBperpPmnc
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: rotBperpZmnc

  !! r and z mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: rMesh, zMesh

  !! poloidal and toroidal current
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: iPol, jTor

  !! precomputed cosinus and sinus values
  REAL(KIND = hp), DIMENSION(:, :, :), ALLOCATABLE :: vCos, vSin, vCosN, vSinN

  !! precomputed integration points
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: uIv, vJv

!> output directory
  CHARACTER(LEN = 255) :: outDir = "."

!> show debug output
  LOGICAL :: debug = .FALSE.

!> show numerical debug output
  LOGICAL :: numDebug = .FALSE.

!> correct VMEC data with spline interpolation
  LOGICAL :: corVMEC = .TRUE.

!> scale magnetic field by sqrtG
  LOGICAL :: useScaledB = .FALSE.

!> include error variance in interpolation
  LOGICAL :: useVar = .TRUE.

!> use last value or linear extrapolation for s > 1
  LOGICAL :: useLastVal = .FALSE.

!> use constant theta path instead of inter-/extrapolate along a line
  LOGICAL :: useConstTheta = .FALSE.

!> use squared flux surface coordinate
  LOGICAL :: useSquared = .TRUE.

!> error intervals
  REAL(KIND = hp), DIMENSION(2) :: errInt = (/.05_hp, .95_hp/)

!> error variance in intervals
  REAL(KIND = hp), DIMENSION(3) :: errVal = (/1E-1_hp, 1E-3_hp, 1E-2_hp/)

!> weight
  REAL(KIND = hp) :: weight = 1E-3_hp

  !! spline coefficients for corrected r, z
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: xr, xz

  !! marker for different spline coefficients
  INTEGER :: ma, mb, mc, md

!> asymptotic behaviour for different modes
  INTEGER, DIMENSION(:), ALLOCATABLE :: asymp

!> asymptotic function for spline approximation
  REAL(KIND = hp), DIMENSION(:, :, :), ALLOCATABLE :: t

  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: argA, argA_nyq
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: sVal, axisM, axisN, theta
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: sValVMEC

  REAL(KIND = hp), DIMENSION(6) :: maxDist = (/-1, 0, 0, 0, 0, 0/)

CONTAINS

! ----------------------------------------------------------------------
!> allocate arrays used for mapping
! ----------------------------------------------------------------------

  SUBROUTINE AllocArrays

    INTEGER :: aStat, aError

    WRITE(*,*)'VMEC ALLOCATE ARRAYS...'
    !! allocate new variables
    aError = 0
    ALLOCATE(lUmnc(mn_mode, stInd:radius), lVmnc(mn_mode, stInd:radius),&
         stat = aStat)
    aError = aError + aStat
    ALLOCATE(kUmnc(mn_mode, radius), kVmnc(mn_mode, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(lU2mns(mn_mode, radius), lV2mns(mn_mode, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(lUVmns(mn_mode, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(cylRmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(cylZmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dRdSmnc(mn_mode, stInd:radius), dZdSmns(mn_mode, stInd:radius),&
         stat = aStat)
    aError = aError + aStat
    ALLOCATE(d2RdS2mnc(mn_mode, stInd:radius), d2ZdS2mns(mn_mode,&
         stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(d2RdUdSmns(mn_mode, stInd:radius), d2ZdUdSmnc(mn_mode,&
         stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(d2RdVdSmns(mn_mode, stInd:radius), d2ZdVdSmnc(mn_mode,&
         stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(d2RdU2mnc(mn_mode, stInd:radius), d2ZdU2mns(mn_mode,&
         stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(d2RdUdVmnc(mn_mode, stInd:radius), d2ZdUdVmns(mn_mode,&
         stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(d2RdV2mnc(mn_mode, stInd:radius), d2ZdV2mns(mn_mode,&
         stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dRdUmns(mn_mode, stInd:radius), dZdUmnc(mn_mode, stInd:radius),&
         stat = aStat)
    aError = aError + aStat
    ALLOCATE(dRdVmns(mn_mode, stInd:radius), dZdVmnc(mn_mode, stInd:radius),&
         stat = aStat)
    aError = aError + aStat
    ALLOCATE(gradSRmnc(mn_mode, radius), gradTRmns(mn_mode, radius), stat =&
         aStat)
    aError = aError + aStat
    ALLOCATE(gradSPmns(mn_mode, radius), gradTPmnc(mn_mode, radius), stat =&
         aStat)
    aError = aError + aStat
    ALLOCATE(gradSZmns(mn_mode, radius), gradTZmnc(mn_mode, radius), stat =&
         aStat)
    aError = aError + aStat
    ALLOCATE(argA(mn_mode, 2), argA_nyq(mn_mode_nyq, 2), stat = aStat)
    aError = aError + aStat
    ALLOCATE(absBmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(absBsmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dBdSmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dBdSsmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bRmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bZmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bPmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dBdUmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dBdVmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dBdUsmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dBdVsmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(gradBRmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(gradBZmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(gradBPmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bgradBRmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bgradBZmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bgradBPmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(divBmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(rotBparAbsmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(rotBperpRmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(rotBperpZmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(rotBperpPmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(rMesh(0:(rGridPoints+1)), zMesh(0:(zGridPoints+1)), stat = aStat)
    aError = aError + aStat
    ALLOCATE(sVal(stInd:radius), theta(radius), sValVMEC(nFluxVMEC), stat = aStat)
    aError = aError + aStat
    ALLOCATE(axisM(nAxis), axisN(nAxis), stat = aStat)
    aError = aError + aStat
    ALLOCATE(iPol(radius), jTor(radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(vCos(mn_mode, intPointsU, intPointsV), stat = aStat)
    aError = aError + aStat
    ALLOCATE(vSin(mn_mode, intPointsU, intPointsV), stat = aStat)
    aError = aError + aStat
    ALLOCATE(vCosN(mn_mode_nyq, intPointsU, intPointsV), stat = aStat)
    aError = aError + aStat
    ALLOCATE(vSinN(mn_mode_nyq, intPointsU, intPointsV), stat = aStat)
    aError = aError + aStat
    ALLOCATE(uIv(intPointsU), vJv(intPointsV), stat = aStat)
    aError = aError + aStat

    IF (aError /= 0) THEN
      PRINT *, "Allocation failure in program VM2MAG!"
      CALL EXIT(1)
    END IF
    WRITE(*,*)'...DONE.'
  END SUBROUTINE AllocArrays

! ----------------------------------------------------------------------
!> read VMEC2000 netcdf output file
! ----------------------------------------------------------------------

  SUBROUTINE ReadVmecOutput(fileName)

    INCLUDE "netcdf.inc"

    CHARACTER(LEN = *), INTENT(IN) :: fileName

    INTEGER :: aStat, aError, ioError, ncid, id

    WRITE(*,*)'VMEC READ WOUT FILE...'
    !! open NetCDF input file
    ioError = NF_OPEN(TRIM(fileName), NF_NOWRITE, ncid)
    IF (ioError /= 0) THEN
      PRINT *, " Cannot open ", TRIM(fileName), " in ReadVmecOutput!"
      CALL EXIT(2)
    END IF

    !! get array dimensions
    !! radial dimension
    ioError = NF_INQ_DIMID(ncid, "radius", id)
    ioError = ioError + NF_INQ_DIMLEN(ncid, id, radius)
    !! number of fourier components of r, z, lambda
    ioError = ioError + NF_INQ_DIMID(ncid, "mn_mode", id)
    ioError = ioError + NF_INQ_DIMLEN(ncid, id, mn_mode)
    !! number of fourier components of b_u, b_v, b_s
    ioError = ioError + NF_INQ_DIMID(ncid, "mn_mode_nyq", id)
    ioError = ioError + NF_INQ_DIMLEN(ncid, id, mn_mode_nyq)
    !! number of fourier components of raxis, zaxis
!    ioError = ioError + NF_INQ_DIMID(ncid, "n-tor", id)
!    ioError = ioError + NF_INQ_DIMLEN(ncid, id, nAxis)

    !! get number of field periods
    ioError = ioError + NF_INQ_VARID(ncid, "nfp", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, nfp)
    !! get dimension of s
    ioError = ioError + NF_INQ_VARID(ncid, "ns", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, nS)
    !! get poloidal mode number
    ioError = ioError + NF_INQ_VARID(ncid, "mpol", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, mPol)
    !! get toroidal mode number
    ioError = ioError + NF_INQ_VARID(ncid, "ntor", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, nTor)
    !! get mnmax
    ioError = ioError + NF_INQ_VARID(ncid, "mnmax", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, mnmax)
    !! get mnmax_nyq
    ioError = ioError + NF_INQ_VARID(ncid, "mnmax_nyq", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, mnmax_nyq)
    !! get B_0
    ioError = ioError + NF_INQ_VARID(ncid, "b0", id)
    ioError = ioError + NF_GET_VAR_DOUBLE(ncid, id, b0)
    !! check the sign of b0
    IF (b0 < 0) THEN
      IF (me_rank == 0) THEN
        PRINT *, "  VMEC run with b0 < 0 !!!"
      END IF
    END IF
    !! absolute value of b0 has to be non-negative
    b0 = ABS(b0)
    !! get signgs
    ioError = ioError + NF_INQ_VARID(ncid, "signgs", id)
    ioError = ioError + NF_GET_VAR_INT(ncid, id, signgs)
    !! get Rmajor_p
    ioError = ioError + NF_INQ_VARID(ncid, "Rmajor_p", id)
    ioError = ioError + NF_GET_VAR_DOUBLE(ncid, id, rMajor)
    !! get volume_p
    ioError = ioError + NF_INQ_VARID(ncid, "volume_p", id)
    ioError = ioError + NF_GET_VAR_DOUBLE(ncid, id, volume)

    IF (ioError /= 0) THEN
      PRINT *, " Cannot read ", TRIM(fileName), "!"
      PRINT *, " Possible wrong file format!"
      CALL EXIT(3)
    END IF

    !! initialize starting index of arrays for polar grid
    stInd = 1

    !! calculate number of flux surfaces
    nFluxVMEC = radius
    IF (mapEuterpe) THEN
      IF (lastExtS < 1) THEN
        nExtraFlux = 0
        lastExtS = 1
      ELSE
        IF (nExtraFlux < 0) THEN
          IF (useSquared) THEN
            nExtraFlux = INT((SQRT(lastExtS) - 1) * (nFluxVMEC - 1)) + 1
          ELSE
            nExtraFlux = INT((lastExtS - 1) * (nFluxVMEC - 1)) + 1
          END IF
          lastExtS = 1 + REAL(nExtraFlux, hp) / REAL(nFluxVMEC - 1, hp)
        END IF
      END IF
      radius = nFluxVMEC + nExtraFlux
      !! calculate real number of flux surfaces for the inner "polar"-like
      !! coordinate system using squared distances
      nInnerS = INT(SQRT(inSVal) * nFluxVMEC) + 1
      !! set starting index of arrays for polar grid
      !! we need extra flux surfaces if linear distances of the flux surfaces
      !! are used
      IF (.NOT.(useSquared)) stInd = -nInnerS
    END IF

    !! allocate memory for fourier arrays
    aError = 0
    ALLOCATE(xm(mn_mode), xn(mn_mode), stat = aStat)
    aError = aError + aStat
    ALLOCATE(xm_nyq(mn_mode_nyq), xn_nyq(mn_mode_nyq), stat = aStat)
    aError = aError + aStat
    ALLOCATE(iotas(stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(presf(radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(phi(radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(phipf(radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(icur(radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(jcur(radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(rmnc(mn_mode, stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(zmns(mn_mode, stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(lmns(mn_mode, stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(kmns(mn_mode, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(gmnc(mn_mode_nyq, stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dgdSmnc(mn_mode_nyq, stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bsupumnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bsupvmnc(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat

    IF (aError /= 0) THEN
      PRINT *, "Allocation failure in subroutine ReadVmecOutput!"
      CALL EXIT(4)
    END IF

    !! initialize some arrays
    iotas(:) = 0
    phipf(:) = 0
    presf(:) = 0

    !! read x_m
    ioError = NF_INQ_VARID(ncid, "xm", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /), (/ mn_mode /),&
         xm(:))
    !! read x_n
    ioError = NF_INQ_VARID(ncid, "xn", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /), (/ mn_mode /),&
         xn(:))
    !! read x_m^nyq
    ioError = NF_INQ_VARID(ncid, "xm_nyq", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ mn_mode_nyq /), xm_nyq(:))
    !! read x_n^nyq
    ioError = NF_INQ_VARID(ncid, "xn_nyq", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ mn_mode_nyq /), xn_nyq(:))
    !! read iotas
    ioError = NF_INQ_VARID(ncid, "iotas", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), iotas(1:))
    !! read presf
    ioError = NF_INQ_VARID(ncid, "presf", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), presf(:))
    !! read phi
    ioError = NF_INQ_VARID(ncid, "phi", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), phi(:))
    !! scale toroidal flux to get internal VMEC Phi
    phi(:nFluxVMEC) = REAL(signgs, hp) * phi(:nFluxVMEC) / TwoPi
    !! read phipf
    ioError = NF_INQ_VARID(ncid, "phipf", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), phipf(:))
    !! scale toroidal flux to get internal VMEC Phi
    phipf(:nFluxVMEC) = REAL(signgs, hp) * phipf(:nFluxVMEC) / TwoPi
    !! read toroidal current
    ioError = NF_INQ_VARID(ncid, "jcuru", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), icur(:))
    !! read poloidal current
    ioError = NF_INQ_VARID(ncid, "jcurv", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), jcur(:))
    !! read R_mn
    ioError = NF_INQ_VARID(ncid, "rmnc", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), rmnc(:, 1:))
    !! read Z_mn
    ioError = NF_INQ_VARID(ncid, "zmns", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), zmns(:, 1:))
    !! read lambda_mn
    ioError = NF_INQ_VARID(ncid, "lmns", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), lmns(:, 1:))
    !! read jacobian_mn
    ioError = NF_INQ_VARID(ncid, "gmnc", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
         mn_mode_nyq, nFluxVMEC /), gmnc(:, 1:))
    !! read b_mn
    ioError = NF_INQ_VARID(ncid, "bmnc", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
         mn_mode_nyq, nFluxVMEC /), bmnc(:, :))
    !! read B^u_mn
    ioError = NF_INQ_VARID(ncid, "bsupumnc", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
         mn_mode_nyq, nFluxVMEC /), bsupumnc(:, :))
    !! read B^v_mn
    ioError = NF_INQ_VARID(ncid, "bsupvmnc", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
         mn_mode_nyq, nFluxVMEC /), bsupvmnc(:, :))

    IF (ioError /= 0) THEN
      PRINT *, " Cannot read variables from ", TRIM(fileName), "!"
      PRINT *, " Possible wrong data!"
      CALL EXIT(5)
    END IF

    ioError = NF_CLOSE(ncid)
    WRITE(*,*)'...DONE.'

  END SUBROUTINE ReadVmecOutput

  SUBROUTINE PreCalcData

    REAL(KIND = hp), DIMENSION(radius) :: hVec
    REAL(KIND = hp) :: help, uI, vJ, rS, zS, lambda, dLambdadU
    REAL(KIND = hp) :: r, dRdS, dRdU, dZdS, dZdU, sqrtG
    REAL(KIND = hp) :: d2RdS2, d2RdUdS, d2ZdS2, d2ZdUdS, dsqrtGdS
    REAL(KIND = hp) :: arg
    INTEGER :: i, j, k, curS

    WRITE(*,*)'VMEC PRECALC DATA ...'
    !! set s values of readed flux surfaces
    DO curS = 1, nFluxVMEC
      sVal(curS) = REAL(curS - 1, hp) / REAL(nFluxVMEC - 1, hp)
      sValVMEC(curS) = sVal(curS)
    END DO

    IF (mapEuterpe) THEN
      DO curS = 1, nExtraFlux
        sVal(curS + nFluxVMEC) = 1 + REAL(curS, hp) *&
             (lastExtS - 1) / REAL(nExtraFlux, hp)
      END DO
      !! calculate extra flux surfaces if linear distances of flux surfaces
      !! are used
      IF (.NOT.(useSquared)) THEN
        DO curS = 1, nInnerS
          sVal(-curS) = REAL(curS - 1, hp) / REAL(nFluxVMEC - 1, hp)
        END DO
      END IF
    END IF

    !! precompute integration points
    DO i = 1, intPointsU
      uIv(i) = TwoPi * REAL(i - 1, hp) / REAL(intPointsU, hp)
    END DO
    DO j = 1, intPointsV
      vJv(j) = TwoPi * REAL(j - 1, hp) / REAL(intPointsV * nfp, hp)
    END DO

    !! precompute cosinus and sinus values for fourier transforms
    DO k = 1, mn_mode
      DO j = 1, intPointsV
        arg = xn(k) * vJv(j)
        DO i = 1, intPointsU
          vCos(k, i, j) = COS(xm(k) * uIv(i) - arg)
          vSin(k, i, j) = SIN(xm(k) * uIv(i) - arg)
        END DO
      END DO
    END DO

    !! precompute cosinus and sinus values for fourier transforms (Nyquist)
    DO k = 1, mn_mode_nyq
      DO j = 1, intPointsV
        arg = xn_nyq(k) * vJv(j)
        DO i = 1, intPointsU
          vCosN(k, i, j) = COS(xm_nyq(k) * uIv(i) - arg)
          vSinN(k, i, j) = SIN(xm_nyq(k) * uIv(i) - arg)
        END DO
      END DO
    END DO

    help = REAL(intPointsU * intPointsV, hp)
    !! interpolate half mesh data to full mesh
    CALL HalfMeshToFullMesh1D(nFluxVMEC, iotas(1:nFluxVMEC))
    CALL HalfMeshToFullMesh2D(mn_mode, nFluxVMEC, lmns(:, 1:nFluxVMEC))
    CALL HalfMeshToFullMesh2D(mn_mode_nyq, nFluxVMEC, gmnc(:, 1:nFluxVMEC))
    CALL HalfMeshToFullMesh2D(mn_mode_nyq, nFluxVMEC, bmnc(:, :nFluxVMEC))
    CALL HalfMeshToFullMesh2D(mn_mode_nyq, nFluxVMEC, bsupumnc(:, :nFluxVMEC))
    CALL HalfMeshToFullMesh2D(mn_mode_nyq, nFluxVMEC, bsupvmnc(:, :nFluxVMEC))

    !! write original R, z values for test purpose
    IF (debug .AND. (me_rank == 0)) THEN
      DO i = 1, mn_mode
        DO j = 1, nFluxVMEC
          WRITE(111, FMT = '(ES13.5, ",", I6, ",", ES13.5, ",", ES13.5, ",",&
               & ES13.5, ",", ES13.5)') sVal(j), i, xm(i), xn(i), rmnc(i, j),&
               zmns(i, j)
        END DO
        WRITE(111, *)
      END DO
      CLOSE(111)
    END IF

    !! smooth cylindrical coordinates
    IF (corVMEC) THEN
      CALL SmoothRZ
      !! recalculate iota on changed grid
      hVec(:) = iotas(1:)
      DO curS = stInd, radius
        IF (curS == 0) CYCLE
        iotas(curS) = IntExtPol3(nFluxVMEC, sValVMEC(:),&
             hVec(:nFluxVMEC), sVal(curS))
      END DO
      !! recalculate psi on changed grid
      hVec(:) = phipf(:)
      DO curS = 1, radius
        phipf(curS) = IntExtPol3(nFluxVMEC, sValVMEC(:),&
             hVec(:nFluxVMEC), sVal(curS))
      END DO
      !! recalculate pressure on changed grid
      hVec(:) = presf(:)
      DO curS = 1, radius
        presf(curS) = IntExtPol3(nFluxVMEC, sValVMEC(:),&
             hVec(:nFluxVMEC), sVal(curS))
      END DO
    END IF

    !! calculate poloidal and toroidal partial derivatives
    DO curS = stInd, radius
      IF (curS == 0) CYCLE
      !! dR/du
      dRdUmns(:, curS) = -xm(:) * rmnc(:, curS)
      !! d^2R/(du^2)
      d2RdU2mnc(:, curS) = xm(:) * dRdUmns(:, curS)
      !! d^2R/(du dv)
      d2RdUdVmnc(:, curS) = -xn(:) * dRdUmns(:, curS)
      !! dR/dv
      dRdVmns(:, curS) = xn(:) * rmnc(:, curS)
      !! d^2R/(dv^2)
      d2RdV2mnc(:, curS) = -xn(:) * dRdVmns(:, curS)
      !! dz/du
      dZdUmnc(:, curS) = xm(:) * zmns(:, curS)
      !! d^2z/(du^2)
      d2ZdU2mns(:, curS) = -xm(:) * dZdUmnc(:, curS)
      !! d^2z/(du dv)
      d2ZdUdVmns(:, curS) = xn(:) * dZdUmnc(:, curS)
      !! dz/dv
      dZdVmnc(:, curS) = -xn(:) * zmns(:, curS)
      !! d^2z/(dv^2)
      d2ZdV2mns(:, curS) = xn(:) * dZdVmnc(:, curS)
    END DO

    !! recompute lambda from interpolated coordinates
    IF (corVMEC) THEN
      !! compute second partial derivatives
      DO curS = stInd, radius
        IF (curS == 0) CYCLE
        !! d^2R/(ds du)
        d2RdUdSmns(:, curS) = -xm(:) * dRdSmnc(:, curS)
        !! d^2R/(ds dv)
        d2RdVdSmns(:, curS) = xn(:) * dRdSmnc(:, curS)
        !! d^2z/(ds du)
        d2ZdUdSmnc(:, curS) = xm(:) * dZdSmns(:, curS)
        !! d^2z/(ds dv)
        d2ZdVdSmnc(:, curS) = -xn(:) * dZdSmns(:, curS)
      END DO
      !! sqrtG
      gmnc(:, :) = 0
      !! dsqrtG/ds
      dgdSmnc(:, :) = 0
      !! numerical integration (Trapez Rule)
      DO j = 1, intPointsV
        vJ = vJv(j)
        DO i = 1, intPointsU
          !! supporting points for u direction
          uI = uIv(i)
          !! radial loop
          DO curS = stInd, radius
            IF (curS == 0) CYCLE
            !! calculate R
            r = CosTransFullMeshF2(mn_mode, rmnc(:, curS), i, j)
            !! dR/ds
            dRdS = CosTransFullMeshF2(mn_mode, dRdSmnc(:, curS), i, j)
            !! dR/du
            dRdU = SinTransFullMeshF2(mn_mode, dRdUmns(:, curS), i, j)
            !! d^2R/ds^2
            d2RdS2 = CosTransFullMeshF2(mn_mode, d2RdS2mnc(:, curS), i, j)
            !! d^2R/(du ds)
            d2RdUdS = SinTransFullMeshF2(mn_mode, d2RdUdSmns(:, curS), i, j)
            !! dz/ds
            dZdS = SinTransFullMeshF2(mn_mode, dZdSmns(:, curS), i, j)
            !! dz/du
            dZdU = CosTransFullMeshF2(mn_mode, dZdUmnc(:, curS), i, j)
            !! d^2z/ds^2
            d2ZdS2 = SinTransFullMeshF2(mn_mode, d2ZdS2mns(:, curS), i, j)
            !! d^2z/(du ds)
            d2ZdUdS = CosTransFullMeshF2(mn_mode, d2ZdUdSmnc(:, curS), i, j)
            !! help arguments
            argA_nyq(1, :) = (/0._hp, 1._hp/)
            argA_nyq(2:, 2) = 2 * vCosN(2:, i, j)
            sqrtG = r * (dRdU * dZdS - dZdU * dRdS)
            dsqrtGdS = dRdS * (dRdU * dZdS - dZdU * dRdS) + r * (d2RdUdS *&
                 dZdS + dRdU * d2ZdS2 - (d2ZdUdS * dRdS + dZdU * d2RdS2))
            gmnc(:, curS) = gmnc(:, curS) + sqrtG * argA_nyq(:, 2)
            dgdSmnc(:, curS) = dgdSmnc(:, curS) + dsqrtGdS * argA_nyq(:, 2)
          END DO
        END DO
      END DO
      gmnc(:, :) = gmnc(:, :) / help
      dgdSmnc(:, :) = dgdSmnc(:, :) / help
      CALL ComputeLambda
    END IF

    !! calculate first and second partial derivatives of lambda
    DO curS = 1, radius
      !! dlambda/du
      lUmnc(:, curS) = xm(:) * lmns(:, curS)
      !! dlambda/dv
      lVmnc(:, curS) = -xn(:) * lmns(:, curS)
      !! d^2lambda/du^2
      lU2mns(:, curS) = -xm(:) * lUmnc(:, curS)
      !! d^2lambda/(du dv)
      lUVmns(:, curS) = xn(:) * lUmnc(:, curS)
      !! d^2lambda/dv^2
      lV2mns(:, curS) = xn(:) * lVmnc(:, curS)
    END DO

    IF (mapEuterpe .AND. .NOT.(useSquared)) THEN
    !! calculate first and second partial derivatives of lambda
      DO curS = stInd, -1
        !! dlambda/du
        lUmnc(:, curS) = xm(:) * lmns(:, curS)
        !! dlambda/dv
        lVmnc(:, curS) = -xn(:) * lmns(:, curS)
      END DO
    END IF
    
    !! recompute magnetic field
    IF (corVMEC) THEN
      !! replace VMEC B field (bmnc)
      CALL ComputeMagnField
    END IF
    !!              2L                                   (    dlambda)
    !! R_mn(s) = ------- Int_0^{2Pi/L} dv Int_0^{2Pi} du (1 + -------) *
    !!           (2Pi)^2                                 (      du   )
    !!
    !!                 * R_uv(s) cos(m(u+lambda)-nv)  if n/=0 or m/=0
    !!
    !!               L                                   (    dlambda)
    !! R_00(s) = ------- Int_0^{2Pi/L} dv Int_0^{2Pi} du (1 + -------) * R_uv(s)
    !!           (2Pi)^2                                 (      du   )
    !!
    !! for Z_mn(s) same equation with R_uv(s) replaced by Z_uv(s) and
    !! cos(..) replaced by sin(..)
    !! Z_00(s) = 0

    !! use the same fourier components for cylindrical R, Z in PEST coordinates
    cylRmnc(:, :) = 0
    cylZmns(:, :) = 0
    !! numerical integration (Trapez Rule)
    DO j = 1, intPointsV
      !! supporting points for v direction
      vJ = vJv(j)
      DO i = 1, intPointsU
        !! supporting points for u direction
        uI = uIv(i)
        !! radial loop
        DO curS = 1, radius
          !! calculate lambda
          lambda = SinTransFullMeshF2(mn_mode, lmns(:, curS), i, j)
          !! calculate 1 + dlambda/du
          dLambdadU = 1._hp + CosTransFullMeshF2(mn_mode, lUmnc(:, curS), i, j)
          !! calculate R
          rS = CosTransFullMeshF2(mn_mode, rmnc(:, curS), i, j)
          !! calculate z
          zS = SinTransFullMeshF2(mn_mode, zmns(:, curS), i, j)
          !! help arguments
          argA_nyq(1, :) = (/0._hp, 1._hp/)
          argA_nyq(2:, 1) = 2 * SIN(xm_nyq(2:) * (uI + lambda) - xn_nyq(2:) * vJ)
          argA_nyq(2:, 2) = 2 * COS(xm_nyq(2:) * (uI + lambda) - xn_nyq(2:) * vJ)
          cylRmnc(:, curS) = cylRmnc(:, curS) + dLambdadU * rS * argA_nyq(:, 2)
          cylZmns(:, curS) = cylZmns(:, curS) + dLambdadU * zS * argA_nyq(:, 1)
        END DO
      END DO
    END DO

    !! scale fourier coefficients
    cylRmnc(:, :) = cylRmnc(:, :) / help
    cylZmns(:, :) = cylZmns(:, :) / help

    !! calculate partial derivatives
    DO curS = 1, radius
      !! dB/du
      dBdUmns(:, curS) = -xm_nyq(:) * bmnc(:, curS)
      !! dB/dv
      dBdVmns(:, curS) = xn_nyq(:) * bmnc(:, curS)
      IF (corVMEC .AND. useScaledB) THEN
        !! scaled dB/du
        dBdUsmns(:, curS) = -xm_nyq(:) * absBsmnc(:, curS)
        !! scaled dB/dv
        dBdVsmns(:, curS) = xn_nyq(:) * absBsmnc(:, curS)
      ELSE
        !! scaled dB/du
        dBdUsmns(:, curS) = dBdUmns(:, curS)
        !! scaled dB/dv
        dBdVsmns(:, curS) = dBdVmns(:, curS)
      END IF
    END DO

    !! store data on r-z-phi mesh
    !! initialize axis mode numbers
    DO i = 1, nAxis
      axisN(i) = i - 1
    END DO
    axisM(:) = 0
    WRITE(*,*)'...DONE.'

  END SUBROUTINE PreCalcData

! ----------------------------------------------------------------------
!> store mapping for Euterpe
! ----------------------------------------------------------------------

  SUBROUTINE EuterpeMapping

    REAL(KIND = hp), DIMENSION(rGridPoints, zGridPoints, radius) :: thRZ
    REAL(KIND = hp), DIMENSION(:, :, :), ALLOCATABLE :: mapRZ
    REAL(KIND = hp), DIMENSION(0:(rGridPoints+1), 0:(zGridPoints+1)) :: bRZ
    REAL(KIND = hp), DIMENSION(0:(rGridPoints+1), 0:(zGridPoints+1)) :: bRRZ
    REAL(KIND = hp), DIMENSION(0:(rGridPoints+1), 0:(zGridPoints+1)) :: bZRZ
    INTEGER, DIMENSION(rGridPoints, zGridPoints) :: minSV, minThV
    REAL(KIND = hp), DIMENSION(3) :: gradS, gradT, gradU
    REAL(KIND = hp), DIMENSION(3) :: gradB, bgradB
    REAL(KIND = hp), DIMENSION(3) :: rotBs, rotB, rotBPerp
    REAL(KIND = hp), DIMENSION(radius) :: dAxis, dRdUV, dRdVV, dZdUV, dZdVV
    REAL(KIND = hp), DIMENSION(radius) :: gradSRV, gradSZV, gradSPV, gradTRV
    REAL(KIND = hp), DIMENSION(radius) :: gradTZV, gradTPV
    REAL(KIND = hp), DIMENSION(radius) :: sqrtGV
    REAL(KIND = hp), DIMENSION(radius) :: rotBparAbsV
    REAL(KIND = hp), DIMENSION(radius) :: dbRdPV, dbZdPV
    REAL(KIND = hp), DIMENSION(radius) :: absBV, dPresdSV
    REAL(KIND = hp), DIMENSION(radius) :: theta, bUV, bVV, lambda, theta2
    REAL(KIND = hp), DIMENSION(radius) :: gB3V, bRV, bZV, thetaSC
    REAL(KIND = hp), DIMENSION(radius) :: lineR, lineZ, minDV
    REAL(KIND = hp), DIMENSION(mn_mode) :: argV, xnv, sinV, cosV
    REAL(KIND = hp), DIMENSION(mn_mode_nyq) :: argVn, xnvn, sinVn, cosVn
    REAL(KIND = hp), DIMENSION(8) :: pMin
    INTEGER, DIMENSION(radius) :: minTV
    CHARACTER(LEN = 255) :: outName
    REAL(KIND = hp) :: h, r, z, rS, zS, rAxis, zAxis, rzS
    REAL(KIND = hp) :: rzTheta, cylPhi, dist, dLambdadU, dLambdadV
    REAL(KIND = hp) :: d2LambdadU2, d2LambdadUdV, d2LambdadV2
    REAL(KIND = hp) :: d2LambdadUdS, d2LambdadVdS
    REAL(KIND = hp) :: dRdS, d2RdUdS, d2RdVdS
    REAL(KIND = hp) :: dRdU, d2RdU2, d2RdUdV, dRdV, d2RdV2
    REAL(KIND = hp) :: dZdS, d2ZdUdS, d2ZdVdS
    REAL(KIND = hp) :: dZdU, d2ZdU2, d2ZdUdV, dZdV, d2ZdV2
    REAL(KIND = hp) :: gradSR, gradSZ, gradSP, gradTR, gradTZ, gradTP
    REAL(KIND = hp) :: absB, bU, bV, bR, bP, bZ
    REAL(KIND = hp) :: dBdSs, dBdS, dBdUs, dBdU, dBdVs, dBdV
    REAL(KIND = hp) :: dbUdSs, dbUdUs, dbUdVs, dbVdSs, dbVdUs, dbVdVs
    REAL(KIND = hp) :: divB, rotBparAbs
    REAL(KIND = hp) :: sqrtG, dsqrtGdS, dsqrtGdU, dsqrtGdV
    REAL(KIND = hp) :: gsu, gsv, guu, guv, gvv, dgsudU, dgsudV, dgsvdU
    REAL(KIND = hp) :: dgsvdV, dguudS, dguudV, dguvdS, dguvdU, dguvdV, dgvvdS
    REAL(KIND = hp) :: dgvvdU, diotadS, rho, sRho
    REAL(KIND = hp) :: gBS, gBU, gBV
    REAL(KIND = hp) :: dbRdS, dbRdU, dbRdV, dbRdR, dbRdZ, dbRdP
    REAL(KIND = hp) :: dbZdS, dbZdU, dbZdV, dbZdR, dbZdZ, dbZdP
    REAL(KIND = hp) :: dMinR, dMaxR, dMinZ, dMaxZ, s, th
    REAL(KIND = hp) :: phiMinR, phiMaxR, phiMinZ, phiMaxZ
    REAL(KIND = hp) :: dLambdadS, iota, psi, dPresdS
    REAL(KIND = hp) :: gB1, gB2, gB3, mDist, rSt, zSt
    REAL(KIND = hp) :: h1, h2, r1, r2, z1, z2, h1m, h2m
    REAL(KIND = hp) :: s_tor
    REAL(KIND = hp), DIMENSION(radius) :: f_pol
    INTEGER :: i, j, k, l, m, curS, ioError, i0, is, ie, is1, mRad
    INTEGER :: minS, maxS, minTh, nequ_b, dir
    LOGICAL :: recalc

    IF (me_rank == 0) THEN
      WRITE(*, FMT = "(A)") "Calculate physical values ..."
    END IF

    !! calculate real number of flux surfaces
    IF ((nFlux < 0) .OR. (nFlux > radius)) nFlux = radius

    IF (me_rank == 0) THEN
      WRITE(*, "(2(A, I0), A)") "  Generating inner polar grid of size ",&
           nInnerS, "x", nInTheta, "!"
    END IF

    !! store iota profile
    IF (me_rank == 0) THEN
      WRITE(outName, FMT = '(A, "/iota.equ")') TRIM(outDir)
      OPEN(1, ACTION = "WRITE", FILE = TRIM(outName), FORM = "UNFORMATTED",&
           IOSTAT = ioError, POSITION = "REWIND", STATUS = "REPLACE")
      IF (ioError /= 0) THEN
        PRINT *, "Could not open output file ", TRIM(outName), "!"
        CALL EXIT(6)
      END IF
      WRITE(1) nFluxVMEC
      DO curS = 1, nFluxVMEC
        IF (useRho) THEN
          WRITE(1) SQRT(sVal(curS)), iotas(curS)
        ELSE
          WRITE(1) sVal(curS), iotas(curS)
        END IF
      END DO
      CLOSE(1)
    END IF

    !! calculate poloidal flux out of iota
    f_pol(1) = 0
    DO curS = 2, radius
      IF (useRho) THEN
        f_pol(curS) = f_pol(curS-1) + 0.5 * (iotas(curS - 1) + iotas(curS)) *&
             (SQRT(sVal(curS)) - SQRT(sVal(curS - 1)))
      ELSE
        f_pol(curS) = f_pol(curS-1) + 0.5 * (iotas(curS - 1) + iotas(curS)) *&
             (sVal(curS) - sVal(curS - 1))
      END IF
    END DO
    !! scale back from VMEC to real value of toroidal flux
    f_pol(:) = f_pol(:) * phipf(:) * REAL(signgs, hp) * TwoPi / REAL(nfp, hp)

    WRITE(outName, FMT = '(A, "/fluxes.equ")') TRIM(outDir)
    OPEN(1, ACTION = "WRITE", FILE = TRIM(outName), FORM = "UNFORMATTED",&
         IOSTAT = ioError, POSITION = "REWIND", STATUS = "REPLACE")
    IF (ioError /= 0) THEN
      PRINT *, "Could not open output file ", TRIM(outName), "!"
      STOP
    END IF
    !! write number of points
    WRITE(1) radius
    !! store poloidal flux
    DO curS = 1, radius
      IF (useRho) THEN
        s_tor = SQRT(sVal(curS))
      ELSE
        s_tor = sVal(curS)
      END IF
      !! scale back from VMEC to real value of toroidal flux
      WRITE(1) s_tor, s_tor * phipf(curS) * REAL(signgs, hp) * TwoPi /&
           REAL(nfp, hp), f_pol(curS)
    END DO
    CLOSE(1)

    !! initialize r-z mesh points
    DO j = 0, rGridPoints + 1
      !! set R value
      rMesh(j) = lowBR + (highBR - lowBR) * REAL(j - 1, hp) /&
           REAL(rGridPoints - 1, hp)
    END DO
    DO k = 0, zGridPoints + 1
      !! set z value
      zMesh(k) = lowBZ + (highBZ - lowBZ) * REAL(k - 1, hp) /&
           REAL(zGridPoints - 1, hp)
    END DO

    !! allocate memory for storing mapping from (s, theta) to (r, z)
    ALLOCATE(mapRZ(nPointsU, nFlux, 2))

    !! calculate minimal distance of points of the outer most flux surface
    !! with the selected box and the corresponding phi cut
    dMinR = highBR - lowBR
    phiMinR = 0
    dMaxR = highBR - lowBR
    phiMaxR = 0
    dMinZ = highBZ - lowBZ
    phiMinZ = 0
    dMaxZ = highBZ - lowBZ
    phiMaxZ = 0

    !! generate mesh
    !! phi cuts
    !! generate only selected cuts
    DO i = startCut, endCut
      IF (MODULO(i - startCut, me_size) == me_rank) THEN
        !! show current cut
        IF (me_size > 1) THEN
          PRINT '(I0, A, I0)', i, ' on ', me_rank
        ELSE
          PRINT '(I0)', i
        END IF
      ELSE
        CYCLE
      END IF
      !! calculate phi angle
      cylPhi = TwoPi * (i - 1) / nPhiCuts

      !! do correction for configuration with multiple field periods
      cylPhi = cylPhi / REAL(nfp, hp)
      xnv(:) = xn(:) * cylPhi
      xnvn(:) = xn_nyq(:) * cylPhi

      !! calculate rAxis and zAxis from inner flux surface
      rAxis = 0
      zAxis = 0
      DO j = 1, nPointsU
        th = TwoPi * REAL(j - 1, hp) / REAL(nPointsU - 1, hp)
        argV(:) = xm(:) * th - xnv(:)
        sinV(:) = SIN(argV(:))
        cosV(:) = COS(argV(:))
        !! calculate axis
        rAxis = rAxis + SUM(rmnc(:, 1) * cosV(:))
        zAxis = zAxis + SUM(zmns(:, 1) * sinV(:))
      END DO
      rAxis = rAxis / nPointsU
      zAxis = zAxis / nPointsU

      !! open output file
      WRITE(outName, FMT = '(A, "/vmec.equ.", I0.4)') TRIM(outDir), i - 1
      OPEN(1, ACTION = "WRITE", FILE = TRIM(outName), FORM = "UNFORMATTED",&
           IOSTAT = ioError, POSITION = "REWIND", STATUS = "REPLACE")
      IF (ioError /= 0) THEN
        PRINT *, "Could not open output file ", TRIM(outName), "!"
        CALL EXIT(7)
      END IF
      !! write version number
      WRITE(1) version
      !! write number of data
      IF (useArcTan) THEN
        nequ_b = 29
      ELSE
        nequ_b = 27
      END IF
      IF (gyroperp) THEN
        nequ_b = nequ_b + 6
      END IF
      WRITE(1) nequ_b
      !! first write some information about the current run
      !! boundaries of selected box
      WRITE(1) lowBR, highBR, lowBZ, highBZ
      !! basic data
      WRITE(1) nfp, rGridPoints, zGridPoints, nPhiCuts
      WRITE(1) nFlux, nTheta
      WRITE(1) bmnc(1, 1), volume
      WRITE(1) useRho
      WRITE(1) useSquared
      WRITE(1) nInnerS, nInTheta

      IF (debug .AND. (me_rank == 0)) THEN
        WRITE(outName, FMT = '("fort.", I0, ".16")') i
        OPEN(116, ACTION = "WRITE", FILE = TRIM(outName), FORM = "FORMATTED",&
             IOSTAT = ioError, POSITION = "REWIND", STATUS = "REPLACE")
        IF (ioError /= 0) THEN
          PRINT *, "Could not open output file ", TRIM(outName), "!"
          CALL EXIT(8)
        END IF
        !! write header
        WRITE(116, FMT = "(A)") "#s, theta(g), phi, r, z, theta"
      END IF
      !! store mesh from PEST coordinates
      !! flux surfaces
      DO curS = 1, nFlux
        s = sVal(curS)
        !! calculate s for writing
        IF (useRho) THEN
          sRho = SQRT(s)
        ELSE
          sRho = s
        END IF
        !! theta
        DO j = 1, nTheta
          th = TwoPi * REAL(j - 1, hp) / REAL(nTheta - 1, hp)
          theta(:) = th
          !! calculate fourier arguments
          argVn(:) = xm_nyq(:) * th - xnvn(:)
          sinVn(:) = SIN(argVn(:))
          cosVn(:) = COS(argVn(:))
          !! R
          r = IntExtPol3(radius, sVal(1:radius),&
               ValsOnFluxSurf(radius, mn_mode_nyq,&
               cylRmnc(:, 1:radius), cosVn(:)), s)
          !! z
          z = IntExtPol3(radius, sVal(1:radius),&
               ValsOnFluxSurf(radius, mn_mode_nyq,&
               cylZmns(:, 1:radius), sinVn(:)), s)
          !! write s, theta, phi, r, z
          WRITE(1) sRho, th, cylPhi, r, z
          IF (debug .AND. (me_rank == 0)) THEN
            IF (curS <= nFluxVMEC) THEN
              WRITE(116, "(6ES13.5)") sRho, th, cylPhi, r, z, SUM(bmnc(:, curS) *&
                   cosVn(:))
            ELSE
              WRITE(116, "(5ES13.5)") sRho, th, cylPhi, r, z
            END IF
          END IF
          !! compute minimal distances
          IF (curS == nFluxVMEC) THEN
            !! distance to lowBR
            IF (r - lowBR < dMinR) THEN
              dMinR = r - lowBR
              phiMinR = cylPhi
            END IF
            !! distance to highBR
            IF (highBR - r < dMaxR) THEN
              dMaxR = highBR - r
              phiMaxR = cylPhi
            END IF
            !! distance to lowBZ
            IF (z - lowBZ < dMinZ) THEN
              dMinZ = z - lowBZ
              phiMinZ = cylPhi
            END IF
            !! distance to highBZ
            IF (highBZ - z < dMaxZ) THEN
              dMaxZ = highBZ - z
              phiMaxZ = cylPhi
            END IF
          END IF
        END DO
      END DO
      IF (debug) CLOSE(116)

      !! initialize mapping
      mapRZ(:, :, :) = 0
      IF (debug .AND. (me_rank == 0)) THEN
        WRITE(outName, FMT = '("fort.", I0, ".19")') i
        OPEN(119, ACTION = "WRITE", FILE = TRIM(outName), FORM = "FORMATTED",&
             IOSTAT = ioError, POSITION = "REWIND", STATUS = "REPLACE")
        IF (ioError /= 0) THEN
          PRINT *, "Could not open output file ", TRIM(outName), "!"
          CALL EXIT(9)
        END IF
      END IF
      !! flux surfaces
      DO curS = 1, nFlux
        !! theta
        DO j = 1, nPointsU
          th = TwoPi * REAL(j - 1, hp) / REAL(nPointsU - 1, hp)
          argV(:) = xm(:) * th - xnv(:)
          sinV(:) = SIN(argV(:))
          cosV(:) = COS(argV(:))
          !! store (s,u)->(r,z) mapping
          mapRZ(j, curS, 1) = SUM(rmnc(:, curS) * cosV(:))
          mapRZ(j, curS, 2) = SUM(zmns(:, curS) * sinV(:))
          IF (debug .AND. (me_rank == 0)) THEN
            WRITE(119, FMT = "(2ES12.4)") mapRZ(j, curS, :)
          END IF
        END DO
        IF (debug .AND. (me_rank == 0)) WRITE(119, FMT = *)
      END DO
      IF (debug .AND. (me_rank == 0)) CLOSE(119)

      !! set output name
      IF (debug .AND. (me_rank == 0)) THEN
        WRITE(outName, FMT = '("fort.", I0, ".17")') i
        OPEN(117, ACTION = "WRITE", FILE = TRIM(outName), FORM = "FORMATTED",&
             IOSTAT = ioError, POSITION = "REWIND", STATUS = "REPLACE")
        IF (ioError /= 0) THEN
          PRINT *, "Could not open output file ", TRIM(outName), "!"
          CALL EXIT(10)
        END IF
      END IF

      !! precalc absB
      !! R
      DO j = 0, rGridPoints + 1
        rS = rMesh(j)
        !! z
        DO k = 0, zGridPoints + 1
          zS = zMesh(k)

          !! get points on the flux surfaces along a ray from the center
          !! through the selected point
          CALL CalcRay(rS, zS)

          !! store theta for performance reasons
          IF ((j > 0) .AND. (j <= rGridPoints) .AND.&
               (k > 0) .AND. (k <= zGridPoints)) THEN
            thRZ(j, k, :) = theta(:)
            !! store nearest point
            minSV(j, k) = minS
            minThV(j, k) = minTh
          END IF

          !! interpolate |B|, b_R, b_z along ray
          CALL IntAlongRay(bRZ(j, k), bRRZ(j, k), bZRZ(j, k))

          IF (debug .AND. (me_rank == 0)) THEN
            WRITE(117, FMT = "(I6, 5ES12.4)") i, rS, zS, bRZ(j, k),&
                 bRRZ(j, k), bZRZ(j, k)
          END IF
        END DO
      END DO
      IF (debug .AND. (me_rank == 0)) CLOSE(117)

      IF (debug .AND. (me_rank == 0)) THEN
        WRITE(outName, FMT = '("fort.", I0, ".13")') i
        OPEN(113, ACTION = "WRITE", FILE = TRIM(outName), FORM = "FORMATTED",&
             IOSTAT = ioError, POSITION = "REWIND", STATUS = "REPLACE")
        IF (ioError /= 0) THEN
          PRINT *, "Could not open output file ", TRIM(outName), "!"
          CALL EXIT(11)
        END IF
      END IF

      IF (debug .AND. (me_rank == 0)) THEN
        WRITE(outName, FMT = '("fort.", I0, ".14")') i
        OPEN(114, ACTION = "WRITE", FILE = TRIM(outName), FORM = "FORMATTED",&
             IOSTAT = ioError, POSITION = "REWIND", STATUS = "REPLACE")
        IF (ioError /= 0) THEN
          PRINT *, "Could not open output file ", TRIM(outName), "!"
          CALL EXIT(12)
        END IF
        !! write header
        IF (useArcTan) THEN
          WRITE(114, FMT = "(A)") "#rS, zS, cylPhi, rzS, rzTheta, absB, bR, bZ,&
               & bP, sqrtG, gradSR, gradSZ, gradSP, gradTR, gradTZ, gradTP,&
               & bgradBR, bgradBZ, bgradBP, divB, rotBparAbs, rotBperpR,&
               & rotBperpZ, rotBperpP, gradBR, gradBZ, gradBP, sinRZTheta,&
               & cosRZTheta, dbRdR, dbRdZ, dbRdP, dbZdR, dbZdZ, dbZdP"
        ELSE
          WRITE(114, FMT = "(A)") "#rS, zS, cylPhi, rzS, rzTheta, absB, bR, bZ,&
               & bP, sqrtG, gradSR, gradSZ, gradSP, gradTR, gradTZ, gradTP,&
               & bgradBR, bgradBZ, bgradBP, divB, rotBparAbs, rotBperpR,&
               & rotBperpZ, rotBperpP, gradBR, gradBZ, gradBP, dbRdS, dbRdU,&
               & dbRdV, dbZdS, dbZdU, dbZdV"
        END IF
        WRITE(outName, FMT = '("fort.", I0, ".141")') i
        OPEN(1141, ACTION = "WRITE", FILE = TRIM(outName), FORM = "FORMATTED",&
             IOSTAT = ioError, POSITION = "REWIND", STATUS = "REPLACE")
        IF (ioError /= 0) THEN
          PRINT *, "Could not open output file ", TRIM(outName), "!"
          CALL EXIT(12)
        END IF
        !! write header
        WRITE(1141, FMT = "(A)") "#r, z, s, theta, gradTR, gradTZ, gradTP, rho^2"
      END IF

      !! generate grid for inner torus using polar-like coordinates
      !! use geometric theta
      !! define direction in arrays
      dir = -1
      IF (useSquared) dir = 1
      DO j = 1, nInTheta
        th = TwoPi * REAL(j - 1, hp) / REAL(nInTheta - 1, hp)
        thetaSC(:nInnerS) = th
        CALL FindInterSectionRZNewton(MIN(dir, dir * nInnerS), MAX(dir, dir *&
             nInnerS), rAxis, zAxis, COS(th), SIN(th), cylPhi,&
             thetaSC(:nInnerS), loopUp = (dir == 1))
        !! calculate lambda
        DO curS = 1, nInnerS
          lambda(curS) = SinTransFullMesh(mn_mode, lmns(:, dir * curS),&
               xm(:), xn(:), thetaSC(curS), cylPhi)
        END DO
        !! inner flux surfaces
        DO curS = 1, nInnerS
          !! calculate s for writing
          IF (useRho) THEN
            sRho = SQRT(sVal(dir * curS))
          ELSE
            sRho = sVal(dir * curS)
          END IF
          !! get interpolated values
          CALL GetInterpolatesPolar(dir * curS, thetaSC(curS))
          !! calculate rho^2
          rho = (r - rAxis) * (r - rAxis) + (z - zAxis) * (z - zAxis)
          !! scale gradT for better interpolation in euterpe
          IF (useRho) THEN
            sqrtG = 2 * sRho * sqrtG / dLambdadU
            gradS(:) = gradS(:) / (2 * sRho)
          END IF
          gradT(:) = gradT(:) * rho
          gradS(:) = gradS(:) * rho
          !! write only neccessary data
          WRITE(1) r, z, sRho, MODULO(thetaSC(curS) + lambda(curS), TwoPi),&
               gradS(:), gradT(:), sqrtG, SQRT(rho)
          !! check if write data for inner polar grid
          IF (debug .AND. (me_rank == 0)) THEN
            WRITE(1141, FMT = "(12ES18.10)") r, z, sRho,&
                 MODULO(thetaSC(curS) + lambda(curS), TwoPi), gradS(:),&
                 gradT(:), sqrtG, SQRT(rho)
          END IF
        END DO
        IF (debug .AND. (me_rank == 0)) THEN
          WRITE(1141, FMT = *)
        END IF
      END DO
      IF (debug .AND. (me_rank == 0)) THEN
        CLOSE(1141)
      END IF

      !! R
      DO j = 1, rGridPoints
        rS = rMesh(j)
        !! z
        DO k = 1, zGridPoints
          zS = zMesh(k)
          !! use precalculated theta values
          theta(:) = thRZ(j, k, :)
          !! ray starting point
          minS = minSV(j, k)
          minTh = minThV(j, k)
          rSt = mapRZ(minTh, minS, 1)
          zSt = mapRZ(minTh, minS, 2)
          !! distance to the ray starting point
          dist = SQRT((rS - rSt)**2 + (zS - zSt)**2)

          !! set inner distance to 0
          dAxis(minS) = 0
          !! set inner lambda to 0
          lambda(minS) = 0
          i0 = minS
          mRad = minS + 10
          !! calculate distance of cutting points to the magnetic axis
          DO curS = minS + 1, minS + 10
            dAxis(curS) = DistanceToPoint(curS, theta(curS), cylPhi, rSt, zSt)
            lambda(curS) = SinTransFullMesh(mn_mode, lmns(:, curS), xm(:),&
                 xn(:), theta(curS), cylPhi)
            !! search position if grid point
            IF (dist > dAxis(curS)) i0 = curS
            !! check if calculated theta values are consistent with the
            !! distance to the magnetic axis
            IF (curS > minS) THEN
              !! check strict monotonic behaviour
              IF (dAxis(curS) < dAxis(curS - 1)) THEN
                mRad = curS - 1
                PRINT "(A, I0, A, I0, A, F0.3, A, F0.3)",&
                     "Maximal number of used flux surfaces is set to ",&
                     mRad - minS, " instead of ", 10, " for r=", rS,&
                     " and z=", zS
                is = MAX(minS + 1, i0 - 3)
                ie = MIN(minS + 10, i0 + 3)
                PRINT "(A, I0, A, I0, A, 10ES11.3)", "theta(", is, ":", ie,&
                     ")= ", theta(is:ie)
                PRINT "(A, I0, A, I0, A, 10ES11.3)", "dist(", is, ":", curS,&
                     ")= ", dAxis(is:curS)
                EXIT
              END IF
            END IF
          END DO
          !! interpolate or extrapolate the data along the line
          !! calculate only nearest neighbouring points along the line
          is = MAX(minS + 1, i0 - 3)
          ie = MIN(mRad, i0 + 3)
          !! help indexes
          is1 = MAX(3, is)
          !! maximal inter- or extrapolated distance of a point
          mDist = constFac * DistanceToPoint(mRad, theta(mRad), cylPhi, rSt,&
               zSt)
          !! interpolate along constant theta
          IF (useConstTheta .AND. (dist > dAxis(minS + 1)) .AND.&
               (dist < mDist)) THEN
            !! interpolate or extrapolate theta along the line
            th = IntExtPol3(ie - is + 1, dAxis(is:ie), theta(is:ie),&
                 dist)
            lineR(minS) = rSt
            lineZ(minS) = zSt
            recalc = .TRUE.
            mRad = minS + 10
            !! recalculate distance along constant theta path
            DO curS = minS + 1, minS + 10
              lambda(curS) = SinTransFullMesh(mn_mode, lmns(:, curS), xm(:),&
                   xn(:), th, cylPhi)
              lineR(curS) = CosTransFullMeshF(mn_mode, rmnc(:, curS), xm(:),&
                   xn(:), th, cylPhi, recalc)
              lineZ(curS) = SinTransFullMeshF(mn_mode, zmns(:, curS), xm(:),&
                   xn(:), th, cylPhi, recalc)
              recalc = .FALSE.
              dAxis(curS) = dAxis(curS - 1) + SQRT((lineR(curS) -&
                   lineR(curS - 1))**2 + (lineZ(curS) - lineZ(curS - 1))**2)
            END DO
            !! maximal inter- or extrapolated distance of a point
            mDist = constFac * dAxis(minS + 10)
            dist = dAxis(i0 - 1) + SQRT((rS - lineR(i0 - 1))**2 + (zS -&
                 lineZ(i0 - 1))**2)
          END IF
          !! s coordinate
          IF (dist < dAxis(minS + 1)) THEN
            !! linear interpolation
            rzS = dist * (sVal(minS + 1) - sVal(minS)) / dAxis(minS + 1) +&
                 sVal(minS)
          ELSE
            rzS = IntExtPol3(ie - is + 1, dAxis(is:ie), sVal(is:ie), dist)
          END IF
          IF (rzS < 0) THEN
            PRINT "(A, ES12.4, 3(A, F0.3))", "adjust ", rzS, " to 0 at ", rS,&
                 ",", zS, ",", cylPhi
            rzS = 0
          END IF
          !! set all values to constant values along the line outside a
          !! specific region
          dist = MIN(dist, mDist)

          !! theta coordinate
          rzTheta = MODULO(IntExtPol3(ie - is + 1, dAxis(is:ie),&
               theta(is:ie) + lambda(is:ie), dist), TwoPi)

          !! calculate values along the line
          DO curS = is, ie
            !! calculate s for writing
            IF (useRho) THEN
              sRho = SQRT(sVal(curS))
              !! get the interpolated values
              CALL GetInterpolates(curS, SQRT(rzS), theta(curS))
            ELSE
              sRho = sVal(curS)
              !! get the interpolated values
              CALL GetInterpolates(curS, rzS, theta(curS))
            END IF

            !! store values
            !! R component of grad s
            gradSRV(curS) = gradS(1)
            !! z component of grad s
            gradSZV(curS) = gradS(2)
            !! phi component of grad s
            gradSPV(curS) = gradS(3)
            !! R component of grad theta
            gradTRV(curS) = gradT(1) * sVal(curS)
            !! z component of grad theta
            gradTZV(curS) = gradT(2) * sVal(curS)
            !! phi component of grad theta
            gradTPV(curS) = gradT(3) * sVal(curS)
            !! |B|
            absBV(curS) = absB
            !! sqrtG
            sqrtGV(curS) = sqrtG / dLambdadU
            !! |rot B parallel|
            rotBparAbsV(curS) = rotBparAbs
            IF (gyroperp) THEN
              !! db_R/dphi
              dbRdPV(curS) = dbRdP
              !! db_z/dphi
              dbZdPV(curS) = dbZdP
            END IF
            !! dpres/ds
            dPresdSV(curS) = GetRadialDerivative(radius, curS, presf(:))
            bUV(curS) = bU
            bVV(curS) = bV
            dRdUV(curS) = dRdU
            dRdVV(curS) = dRdV
            dZdUV(curS) = dZdU
            dZdVV(curS) = dZdV
            gB3V(curS) = gBV
            IF (debug .AND. (me_rank == 0) .AND. (curS == i0)) THEN
              IF (gyroperp) THEN
                WRITE(113, FMT = "(21ES13.5)") r, z, cylPhi, sRho,&
                     theta(curS) + lambda(curS), h, absB, dBdS, dBdU, dRdS,&
                     dRdU, dZdS, dZdU, gBS, gBU, gBV, dbRdP, dbZdP, bR, bZ,&
                     rotBparAbs
              ELSE
                WRITE(113, FMT = "(23ES13.5)") r, z, cylPhi, sRho,&
                     theta(curS) + lambda(curS), h, absB, dBdS, dBdU, dRdS,&
                     dRdU, dZdS, dZdU, gBS, gBU, gBV, bR, bZ, bP, rotBparAbs,&
                     rotB(:)
              END IF
            END IF
          END DO

          !! interpolate all values
          !! |B|
          absB = IntExtPol3(ie - is + 1, dAxis(is:ie), absBV(is:ie), dist)
          !! consistency check for |B|
          IF (absB < 0) THEN
            PRINT "(A, 3ES12.4, A)", "|B| less than zero at ", rS, zS,&
                 cylPhi, "!!"
          END IF
          !! sqrtG
          sqrtG = IntExtPol3(ie - is + 1, dAxis(is:ie), sqrtGV(is:ie), dist)
          !! R component of grad s
          gradSR = IntExtPol3(ie - is + 1, dAxis(is:ie), gradSRV(is:ie), dist)
          !! z component of grad s
          gradSZ = IntExtPol3(ie - is + 1, dAxis(is:ie), gradSZV(is:ie), dist)
          !! phi component of grad s
          gradSP = IntExtPol3(ie - is + 1, dAxis(is:ie), gradSPV(is:ie), dist)
          !! R component of grad theta
          gradTR = IntExtPol3(ie - is + 1, dAxis(is:ie), gradTRV(is:ie),&
               dist) / rzS
          !! z component of grad theta
          gradTZ = IntExtPol3(ie - is + 1, dAxis(is:ie), gradTZV(is:ie),&
               dist) / rzS
          !! phi component of grad theta
          gradTP = IntExtPol3(ie - is + 1, dAxis(is:ie), gradTPV(is:ie),&
               dist) / rzS

          bU = IntExtPol3(ie - is + 1, dAxis(is:ie), bUV(is:ie), dist)
          bV = IntExtPol3(ie - is + 1, dAxis(is:ie), bVV(is:ie), dist)

          dRdU = IntExtPol3(ie - is + 1, dAxis(is:ie), dRdUV(is:ie), dist)
          dRdV = IntExtPol3(ie - is + 1, dAxis(is:ie), dRdVV(is:ie), dist)
          dZdU = IntExtPol3(ie - is + 1, dAxis(is:ie), dZdUV(is:ie), dist)
          dZdV = IntExtPol3(ie - is + 1, dAxis(is:ie), dZdVV(is:ie), dist)

          !! calculate B_R = dR/du * B^u + dR/dv * B^v
          bR = (bU * dRdU + bV * dRdV) / absB
          !! calculate B_z = dz/du * B^u + dz/dv * B^v
          bZ = (bU * dZdU + bV * dZdV) / absB
          !! calculate B_phi
          bP = SQRT(1._hp - bR * bR - bZ * bZ) * SIGN(1._hp, bV * rS)

          !! use precalculated magnetic field for the derivative of B
          gB1 = .125_hp * REAL(rGridPoints - 1, hp) * (bRZ(j + 1, k + 1) -&
               bRZ(j - 1, k + 1) + 2 * (bRZ(j + 1, k) - bRZ(j - 1, k)) +&
               bRZ(j + 1, k - 1) - bRZ(j - 1, k - 1)) / (highBR - lowBR)
          gB2 = .125_hp * REAL(zGridPoints - 1, hp) * (bRZ(j + 1, k + 1) -&
               bRZ(j + 1, k - 1) + 2 * (bRZ(j, k + 1) - bRZ(j, k - 1)) +&
               bRZ(j - 1, k + 1) - bRZ(j - 1, k - 1)) / (highBZ - lowBZ)
          gB3 = IntExtPol3(ie - is1 + 1, dAxis(is1:ie), gB3V(is1:ie), dist)
          gradB(1) = gB1
          gradB(2) = gB2
          gradB(3) = gB3
          !! calculate b x grad B / |B|
          bgradB(:) = CrossL((/bR, bZ, bP/), gradB(:)) / absB
          !! calculate divergence of B
          divB = -DOT_PRODUCT((/bR, bZ, bP/), gradB(:)) / absB
          !! |rot B parallel|
          rotBparAbs = IntExtPol3(ie - is + 1, dAxis(is:ie),&
               rotBparAbsV(is:ie), dist)
          IF (gyroperp) THEN
            !! use precalculated magnetic field for the derivative of b_R
            dbRdR = .125_hp * REAL(rGridPoints - 1, hp) * (bRRZ(j + 1, k + 1) -&
                 bRRZ(j - 1, k + 1) + 2 * (bRRZ(j + 1, k) - bRRZ(j - 1, k)) +&
                 bRRZ(j + 1, k - 1) - bRRZ(j - 1, k - 1)) / (highBR - lowBR)
            dbRdZ = .125_hp * REAL(zGridPoints - 1, hp) * (bRRZ(j + 1, k + 1) -&
                 bRRZ(j + 1, k - 1) + 2 * (bRRZ(j, k + 1) - bRRZ(j, k - 1)) +&
                 bRRZ(j - 1, k + 1) - bRRZ(j - 1, k - 1)) / (highBZ - lowBZ)
            !! db_R/dphi
            dbRdP = IntExtPol3(ie - is + 1, dAxis(is:ie), dbRdPV(is:ie), dist)
            !! use precalculated magnetic field for the derivative of b_Z
            dbZdR = .125_hp * REAL(rGridPoints - 1, hp) * (bZRZ(j + 1, k + 1) -&
                 bZRZ(j - 1, k + 1) + 2 * (bZRZ(j + 1, k) - bZRZ(j - 1, k)) +&
                 bZRZ(j + 1, k - 1) - bZRZ(j - 1, k - 1)) / (highBR - lowBR)
            dbZdZ = .125_hp * REAL(zGridPoints - 1, hp) * (bZRZ(j + 1, k + 1) -&
                 bZRZ(j + 1, k - 1) + 2 * (bZRZ(j, k + 1) - bZRZ(j, k - 1)) +&
                 bZRZ(j - 1, k + 1) - bZRZ(j - 1, k - 1)) / (highBZ - lowBZ)
            !! db_z/dphi
            dbZdP = IntExtPol3(ie - is + 1, dAxis(is:ie), dbZdPV(is:ie), dist)
          END IF
          !! dpres/ds
          dPresdS = IntExtPol3(ie - is + 1, dAxis(is:ie), dPresdSV(is:ie),&
               dist)
          !! calc rotB perpendicular from pressure gradient because it is
          !! stable over the whole domain
          !! rotBperp = mu0 (b x (\nabla s)) p' / |B|
          rotBPerp(:) = mu0 * CrossL((/bR, bZ, bP/),&
               (/gradSR, gradSZ, gradSP / rS/)) * dPresdS / absB

          !! scale only s dependent data for representating them in rho
          IF (useRho) THEN
            rzS = SQRT(rzS)
            sqrtG = 2 * rzS * sqrtG
            IF (rzS > 1E-9) THEN
              gradSR = gradSR / (2 * rzS)
              gradSZ = gradSZ / (2 * rzS)
              gradSP = gradSP / (2 * rzS)
            ELSE
              gradSR = 0
              gradSZ = 0
              gradSP = 0
            END IF
          END IF

          !! write all data
          IF (useArcTan) THEN
            !! check if write data for perpendicular gyro ring
            IF (gyroperp) THEN
              WRITE(1) rS, zS, cylPhi, rzS, rzTheta, absB, bR, bZ, bP, sqrtG,&
                   gradSR, gradSZ, gradSP, gradTR, gradTZ, gradTP, bgradB(:),&
                   divB, rotBparAbs, rotBperp(:), gradB(:), SQRT(rzS) *&
                   SIN(rzTheta), SQRT(rzS) * COS(rzTheta), dbRdR, dbRdZ,&
                   dbRdP, dbZdR, dbZdZ, dbZdP
            ELSE
              WRITE(1) rS, zS, cylPhi, rzS, rzTheta, absB, bR, bZ, bP, sqrtG,&
                   gradSR, gradSZ, gradSP, gradTR, gradTZ, gradTP, bgradB(:),&
                   divB, rotBparAbs, rotBperp(:), gradB(:), SQRT(rzS) *&
                   SIN(rzTheta), SQRT(rzS) * COS(rzTheta)
            END IF
          ELSE
            !! check if write data for perpendicular gyro ring
            IF (gyroperp) THEN
              WRITE(1) rS, zS, cylPhi, rzS, rzTheta, absB, bR, bZ, bP, sqrtG,&
                   gradSR, gradSZ, gradSP, gradTR, gradTZ, gradTP, bgradB(:),&
                   divB, rotBparAbs, rotBperp(:), gradB(:), dbRdR, dbRdZ,&
                   dbRdP, dbZdR, dbZdZ, dbZdP
            ELSE
              WRITE(1) rS, zS, cylPhi, rzS, rzTheta, absB, bR, bZ, bP, sqrtG,&
                   gradSR, gradSZ, gradSP, gradTR, gradTZ, gradTP, bgradB(:),&
                   divB, rotBparAbs, rotBperp(:), gradB(:)
            END IF
          END IF
          IF (debug .AND. (me_rank == 0)) THEN
            IF (useArcTan) THEN
              !! check if write data for perpendicular gyro ring
              IF (gyroperp) THEN
                WRITE(114, FMT = "(35ES18.10)") rS, zS, cylPhi, rzS, rzTheta,&
                     absB, bR, bZ, bP, sqrtG, gradSR, gradSZ, gradSP, gradTR,&
                     gradTZ, gradTP, bgradB(:), divB, rotBparAbs, rotBperp(:),&
                     gradB(:), SQRT(rzS) * SIN(rzTheta), SQRT(rzS) *&
                     COS(rzTheta), dbRdR, dbRdZ, dbRdP, dbZdR, dbZdZ, dbZdP
              ELSE
                WRITE(114, FMT = "(29ES18.10)") rS, zS, cylPhi, rzS, rzTheta,&
                     absB, bR, bZ, bP, sqrtG, gradSR, gradSZ, gradSP, gradTR,&
                     gradTZ, gradTP, bgradB(:), divB, rotBparAbs, rotBperp(:),&
                     gradB(:), SQRT(rzS) * SIN(rzTheta), SQRT(rzS) *&
                     COS(rzTheta)
              END IF
            ELSE
              !! check if write data for perpendicular gyro ring
              IF (gyroperp) THEN
                WRITE(114, FMT = "(33ES18.10)") rS, zS, cylPhi, rzS, rzTheta,&
                     absB, bR, bZ, bP, sqrtG, gradSR, gradSZ, gradSP, gradTR,&
                     gradTZ, gradTP, bgradB(:), divB, rotBparAbs, rotBperp(:),&
                     gradB(:), dbRdR, dbRdZ, dbRdP, dbZdR, dbZdZ, dbZdP
              ELSE
                WRITE(114, FMT = "(27ES18.10)") rS, zS, cylPhi, rzS, rzTheta,&
                     absB, bR, bZ, bP, sqrtG, gradSR, gradSZ, gradSP, gradTR,&
                     gradTZ, gradTP, bgradB(:), divB, rotBparAbs, rotBperp(:),&
                     gradB(:)
              END IF
            END IF
          END IF
        END DO
      END DO
      IF (debug .AND. (me_rank == 0)) THEN
        CLOSE(113)
        CLOSE(114)
      END IF
      CLOSE(1)
    END DO

    pMin(1) = dMinR
    pMin(2) = phiMinR
    pMin(3) = dMaxR
    pMin(4) = phiMaxR
    pMin(1) = dMinZ
    pMin(2) = phiMinZ
    pMin(3) = dMaxZ
    pMin(4) = phiMaxZ

    IF (me_rank == 0) THEN
      WRITE(*, *) " done!"
      PRINT *, "minimal distances to the box boundaries"
      PRINT '(A, ES11.4, " ( phi =", ES11.4, ")")', "to minR: ", pMin(1:2)
      PRINT '(A, ES11.4, " ( phi =", ES11.4, ")")', "to maxR: ", pMin(3:4)
      PRINT '(A, ES11.4, " ( phi =", ES11.4, ")")', "to minZ: ", pMin(5:6)
      PRINT '(A, ES11.4, " ( phi =", ES11.4, ")")', "to maxZ: ", pMin(7:8)
    END IF

  CONTAINS

    SUBROUTINE CalcRay(rS, zS)

!> cylinder coordinates
      REAL(KIND = hp), INTENT(IN) :: rS, zS

      INTEGER :: curS

      !! look for nearest point
      !! initialization
      minDV(1) = 1D39
      minTV(1) = -1
      minTh = 1
      !! all points from the (s,u) -> (r,z) mapping
      DO l = 1, nFlux
        minDV(l) = (mapRZ(minTh, l, 1) - rS) *&
             (mapRZ(minTh, l, 1) - rS) + (mapRZ(minTh, l, 2) - zS) *&
             (mapRZ(minTh, l, 2) - zS)
        minTV(l) = minTh
        DO m = 1, nPointsU
          !! L2-distance to current grid point
          h1 = (mapRZ(m, l, 1) - rS) * (mapRZ(m, l, 1) - rS)
          IF (h1 < minDV(l)) THEN
            h2 = (mapRZ(m, l, 2) - zS) * (mapRZ(m, l, 2) - zS)
            !! check for minimum
            IF (h1 + h2 < minDV(l)) THEN
              !! store minimal position
              minTV(l) = m
              minTh = m
              minDV(l) = h1 + h2
            END IF
          END IF
        END DO
      END DO
      !! use 10 points around the nearest point for interpolation
      minS = MINLOC(minDV(:), 1)
      minTh = minTV(minS)
      maxS = MIN(radius, minS + 5)
      minS = MAX(maxS - 10, 2)
      !! set the starting point of the ray
      IF (minS > 3) THEN
        rSt = mapRZ(minTh, minS, 1)
        zSt = mapRZ(minTh, minS, 2)
      ELSE
        minS = 1
        maxS = MIN(radius, minS + 10)
        rSt = rAxis
        zSt = zAxis
      END IF

      !! distance to the ray starting point
      dist = SQRT((rS - rSt)**2 + (zS - zSt)**2)
      !! look for all cuts of the flux surfaces with a the ray through
      !! the current grid points for theta
      DO curS = minS + 1, minS + 10
        !! search segment of the flux surface which is cutted by the ray
        m = 0
        h2m = 1E39_hp
        !! all points of the approximation of the flux surface
        DO l = 1, nPointsU - 1
          !! first point of segment
          r1 = mapRZ(l, curS, 1)
          z1 = mapRZ(l, curS, 2)
          !! second point of segment
          r2 = mapRZ(l + 1, curS, 1)
          z2 = mapRZ(l + 1, curS, 2)
          !! calculate crossing point by
          !!  (rSt)        (rS - rSt)   (r1)        (r2 - r1)
          !!  (   ) + h1 * (        ) = (  ) + h2 * (       )
          !!  (zSt)        (zS - zSt)   (z1)        (z2 - z1)
          h = (rS - rSt) * (z2 - z1) - (zS - zSt) * (r2 - r1)
          h1 = (rSt * (z1 - zS) + r1 * (zS - zSt) + rS * (zSt - z1)) / h
          !! check if crossing point inside the segment and on the
          !! correct side
          IF ((h1 >= 0) .AND. (h1 <= 1)) THEN
            h2 = (rSt * (z1 - z2) + r1 * (z2 - zSt) + r2 * (zSt - z1)) / h
            IF ((h2 > 0) .AND. (h2 < h2m)) THEN
              m = l
              h1m = h1
              h2m = h2
            END IF
          END IF
        END DO
        !! calculate approximated value of u
        h = TwoPi * REAL(m - 1, hp) / REAL(nPointsU - 1, hp)
        theta(curS) = h + TwoPi * h1m / REAL(nPointsU - 1, hp)
        !! store a second point if secant method is used
        IF (useSecant) THEN
          !! check if the approximated point is a vertex
          IF ((h1m < .5_hp) .AND. (h1m > 0)) THEN
            theta2(curS) = h
          ELSE
            theta2(curS) = TwoPi * REAL(m, hp) / REAL(nPointsU - 1, hp)
          END IF
        END IF
      END DO
      IF (useSecant) THEN
        CALL FindInterSectionRZSec(minS + 1, minS + 10, rSt, zSt, rS - rSt,&
             zS - zSt, cylPhi, theta(:), theta2(:))
      ELSE
        CALL FindInterSectionRZNewton(minS + 1, minS + 10, rSt, zSt,&
             rS - rSt, zS - zSt, cylPhi, theta(:))
      END IF

    END SUBROUTINE CalcRay

    SUBROUTINE IntAlongRay(b, bR, bZ)

!> resulting values for |B|, b_R, b_z
      REAL(KIND = hp), INTENT(out) :: b, bR, bZ

      INTEGER :: curS
      REAL(KIND = hp) :: iota, psi, sqrtG, dLambdadU, dLambdadV
      REAL(KIND = hp) :: dRdU, dRdV, dZdU, dZdV, bU, bV

      !! set inner distance to 0
      dAxis(minS) = 0
      i0 = minS
      mRad = minS + 10
      !! calculate distance of cutting points to the magnetic axis
      DO curS = minS + 1, minS + 10
        dAxis(curS) = DistanceToPoint(curS, theta(curS), cylPhi, rSt, zSt)
        !! search position if grid point
        IF (dist > dAxis(curS)) i0 = curS
        !! check if calculated theta values are consistent with the
        !! distance to the ray starting point
        IF (curS > minS) THEN
          !! check strict monotonic behaviour
          IF (dAxis(curS) < dAxis(curS - 1)) THEN
            mRad = curS - 1
            EXIT
          END IF
        END IF
      END DO
      !! interpolate or extrapolate the data along the line
      !! calculate only nearest neighbouring points along the line
      is = MAX(minS + 1, i0 - 3)
      ie = MIN(mRad, i0 + 3)
      !! maximal inter- or extrapolated distance of a point
      mDist = constFac * DistanceToPoint(mRad, theta(mRad), cylPhi, rSt,&
           zSt)
      !! interpolate along constant theta
      IF (useConstTheta .AND. (dist > dAxis(minS + 1)) .AND.&
           (dist < mDist)) THEN
        !! interpolate or extrapolate theta along the line
        theta(:) = IntExtPol3(ie - is + 1, dAxis(is:ie), theta(is:ie),&
             dist)
        lineR(minS) = rSt
        lineZ(minS) = zSt
        recalc = .TRUE.
        mRad = minS + 10
        !! recalculate distance along constant theta path
        DO curS = minS + 1, minS + 10
          lineR(curS) = CosTransFullMeshF(mn_mode, rmnc(:, curS), xm(:),&
               xn(:), theta(1), cylPhi, recalc)
          lineZ(curS) = SinTransFullMeshF(mn_mode, zmns(:, curS), xm(:),&
               xn(:), theta(1), cylPhi, recalc)
          recalc = .FALSE.
          dAxis(curS) = dAxis(curS - 1) + SQRT((lineR(curS) -&
               lineR(curS - 1))**2 + (lineZ(curS) - lineZ(curS - 1))**2)
        END DO
        !! maximal inter- or extrapolated distance of a point
        mDist = constFac * dAxis(minS + 10)
        dist = dAxis(i0 - 1) + SQRT((rS - lineR(i0 - 1))**2 + (zS -&
             lineZ(i0 - 1))**2)
      END IF

      !! calculate values along the line
      DO curS = is, ie
        argV(:) = xm(:) * theta(curS) - xnv(:)
        sinV(:) = SIN(argV(:))
        cosV(:) = COS(argV(:))
        argVn(:) = xm_nyq(:) * theta(curS) - xnvn(:)
        cosVn(:) = COS(argVn(:))
        !! iota
        iota = iotas(curS)
        !! psi
        psi = phipf(curS)
        !! sqrtG from VMEC data
        sqrtG = SUM(gmnc(:, curS) * cosVn(:))
        !! 1 + dlambda/du
        dLambdadU = 1 + SUM(lUmnc(:, curS) * cosV(:))
        !! dlambda/dv
        dLambdadV = SUM(lVmnc(:, curS) * cosV(:))
        !! dR/du
        dRdU = SUM(dRdUmns(:, curS) * sinV(:))
        !! dR/dv
        dRdV = SUM(dRdVmns(:, curS) * sinV(:))
        !! dz/du
        dZdU = SUM(dZdUmnc(:, curS) * cosV(:))
        !! dz/dv
        dZdV = SUM(dZdVmnc(:, curS) * cosV(:))
        IF (corVMEC) THEN
          !! B^u = (iota - dlambda/dv) * psi / sqrtG
          bU = (iota - dLambdadV) * psi / sqrtG
          !! B^v = (1 + dlambda/du) * psi / sqrtG
          bV = dLambdadU * psi / sqrtG
        ELSE
          bU = SUM(bsupumnc(:, curS) * cosVn(:))
          bV = SUM(bsupvmnc(:, curS) * cosVn(:))
        END IF
        !! calculate B_R = dR/du * B^u + dR/dv * B^v
        bR = bU * dRdU + bV * dRdV
        !! calculate B_z = dz/du * B^u + dz/dv * B^v
        bZ = bU * dZdU + bV * dZdV
        !! |B|
        absBV(curS) = CosTransFullMesh(mn_mode_nyq, bmnc(:, curS),&
             xm_nyq(:), xn_nyq(:), theta(curS), cylPhi)
        !! calculate b_R
        bRV(curS) = bR / absBV(curS)
        !! calculate b_z
        bZV(curS) = bZ / absBV(curS)
      END DO
      !! |B|
      b = IntExtPol3(ie - is + 1, dAxis(is:ie), absBV(is:ie),&
           MIN(dist, mDist))
      !! b_R
      bR = IntExtPol3(ie - is + 1, dAxis(is:ie), bRV(is:ie),&
           dist)
      !! b_z
      bZ = IntExtPol3(ie - is + 1, dAxis(is:ie), bZV(is:ie),&
           MIN(dist, mDist))

    END SUBROUTINE IntAlongRay

    SUBROUTINE GetInterpolates(curS, rzS, theta)

      !! number of flux surface
      INTEGER, INTENT(IN) :: curS

      !! s coordinate
      REAL(KIND = hp), INTENT(IN) :: rzS

      !! value of theta
      REAL(KIND = hp), INTENT(IN) :: theta

      argV(:) = xm(:) * theta - xnv(:)
      sinV(:) = SIN(argV(:))
      cosV(:) = COS(argV(:))
      argVn(:) = xm_nyq(:) * theta - xnvn(:)
      sinVn(:) = SIN(argVn(:))
      cosVn(:) = COS(argVn(:))
      !! R
      r = SUM(rmnc(:, curS) * cosV(:))
      !! z
      z = SUM(zmns(:, curS) * sinV(:))
      !! dlambda/ds
      dLambdadS = GetRadialDerivativeF(radius, curS, mn_mode,&
           lmns(:, 1:), sinV(:))
      !! 1 + dlambda/du
      dLambdadU = 1 + SUM(lUmnc(:, curS) * cosV(:))
      !! dlambda/dv
      dLambdadV = SUM(lVmnc(:, curS) * cosV(:))
      !! d^2lambda / (ds du)
      d2LambdadUdS = GetRadialDerivativeF(radius, curS, mn_mode,&
           lUmnc(:, 1:), cosV(:))
      !! d^2lambda / (ds dv)
      d2LambdadVdS = GetRadialDerivativeF(radius, curS, mn_mode,&
           lVmnc(:, 1:), cosV(:))
      !! d^2lambda / (du^2)
      d2LambdadU2 = SUM(lU2mns(:, curS) * sinV(:))
      !! d^2lambda / (du dv)
      d2LambdadUdV = SUM(lUVmns(:, curS) * sinV(:))
      !! d^2lambda / (dv^2)
      d2LambdadV2 = SUM(lV2mns(:, curS) * sinV(:))
      !! dR/dS, dz/dS and second partial derivates from fourier series if
      !! available
      IF (corVMEC) THEN
        !! dR/ds
        dRdS = SUM(dRdSmnc(:, curS) * cosV(:))
        !! d^2R/(du ds)
        d2RdUdS = SUM(d2RdUdSmns(:, curS) * sinV(:))
        !! d^2R/(dv ds)
        d2RdVdS = SUM(d2RdVdSmns(:, curS) * sinV(:))
        !! dz/ds
        dZdS = SUM(dZdSmns(:, curS) * sinV(:))
        !! d^2z/(du ds)
        d2ZdUdS = SUM(d2ZdUdSmnc(:, curS) * cosV(:))
        !! d^2z/(dv ds)
        d2ZdVdS = SUM(d2ZdVdSmnc(:, curS) * cosV(:))
      ELSE
        !! dR/ds from VMEC data
        dRdS = GetRadialDerivativeF(radius, curS, mn_mode,&
             rmnc(:, 1:), cosV(:))
        !! d^2R/(du ds)
        d2RdUdS = GetRadialDerivativeF(radius, curS, mn_mode,&
             dRdUmns(:, 1:), sinV(:))
        !! d^2R/(dv ds)
        d2RdVdS = GetRadialDerivativeF(radius, curS, mn_mode,&
             dRdVmns(:, 1:), sinV(:))
        !! dz/ds from VMEC data
        dZdS = GetRadialDerivativeF(radius, curS, mn_mode,&
             zmns(:, 1:), sinV(:))
        !! d^2z/(du ds)
        d2ZdUdS = GetRadialDerivativeF(radius, curS, mn_mode,&
             dZdUmnc(:, 1:), cosV(:))
        !! d^2z/(dv ds)
        d2ZdVdS = GetRadialDerivativeF(radius, curS, mn_mode,&
             dZdVmnc(:, 1:), cosV(:))
      END IF
      !! dR/du
      dRdU = SUM(dRdUmns(:, curS) * sinV(:))
      !! d2R/du2
      d2RdU2 = SUM(d2RdU2mnc(:, curS) * cosV(:))
      !! d2R/(du dv)
      d2RdUdV = SUM(d2RdUdVmnc(:, curS) * cosV(:))
      !! dR/dv
      dRdV = SUM(dRdVmns(:, curS) * sinV(:))
      !! d2R/dv2
      d2RdV2 = SUM(d2RdV2mnc(:, curS) * cosV(:))
      !! dz/du
      dZdU = SUM(dZdUmnc(:, curS) * cosV(:))
      !! d2z/du2
      d2ZdU2 = SUM(d2ZdU2mns(:, curS) * sinV(:))
      !! d2z/(du dv)
      d2ZdUdV = SUM(d2ZdUdVmns(:, curS) * sinV(:))
      !! dz/dv
      dZdV = SUM(dZdVmnc(:, curS) * cosV(:))
      !! d2z/du2
      d2ZdV2 = SUM(d2ZdV2mns(:, curS) * sinV(:))
      !! calculate jacobian and partial derivatives from following
      !! expression
      !!             (dR dz   dR dz)
      !!   sqrtG = R (-- -- - -- --)
      !!             (du ds   ds du)
      !! sqrtG
      h = dRdU * dZdS - dZdU * dRdS
      !! sqrtG from VMEC data
      sqrtG = SUM(gmnc(:, curS) * cosVn(:))
      !! correction of h if curS == 1
      IF (h == 0.0) h = sqrtG / r
      !! dsqrtG/ds
      IF (corVMEC) THEN
        dsqrtGdS = SUM(dgdSmnc(:, curS) * cosVn(:))
      ELSE
        dsqrtGdS = GetRadialDerivativeF(radius, curS,&
             mn_mode_nyq, gmnc(:, 1:), cosVn(:))
      END IF
      !! dsqrtG/du
      dsqrtGdU = SUM(-xm_nyq(:) * gmnc(:, curS) * sinVn(:))
      !! dsqrtG/dv
      dsqrtGdV = SUM(xn_nyq(:) * gmnc(:, curS) * sinVn(:))
      !! |B|
      absB = SUM(bmnc(:, curS) * cosVn(:))
      IF (corVMEC) THEN
        IF (useScaledB .AND. (rzS <= errInt(1))) THEN
          dBdSs = SUM(dBdSsmnc(:, curS) * cosVn(:))
          !! d|B|/ds = d(|Bs|/sqrtG)/ds
          !!         = (d|Bs|/ds * sqrtG - |Bs| * dsqrtG/ds) / sqrtG^2
          dBdS = (signgs * dBdSs - absB * dsqrtGdS) / sqrtG
        ELSE
          !! dB/ds
          dBdS = SUM(dBdSmnc(:, curS) * cosVn(:))
        END IF
      ELSE
        dBdS = GetRadialDerivativeF(radius, curS, mn_mode_nyq,&
             bmnc(:, 1:), cosVn(:))
      END IF

      !! dB/du
      dBdU = SUM(dBdUmns(:, curS) * sinVn(:))
      IF (corVMEC .AND. useScaledB) THEN
        dBdUs = SUM(dBdUsmns(:, curS) * sinVn(:))
        !! d|B|/du = d(|Bs|/sqrtG)/du
        !!         = (d|Bs|/du * sqrtG - |Bs| * dsqrtG/du) / sqrtG^2
        dBdU = (signgs * dBdUs - dsqrtGdU * absB) / sqrtG
      END IF

      !! dB/dv
      dBdV = SUM(dBdVmns(:, curS) * sinVn(:))
      IF (corVMEC .AND. useScaledB) THEN
        dBdVs = SUM(dBdVsmns(:, curS) * sinVn(:))
        !! d|B|/dv = d(|Bs|/sqrtG)/dv
        !!         = (d|Bs|/dv * sqrtG - |Bs| * dsqrtG/dv) / sqrtG^2
        dBdV = (signgs * dBdVs - dsqrtGdV * absB) / sqrtG
      END IF

      !! diota/ds
      dIotadS = GetRadialDerivative(radius, curS, iotas(1:))
      !! iota
      iota = iotas(curS)
      !! psi
      psi = phipf(curS)

      !! calculate grad s and grad u from inverse of
      !!        (dR/ds dR/du dR/dv)
      !!        (dz/ds dz/du dz/dv)
      !!        (  0     0     1  )
      !!
      !!          (         -dz/du          )
      !! grad s = (          dR/du          )/(dR/du dz/ds - dz/du dR/ds)
      !!          (dR/dv dz/du - dR/du dz/dv)
      gradS(1) = -dZdU / h
      gradS(2) = dRdU / h
      gradS(3) = (dRdV * dZdU - dRdU * dZdV) / h
      !!          (          dz/ds          )
      !! grad u = (         -dR/ds          )/(dR/du dz/ds - dz/du dR/ds)
      !!          (dR/ds dz/dv - dR/dv dz/ds)
      gradU(1) = dZdS / h
      gradU(2) = -dRdS / h
      gradU(3) = (dRdS * dZdV - dRdV * dZdS) / h
      !! calculate grad theta = (1 + dlambda/du) * grad u +
      !!                   dlambda/ds * grad s + dlambda/dv * grad v
      gradT(:) = dLambdadU * gradU(:) + dLambdadS * gradS(:)
      gradT(3) = gradT(3) + dLambdadV

      IF (corVMEC) THEN
        !! B^u = (iota - dlambda/dv) * psi / sqrtG
        bU = (iota - dLambdadV) * psi / sqrtG
        !! B^v = (1 + dlambda/du) * psi / sqrtG
        bV = dLambdadU * psi / sqrtG
      ELSE
        bU = SUM(bsupumnc(:, curS) * cosVn(:))
        bV = SUM(bsupvmnc(:, curS) * cosVn(:))
      END IF

      !! calculate B_R = dR/du * B^u + dR/dv * B^v
      bR = (bU * dRdU + bV * dRdV) / absB
      !! calculate B_z = dz/du * B^u + dz/dv * B^v
      bZ = (bU * dZdU + bV * dZdV) / absB
      !! calculate B_phi
      !! bP = bV * r
      !! for consistency we use B_phi = SQRT(|B|^2-B_R^2-B_Z^2) instead of
      !! the formula B_phi = B^v * R since the later gives slightly
      !! inconsistent results near the axis
      bP = SQRT(1._hp - bR * bR - bZ * bZ) * SIGN(1._hp, bV * r)

      !! calculate rot B
      !! metric
      !!          dR dR   dz dz
      !! g_{su} = -- -- + -- --
      !!          ds du   ds du
      gsu = dRdS * dRdU + dZdS * dZdU
      !!          dR dR   dz dz
      !! g_{sv} = -- -- + -- --
      !!          ds dv   ds dv
      gsv = dRdS * dRdV + dZdS * dZdV
      !!          dR dR   dz dz
      !! g_{uu} = -- -- + -- --
      !!          du du   du du
      guu = dRdU * dRdU + dZdU * dZdU
      !!          dR dR   dz dz
      !! g_{uv} = -- -- + -- --
      !!          du dv   du dv
      guv = dRdU * dRdV + dZdU * dZdV
      !!          dR dR   dz dz
      !! g_{vv} = -- -- + -- -- + r^2
      !!          dv dv   dv dv
      gvv = dRdV * dRdV + dZdV * dZdV + r * r
      !! for consistency one could set g_{vv} with
      !! gvv = (gsv * gsv * guu - 2 * gsu * gsv * guv + gss * guv * guv +&
      !!     sqrtG * sqrtG) / (gss * guu - gsu * gsu)
      !! christoffel symbols
      !! dg_{su}/du
      dgsudU = d2RdUdS * dRdU + dRdS * d2RdU2 + d2ZdUdS * dZdU + dZdS *&
           d2ZdU2
      !! dg_{su}/dv
      dgsudV = d2RdVdS * dRdU + dRdS * d2RdUdV + d2ZdVdS * dZdU + dZdS *&
           d2ZdUdV
      !! dg_{sv}/du
      dgsvdU = d2RdUdS * dRdV + dRdS * d2RdUdV + d2ZdUdS * dZdV + dZdS *&
           d2ZdUdV
      !! dg_{sv}/dv
      dgsvdV = d2RdVdS * dRdV + dRdS * d2RdV2 + d2ZdVdS * dZdV + dZdS *&
           d2ZdV2
      !! dg_{uu}/ds
      dguudS = 2 * (dRdU * d2RdUdS + dZdU * d2ZdUdS)
      !! dg_{uu}/dv
      dguudV = 2 * (dRdU * d2RdUdV + dZdU * d2ZdUdV)
      !! dg_{uv}/ds
      dguvdS = d2RdUdS * dRdV + dRdU * d2RdVdS + d2ZdUdS * dZdV + dZdU *&
           d2ZdVdS
      !! dg_{uv}/du
      dguvdU = d2RdU2 * dRdV + dRdU * d2RdUdV + d2ZdU2 * dZdV + dZdU *&
           d2ZdUdV
      !! dg_{uv}/dv
      dguvdV = d2RdUdV * dRdV + dRdU * d2RdV2 + d2ZdUdV * dZdV + dZdU *&
           d2ZdV2
      !! dg_{vv}/ds
      dgvvdS = 2 * (dRdV * d2RdVdS + dZdV * d2ZdVdS + r * dRdS)
      !! dg_{vv}/du
      dgvvdU = 2 * (dRdV * d2RdUdV + dZdV * d2ZdUdV + r * dRdU)
      !!       (diota   d^2lambda)   (       dLambda) dsqrtG
      !!       (----- - ---------) - (iota - -------) ------ / sqrtG
      !! dB^u  ( ds       ds dv  )   (         dv   )   ds
      !! ---- =----------------------------------------------------- psi
      !!  ds                        sqrtG
      !!
      dbUdSs = (dIotadS - d2LambdadVdS - (iota - dLambdadV) *&
           dsqrtGdS / sqrtG) * psi
      !!          d^2lambda   (       dlambda) dsqrtG
      !!        - --------- - (iota - -------) ------ / sqrtG
      !! dB^u      du dv      (         dv   )   du
      !! ---- = --------------------------------------------- psi
      !!  du                       sqrtG
      !!
      dbUdUs = -(d2LambdadUdV + (iota - dLambdadV) * dsqrtGdU&
           / sqrtG) * psi
      !!          d^2lambda   (       dlambda) dsqrtG
      !!        - --------- - (iota - -------) ------ / sqrtG
      !! dB^u       dv^2      (         dv   )   dv
      !! ---- = --------------------------------------------- psi
      !!  dv                       sqrtG
      !!
      dbUdVs = -(d2LambdadV2 + (iota - dLambdadV) * dsqrtGdV / sqrtG) *&
           psi
      !!       d^2lambda   (    dlambda) dsqrtG
      !!       --------- - (1 + -------) ------ / sqrtG
      !! dB^v    ds du     (      du   )   ds
      !! ---- =---------------------------------------- psi
      !!  ds                    sqrtG
      !!
      dbVdSs = (d2LambdadUdS - dLambdadU * dsqrtGdS / sqrtG) * psi
      !!        d^2lambda   (    dlambda) dsqrtG
      !!        --------- - (1 + -------) ------ / sqrtG
      !! dB^v     du^2      (      du   )   du
      !! ---- = ---------------------------------------- psi
      !!  du                    sqrtG
      !!
      dbVdUs = (d2LambdadU2 - dLambdadU * dsqrtGdU / sqrtG) * psi
      !!        d^2lambda   (    dlambda) dsqrtG
      !!        --------- - (1 + -------) ------ / sqrtG
      !! dB^v     du dv     (      du   )   dv
      !! ---- = ---------------------------------------- psi
      !!  dv                     sqrtG
      !!
      dbVdVs = (d2LambdadUdV - dLambdadU * dsqrtGdV / sqrtG) * psi
      !!               (dB_v   dB_u)
      !!               (---- - ----)
      !!               ( du     dv )
      !!           1   (dB_s   dB_v)
      !! rot B = ----- (---- - ----)
      !!         sqrtG ( dv     ds )
      !!               (dB_u   dB_s)
      !!               (---- - ----)
      !!               ( ds     du )
      !! where
      !!  B_s = g_{su} B^u + g_{sv} B^v
      !!  B_u = g_{uu} B^u + g_{uv} B^v
      !!  B_v = g_{uv} B^u + g_{vv} B^v
      rotBs(1) = ((dguvdU - dguudV) * bU + (dgvvdU - dguvdV) * bV +&
           (guv * dbUdUs - guu * dbUdVs + gvv * dbVdUs - guv * dbVdVs) &
           / sqrtG) / sqrtG
      rotBs(2) = ((dgsudV - dguvdS) * bU + (dgsvdV - dgvvdS) * bV +&
           (gsu * dbUdVs - guv * dbUdSs + gsv * dbVdVs - gvv * dbVdSs) &
           / sqrtG) / sqrtG
      rotBs(3) = ((dguudS - dgsudU) * bU + (dguvdS - dgsvdU) * bV +&
           (guu * dbUdSs - gsu * dbUdUs + guv * dbVdSs - gsv * dbVdUs) &
           / sqrtG) / sqrtG
      !! convert to cylindrical coordinates
      rotB(1) = dRdS * rotBs(1) + dRdU * rotBs(2) + dRdV * rotBs(3)
      rotB(2) = dZdS * rotBs(1) + dZdU * rotBs(2) + dZdV * rotBs(3)
      rotB(3) = rotBs(3) * r
      !! get rotB parallel
      !! |rotBpar| = (b . (\nabla x B))
      rotBParAbs = DOT_PRODUCT((/bR, bZ, bP/), rotB(:))
      gBS = (dBdU * dZdS - dBdS * dZdU) / h
      gBU = (dBdS * dRdU - dBdU * dRdS) / h
      gBV = ((dBdS * (dRdV * dZdU - dRdU * dZdV) +&
           dBdU * (dRdS * dZdV - dRdV * dZdS)) / h + dBdV) / r

      IF (gyroperp) THEN
        !! dB_R   db^u dR   db^v dR       d^2 R       d^2 R
        !! ---- = ---- -- + ---- -- + b^u ----- + b^v -----
        !!  ds     ds  du    ds  dv       ds du       ds dv
        dbRdS = (dbUdSs * dRdU + dbVdSs * dRdV) / sqrtG + bU * d2RdUdS +&
             bV * d2RdVdS
        !! db_R   (dB_R       d|b|)
        !! ---- = (---- - b_R ----) / |b|
        !!  ds    ( ds         ds )
        dbRdS = (dbRdS - bR * dBdS) / absB

        !! dB_R   db^u dR   db^v dR       d^2R       d^2 R
        !! ---- = ---- -- + ---- -- + b^u ---- + b^v -----
        !!  du     du  du    du  dv       du^2       du dv
        dbRdU = (dbUdUs * dRdU + dbVdUs * dRdV) / sqrtG + bU * d2RdU2 +&
             bV * d2RdUdV
        !! db_R   (dB_R       d|b|)
        !! ---- = (---- - b_R ----) / |b|
        !!  du    ( du         du )
        dbRdU = (dbRdU - bR * dBdU) / absB

        !! dB_R   db^u dR   db^v dR       d^2 R       d^2R
        !! ---- = ---- -- + ---- -- + b^u ----- + b^v ----
        !!  dv     dv  du    dv  dv       du dv       dv^2
        dbRdV = (dbUdVs * dRdU + dbVdVs * dRdV) / sqrtG + bU * d2RdUdV +&
             bV * d2RdV2
        !! db_R   (dB_R       d|b|)
        !! ---- = (---- - b_R ----) / |b|
        !!  dv    ( dv         dv )
        dbRdV = (dbRdV - bR * dBdV) / absB

        !! dB_z   db^u dz   db^v dz       d^2 z       d^2 z
        !! ---- = ---- -- + ---- -- + b^u ----- + b^v -----
        !!  ds     ds  du    ds  dv       ds du       ds dv
        dbZdS = (dbUdSs * dZdU + dbVdSs * dZdV) / sqrtG + bU * d2ZdUdS +&
             bV * d2ZdVdS
        !! db_z   (dB_z       d|b|)
        !! ---- = (---- - b_z ----) / |b|
        !!  ds    ( ds         ds )
        dbZdS = (dbZdS - bZ * dBdS) / absB

        !! dB_z   db^u dz   db^v dz       d^2z       d^2 z
        !! ---- = ---- -- + ---- -- + b^u ---- + b^v -----
        !!  du     du  du    du  dv       du^2       du dv
        dbZdU = (dbUdUs * dZdU + dbVdUs * dZdV) / sqrtG + bU * d2ZdU2 +&
             bV * d2ZdUdV
        !! db_z   (dB_z       d|b|)
        !! ---- = (---- - b_z ----) / |b|
        !!  du    ( du         du )
        dbZdU = (dbZdU - bZ * dBdU) / absB

        !! dB_z   db^u dz   db^v dz       d^2 z       d^2z
        !! ---- = ---- -- + ---- -- + b^u ----- + b^v ----
        !!  dv     dv  du    dv  dv       du dv       dv^2
        dbZdV = (dbUdVs * dZdU + dbVdVs * dZdV) / sqrtG + bU * d2ZdUdV +&
             bV * d2ZdV2
        !! db_z   (dB_z       d|b|)
        !! ---- = (---- - b_z ----) / |b|
        !!  dv    ( dv         dv )
        dbZdV = (dbZdV - bZ * dBdV) / absB

        !! db_R   db_R  ds    db_R  du    db_R
        !! ---- = ---- ---- + ---- ---- + ----
        !! dphi    ds  dphi    du  dphi    dv
        dbRdP = dbRdS * gradS(3) + dbRdU * gradU(3) + dbRdV
        !! db_z   db_z  ds    db_z  du    db_z
        !! ---- = ---- ---- + ---- ---- + ----
        !! dphi    ds  dphi    du  dphi    dv
        dbZdP = dbZdS * gradS(3) + dbZdU * gradU(3) + dbZdV
      END IF

    END SUBROUTINE GetInterpolates

    SUBROUTINE GetInterpolatesPolar(curS, theta)

      !! number of flux surface
      INTEGER, INTENT(IN) :: curS

      !! value of theta
      REAL(KIND = hp), INTENT(IN) :: theta

      argV(:) = xm(:) * theta - xnv(:)
      sinV(:) = SIN(argV(:))
      cosV(:) = COS(argV(:))
      argVn(:) = xm_nyq(:) * theta - xnvn(:)
      cosVn(:) = COS(argVn(:))
      !! R
      r = SUM(rmnc(:, curS) * cosV(:))
      !! z
      z = SUM(zmns(:, curS) * sinV(:))
      !! dlambda/ds
      dLambdadS = GetRadialDerivativeFPolar(nInnerS, curS, mn_mode,&
           lmns(:, stInd:(stInd + nInnerS - 1)), sinV(:))
      !! 1 + dlambda/du
      dLambdadU = 1 + SUM(lUmnc(:, curS) * cosV(:))
      !! dlambda/dv
      dLambdadV = SUM(lVmnc(:, curS) * cosV(:))
      !! dR/dS, dz/dS
      IF (corVMEC) THEN
        !! dR/ds
        dRdS = SUM(dRdSmnc(:, curS) * cosV(:))
        !! dz/ds
        dZdS = SUM(dZdSmns(:, curS) * sinV(:))
      ELSE
        !! dR/ds from VMEC data
        dRdS = dir * GetRadialDerivativeFPolar(nInnerS, dir * curS, mn_mode,&
             rmnc(:, stInd:(stInd + nInnerS - 1)), cosV(:))
        !! dz/ds from VMEC data
        dZdS = dir * GetRadialDerivativeFPolar(nInnerS, dir * curS, mn_mode,&
             zmns(:, stInd:(stInd + nInnerS - 1)), sinV(:))
      END IF
      !! dR/du
      dRdU = SUM(dRdUmns(:, curS) * sinV(:))
      !! dR/dv
      dRdV = SUM(dRdVmns(:, curS) * sinV(:))
      !! dz/du
      dZdU = SUM(dZdUmnc(:, curS) * cosV(:))
      !! dz/dv
      dZdV = SUM(dZdVmnc(:, curS) * cosV(:))
      !! calculate jacobian and partial derivatives from following
      !! expression
      !!             (dR dz   dR dz)
      !!   sqrtG = R (-- -- - -- --)
      !!             (du ds   ds du)
      !! sqrtG
      h = dRdU * dZdS - dZdU * dRdS
      !! sqrtG from VMEC data
      sqrtG = SUM(gmnc(:, curS) * cosVn(:))

      !! calculate grad s and grad u from inverse of
      !!        (dR/ds dR/du dR/dv)
      !!        (dz/ds dz/du dz/dv)
      !!        (  0     0     1  )
      !!
      !!          (         -dz/du          )
      !! grad s = (          dR/du          )/(dR/du dz/ds - dz/du dR/ds)
      !!          (dR/dv dz/du - dR/du dz/dv)
      gradS(1) = -dZdU / h
      gradS(2) = dRdU / h
      gradS(3) = (dRdV * dZdU - dRdU * dZdV) / h
      !!          (          dz/ds          )
      !! grad u = (         -dR/ds          )/(dR/du dz/ds - dz/du dR/ds)
      !!          (dR/ds dz/dv - dR/dv dz/ds)
      gradU(1) = dZdS / h
      gradU(2) = -dRdS / h
      gradU(3) = (dRdS * dZdV - dRdV * dZdS) / h
      !! calculate grad theta = (1 + dlambda/du) * grad u +
      !!                   dlambda/ds * grad s + dlambda/dv * grad v
      gradT(:) = dLambdadU * gradU(:) + dLambdadS * gradS(:)
      gradT(3) = gradT(3) + dLambdadV

    END SUBROUTINE GetInterpolatesPolar

  END SUBROUTINE EuterpeMapping

! ----------------------------------------------------------------------
!> sort array started from largest value (only indices will be changed)
! ----------------------------------------------------------------------

  SUBROUTINE Sort(indices, values)

!> index array
    INTEGER, DIMENSION(:), INTENT(INOUT) :: indices

!> value array
    REAL(KIND = hp), DIMENSION(:), INTENT(IN) :: values

    INTEGER :: count, i, j

    count = SIZE(indices)

    DO WHILE (count > 1)
      count = count - 1
      DO i = 1, count
        IF (values(indices(i + 1)) - values(indices(i)) > 1E-12_hp) THEN
          j = indices(i)
          indices(i) = indices(i + 1)
          indices(i + 1) = j
          count = i
        END IF
      END DO
    END DO

  END SUBROUTINE Sort

! ----------------------------------------------------------------------
!> Smooth fourier components of r and z from VMEC output using a test
!> function and smoothing splines:
!>    spl(s) = t(s) (a + b s + c s^2 + d s^3),  with t(s) = s^(p/2)
!>
!> We use p = m if mode m < 2 and p = 2 for all other modes.
!> For the given "approximate" data one can apply an estimate of the
!> standard deviation (errInt - interval borders, errVal - values).
!> With a "weight" one can decide if the smoothed function should be
!> more smooth or more exact.
!>
!> To get the smoothing splines we minimize a function which is a sum of
!> three terms:
!> 1. approximation term
!>
!> \f$ \sum_{i=1}^{n-1} \big(\frac{t(s_{i+1})(a_i + b_i h_i + c_i h_i^2 +
!>        d_i h_i^3) - y_i}{errVar_i}\big)^2 \f$,
!>
!> 2. curvature term
!>
!> \f$ \sum_{i=1}^{n-1} d_i^2 \f$,
!>
!> 3. smoothing term
!>
!> \f$ \sum_{i=1}^{n-2} \alpha_i (a_i + b_i h_i + c_i h_i^2 + d_i h_i^3 -
!>        a_{i+1}) + \beta_i (b_i h_i + 2 c_i h_i + 3 d_i h_i^2 - b_{i+1}) +
!>        \gamma_i (c_i + 3 d_i h_i - c_{i+1}) \f$.
!>
!> The approximation term and the smoothing term will be weighted by
!> (1-weight) and weight.
!>
! ----------------------------------------------------------------------

  SUBROUTINE SmoothRZ

    REAL(KIND = hp), DIMENSION(:, :, :), ALLOCATABLE :: a
    REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: br, bz
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: perm
    REAL(KIND = hp), DIMENSION(nFluxVMEC) :: errVar
    REAL(KIND = hp), DIMENSION(stInd:radius - 1) :: h
    INTEGER, DIMENSION(MIN(0,stInd+1):radius - 1) :: map
    REAL(KIND = hp) :: lr, lz, l, h1, h2, h3
    INTEGER :: n, mode, i, j, k, aStat, aError, func
    INTEGER :: malpha, mbeta, mgamma

    IF (me_rank == 0) THEN
      WRITE(*, FMT = "(A)", ADVANCE = "no") "Smooth R and z ..."
    END IF

    !! number of parameters
    n = 7 * nFluxVMEC - 10

    !! allocate memory for arrays
    aError = 0
    ALLOCATE(a(0:4, n, n), perm(0:4, n), stat = aStat)
    aError = aError + aStat
    ALLOCATE(xr(mn_mode, n), xz(mn_mode, n), t(stInd:radius, 0:4, 0:3), stat = aStat)
    aError = aError + aStat
    ALLOCATE(br(n), bz(n), asymp(mn_mode), stat = aStat)
    aError = aError + aStat
    IF (aError /= 0) THEN
      PRINT *, "Allocation failure in subroutine SmoothRZ!"
      CALL EXIT(17)
    END IF

    !! marker to column number
    ma = 0
    mb = ma + nFluxVMEC - 1
    mc = mb + nFluxVMEC - 1
    md = mc + nFluxVMEC - 1
    malpha = md + nFluxVMEC - 1
    mbeta = malpha + nFluxVMEC - 2
    mgamma = mbeta + nFluxVMEC - 2

    !! help vector
    DO i = 1, nFluxVMEC - 1
      h(i) = sVal(i + 1) - sVal(i)
    END DO

    !! scale vector to smooth values on border
    errVar(:) = 1
    IF (useVar) THEN
      DO i = 1, nFluxVMEC
        IF (sVal(i) < errInt(1)) THEN
          errVar(i) = errVal(1)
        ELSE IF (sVal(i) < errInt(2)) THEN
          errVar(i) = errVal(2)
        ELSE
          errVar(i) = errVal(3)
        END IF
      END DO
    END IF

    !! initialize test function
    t(:, :, :) = 0
    t(:, 0, 0) = 1
    DO i = 2, nFluxVMEC
      !! t(s) = s^(1/2)
      t(i, 1, 0) = SQRT(sVal(i))
      !! t'(s) = 1/2 s^(-1/2)
      t(i, 1, 1) = .5_hp / SQRT(sVal(i))
      !! t''(s) = -1/4 s^(-3/2)
      t(i, 1, 2) = -.25_hp / (sVal(i) * SQRT(sVal(i)))
      !! t'''(s) = 3/8 s^(-5/2)
      t(i, 1, 3) = .375_hp / (sVal(i) * sVal(i) * SQRT(sVal(i)))

      !! t(s) = s
      t(i, 2, 0) = sVal(i)
      !! t'(s) = 1
      t(i, 2, 1) = 1
      !! t''(s) = 0
      t(i, 2, 2) = 0
      !! t'''(s) = 0
      t(i, 2, 3) = 0

      !! t(s) = s^(3/2)
      t(i, 3, 0) = SQRT(sVal(i)) * sVal(i)
      !! t'(s) = 3/2 s^(1/2)
      t(i, 3, 1) = 1.5_hp * SQRT(sVal(i))
      !! t''(s) = 3/4 s^(-1/2)
      t(i, 3, 2) =  .75_hp / SQRT(sVal(i))
      !! t'''(s) = -3/8 s^(-3/2)
      t(i, 3, 3) = -.375_hp / (sVal(i) * SQRT(sVal(i)))

      !! t(s) = s^2
      t(i, 4, 0) = sVal(i) * sVal(i)
      !! t'(s) = 2 s
      t(i, 4, 1) = 2 * sVal(i)
      !! t''(s) = 2
      t(i, 4, 2) = 2
      !! t'''(s) = 0
      t(i, 4, 3) = 0
    END DO

    !! initialize matrices for fitting
    a(:, :, :) = 0
    DO func = 0, 2
      DO j = 1, nFluxVMEC - 1
        !! \sum_{i=1}^{n-1} \big(\frac{t(s_{i+1})(a_i + b_i h_i +
        !!   c_i h_i^2 + d_i h_i^3) - y_i}{errVar_i}\big)^2 +
        h1 = 2 * (1 - weight) * t(j + 1, func, 0) /&
             (errVar(j + 1) * errVar(j + 1))
        a(func, ma + j, ma + j) = h1 * t(j + 1, func, 0)
        a(func, ma + j, mb + j) = a(func, ma + j, ma + j) * h(j)
        a(func, ma + j, mc + j) = a(func, ma + j, mb + j) * h(j)
        a(func, ma + j, md + j) = a(func, ma + j, mc + j) * h(j)

        a(func, mb + j, ma + j) = a(func, ma + j, mb + j)
        a(func, mb + j, mb + j) = a(func, ma + j, mc + j)
        a(func, mb + j, mc + j) = a(func, ma + j, md + j)
        a(func, mb + j, md + j) = a(func, ma + j, md + j) * h(j)

        a(func, mc + j, ma + j) = a(func, ma + j, mc + j)
        a(func, mc + j, mb + j) = a(func, mb + j, mc + j)
        a(func, mc + j, mc + j) = a(func, mb + j, md + j)
        a(func, mc + j, md + j) = a(func, mb + j, md + j) * h(j)

        !! curvature term
        !!
        !! \sum_{i=i+1}^{n-1} d_i^2
        a(func, md + j, ma + j) = a(func, ma + j, md + j)
        a(func, md + j, mb + j) = a(func, mb + j, md + j)
        a(func, md + j, mc + j) = a(func, mc + j, md + j)
        a(func, md + j, md + j) = a(func, mc + j, md + j) * h(j) +&
             2 * weight

        !! smoothing term
        !!
        !! \sum_{i=1}^{n-2} \alpha_i (a_i + b_i h_i + c_i h_i^2 +
        !!   d_i h_i^3 - a_{i+1}) + \beta_i (b_i + 2 c_i h_i +
        !!   3 d_i h_i^2 - b_{i+1}) + \gamma_i (c_i + 3 d_i h_i -
        !!   c_{i+1})
        IF (j < nFluxVMEC - 1) THEN
          a(func, ma + j, malpha + j) = 1

          a(func, mb + j, malpha + j) = h(j)
          a(func, mb + j, mbeta + j) = 1

          a(func, mc + j, malpha + j) = h(j) * h(j)
          a(func, mc + j, mbeta + j) = 2 * h(j)
          a(func, mc + j, mgamma + j) = 1

          a(func, md + j, malpha + j) = h(j) * h(j) * h(j)
          a(func, md + j, mbeta + j) = 3 * h(j) * h(j)
          a(func, md + j, mgamma + j) = 3 * h(j)

          a(func, malpha + j, ma + j) = 1
          a(func, malpha + j, ma + j + 1) = -1
          a(func, malpha + j, mb + j) = h(j)
          a(func, malpha + j, mc + j) = h(j) * h(j)
          a(func, malpha + j, md + j) = h(j) * h(j) * h(j)

          a(func, mbeta + j, mb + j) = 1
          a(func, mbeta + j, mb + j + 1) = -1
          a(func, mbeta + j, mc + j) = 2 * h(j)
          a(func, mbeta + j, md + j) = 3 * h(j) * h(j)

          a(func, mgamma + j, mc + j) = 1
          a(func, mgamma + j, mc + j + 1) = -1
          a(func, mgamma + j, md + j) = 3 * h(j)
        END IF
        IF (j > 1) THEN
          a(func, ma + j, malpha + j - 1) = -1
          a(func, mb + j, mbeta + j - 1) = -1
          a(func, mc + j, mgamma + j - 1) = -1
        END IF
      END DO
      !! initialize solution of system of equations
      CALL LUDecompose(n, a(func, :, :), perm(func, :))
    END DO

    !! do smoothing for every mode
    DO mode = 1, mn_mode
      IF (debug .AND. (me_rank == 0)) THEN
        WRITE(*, FMT = '(2I5)', ADVANCE = "NO") NINT(xm(mode)), NINT(xn(mode))
      END IF
      !! select test function
      IF (xm(mode) < 2) THEN
        func = INT(xm(mode))
      ELSE
        !! select 2 for m > 1
        func = 2
      END IF

      !! calculate right-hand side
      br(:) = 0
      bz(:) = 0
      DO j = 1, nFluxVMEC - 1
        !! approximation term
        h1 = 2 * (1 - weight) * t(j + 1, func, 0) /&
             (errVar(j + 1) * errVar(j + 1))

        !! 2 * (1 - w) * t(s_{j+1}) * y_{i+1} / rho_{j+1}^2
        br(ma + j) = h1 * rmnc(mode, j + 1)
        bz(ma + j) = h1 * zmns(mode, j + 1)

        !! 2 * (1 - w) * t(s_{j+1}) * h_j * y_{i+1} / rho_{j+1}^2
        br(mb + j) = (h1 * h(j)) * rmnc(mode, j + 1)
        bz(mb + j) = (h1 * h(j)) * zmns(mode, j + 1)

        !! 2 * (1 - w) * t(s_{j+1}) * h_j^2 * y_{i+1} / rho_{j+1}^2
        br(mc + j) = (h1 * h(j) * h(j)) * rmnc(mode, j + 1)
        bz(mc + j) = (h1 * h(j) * h(j)) * zmns(mode, j + 1)

        !! 2 * (1 - w) * t(s_{j+1}) * h_j^3 * y_{i+1} / rho_{j+1}^2
        br(md + j) = (h1 * h(j) * h(j) * h(j)) * rmnc(mode, j + 1)
        bz(md + j) = (h1 * h(j) * h(j) * h(j)) * zmns(mode, j + 1)
      END DO
      !! solve system of equation
      CALL LUBacksolve(n, a(func, :, :), perm(func, :), br(:))
      CALL LUBacksolve(n, a(func, :, :), perm(func, :), bz(:))
      !! store test function
      asymp(mode) = func
      !! store solution
      xr(mode, :) = br(:)
      xz(mode, :) = bz(:)
      IF (debug .AND. (me_rank == 0)) PRINT '(A, I0)', " => ", asymp(mode)
    END DO

    !! set new flux surfaces for better central resolution
    IF (mapEuterpe) THEN
      IF (useSquared) THEN
        sVal(:) = sVal(:) * sVal(:)
      ELSE
        sVal(stInd:-1) = sVal(stInd:-1) * sVal(stInd:-1)
      END IF
    END IF
    sVal(1) = 1E-9
    IF (mapEuterpe .AND. .NOT.(useSquared)) sVal(-1) = 1E-9
    IF (useRho .AND. mapEuterpe) THEN
      IF (.NOT.(useSquared)) sVal(0) = 0
      sVal(:) = sVal(:) * sVal(:)
    END IF

    !! calculate new distances to original flux surfaces and store flux label
    j = 1
    map(0) = 1
    DO i = 1, radius - 1
      !! check if last original flux surface is reached and test if inside a
      !! given original flux surface interval
      DO WHILE ((j < nFluxVMEC - 1) .AND. (ABS((nFluxVMEC - 1) * (sVal(i + 1) -&
           sValVMEC(j)) - .5) - .5 > 1E-8))
        j = j + 1
      END DO
      !! calculate distance
      h(i) = (sVal(i + 1) - sValVMEC(j))
      !! store flux label
      map(i) = j
    END DO

    IF (mapEuterpe .AND. .NOT.(useSquared)) THEN
      j = 1
      DO i = -1, stInd + 1, -1
        !! check if last original flux surface is reached and test if inside a
        !! given original flux surface interval
        DO WHILE ((-j > stInd + 1) .AND. (ABS((nFluxVMEC - 1) * (sVal(i - 1) -&
             sValVMEC(j)) - .5) - .5 > 1E-8))
          j = j + 1
        END DO
        !! calculate distance
        h(i) = (sVal(i - 1) - sValVMEC(j))
        !! store flux label
        map(i) = j
      END DO
    END IF

    !! recalculate test function
    t(:, :, :) = 0
    t(:, 0, 0) = 1
    DO i = stInd, radius
      IF (i == 0) CYCLE
      !! t(s) = s^(1/2)
      t(i, 1, 0) = SQRT(sVal(i))
      !! t'(s) = 1/2 s^(-1/2)
      t(i, 1, 1) = .5_hp / SQRT(sVal(i))
      !! t''(s) = -1/4 s^(-3/2)
      t(i, 1, 2) = -.25_hp / (sVal(i) * SQRT(sVal(i)))
      !! t'''(s) = 3/8 s^(-5/2)
      t(i, 1, 3) = .375_hp / (sVal(i) * sVal(i) * SQRT(sVal(i)))

      !! t(s) = s
      t(i, 2, 0) = sVal(i)
      !! t'(s) = 1
      t(i, 2, 1) = 1
      !! t''(s) = 0
      t(i, 2, 2) = 0
      !! t'''(s) = 0
      t(i, 2, 3) = 0

      !! t(s) = s^(3/2)
      t(i, 3, 0) = SQRT(sVal(i)) * sVal(i)
      !! t'(s) = 3/2 s^(1/2)
      t(i, 3, 1) = 1.5_hp * SQRT(sVal(i))
      !! t''(s) = 3/4 s^(-1/2)
      t(i, 3, 2) =  .75_hp / SQRT(sVal(i))
      !! t'''(s) = -3/8 s^(-3/2)
      t(i, 3, 3) = -.375_hp / (sVal(i) * SQRT(sVal(i)))

      !! t(s) = s^2
      t(i, 4, 0) = sVal(i) * sVal(i)
      !! t'(s) = 2 s
      t(i, 4, 1) = 2 * sVal(i)
      !! t''(s) = 2
      t(i, 4, 2) = 2
      !! t'''(s) = 0
      t(i, 4, 3) = 0
    END DO

    !! set new fourier coefficients
    DO i = 1, mn_mode
      !! calculate all other coefficients with h = h(j)
      DO j = 0, radius - 1
        !! t(s_{j+1})
        h1 = t(j + 1, asymp(i), 0)
        !! t'(s_{j+1})
        h2 = t(j + 1, asymp(i), 1)
        !! t''(s_{j+1})
        h3 = t(j + 1, asymp(i), 2)
        !! h(j) = s_{j+1} - s_j
        IF (j > 0) THEN
          l = h(j)
        ELSE
          l = sVal(1)
        END IF
        !! r_{j+1} = t(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j * h^3)
        rmnc(i, j + 1) = h1 * (xr(i, ma + map(j)) + l * (xr(i, mb + map(j)) +&
             l * (xr(i, mc + map(j)) + l * xr(i, md + map(j)))))
        !! dr_{j+1}/ds = t'(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j *
        !!                 h^3) + t(s_{j+1}) * (b_j + 2 * c_j * h + 3 * d_j *
        !!                 h^2)
        dRdSmnc(i, j + 1) = h2 * xr(i, ma + map(j)) + (h1 + h2 * l) *&
             xr(i, mb + map(j)) + l * ((2 * h1 + h2 * l) * xr(i, mc + map(j)) +&
             l * (3 * h1 + h2 * l) * xr(i, md + map(j)))
        !! d^2r_{j+1}/ds^2 = t''(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j *
        !!                     h^3) + 2 * t'(s_{j+1}) * (b_j + 2 * c_j * h +
        !!                     3 * d_j * h^2) + t(s_{j+1}) * (2 * c_j + 6 *
        !!                     d_j * h)
        d2RdS2mnc(i, j + 1) = h3 * xr(i, ma + map(j)) + (2 * h2 + h3 * l) *&
             xr(i, mb + map(j)) + (2 * h1 + (4 * h2 + h3 * l) * l) *&
             xr(i, mc + map(j)) + (6 * h1 + (6 * h2 + h3 * l) * l) * l *&
             xr(i, md + map(j))
        !! z_{j+1} = t(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j * h^3)
        zmns(i, j + 1) = h1 * (xz(i, ma + map(j)) + l * (xz(i, mb + map(j)) +&
             l * (xz(i, mc + map(j)) + l * xz(i, md + map(j)))))
        !! dz_{j+1}/ds = t'(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j *
        !!                 h^3) + t(s_{j+1}) * (b_j + 2 * c_j * h + 3 * d_j *
        !!                 h^2)
        dZdSmns(i, j + 1) = h2 * xz(i, ma + map(j)) + (h1 + h2 * l) *&
             xz(i, mb + map(j)) + l * ((2 * h1 + h2 * l) * xz(i, mc + map(j)) +&
             l * (3 * h1 + h2 * l) * xz(i, md + map(j)))
        !! d^2z_{j+1}/ds^2 = t''(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j *
        !!                     h^3) + 2 * t'(s_{j+1}) * (b_j + 2 * c_j * h +
        !!                     3 * d_j * h^2) + t(s_{j+1}) * (2 * c_j + 6 *
        !!                     d_j * h)
        d2ZdS2mns(i, j + 1) = h3 * xz(i, ma + map(j)) + (2 * h2 + h3 * l) *&
             xz(i, mb + map(j)) + (2 * h1 + (4 * h2 + h3 * l) * l) *&
             xz(i, mc + map(j)) + (6 * h1 + (6 * h2 + h3 * l) * l) * l *&
             xz(i, md + map(j))
      END DO
      IF (mapEuterpe .AND. (.NOT.(useSquared))) THEN
        DO j = 0, stInd + 1, -1
          !! t(s_{j+1})
          h1 = t(j - 1, asymp(i), 0)
          !! t'(s_{j+1})
          h2 = t(j - 1, asymp(i), 1)
          !! t''(s_{j+1})
          h3 = t(j - 1, asymp(i), 2)
          !! h(j) = s_{j+1} - s_j
          IF (j < 0) THEN
            l = h(j)
          ELSE
            l = sVal(-1)
          END IF
          !! r_{j+1} = t(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j * h^3)
          rmnc(i, j - 1) = h1 * (xr(i, ma + map(j)) + l * (xr(i, mb + map(j)) +&
               l * (xr(i, mc + map(j)) + l * xr(i, md + map(j)))))
          !! dr_{j+1}/ds = t'(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j *
          !!                 h^3) + t(s_{j+1}) * (b_j + 2 * c_j * h + 3 * d_j *
          !!                 h^2)
          dRdSmnc(i, j - 1) = h2 * xr(i, ma + map(j)) + (h1 + h2 * l) *&
               xr(i, mb + map(j)) + l * ((2 * h1 + h2 * l) * xr(i, mc + map(j)) +&
               l * (3 * h1 + h2 * l) * xr(i, md + map(j)))
          !! d^2r_{j+1}/ds^2 = t''(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j *
          !!                     h^3) + 2 * t'(s_{j+1}) * (b_j + 2 * c_j * h +
          !!                     3 * d_j * h^2) + t(s_{j+1}) * (2 * c_j + 6 *
          !!                     d_j * h)
          d2RdS2mnc(i, j - 1) = h3 * xr(i, ma + map(j)) + (2 * h2 + h3 * l) *&
               xr(i, mb + map(j)) + (2 * h1 + (4 * h2 + h3 * l) * l) *&
               xr(i, mc + map(j)) + (6 * h1 + (6 * h2 + h3 * l) * l) * l *&
               xr(i, md + map(j))
          !! z_{j+1} = t(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j * h^3)
          zmns(i, j - 1) = h1 * (xz(i, ma + map(j)) + l * (xz(i, mb + map(j)) +&
               l * (xz(i, mc + map(j)) + l * xz(i, md + map(j)))))
          !! dz_{j+1}/ds = t'(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j *
          !!                 h^3) + t(s_{j+1}) * (b_j + 2 * c_j * h + 3 * d_j *
          !!                 h^2)
          dZdSmns(i, j - 1) = h2 * xz(i, ma + map(j)) + (h1 + h2 * l) *&
               xz(i, mb + map(j)) + l * ((2 * h1 + h2 * l) * xz(i, mc + map(j)) +&
               l * (3 * h1 + h2 * l) * xz(i, md + map(j)))
          !! d^2z_{j+1}/ds^2 = t''(s_{j+1}) * (a_j + b_j * h + c_j * h^2 + d_j *
          !!                     h^3) + 2 * t'(s_{j+1}) * (b_j + 2 * c_j * h +
          !!                     3 * d_j * h^2) + t(s_{j+1}) * (2 * c_j + 6 *
          !!                     d_j * h)
          d2ZdS2mns(i, j - 1) = h3 * xz(i, ma + map(j)) + (2 * h2 + h3 * l) *&
               xz(i, mb + map(j)) + (2 * h1 + (4 * h2 + h3 * l) * l) *&
               xz(i, mc + map(j)) + (6 * h1 + (6 * h2 + h3 * l) * l) * l *&
               xz(i, md + map(j))
        END DO
      END IF
    END DO
    
    !! write interpolated R, z values for test purpose
    IF (debug .AND. (me_rank == 0)) THEN
      DO i = 1, mn_mode
        DO j = 1, nFluxVMEC - 1
          DO k = 0, 10
            l = REAL(k, hp) * h(j) / 10._hp
            lr = xr(i, ma + map(j)) + l * (xr(i, mb + map(j)) + l *&
                 (xr(i, mc + map(j)) + l * xr(i, md + map(j))))
            lz = xz(i, ma + map(j)) + l * (xz(i, mb + map(j)) + l *&
                 (xz(i, mc + map(j)) + l * xz(i, md + map(j))))
            h1 = (sVal(j) + l) ** (REAL(asymp(i), hp) / 2)
            WRITE(112, FMT = '(ES13.5, ",", I6, ",", ES13.5, ",", ES13.5, ",",&
                 & ES13.5, ",", ES13.5, ",", I6)') sVal(j) + l, i,&
                 xm(i), xn(i), lr * h1, lz * h1, asymp(i)
          END DO
        END DO
        WRITE(112, *)
      END DO
      CLOSE(112)
    END IF

    !! deallocate memory
    DEALLOCATE(a, perm, br, bz, stat = aStat)
    IF (me_rank == 0) THEN
      WRITE(*, *) " done!"
    END IF

  END SUBROUTINE SmoothRZ

! ----------------------------------------------------------------------
!> Compute the lambda function from the VMEC code with an algorithm
!> suggested by G. Leitold (G. Leitold, Magnetische Koordinatensysteme
!> in Stellaratoren, TU Graz, 2000)
!>
!>       dlambda     dlambda     d^2lambda       d^2lambda      d^2lambda
!> H = A ------- + B ------- + C --------- + 2D ----------- + E ---------
!>       dtheta       dphi       dtheta^2       dtheta dphi      dphi^2
!>
! ----------------------------------------------------------------------

  SUBROUTINE ComputeLambda

    REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: cmnc, dmnc, emnc
    REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: amns, bmns, hmns
    REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: aSmooth, lambda, argA
    REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: bSmooth, res
    REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: xm2, xn2, cxm, cxn
    INTEGER, DIMENSION(:), ALLOCATABLE :: perm
    REAL(KIND = hp) :: a, b, c, d, e, help, uI, vJ
    REAL(KIND = hp) :: a0, b0, c0, d0, e0, a1, b1, c1, d1, e1
    REAL(KIND = hp) :: a2, b2, c2, d2, e2, a3, b3, c3, d3, e3
    REAL(KIND = hp) :: r, dRdU, dRdV, dZdU, dZdV, sqrtG
    REAL(KIND = hp) :: guu, guv, gvv
    REAL(KIND = hp) :: d2RdU2, d2RdUdV, d2RdV2
    REAL(KIND = hp) :: d2ZdU2, d2ZdUdV, d2ZdV2
    REAL(KIND = hp) :: dsqrtGdU, dsqrtGdV
    INTEGER :: i, j, curS, mode, mode1, mode2, aStat, aError, n, dn
    INTEGER :: mn_mode2, n_start, cmmode, cnmode, cmn_mode, startS, endS

    IF (me_rank == 0) THEN
      WRITE(*, FMT = "(A)", ADVANCE = "no") "Compute lambda function ..."
    END IF
    !! calculate flux surfaces for current node
    !! number of flux surfaces per node
    n = INT(REAL(radius + 1 - stInd) / REAL(me_size))
    !! distribute the remaining surfaces
    dn = radius + 1 - stInd - n * me_size
    !! calculate the number of the first surface of the current node
    IF (dn > 0) THEN
      !! include the remaining surfaces
      IF (me_rank < dn) THEN
        n = n + 1
        startS = stInd + me_rank * n
      ELSE
        startS = stInd + dn * (n + 1) + (me_rank - dn) * n
      END IF
    ELSE
      startS = stInd + me_rank * n
    END IF
    !! calculate the number of the last surface of the current node
    endS = startS + n - 1

    !! we could not use factor 4 of the number of modes, because the
    !! the resulting system of equations is problematic
    !! (at least for one W7-X configuration)
    cmmode = 2 * INT(MAXVAL(xm(:)))
    cnmode = 2 * INT(MAXVAL(xn(:)) / nfp)
    cmn_mode = (2 * cnmode + 1) * cmmode + cnmode + 1
    mn_mode2 = (4 * cnmode + 1) * 2 * cmmode + 2 * cnmode + 1
    !! allocate memory for fourier series of help variables
    !! and for solving linear system of equations
    aError = 0
    ALLOCATE(amns(mn_mode2, startS:endS), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bmns(mn_mode2, startS:endS), stat = aStat)
    aError = aError + aStat
    ALLOCATE(cmnc(mn_mode2, startS:endS), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dmnc(mn_mode2, startS:endS), stat = aStat)
    aError = aError + aStat
    ALLOCATE(emnc(mn_mode2, startS:endS), stat = aStat)
    aError = aError + aStat
    ALLOCATE(hmns(mn_mode2, startS:endS), stat = aStat)
    aError = aError + aStat
    ALLOCATE(argA(mn_mode2, 2), stat = aStat)
    aError = aError + aStat
    ALLOCATE(aSmooth(cmn_mode, cmn_mode), stat = aStat)
    aError = aError + aStat
    ALLOCATE(lambda(cmn_mode, stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(bSmooth(cmn_mode), res(cmn_mode), perm(cmn_mode), stat = aStat)
    aError = aError + aStat
    ALLOCATE(cxm(cmn_mode), cxn(cmn_mode), stat = aStat)
    aError = aError + aStat
    ALLOCATE(xm2(mn_mode2), xn2(mn_mode2), stat = aStat)
    aError = aError + aStat

    !! initiate fourier components (doubled resolution)
    i = 0
    DO mode1 = 0, cmmode
      n_start = -cnmode
      IF (mode1 == 0) n_start = 0
      DO mode2 = n_start, cnmode
        i = i + 1
        cxm(i) = mode1
        cxn(i) = nfp * mode2
      END DO
    END DO
    i = 0
    DO mode1 = 0, 2 * cmmode
      n_start = -2 * cnmode
      IF (mode1 == 0) n_start = 0
      DO mode2 = n_start, 2 * cnmode
        i = i + 1
        xm2(i) = mode1
        xn2(i) = nfp * mode2
      END DO
    END DO

    !! get fourier representation (doubled resolution) of help variables
    !!
    !!       g_{phi,phi}      g_{theta,phi}        g_{theta,theta}
    !! c = - -----------, d = -------------, e = - ---------------
    !!          sqrtg             sqrtg                 sqrtg
    !!
    !!
    !!      dd      dc          dd      de
    !! a = ---- + ------, b = ------ + ----
    !!     dphi   dtheta      dtheta   dphi
    !!
    amns(:, :) = 0
    bmns(:, :) = 0
    cmnc(:, :) = 0
    dmnc(:, :) = 0
    emnc(:, :) = 0
    !! numerical integration (Trapez Rule)
    DO j = 1, intPointsV
      !! supporting points for v direction
      vJ = vJv(j)
      DO i = 1, intPointsU
        !! supporting points for u direction
        uI = uIv(i)
        !! argument vector
        argA(1, :) = (/0._hp, 1._hp/)
        argA(2:, 1) = 2 * SIN(xm2(2:) * uI - xn2(2:) * vJ)
        argA(2:, 2) = 2 * COS(xm2(2:) * uI - xn2(2:) * vJ)
        !! radial loop
        DO curS = startS, MIN(endS, radius)
          IF (curS == 0) CYCLE
          !! calculate R
          r = CosTransFullMeshF2(mn_mode, rmnc(:, curS), i, j)
          !! calculate dR/du
          dRdU = SinTransFullMeshF2(mn_mode, dRdUmns(:, curS), i, j)
          !! calculate dR/dv
          dRdV = SinTransFullMeshF2(mn_mode, dRdVmns(:, curS), i, j)
          !! calculate dz/du
          dZdU = CosTransFullMeshF2(mn_mode, dZdUmnc(:, curS), i, j)
          !! calculate dz/dv
          dZdV = CosTransFullMeshF2(mn_mode, dZdVmnc(:, curS), i, j)
          !! d2R/du2
          d2RdU2 = CosTransFullMeshF2(mn_mode, d2RdU2mnc(:, curS), i, j)
          !! d2R/(du dv)
          d2RdUdV = CosTransFullMeshF2(mn_mode, d2RdUdVmnc(:, curS), i, j)
          !! d2R/dv2
          d2RdV2 = CosTransFullMeshF2(mn_mode, d2RdV2mnc(:, curS), i, j)
          !! d2z/du2
          d2ZdU2 = SinTransFullMeshF2(mn_mode, d2ZdU2mns(:, curS), i, j)
          !! d2z/(du dv)
          d2ZdUdV = SinTransFullMeshF2(mn_mode, d2ZdUdVmns(:, curS), i, j)
          !! d2z/du2
          d2ZdV2 = SinTransFullMeshF2(mn_mode, d2ZdV2mns(:, curS), i, j)
          sqrtG = CosTransFullMeshF2(mn_mode_nyq, gmnc(:, curS), i, j)
          dsqrtGdU = SinTransFullMeshF2(mn_mode_nyq,&
               -xm_nyq(:) * gmnc(:, curS), i, j)
          dsqrtGdV = SinTransFullMeshF2(mn_mode_nyq, xn_nyq(:) * gmnc(:, curS),&
               i, j)
          !!          dR dR   dz dz
          !! g_{uu} = -- -- + -- --
          !!          du du   du du
          guu = dRdU * dRdU + dZdU * dZdU
          !!          dR dR   dz dz
          !! g_{uv} = -- -- + -- --
          !!          du dv   du dv
          guv = dRdU * dRdV + dZdU * dZdV
          !!          dR dR   dz dz
          !! g_{vv} = -- -- + -- -- + r^2
          !!          dv dv   dv dv
          gvv = dRdV * dRdV + dZdV * dZdV + r * r
          a = (d2RdV2 * dRdU - d2RdUdV * dRdV + d2ZdV2 * dZdU - d2ZdUdV *&
               dZdV - 2. * dRdU * r) * sqrtG - guv * dsqrtGdV + gvv * dsqrtGdU
          b = (d2RdU2 * dRdV - d2RdUdV * dRdU + d2ZdU2 * dZdV - d2ZdUdV *&
               dZdU) * sqrtG - guv * dsqrtGdU + guu * dsqrtGdV
          c = -gvv * sqrtG
          d = guv * sqrtG
          e = -guu * sqrtG
          !! generate fourier coefficients
          !! sum up integrands
          amns(:, curS) = amns(:, curS) + a * argA(:, 1)
          bmns(:, curS) = bmns(:, curS) + b * argA(:, 1)
          cmnc(:, curS) = cmnc(:, curS) + c * argA(:, 2)
          dmnc(:, curS) = dmnc(:, curS) + d * argA(:, 2)
          emnc(:, curS) = emnc(:, curS) + e * argA(:, 2)
        END DO
      END DO
    END DO
    !! normalize vectors
    help = REAL(intPointsU * intPointsV, hp)
    amns(:, :) = amns(:, :) / help
    bmns(:, :) = bmns(:, :) / help
    cmnc(:, :) = cmnc(:, :) / help
    dmnc(:, :) = dmnc(:, :) / help
    emnc(:, :) = emnc(:, :) / help

    !! generate linear system of equation for lambda for every flux surface
    lambda(:, :) = 0
    DO curS = startS, MIN(endS, radius)
      IF (curS == 0) CYCLE
      IF (debug) THEN
        IF (me_size > 1) THEN
          PRINT *, curS, " on ", me_rank
        ELSE
          PRINT *, curS
        END IF
      END IF
      !! generate remaining help variables
      !!
      !! H = chi' B - psi' A
      !!
      hmns(:, curS) = iotas(curS) * bmns(:, curS) - amns(:, curS)
      !! calculate matrix for following system of equations
      !!
      !! 2 h_{o,p} = \sum_{m,n} \lambda_{m,n}(m(a_{o-m,p-n}+a_{o+m,p+n}-
      !!               a_{m-o,n-p})-n(b_{o-m,p-n}+b_{o+m,p+n}-b_{m-o,n-p})-
      !!               m^2(c_{o-m,p-n}-c_{o+m,p+n}+c_{m-o,n-p})+
      !!               2mn(d_{o-m,p-n}-d_{o+m,p+n}+d_{m-o,n-p})-
      !!               n^2(e_{o-m,p-n}-e_{o+m,p+n}+e_{m-o,n-p}))
      !!
      aSmooth(:, :) = 0
      !! loop over (o,p)
      DO mode1 = 2, cmn_mode
        !! map smaller resolution to higher
        mode = INT(cxm(mode1) * (4 * cnmode + 1) + cxn(mode1) / nfp + 1)
        !! right hand side
        bSmooth(mode1) = 2 * hmns(mode, curS)
        !! loop over (m,n)
        DO mode2 = 1, cmn_mode
          !! o - m
          i = INT(cxm(mode1) - cxm(mode2))
          !! p - n
          j = INT((cxn(mode1) - cxn(mode2)) / nfp)
          !! set value only if (i,j) are valid modes
          IF ((i < 0) .OR. (ABS(i) > 2 * cmmode) .OR.&
               (ABS(j) > 2 * cnmode) .OR. ((i == 0) .AND. (j < 0))) THEN
            a0 = 0
            b0 = 0
            c0 = 0
            d0 = 0
            e0 = 0
          ELSE
            !! get mode number
            mode = i * (4 * cnmode + 1) + j + 1
            a0 = amns(mode, curS)
            b0 = bmns(mode, curS)
            c0 = cmnc(mode, curS)
            d0 = dmnc(mode, curS)
            e0 = emnc(mode, curS)
          END IF
          !! -(o + m)
          i = -INT(cxm(mode1) + cxm(mode2))
          !! -(p + n)
          j = -INT((cxn(mode1) + cxn(mode2)) / nfp)
          !! set value only if (i,j) are valid modes
          IF ((i < 0) .OR. (ABS(i) > 2 * cmmode) .OR.&
               (ABS(j) > 2 * cnmode) .OR. ((i == 0) .AND. (j < 0))) THEN
            a1 = 0
            b1 = 0
            c1 = 0
            d1 = 0
            e1 = 0
          ELSE
            !! get mode number
            mode = i * (4 * cnmode + 1) + j + 1
            a1 = amns(mode, curS)
            b1 = bmns(mode, curS)
            c1 = cmnc(mode, curS)
            d1 = dmnc(mode, curS)
            e1 = emnc(mode, curS)
          END IF
          !! o + m
          i = INT(cxm(mode1) + cxm(mode2))
          !! p + n
          j = INT((cxn(mode1) + cxn(mode2)) / nfp)
          !! set value only if (i,j) are valid modes
          IF ((i < 0) .OR. (ABS(i) > 2 * cmmode) .OR.&
               (ABS(j) > 2 * cnmode) .OR. ((i == 0) .AND. (j < 0))) THEN
            a2 = 0
            b2 = 0
            c2 = 0
            d2 = 0
            e2 = 0
          ELSE
            !! get mode number
            mode = i * (4 * cnmode + 1) + j + 1
            a2 = amns(mode, curS)
            b2 = bmns(mode, curS)
            c2 = cmnc(mode, curS)
            d2 = dmnc(mode, curS)
            e2 = emnc(mode, curS)
          END IF
          !! -(o - m)
          i = -INT(cxm(mode1) - cxm(mode2))
          !! -(p - n)
          j = -INT((cxn(mode1) - cxn(mode2)) / nfp)
          !! set value only if (i,j) are valid modes
          IF ((i < 0) .OR. (ABS(i) > 2 * cmmode) .OR.&
               (ABS(j) > 2 * cnmode) .OR. ((i == 0) .AND. (j < 0))) THEN
            a3 = 0
            b3 = 0
            c3 = 0
            d3 = 0
            e3 = 0
          ELSE
            !! get mode number
            mode = i * (4 * cnmode + 1) + j + 1
            a3 = amns(mode, curS)
            b3 = bmns(mode, curS)
            c3 = cmnc(mode, curS)
            d3 = dmnc(mode, curS)
            e3 = emnc(mode, curS)
          END IF
          !! update matrix element
          aSmooth(mode1, mode2) = aSmooth(mode1, mode2) +&
               cxm(mode2) * (a0 - a1 + a2 - a3) -&
               cxn(mode2) * (b0 - b1 + b2 - b3) -&
               cxm(mode2) * cxm(mode2) * (c0 - c1 - c2 + c3) +&
               2 * cxm(mode2) * cxn(mode2) * (d0 - d1 - d2 + d3) -&
               cxn(mode2) * cxn(mode2) * (e0 - e1 - e2 + e3)
        END DO
      END DO

      !! solve linear system of equations
      !! one could get the (0,0) component, because lambda is represented by
      !! a fourier sine transform
      CALL GMRES_Solve(800, cmn_mode - 1, aSmooth(2:, 2:), bSmooth(2:),&
           res(2:), 1E-12_hp)
      res(1) = 0
      !! set the newly calculated lambda components
      lambda(:, curS) = res(:)
    END DO


    !! set new lambda
    DO mode1 = 1, mn_mode
      mode = INT(xm(mode1) * (2 * cnmode + 1) + xn(mode1) / nfp + 1)
      DO curS = stInd, radius
        IF (ABS(xn(mode1) * iotas(curS) - xm(mode1)) > 1E-6) THEN
          lmns(mode1, curS) = lambda(mode, curS)
        ELSE
          lmns(mode1, curS) = 0
        END IF
      END DO
    END DO

    DEALLOCATE(xn2, xm2, cxn, cxm, perm, res, bSmooth)
    DEALLOCATE(aSmooth, lambda, argA, amns, bmns, cmnc, dmnc, emnc, hmns)
    IF (me_rank == 0) THEN
      WRITE(*, *) " done!"
    END IF

  END SUBROUTINE ComputeLambda

  SUBROUTINE ComputeMagnField

    REAL(KIND = hp), DIMENSION(mn_mode_nyq, radius) :: absBmnc
    REAL(KIND = hp), DIMENSION(mn_mode_nyq) :: argA
    REAL(KIND = hp), DIMENSION(radius) :: absB, absBs
    REAL(KIND = hp) :: r, dRdU, dRdV, dZdU, dZdV
    REAL(KIND = hp) :: dLambdadU, dLambdadV, bR, bP, bZ, sqrtG
    !! numerical integration (Trapez Rule)
    INTEGER :: i, j, curS

    IF (me_rank == 0) THEN
      WRITE(*, FMT = "(A)", ADVANCE = "no") "Compute magnetic field ..."
    END IF

    !! initialization
    absBmnc(:, :) = 0
    absBsmnc(:, :) = 0
    dBdSmnc(:, :) = 0
    dBdSsmnc(:, :) = 0
    !! numerical integration (Trapez Rule)
    DO j = 1, intPointsV
      DO i = 1, intPointsU
        !! argument vector
        argA(1) = 1._hp
        argA(2:) = 2 * vCosN(2:, i, j)
        !! radial loop
        DO curS = 1, radius
          !! calculate R
          r = CosTransFullMeshF2(mn_mode, rmnc(:, curS), i, j)
          !! calculate dR/du
          dRdU = SinTransFullMeshF2(mn_mode, dRdUmns(:, curS), i, j)
          !! calculate dR/dv
          dRdV = SinTransFullMeshF2(mn_mode, dRdVmns(:, curS), i, j)
          !! calculate dz/du
          dZdU = CosTransFullMeshF2(mn_mode, dZdUmnc(:, curS), i, j)
          !! calculate dz/dv
          dZdV = CosTransFullMeshF2(mn_mode, dZdVmnc(:, curS), i, j)
          !! calculate dlambda/du
          dLambdadU = CosTransFullMeshF2(mn_mode, lUmnc(:, curS), i, j)
          !! calculate dlambda/dv
          dLambdadV = CosTransFullMeshF2(mn_mode, lVmnc(:, curS), i, j)
          !! calculate sqrtG
          sqrtG = CosTransFullMeshF2(mn_mode_nyq, gmnc(:, curS), i, j)
          !! calculate B_R = dR/du B^u + dR/dv B^v
          bR = dRdU * (iotas(curS) - dLambdadV) + dRdV * (1._hp + dLambdadU)
          !! calculate B_phi = r B^v
          bP = r * (1._hp + dLambdadU)
          !! calculate B_z = dz/du B^u + dz/dv B^v
          bZ = dZdU * (iotas(curS) - dLambdadV) + dZdV * (1._hp + dLambdadU)
          !! calculate absB
          !! scaled magnetic field
          absBs(curS) = SQRT(bR * bR + bP * bP + bZ * bZ) * ABS(phipf(curS))
          !! magnetic field
          absB(curS) = absBs(curS) / ABS(sqrtG)
          !! generate fourier coefficients for both fields
          !! sum up integrands
          absBsmnc(:, curS) = absBsmnc(:, curS) + absBs(curS) * argA(:)
          absBmnc(:, curS) = absBmnc(:, curS) + absB(curS) * argA(:)
        END DO
      END DO
    END DO
    !! normalize vectors
    absBsmnc(:, :) = absBsmnc(:, :) / REAL(intPointsU * intPointsV, hp)
    bmnc(:, :) = absBmnc(:, :) / REAL(intPointsU * intPointsV, hp)
    !! get fourier coefficients of the radial derivative of |B|
    DO curS = 1, radius
      DO i = 1, mn_mode_nyq
        !! scaled |B|
        dBdSsmnc(i, curS) = GetRadialDerivative(radius, curS, absBsmnc(i, :))
        !! |B|
        dBdSmnc(i, curS) = GetRadialDerivative(radius, curS, bmnc(i, :))
      END DO
    END DO

    IF (me_rank == 0) THEN
      WRITE(*, *) " done!"
    END IF

  END SUBROUTINE ComputeMagnField

! ----------------------------------------------------------------------
!> calculate distance between point on a flux surface and a selected
!> point
! ----------------------------------------------------------------------

  FUNCTION DistanceToPoint(curS, th, phi, r, z) RESULT(dist)

!> selected flux surface
    INTEGER, INTENT(IN) :: curS

!> theta coordinate of point on flux surface
    REAL(KIND = hp), INTENT(IN) :: th

!> selected phi cut
    REAL(KIND = hp), INTENT(IN) :: phi

!> r coordinate of selected point
    REAL(KIND = hp), INTENT(IN) :: r

!> z coordinate of selected point
    REAL(KIND = hp), INTENT(IN) :: z

!> distance between the two points
    REAL(KIND = hp) :: dist

    REAL(KIND = hp), DIMENSION(mn_mode) :: argV
    REAL(KIND = hp) :: dR, dZ

    argV(:) = th * xm(:) - phi * xn(:)
    dR = SUM(rmnc(:, curS) * COS(argV(:))) - r
    dZ = SUM(zmns(:, curS) * SIN(argV(:))) - z
    dist = SQRT(dR * dR + dZ * dZ)

  END FUNCTION DistanceToPoint

! ----------------------------------------------------------------------
!> interpolate half mesh data to full mesh
!>
!> half mesh data is know at radial grid points
!>             (i - 1.5)ds   i = 2,...,n; ds = sMax/(n - 1)
!> half mesh data is zero on first radial point, which is used as value on
!> magnetic axis
!>
!> full mesh data should be known at radial grid points
!>             (i - 1) ds          i = 1,...¸n; ds = sMax/(n - 1)
!>
!> linear interpolation
! ----------------------------------------------------------------------

  SUBROUTINE HalfMeshToFullMesh1D(n, a)

!> number of radial grid point
    INTEGER, INTENT(IN) :: n

!> fourier components
    REAL(KIND = hp), DIMENSION(n), INTENT(INOUT) :: a

    REAL(KIND = hp) :: extPol

    !! extrapolate point on first surface
    extPol = .5_hp * (3 * a(2) - a(3))
    !! set first value
    a(1) = extPol
    !! extrapolate point on last surface
    extPol = .5_hp * (3 * a(n) - a(n - 1))
    !! loop over radial points
    a(2:(n - 1)) = .5_hp * (a(2:(n - 1)) + a(3:n))
    !! set last value
    a(n) = extPol

  END SUBROUTINE HalfMeshToFullMesh1D

  SUBROUTINE HalfMeshToFullMesh2D(m, n, a)

!> number of fourier components
    INTEGER, INTENT(IN) :: m

!> number of radial grid point
    INTEGER, INTENT(IN) :: n

!> fourier components
    REAL(KIND = hp), DIMENSION(m, n), INTENT(INOUT) :: a

    INTEGER :: i
    REAL(KIND = hp) :: extPol

    !! loop over all fourier components
    DO i = 1, m
      !! extrapolate point on first surface
      extPol = .5_hp * (3 * a(i, 2) - a(i, 3))
      !! set first value
      a(i, 1) = extPol
      !! extrapolate point on last surface
      extPol = .5_hp * (3 * a(i, n) - a(i, n - 1))
      !! loop over radial points
      a(i, 2:(n - 1)) = .5_hp * (a(i, 2:(n - 1)) + a(i, 3:n))
      !! set last value
      a(i, n) = extPol
    END DO

  END SUBROUTINE HalfMeshToFullMesh2D

! ----------------------------------------------------------------------
!> get the radial partial derivative of a function using central 
!> difference 
! ----------------------------------------------------------------------

  FUNCTION GetRadialDerivativeFourier(radius, curS, mn_mode, coeff, xm, xn,&
       uI, vJ, mode) RESULT(dVal)

!> number of flux surfaces
    INTEGER, INTENT(IN) :: radius

!> s position
    INTEGER, INTENT(IN) :: curS

!> number of fourier modes
    INTEGER, INTENT(IN) :: mn_mode

!> fourier coefficients
    REAL(KIND = hp), DIMENSION(mn_mode, radius), INTENT(IN) :: coeff

!> poloidal mode numbers
    REAL(KIND = hp), DIMENSION(mn_mode), INTENT(IN) :: xm

!> toroidal mode numbers
    REAL(KIND = hp), DIMENSION(mn_mode), INTENT(IN) :: xn

!> current poloidal value
    REAL(KIND = hp), INTENT(IN) :: uI

!> current toroidal value
    REAL(KIND = hp), INTENT(IN) :: vJ

!> fourier mode (1 - cosine, 2 - sine)
    INTEGER, INTENT(IN) :: mode

!> resulting derivative
    REAL(KIND = hp) :: dVal

    REAL(KIND = hp) :: h1, h2, h3
    REAL(KIND = hp) :: w1, w2, w3

    SELECT CASE (mode)
    CASE(1)
      !! cosine fourier transform
      IF (curS <= 2) THEN
        !! get supporting points
        h1 = CosTransFullMesh(mn_mode, coeff(:, curS), xm(:), xn(:),&
             uI, vJ)
        h2 = CosTransFullMesh(mn_mode, coeff(:, curS + 1), xm(:), xn(:),&
             uI, vJ)
        h3 = CosTransFullMesh(mn_mode, coeff(:, curS + 2), xm(:), xn(:),&
             uI, vJ)
        !! weights for central differencing
        CALL GetWeights(1, curS, w1, w2, w3)
      ELSE IF (curS == radius) THEN
        !! get supporting points
        h1 = CosTransFullMesh(mn_mode, coeff(:, curS), xm(:), xn(:),&
             uI, vJ)
        h2 = CosTransFullMesh(mn_mode, coeff(:, curS - 1), xm(:), xn(:),&
             uI, vJ)
        h3 = CosTransFullMesh(mn_mode, coeff(:, curS - 2), xm(:), xn(:),&
             uI, vJ)
        !! weights for central differencing
        CALL GetWeights(2, curS, w1, w2, w3)
      ELSE
        !! get supporting points
        h1 = CosTransFullMesh(mn_mode, coeff(:, curS), xm(:), xn(:),&
             uI, vJ)
        h2 = CosTransFullMesh(mn_mode, coeff(:, curS - 1), xm(:), xn(:),&
             uI, vJ)
        h3 = CosTransFullMesh(mn_mode, coeff(:, curS + 1), xm(:), xn(:),&
             uI, vJ)
        !! weights for central differencing
        CALL GetWeights(3, curS, w1, w2, w3)
      END IF
    CASE(2)
      !! sine fourier transform
      IF (curS <= 2) THEN
        !! get supporting points
        h1 = SinTransFullMesh(mn_mode, coeff(:, curS), xm(:), xn(:),&
             uI, vJ)
        h2 = SinTransFullMesh(mn_mode, coeff(:, curS + 1), xm(:), xn(:),&
             uI, vJ)
        h3 = SinTransFullMesh(mn_mode, coeff(:, curS + 2), xm(:), xn(:),&
             uI, vJ)
        !! weights for central differencing
        CALL GetWeights(1, curS, w1, w2, w3)
      ELSE IF (curS == radius) THEN
        !! get supporting points
        h1 = SinTransFullMesh(mn_mode, coeff(:, curS), xm(:), xn(:),&
             uI, vJ)
        h2 = SinTransFullMesh(mn_mode, coeff(:, curS - 1), xm(:), xn(:),&
             uI, vJ)
        h3 = SinTransFullMesh(mn_mode, coeff(:, curS - 2), xm(:), xn(:),&
             uI, vJ)
        !! weights for central differencing
        CALL GetWeights(2, curS, w1, w2, w3)
      ELSE
        !! get supporting points
        h1 = SinTransFullMesh(mn_mode, coeff(:, curS), xm(:), xn(:),&
             uI, vJ)
        h2 = SinTransFullMesh(mn_mode, coeff(:, curS - 1), xm(:), xn(:),&
             uI, vJ)
        h3 = SinTransFullMesh(mn_mode, coeff(:, curS + 1), xm(:), xn(:),&
             uI, vJ)
        !! weights for central differencing
        CALL GetWeights(3, curS, w1, w2, w3)
      END IF
    CASE DEFAULT
      PRINT *, "Unknown fourier transform!"
      CALL EXIT(18)
    END SELECT
    !! calculate central differencing
    dVal = w1 * h1 + w2 * h2 + w3 * h3

  END FUNCTION GetRadialDerivativeFourier
 
! ----------------------------------------------------------------------
!> get the radial partial derivative of a function using central 
!> difference 
! ----------------------------------------------------------------------

  FUNCTION GetRadialDerivativeF(radius, curS, mn_mode, coeff, v) RESULT(dVal)

!> number of flux surfaces
    INTEGER, INTENT(IN) :: radius

!> s position
    INTEGER, INTENT(IN) :: curS

!> number of fourier modes
    INTEGER, INTENT(IN) :: mn_mode

!> fourier coefficients
    REAL(KIND = hp), DIMENSION(mn_mode, radius), INTENT(IN) :: coeff

!> poloidal mode numbers
    REAL(KIND = hp), DIMENSION(mn_mode), INTENT(IN) :: v

!> resulting derivative
    REAL(KIND = hp) :: dVal

    REAL(KIND = hp) :: h1, h2, h3
    REAL(KIND = hp) :: w1, w2, w3

    IF (curS <= 2) THEN
      !! get supporting points
      h1 = SUM(coeff(:, curS) * v(:))
      h2 = SUM(coeff(:, curS + 1) * v(:))
      h3 = SUM(coeff(:, curS + 2) * v(:))
      !! weights for central differencing
      CALL GetWeights(1, curS, w1, w2, w3)
    ELSE IF (curS == radius) THEN
      !! get supporting points
      h1 = SUM(coeff(:, curS) * v(:))
      h2 = SUM(coeff(:, curS - 1) * v(:))
      h3 = SUM(coeff(:, curS - 2) * v(:))
      !! weights for central differencing
      CALL GetWeights(2, curS, w1, w2, w3)
    ELSE
      !! get supporting points
      h1 = SUM(coeff(:, curS) * v(:))
      h2 = SUM(coeff(:, curS - 1) * v(:))
      h3 = SUM(coeff(:, curS + 1) * v(:))
      !! weights for central differencing
      CALL GetWeights(3, curS, w1, w2, w3)
    END IF
    !! calculate central differencing
    dVal = w1 * h1 + w2 * h2 + w3 * h3

  END FUNCTION GetRadialDerivativeF

  FUNCTION GetRadialDerivativeFPolar(radius, curS, mn_mode, coeff, v) RESULT(dVal)

!> number of flux surfaces
    INTEGER, INTENT(IN) :: radius

!> s position
    INTEGER, INTENT(IN) :: curS

!> number of fourier modes
    INTEGER, INTENT(IN) :: mn_mode

!> fourier coefficients
    REAL(KIND = hp), DIMENSION(mn_mode, stInd:(stInd + radius - 1)),&
         INTENT(IN) :: coeff

!> poloidal mode numbers
    REAL(KIND = hp), DIMENSION(mn_mode), INTENT(IN) :: v

!> resulting derivative
    REAL(KIND = hp) :: dVal

    REAL(KIND = hp) :: h1, h2, h3
    REAL(KIND = hp) :: w1, w2, w3

    IF (curS > 0) THEN
      IF (curS <= 2) THEN
        !! get supporting points
        h1 = SUM(coeff(:, curS) * v(:))
        h2 = SUM(coeff(:, curS + 1) * v(:))
        h3 = SUM(coeff(:, curS + 2) * v(:))
        !! weights for central differencing
        CALL GetWeights(1, curS, w1, w2, w3)
      ELSE IF (curS == radius) THEN
        !! get supporting points
        h1 = SUM(coeff(:, curS) * v(:))
        h2 = SUM(coeff(:, curS - 1) * v(:))
        h3 = SUM(coeff(:, curS - 2) * v(:))
        !! weights for central differencing
        CALL GetWeights(2, curS, w1, w2, w3)
      ELSE
        !! get supporting points
        h1 = SUM(coeff(:, curS) * v(:))
        h2 = SUM(coeff(:, curS - 1) * v(:))
        h3 = SUM(coeff(:, curS + 1) * v(:))
        !! weights for central differencing
        CALL GetWeights(3, curS, w1, w2, w3)
      END IF
    ELSE
      IF (curS >= -2) THEN
        !! get supporting points
        h1 = SUM(coeff(:, curS) * v(:))
        h2 = SUM(coeff(:, curS - 1) * v(:))
        h3 = SUM(coeff(:, curS - 2) * v(:))
        !! weights for central differencing
        CALL GetWeights(2, curS, w1, w2, w3)
      ELSE IF (curS == -radius) THEN
        !! get supporting points
        h1 = SUM(coeff(:, curS) * v(:))
        h2 = SUM(coeff(:, curS + 1) * v(:))
        h3 = SUM(coeff(:, curS + 2) * v(:))
        !! weights for central differencing
        CALL GetWeights(1, curS, w1, w2, w3)
      ELSE
        !! get supporting points
        h1 = SUM(coeff(:, curS) * v(:))
        h2 = SUM(coeff(:, curS - 1) * v(:))
        h3 = SUM(coeff(:, curS + 1) * v(:))
        !! weights for central differencing
        CALL GetWeights(3, curS, w1, w2, w3)
      END IF
    END IF
    !! calculate central differencing
    dVal = w1 * h1 + w2 * h2 + w3 * h3

  END FUNCTION GetRadialDerivativeFPolar

! ----------------------------------------------------------------------
!> get the radial partial derivative of a function using central 
!> difference 
! ----------------------------------------------------------------------

  FUNCTION GetRadialDerivative(radius, curS, fVal) RESULT(dVal)

!> number of flux surfaces
    INTEGER, INTENT(IN) :: radius

!> s position
    INTEGER, INTENT(IN) :: curS

!> function values
    REAL(KIND = hp), DIMENSION(radius), INTENT(IN) :: fVal

!> resulting derivative
    REAL(KIND = hp) :: dVal

    REAL(KIND = hp) :: h1, h2, h3
    REAL(KIND = hp) :: w1, w2, w3

    IF (curS <= 2) THEN
      !! get supporting points
      h1 = fVal(curS)
      h2 = fVal(curS + 1)
      h3 = fVal(curS + 2)
      !! weights for central differencing
      CALL GetWeights(1, curS, w1, w2, w3)
    ELSE IF (curS == radius) THEN
      !! get supporting points
      h1 = fVal(curS)
      h2 = fVal(curS - 1)
      h3 = fVal(curS - 2)
      !! weights for central differencing
      CALL GetWeights(2, curS, w1, w2, w3)
    ELSE
      !! get supporting points
      h1 = fVal(curS)
      h2 = fVal(curS - 1)
      h3 = fVal(curS + 1)
      !! weights for central differencing
      CALL GetWeights(3, curS, w1, w2, w3)
    END IF
    !! calculate central differencing
    dVal = w1 * h1 + w2 * h2 + w3 * h3

  END FUNCTION GetRadialDerivative

! ----------------------------------------------------------------------
!> get the radial partial derivative of a function using central 
!> difference 
! ----------------------------------------------------------------------

  FUNCTION GetRadialDerivative1(radius, curS, fVal) RESULT(dVal)

!> number of flux surfaces
    INTEGER, INTENT(IN) :: radius

!> s position
    INTEGER, INTENT(IN) :: curS

!> function values
    REAL(KIND = hp), DIMENSION(radius), INTENT(IN) :: fVal

!> resulting derivative
    REAL(KIND = hp) :: dVal

    REAL(KIND = hp) :: h1, h2, h3
    REAL(KIND = hp) :: w1, w2, w3

    IF (curS == 1) THEN
      !! get supporting points
      h1 = fVal(curS)
      h2 = fVal(curS + 1)
      h3 = fVal(curS + 2)
      !! weights for central differencing
      CALL GetWeights(1, curS, w1, w2, w3)
    ELSE IF (curS == radius) THEN
      !! get supporting points
      h1 = fVal(curS)
      h2 = fVal(curS - 1)
      h3 = fVal(curS - 2)
      !! weights for central differencing
      CALL GetWeights(2, curS, w1, w2, w3)
    ELSE
      !! get supporting points
      h1 = fVal(curS)
      h2 = fVal(curS - 1)
      h3 = fVal(curS + 1)
      !! weights for central differencing
      CALL GetWeights(3, curS, w1, w2, w3)
    END IF
    !! calculate central differencing
    dVal = w1 * h1 + w2 * h2 + w3 * h3

  END FUNCTION GetRadialDerivative1

! ----------------------------------------------------------------------
!> calculate weights for central differencing
! g1 = f(x)
! g2 = f(x + h1) = f(x) + h1 f'(x) + h1^2/2 f''(x)
! g3 = f(x + h2) = f(x) + h2 f'(x) + h2^2/2 f''(x)
!
! solve:       a g1 + b g2 + c g3 = f'(x)
!
! (a + b + c) f(x) + (b h1 + c h2) f'(x) + .5 (b h1^2 + c h2^2) f''(x) = f'(x)
!
!  a -> -((h1 + h2) / (h1 h2))
!  b -> h2 / (h1 (h2 - h1)) 
!  c -> -h1 / (h2 (h2 - h1))
! ----------------------------------------------------------------------

  SUBROUTINE GetWeights(mode, curS, w1, w2, w3)

!> central differencing mode (1- left boundary, 2- right boundary, 3- normal)
    INTEGER, INTENT(IN) :: mode

!> current flux surface
    INTEGER, INTENT(IN) :: curS

!> weights
    REAL(KIND = hp), INTENT(out) :: w1, w2, w3

    REAL(KIND = hp) :: h1, h2

    SELECT CASE (mode)
    CASE (1)
      h1 = sVal(curS + 1) - sVal(curS)
      h2 = sVal(curS + 2) - sVal(curS)
    CASE (2)
      h1 = sVal(curS - 1) - sVal(curS)
      h2 = sVal(curS - 2) - sVal(curS)
    CASE (3)
      h1 = sVal(curS - 1) - sVal(curS)
      h2 = sVal(curS + 1) - sVal(curS)
    END SELECT
    w1 = -(h1 + h2) / (h1 * h2)
    w2 = h2 / (h1 * (h2 - h1))
    w3 = -h1 / (h2 * (h2 - h1))

  END SUBROUTINE GetWeights

! ----------------------------------------------------------------------
!> calculate backward cosine fourier transform with reusing cosine
!> vector
! ----------------------------------------------------------------------

  FUNCTION CosTransFullMeshF2(count, aij, i, j)

!> number of fourier elements (mn_mode or mn_mode_nyq)
    INTEGER, INTENT(IN) :: count

!> cosine fourier elements
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: aij

!> indices for u and v
    INTEGER, INTENT(IN) :: i, j

!> resulting function value
    REAL(KIND = hp) :: CosTransFullMeshF2

    IF (count == mn_mode) THEN
      CosTransFullMeshF2 = SUM(aij(:) * vCos(:, i, j))
    ELSE IF (count == mn_mode_nyq) THEN
      CosTransFullMeshF2 = SUM(aij(:) * vCosN(:, i, j))
    ELSE
      PRINT *, "Wrong number of fourier modes in CosTransFullMeshF2!"
      CALL EXIT(19)
    END IF

  END FUNCTION CosTransFullMeshF2

! ----------------------------------------------------------------------
!> calculate backward cosine fourier transform with reusing cosine
!> vector
! ----------------------------------------------------------------------

  FUNCTION CosTransFullMeshF(count, aij, indU, indV, u, v, recalc)

!> number of fourier elements
    INTEGER, INTENT(IN) :: count

!> cosine fourier elements
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: aij

!> indices for u and v
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: indU, indV

!> arguments
    REAL(KIND = hp), INTENT(IN) :: u, v

!> recalc cosine vector
    LOGICAL, INTENT(IN) :: recalc

!> resulting function value
    REAL(KIND = hp) :: CosTransFullMeshF

    REAL(KIND = hp), DIMENSION(:), ALLOCATABLE, SAVE :: cosV

    IF (recalc) THEN
      IF (ALLOCATED(cosV)) DEALLOCATE(cosV)
      ALLOCATE(cosV(count))
      cosV(:) = COS(indU(:) * u - indV(:) * v)
    END IF

    CosTransFullMeshF = SUM(aij(:) * cosV(:))

  END FUNCTION CosTransFullMeshF

! ----------------------------------------------------------------------
!> calculate backward cosine fourier transform
! ----------------------------------------------------------------------

  FUNCTION CosTransFullMesh(count, aij, indU, indV, u, v)

!> number of fourier elements
    INTEGER, INTENT(IN) :: count

!> cosine fourier elements
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: aij

!> indices for u and v
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: indU, indV

!> arguments
    REAL(KIND = hp), INTENT(IN) :: u, v

!> resulting function value
    REAL(KIND = hp) :: CosTransFullMesh

    CosTransFullMesh = SUM(aij(:) * COS(indU(:) * u - indV(:) * v))

  END FUNCTION CosTransFullMesh

! ----------------------------------------------------------------------
!> calculate cosine transform at all flux surfaces for given theta
!> and phi
! ----------------------------------------------------------------------

  FUNCTION CosValsOnFluxSurf(n, u, count, fmn, indU, indV, v) RESULT(fVal)

!> number of flux surfaces
    INTEGER, INTENT(IN) :: n

!> u values on surfaces
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: u

!> number of fourier elements
    INTEGER, INTENT(IN) :: count

!> cosine fourier elements
    REAL(KIND = hp), DIMENSION(count, n), INTENT(IN) :: fmn

!> indices for u and v
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: indU, indV

!> selected v cut
    REAL(KIND = hp), INTENT(IN) :: v

!> resulting function values
    REAL(KIND = hp), DIMENSION(n) :: fVal

    INTEGER :: curS

    DO curS = 1, n
      fVal(curS) = CosTransFullMesh(count, fmn(:, curS), indU, indV,&
           u(curS), v)
    END DO

  END FUNCTION CosValsOnFluxSurf

! ----------------------------------------------------------------------
!> calculate backward sine fourier transform with reusing cosine
!> vector
! ----------------------------------------------------------------------

  FUNCTION SinTransFullMeshF2(count, aij, i, j)

!> number of fourier elements (mn_mode or mn_mode_nyq)
    INTEGER, INTENT(IN) :: count

!> cosine fourier elements
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: aij

!> indices for u and v
    INTEGER, INTENT(IN) :: i, j

!> resulting function value
    REAL(KIND = hp) :: SinTransFullMeshF2

    IF (count == mn_mode) THEN
      SinTransFullMeshF2 = SUM(aij(:) * vSin(:, i, j))
    ELSE IF (count == mn_mode_nyq) THEN
      SinTransFullMeshF2 = SUM(aij(:) * vSinN(:, i, j))
    ELSE
      PRINT *, "Wrong number of fourier modes in SinTransFullMeshF2!"
      CALL EXIT(20)
    END IF

  END FUNCTION SinTransFullMeshF2

! ----------------------------------------------------------------------
!> calculate backward sine fourier transform with reusing cosine vector
! ----------------------------------------------------------------------

  FUNCTION SinTransFullMeshF(count, aij, indU, indV, u, v, recalc)

!> number of fourier elements
    INTEGER, INTENT(IN) :: count

!> cosine fourier elements
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: aij

!> indices for u and v
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: indU, indV

!> arguments
    REAL(KIND = hp), INTENT(IN) :: u, v

!> recalc sine vector
    LOGICAL, INTENT(IN) :: recalc

!> resulting function value
    REAL(KIND = hp) :: SinTransFullMeshF

    REAL(KIND = hp), DIMENSION(:), ALLOCATABLE, SAVE :: sinV

    IF (recalc) THEN
      IF (ALLOCATED(sinV)) DEALLOCATE(sinV)
      ALLOCATE(sinV(count))
      sinV(:) = SIN(indU(:) * u - indV(:) * v)
    END IF

    SinTransFullMeshF = SUM(aij(:) * sinV(:))

  END FUNCTION SinTransFullMeshF

! ----------------------------------------------------------------------
!> calculate backward sine fourier transform
! ----------------------------------------------------------------------

  FUNCTION SinTransFullMesh(count, aij, indU, indV, u, v)

!> number of fourier elements
    INTEGER, INTENT(IN) :: count

!> cosine fourier elements
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: aij

!> indices for u and v
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: indU, indV

!> arguments
    REAL(KIND = hp), INTENT(IN) :: u, v

!> resulting function value
    REAL(KIND = hp) :: SinTransFullMesh

    SinTransFullMesh = SUM(aij(:) * SIN(indU(:) * u - indV(:) * v))

  END FUNCTION SinTransFullMesh

! ----------------------------------------------------------------------
!> calculate sine transform at all flux surfaces for given theta and phi
! ----------------------------------------------------------------------

  FUNCTION SinValsOnFluxSurf(n, u, count, fmn, indU, indV, v) RESULT(fVal)

!> number of flux surfaces
    INTEGER, INTENT(IN) :: n

!> u values on surfaces
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: u

!> number of fourier elements
    INTEGER, INTENT(IN) :: count

!> cosine fourier elements
    REAL(KIND = hp), DIMENSION(count, n), INTENT(IN) :: fmn

!> indices for u and v
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: indU, indV

!> selected v cut
    REAL(KIND = hp), INTENT(IN) :: v

!> resulting function values
    REAL(KIND = hp), DIMENSION(n) :: fVal

    INTEGER :: curS

    DO curS = 1, n
      fVal(curS) = SinTransFullMesh(count, fmn(:, curS), indU, indV,&
           u(curS), v)
    END DO

  END FUNCTION SinValsOnFluxSurf

! ----------------------------------------------------------------------
!> calculate (co)sine transform at all flux surfaces for given theta
!> and phi
! ----------------------------------------------------------------------

  FUNCTION ValsOnFluxSurf(n, count, fmn, v) RESULT(fVal)

!> number of flux surfaces
    INTEGER, INTENT(IN) :: n

!> number of fourier elements
    INTEGER, INTENT(IN) :: count

!> cosine fourier elements
    REAL(KIND = hp), DIMENSION(count, n), INTENT(IN) :: fmn

!> indices for u and v
    REAL(KIND = hp), DIMENSION(count), INTENT(IN) :: v

!> resulting function values
    REAL(KIND = hp), DIMENSION(n) :: fVal

    INTEGER :: curS

    DO curS = 1, n
      fVal(curS) = SUM(fmn(:, curS) * v(:))
    END DO

  END FUNCTION ValsOnFluxSurf

! ----------------------------------------------------------------------
!> Lagrange interpolation/extrapolation
!> taken from Numerical Recipes
!> with a given number of supporting points one gets an interpolation or
!> extrapolation using following expression
!>
!>         (x-x2)(x-x3)...(x-xn)             (x-x1)(x-x2)...(x-xn-1)
!> f(x) = ------------------------f(x1)+...+--------------------------f(xn)
!>        (x1-x2)(x1-x3)...(x1-xn)          (xn-x1)(xn-x2)...(xn-xn-1)
! ----------------------------------------------------------------------

  FUNCTION IntExtPol(n, xArr, fArr, x, errEst) RESULT(f)

!> number of supporting points
    INTEGER, INTENT(IN) :: n

!> supporting points
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: xArr, fArr

!> point to interpolate
    REAL(KIND = hp), INTENT(IN) :: x

!> error estimate
    REAL(KIND = hp), INTENT(OUT), OPTIONAL :: errEst

!> resulting function value
    REAL(KIND = hp) :: f

    REAL(KIND = hp), DIMENSION(n) :: c, d
    REAL(KIND = hp) :: h, h1, h2, w, df
    INTEGER :: minI, i, j

    !! look for the nearest supporting point
    minI = MINLOC(ABS(x - xArr(:)), 1)

    !! initialize tableau of c's and d's
    c(:) = fArr(:)
    d(:) = fArr(:)

    !! initial approximation of function value
    f = fArr(minI)
    minI = minI - 1

    !! for each column of the tableau
    DO j = 1, n - 1
      !! loop over the current c's and d's and update them
      DO i = 1, n - j
        h1 = xArr(i) - x
        h2 = xArr(i + j) - x
        w = c(i + 1) - d(i)
        h = xArr(i) - xArr(i + j)
        !! stop if two supporting points are identical
        IF (ABS(h) < tiny) THEN
          PRINT *, "Interpolation was stopped, because 2 supporting points are&
               & identical!"
          RETURN
        END IF
        h = w / h
        !! update the c's and d's
        c(i) = h1 * h
        d(i) = h2 * h
      END DO

      !! After each column in the tablaeu is completed, we decide which
      !! correction, c or d, we want to add to our accumulating value of f,
      !! i.e. which path to take through the tablaeu - forking up or down. We
      !! do this in such a way as to take the most "straight line" route
      !! through the tableau to its apex, updating minI accordingly to keep
      !! track of where we are. This route keeps the partial approximations
      !! centered (insofar as possible) on the target x. The last df added is
      !! thus the error estimate.
      IF (2 * minI < n - j) THEN
        df = c(minI + 1)
      ELSE
        df = d(minI)
        minI = minI - 1
      END IF
      f = f + df
    END DO
    IF (PRESENT(errEst)) errEst = df

  END FUNCTION IntExtPol

! ----------------------------------------------------------------------
!> interpolate or extrapolate with maximal 3rd order, use only the 4 nearest
!> points
! ----------------------------------------------------------------------

  FUNCTION IntExtPol3(n, xArr, fArr, x, f0, errEst) RESULT(f)

!> number of supporting points
    INTEGER, INTENT(IN) :: n

!> supporting points
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: xArr, fArr

!> point to interpolate
    REAL(KIND = hp), INTENT(IN) :: x

!> extra value on axis
    REAL(KIND = hp), INTENT(IN), OPTIONAL :: f0

!> error estimate
    REAL(KIND = hp), INTENT(OUT), OPTIONAL :: errEst

!> resulting function value
    REAL(KIND = hp) :: f

    REAL(KIND = hp), DIMENSION(4) :: xVal, fVal
    INTEGER :: i, s, e, l

    !! look for last point which has smaller x value
    DO i = 1, n
      IF (x < xArr(i)) EXIT
    END DO
    !! set starting point
    IF (PRESENT(f0)) THEN
      s = MAX(0, i - 3)
    ELSE
      s = MAX(1, i - 3)
    END IF
    !! set end point
    e = MIN(n, 3 + s)
    !! use last value or linear extrapolation for outer values
    IF (useLastVal .AND. (i > n - 1)) s = e
    !! calculate number of points for interpolation
    l = e - s + 1
    !! store function values
    IF (s < 1) THEN
      !! if extra value given
      xVal(1) = 0._hp
      xVal(2:l) = xArr(1:e)
      fVal(1) = f0
      fVal(2:l) = fArr(1:e)
    ELSE
      !! no extra value given
      xVal(1:l) = xArr(s:e)
      fVal(1:l) = fArr(s:e)
    END IF

    !! get interpolation value
    IF (PRESENT(errEst)) THEN
      f = IntExtPol(l, xVal(1:l), fVal(1:l), x, errEst)
    ELSE
      f = IntExtPol(l, xVal(1:l), fVal(1:l), x)
    END IF

  END FUNCTION IntExtPol3

! ----------------------------------------------------------------------
!> look for intersection of flux surface and connection line between axis
!> and selected point
!> use theorem on intersecting lines with a secant method
! ----------------------------------------------------------------------

  SUBROUTINE FindInterSectionRZSec(minS, maxS, rAxis, zAxis, delR, delZ, phi,&
       theta, theta2, loopUp)

!> starting flux surface
    INTEGER, INTENT(IN) :: minS

!> last flux surface
    INTEGER, INTENT(IN) :: maxS

!> coordinates of the axis
    REAL(KIND = hp), INTENT(IN) :: rAxis, zAxis

!> selected point
!> distance to axis
    REAL(KIND = hp), INTENT(IN) :: delR, delZ
    
!> phi coordinate
    REAL(KIND = hp), INTENT(IN) :: phi

!> starting angle and resulting angle theta
    REAL(KIND = hp), DIMENSION(:), INTENT(INOUT) :: theta

!> second start angle for secant method
    REAL(KIND = hp), DIMENSION(:), INTENT(IN) :: theta2

!> direction of do loop
    LOGICAL, INTENT(IN), OPTIONAL :: loopUp


    REAL(KIND = hp), DIMENSION(mn_mode) :: xnv, argV
    REAL(KIND = hp), DIMENSION(3) :: h
    REAL(KIND = hp) :: th, f, delta, cR, cZ
    REAL(KIND = hp) :: initTh, th2, f2, cR2, cZ2, thN
    INTEGER :: curS, curIt, startS, endS, dir
    LOGICAL :: loopDir

    !! find solution for theta with a given phi for every flux surfaces with
    !! theorem on intersecting lines
    !!
    !! z - zAxis   SinTransFullMeshF(zmns) - zAxis
    !! --------- = -------------------------------
    !! r - rAxis   CosTransFullMeshF(rmnc) - rAxis
    !!
    !! with the secant method for
    !!
    !! f(theta) = (z - zAxis) * CosTransFullMeshF(rmnc) -
    !!             (r - rAxis) * SinTransFullMeshF(zmns) +
    !!             (r - rAxis) * zAxis - (z - zAxis) * rAxis

    loopDir = .TRUE.
    IF (PRESENT(loopUp)) loopDir = loopUp

    IF (loopDir) THEN
      startS = minS
      endS = maxS
      dir = 1
    ELSE
      startS = maxS
      endS = minS
      dir = -1
    END IF

    !! store theta independent part of fourier coefficients
    xnv(:) = xn(:) * phi

    delta = delR * zAxis - delZ * rAxis

    !! loop over all flux surfaces
    DO curS = startS, endS, dir
      !! actual theta
      th = theta(curS)
      initTh = th
      th2 = theta2(curS)
      !! check if angles do not differ to much
      IF (ABS(th2 - th) > Pi) THEN
        !! check against last flux surface
        h(1) = ABS(MODULO(th2, TwoPi) - TwoPi - th)
        h(2) = ABS(MODULO(th2, TwoPi) - th)
        h(3) = ABS(MODULO(th2, TwoPi) + TwoPi - th)
        th2 = MODULO(th2, TwoPi) + (MINLOC(h(:), 1) - 2) * TwoPi
      END IF
      curIt = 0
      !! run until solution fulfil the convergence criterion
      !! f(theta) = (z - zAxis) * CosTransFullMeshF(rmnc) -
      !!             (r - rAxis) * SinTransFullMeshF(zmns) +
      !!             (r - rAxis) * zAxis - (z - zAxis) * rAxis
      !! argument for the first point
      argV(:) = xm(:) * th - xnv(:)
      !! get the point on flux surface
      cR = SUM(rmnc(:, curS) * COS(argV(:)))
      cZ = SUM(zmns(:, curS) * SIN(argV(:)))
      !! function value for the first point
      f = delZ * cR - delR * cZ + delta
      !! argument for the second point
      argV(:) = xm(:) * th2 - xnv(:)
      !! get the point on flux surface
      cR2 = SUM(rmnc(:, curS) * COS(argV(:)))
      cZ2 = SUM(zmns(:, curS) * SIN(argV(:)))
      !! function value for the second point
      f2 = delZ * cR2 - delR * cZ2 + delta
      !! run until convergence or maximum iterations
      DO WHILE ((ABS(f) > newtonEps) .AND. (curIt <= newtonMax))
        !! count number of iterations
        curIt = curIt + 1
        !! secant step
        IF (ABS(f - f2) < 1D-18) THEN
          thN = (th + th2) / 2
        ELSE
          thN = th - (th - th2) * f / (f - f2)
        END IF
        !! replace second point
        th2 = th
        f2 = f
        !! set new first point
        th = thN
        argV(:) = xm(:) * th - xnv(:)
        !! f(theta) = (z - zAxis) * CosTransFullMeshF(rmnc) -
        !!             (r - rAxis) * SinTransFullMeshF(zmns) +
        !!             (r - rAxis) * zAxis - (z - zAxis) * rAxis
        cR = SUM(rmnc(:, curS) * COS(argV(:)))
        cZ = SUM(zmns(:, curS) * SIN(argV(:)))
        f = delZ * cR - delR * cZ + delta
      END DO
      !! check if convergence is reached
      IF (curIt > newtonMax) THEN
        PRINT '(A, I6, 2ES12.4)', "FindInterSectionRZSec:", curS, rAxis +&
             delR, zAxis + delZ
        th = initTh
      ELSE IF (curS /= startS) THEN
        th = MODULO(th, TwoPi)
        !! check if angles do not differ to much
        IF (ABS(th - theta(curS - dir)) > Pi) THEN
          !! check against last flux surface
          h(1) = ABS(MODULO(th, TwoPi) - TwoPi - theta(curS - dir))
          h(2) = ABS(MODULO(th, TwoPi) - theta(curS - dir))
          h(3) = ABS(MODULO(th, TwoPi) + TwoPi - theta(curS - dir))
          th = MODULO(th, TwoPi) + (MINLOC(h(:), 1) - 2) * TwoPi
        END IF
      ELSE IF (ABS(th) > Pi) THEN
        !! start with an angle near zero
        h(1) = ABS(MODULO(th, TwoPi) - TwoPi)
        h(2) = ABS(MODULO(th, TwoPi))
        h(3) = ABS(MODULO(th, TwoPi) + TwoPi)
        th = MODULO(th, TwoPi) + (MINLOC(h(:), 1) - 2) * TwoPi
      END IF
      !! store solution
      theta(curS) = th
    END DO

  END SUBROUTINE FindInterSectionRZSec

! ----------------------------------------------------------------------
!> look for intersection of flux surface and connection line between axis
!> and selected point
!> use theorem on intersecting lines with a Newton method
! ----------------------------------------------------------------------

  SUBROUTINE FindInterSectionRZNewton(minS, maxS, rAxis, zAxis, delR, delZ, phi,&
       theta, useLast, loopUp)

!> starting flux surface
    INTEGER, INTENT(IN) :: minS

!> last flux surface
    INTEGER, INTENT(IN) :: maxS

!> coordinates of the axis
    REAL(KIND = hp), INTENT(IN) :: rAxis, zAxis

!> selected point
!> distance to axis
    REAL(KIND = hp), INTENT(IN) :: delR, delZ
    
!> phi coordinate
    REAL(KIND = hp), INTENT(IN) :: phi

!> resulting angle theta
    REAL(KIND = hp), DIMENSION(:), INTENT(INOUT) :: theta

!> use theta value from previous flux surface as initial value for theta
    LOGICAL, INTENT(IN), OPTIONAL :: useLast

!> direction of do loop
    LOGICAL, INTENT(IN), OPTIONAL :: loopUp


    REAL(KIND = hp), DIMENSION(mn_mode) :: xnv, argV, sinV, cosV
    REAL(KIND = hp), DIMENSION(3) :: h
    REAL(KIND = hp) :: th, angle, f, fD, delta, newRelax, cR, cZ
    REAL(KIND = hp) :: initTh
    INTEGER :: curS, curIt, curRec, startS, endS, dir
    LOGICAL :: restart, loopDir, useLastVal

    !! find solution for theta with a given phi for every flux surfaces with
    !! theorem on intersecting lines
    !!
    !! z - zAxis   SinTransFullMeshF(zmns) - zAxis
    !! --------- = -------------------------------
    !! r - rAxis   CosTransFullMeshF(rmnc) - rAxis
    !!
    !! with the newton method for
    !!
    !! f(theta) = (z - zAxis) * CosTransFullMeshF(rmnc) -
    !!             (r - rAxis) * SinTransFullMeshF(zmns) +
    !!             (r - rAxis) * zAxis - (z - zAxis) * rAxis

    useLastVal = .FALSE.
    IF (PRESENT(useLast)) useLastVal = useLast
    loopDir = .TRUE.
    IF (PRESENT(loopUp)) loopDir = loopUp

    IF (loopDir) THEN
      startS = minS
      endS = maxS
      dir = 1
    ELSE
      startS = maxS
      endS = minS
      dir = -1
    END IF

    !! store theta independent part of fourier coefficients
    xnv(:) = xn(:) * phi

    delta = delR * zAxis - delZ * rAxis
    !! initialize theta with normal angle
    angle = MODULO(ATAN2(delZ, delR), TwoPi)

    !! loop over all flux surfaces
    DO curS = startS, endS, dir
      !! check if last value should be used
      IF (useLastVal) THEN
        IF (curS == startS) THEN
          initTh = theta(ABS(curS))
        ELSE
          initTh = theta(ABS(curS - dir))
        END IF
      ELSE
        !! original theta
        initTh = theta(ABS(curS))
      END IF
      !! actual theta
      th = initTh
      curIt = 0
      curRec = 1
      restart = .FALSE.
      newRelax = 0.8_hp
      !! run until solution fulfil the convergence criterion
      !! f(theta) = (z - zAxis) * CosTransFullMeshF(rmnc) -
      !!             (r - rAxis) * SinTransFullMeshF(zmns) +
      !!             (r - rAxis) * zAxis - (z - zAxis) * rAxis
      argV(:) = xm(:) * th - xnv(:)
      sinV(:) = SIN(argV(:))
      cosV(:) = COS(argV(:))
      cR = SUM(rmnc(:, curS) * cosV(:))
      cZ = SUM(zmns(:, curS) * sinV(:))
      f = delZ * cR - delR * cZ + delta
      !! run until convergence or maximum iterations
      DO WHILE ((ABS(f) > newtonEps) .AND. (curIt <= newtonMax))
        curIt = curIt + 1
        !! calculate df/du =
        !!  (z - zAxis) * SinTransFullMeshF(dRdUmns) -
        !!   (r - rAxis) * CosTransFullMeshF(dZdUmnc)
        fD = delZ * SUM(dRdUmns(:, curS) * sinV(:)) - delR *&
             SUM(dZdUmnc(:, curS) * cosV(:)) + 1E-28_hp
        !! Newton step
        th = th - newRelax * f / fD
        argV(:) = xm(:) * th - xnv(:)
        sinV(:) = SIN(argV(:))
        cosV(:) = COS(argV(:))
        !! check for divergence
        IF (ABS(th) > 1E2_hp) THEN
          curIt = newtonMax + 1
        ELSE
          !! f(theta) = (z - zAxis) * CosTransFullMeshF(rmnc) -
          !!             (r - rAxis) * SinTransFullMeshF(zmns) +
          !!             (r - rAxis) * zAxis - (z - zAxis) * rAxis
          cR = SUM(rmnc(:, curS) * cosV(:))
          cZ = SUM(zmns(:, curS) * sinV(:))
          f = delZ * cR - delR * cZ + delta
        END IF
        !! check if convergence is reached
        IF (curIt > newtonMax) THEN
          IF (curS /= 1) THEN
            curRec = curRec + 1
            IF (curRec <= maxRec) THEN
              restart = .TRUE.
            END IF
          END IF
        ELSE IF (ABS(f) < newtonEps) THEN
          !! check if angles on the same side
          IF (ABS(MODULO(ATAN2(cZ - zAxis, cR - rAxis), TwoPi) -&
               angle) > 1) THEN
            curRec = curRec + 1
            IF (curRec <= maxRec) THEN
              restart = .TRUE.
            END IF
          END IF
        END IF
        !! restart Newton iteration with new relaxation
        IF (restart) THEN
          curIt = 0
          th = initTh
          newRelax = newRelax * relax
          argV(:) = xm(:) * th - xnv(:)
          sinV(:) = SIN(argV(:))
          cosV(:) = COS(argV(:))
          cR = SUM(rmnc(:, curS) * cosV(:))
          cZ = SUM(zmns(:, curS) * sinV(:))
          f = delZ * cR - delR * cZ + delta
          restart = .FALSE.
        END IF
      END DO
      !! check if convergence is reached
      IF (curIt > newtonMax) THEN
        PRINT '(A, I6, 2ES12.4)', "FindInterSectionRZNewton:", curS, rAxis +&
             delR, zAxis + delZ
        th = initTh
      ELSE IF (curS /= startS) THEN
        th = MODULO(th, TwoPi)
        !! check if angles do not differ to much
        IF (ABS(th - theta(ABS(curS - dir))) > Pi) THEN
          !! check against last flux surface
          h(1) = ABS(MODULO(th, TwoPi) - TwoPi - theta(ABS(curS - dir)))
          h(2) = ABS(MODULO(th, TwoPi) - theta(ABS(curS - dir)))
          h(3) = ABS(MODULO(th, TwoPi) + TwoPi - theta(ABS(curS - dir)))
          th = MODULO(th, TwoPi) + (MINLOC(h(:), 1) - 2) * TwoPi
        END IF
      ELSE IF (ABS(th) > Pi) THEN
        !! start with an angle near zero
        h(1) = ABS(MODULO(th, TwoPi) - TwoPi)
        h(2) = ABS(MODULO(th, TwoPi))
        h(3) = ABS(MODULO(th, TwoPi) + TwoPi)
        th = MODULO(th, TwoPi) + (MINLOC(h(:), 1) - 2) * TwoPi
      END IF
      !! store solution
      theta(ABS(curS)) = th
    END DO

  END SUBROUTINE FindInterSectionRZNewton

! ----------------------------------------------------------------------
!> fit data to
!> f(s) = a s^2 + b s + c s^(1/2)
! ----------------------------------------------------------------------

  SUBROUTINE CalcLinFitFromDerivative(n, x, f, a, b)

!> number of points
    INTEGER, INTENT(IN) :: n

!> points
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: x

!> function values
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: f

!> fitting parameter (ax+b)
    REAL(KIND = hp), INTENT(OUT) :: a, b

    REAL(KIND = hp), DIMENSION(n) :: sqrtX
    REAL(KIND = hp) :: s1, sy, sx12, sx12y, sx, sx32, sx32y, sx2, sx3, d

    sqrtX(:) = SQRT(x(:))
    s1 = REAL(n, hp)
    sy = SUM(f(:))
    sx12 = SUM(sqrtX(:))
    sx12y = SUM(sqrtX(:) * f(:))
    sx = SUM(x(:))
    sx32 = SUM(x(:) * sqrtX(:))
    sx32y = SUM(x(:) * sqrtX(:) * f(:))
    sx2 = SUM(x(:) * x(:))
    sx3 = SUM(x(:) * x(:) * x(:))
    d = sx32 * (sx12 * sx2 - sx * sx32) + sx12 * (sx2 * sx32 -  sx12 * sx3) +&
         s1 * (sx * sx3 - sx2 * sx2)

    a = 1.5_hp * (sx32 * (sx12 * sx12y - sx * sy) + sx12 *&
         (sx2 * sy - sx12 * sx32y) + s1 * (sx * sx32y - sx12y * sx2)) / d

    b = 0.5_hp * (sx32 * (sx2 * sy - sx12y * sx32) + sx12 *&
         (sx32 * sx32y - sx3 * sy) + s1 * (sx12y * sx3 - sx2 * sx32y)) / d

  END SUBROUTINE CalcLinFitFromDerivative

! ----------------------------------------------------------------------
!> fit data to linear function
! ----------------------------------------------------------------------

  SUBROUTINE CalcLinFit(n, x, f, a, b)

!> number of points
    INTEGER, INTENT(IN) :: n

!> points
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: x

!> function values
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: f

!> fitting parameter (ax+b)
    REAL(KIND = hp), INTENT(OUT) :: a, b

    REAL(KIND = hp) :: sx, sy, sxy, sx2, d

    sx = SUM(x(:))
    sy = SUM(f(:))
    sxy = SUM(x(:) * f(:))
    sx2 = SUM(x(:) * x(:))
    d = REAL(n, hp) * sx2 - sx * sx
    a = (REAL(n, hp) * sxy - sx * sy) / d
    b = (sx2 * sy - sx * sxy) / d

  END SUBROUTINE CalcLinFit

! ----------------------------------------------------------------------
!> compute "normal" right-handed cross product of two vectors
! ----------------------------------------------------------------------

  FUNCTION Cross(vec1, vec2) RESULT(res)

!> two vectors
    REAL(KIND = hp), DIMENSION(3), INTENT(IN) :: vec1, vec2

!> cross product
    REAL(KIND = hp), DIMENSION(3) :: res

    res(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    res(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    res(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

  END FUNCTION Cross

! ----------------------------------------------------------------------
!> compute left-handed cross product of two vectors
! ----------------------------------------------------------------------

  FUNCTION CrossL(vec1, vec2) RESULT(res)

!> two vectors
    REAL(KIND = hp), DIMENSION(3), INTENT(IN) :: vec1, vec2

!> cross product
    REAL(KIND = hp), DIMENSION(3) :: res

    res(1) = vec1(3) * vec2(2) - vec1(2) * vec2(3)
    res(2) = vec1(1) * vec2(3) - vec1(3) * vec2(1)
    res(3) = vec1(2) * vec2(1) - vec1(1) * vec2(2)

  END FUNCTION CrossL

  SUBROUTINE LUDecompose(n, a, ind)

!> size of matrix
    INTEGER, INTENT(IN) :: n

!> matrix (in) and lu decomposition (out)
    REAL(KIND = hp), DIMENSION(n, n), INTENT(INOUT) :: a

!> permutation matrix
    INTEGER, DIMENSION(n), INTENT(OUT) :: ind

    REAL(KIND = hp) :: aMax, h
    REAL(KIND = hp), DIMENSION(n) :: vv
    INTEGER :: i, j, k, iMax

    DO i = 1, n
      vv(i) = 1._hp / MAXVAL(ABS(a(i, :)))
    END DO

    DO j = 1, n
      DO i = 1, j - 1
        h = a(i, j)
        DO k = 1, i - 1
          h = h - a(i, k) * a(k, j)
        END DO
        a(i, j) = h
      END DO
      aMax = 0
      DO i = j, n
        h = a(i, j)
        DO k = 1, j - 1
          h = h - a(i, k) * a(k, j)
        END DO
        a(i, j) = h
        h = vv(i) * ABS(a(i, j))
        IF (h >= aMax) THEN
          iMax = i
          aMax = h
        END IF
      END DO

      IF (j /= iMax) THEN
        DO k = 1, n
          h = a(iMax, k)
          a(iMax, k) = a(j, k)
          a(j, k) = h
        END DO
        vv(iMax) = vv(j)
      END IF

      ind(j) = iMax
      IF (ABS(a(j, j)) < tiny) a(j, j) = tiny

      IF (j /= n) THEN
        a((j + 1):n, j) = a((j + 1):n, j) / a(j, j)
      END IF
    END DO

  END SUBROUTINE LUDecompose

  SUBROUTINE LUBacksolve(n, a, ind, b)

!> size of matrix
    INTEGER, INTENT(IN) :: n

!> lu decomposition
    REAL(KIND = hp), DIMENSION(n, n), INTENT(IN) :: a

!> permutation matrix
    INTEGER, DIMENSION(n), INTENT(IN) :: ind

!> right-hand side
    REAL(KIND = hp), DIMENSION(n), INTENT(INOUT) :: b

    INTEGER :: i, j
    REAL(KIND = hp) :: h

    j = 0
    !! forward substitution
    DO i = 1, n
      h = b(ind(i))
      b(ind(i)) = b(i)
      IF (j /= 0) THEN
        IF (i > j) h = h - DOT_PRODUCT(a(i, j:(i - 1)), b(j:(i - 1)))
      ELSE IF (ABS(h) > tiny) THEN
        j = i
      END IF
      b(i) = h
    END DO

    !! backward substitution
    DO i = n, 1, -1
      h = b(i)
      IF (i < n) THEN
        h = h - DOT_PRODUCT(a(i, (i + 1):n), b((i + 1):n))
      END IF
      b(i) = h / a(i, i)
    END DO

  END SUBROUTINE LUBacksolve

  SUBROUTINE GMRES_Solve(order, n, a, b, x, stopVal, startVec, resOut, err)

!> truncation order of GMRES(order)
    INTEGER, INTENT(IN) :: order

!> order of matrix
    INTEGER, INTENT(IN) :: n

!> non-zero values of the matrix
    REAL(KIND = hp), DIMENSION(n, n), INTENT(IN) :: a

!> right-hand side of the equation
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: b

!> solution of the equations
    REAL(KIND = hp), DIMENSION(n), INTENT(OUT) :: x

!> stopping residuum
    REAL(KIND = hp), INTENT(IN) :: stopVal

!> starting vector for the iteration
    REAL(KIND = hp), DIMENSION(n), INTENT(IN), OPTIONAL :: startVec

!> residuum after iteration
    REAL(KIND = hp), INTENT(OUT), OPTIONAL :: resOut

!> error during calculations
    INTEGER, INTENT(OUT), OPTIONAL :: err

    REAL(KIND = hp), DIMENSION(n) :: r, rP, w, wP, diag
    REAL(KIND = hp), DIMENSION(n, order + 1) :: v
    REAL(KIND = hp), DIMENSION(order + 1, order) :: h
    REAL(KIND = hp), DIMENSION(order + 1) :: s
    REAL(KIND = hp), DIMENSION(order) :: co, si, y
    REAL(KIND = hp) :: alpha, h1, resi, resiOld
    INTEGER :: curIter, i, j, k

    !! get diagonal
    DO i = 1, n
      IF (ABS(a(i, i)) > tiny) THEN
        diag(i) = a(i, i)
      ELSE
        diag(i) = tiny
      END IF
    END DO

    !! initialization
    IF (PRESENT(startVec)) THEN
      x(:) = startVec(:)
    ELSE
      x(:) = b(:) / diag(:)
    END IF

    curIter = 0
    IF (PRESENT(err)) err = -1

    DO WHILE (curIter < itSolMaxIter)
      !! rP = A x (distributed)
      rP(:) = MATMUL(a(:, :), x(:))
      !! r = b - A x
      r(:) = b(:) - rP(:)
      !! s = ||r|| e_1
      s(1) = SQRT(DOT_PRODUCT(r(:), r(:)))
      resi = s(1)
      !! v(1) = r / ||r||
      v(:, 1) = r(:) / (s(1) + tiny)
      DO i = 1, MIN(itSolMaxIter, order)
        curIter = curIter + 1
        wP(:) = SGS(n, a, diag(:), v(:, i))
        !! w = A w (distributed)
        w(:) = MATMUL(a(:, :), wP(:))
        DO k = 1, i
          !! h_(k,i) = (w, v(k))
          h(k, i) = DOT_PRODUCT(w(:), v(:, k))
          !! w = w - h_(k,i) v(k)
          w(:) = w(:) - h(k, i) * v(:, k)
        END DO
        !! h_(i+1,i) = ||w||
        h(i + 1, i) = SQRT(DOT_PRODUCT(w(:), w(:)))
        !! v(i+1) = w / ||w||
        v(:, i + 1) = w(:) / (h(i + 1, i) + tiny)
        !! J_(i - 1) ... J_1 (h_(1, i), ..., h_(i+1,i))
        !!       ( co_1 si_1 0  .. )
        !!       (-si_1 co_1 0  .. )
        !! J_1 = (  0    0   1  .. ), ...
        !!       ( ..    ..  .. .. )
        !!       (  0    0   .. 1  )
        DO k = 1, i - 1
          h1 = h(k, i)
          h(k, i)     =  co(k) * h1 + si(k) * h(k + 1, i)
          h(k + 1, i) = -si(k) * h1 + co(k) * h(k + 1, i)
        END DO
        !! co_i = h_(i,i) / sqrt(h_(i, i)^2+h_(i+1,i)^2)
        !! si_i = h_(i+1,i) / sqrt(h_(i, i)^2+h_(i+1,i)^2)
        alpha = SQRT(h(i, i) * h(i, i) + h(i + 1, i) * h(i + 1, i))
        co(i) = h(i, i) / (alpha + tiny)
        si(i) = h(i + 1, i) / (alpha + tiny)
        !! s = J_i s
        s(i + 1) = -si(i) * s(i)
        s(i) = co(i) * s(i)
        !! update diagonal element
        h(i, i) = co(i) * h(i, i) + si(i) * h(i + 1, i)
        !! truncation test
        resiOld = resi
        resi = ABS(s(i + 1))
        IF (numDebug .AND. (Me_rank == 0)) THEN
          PRINT '(I6, ES25.15)', curIter, resi
        END IF
        IF ((curIter > minIter) .AND. (ABS(resi - resiOld) < stopVal)) EXIT
      END DO
      !! correct i if previous loop ends normally
      IF (i > MIN(itSolMaxIter, order)) THEN
        i = MIN(itSolMaxIter, order)
      END IF
      !! Solve H y = s'
      y(i) = s(i) / (h(i, i) + tiny)
      DO j = i - 1, 1, -1
        y(j) = (s(j) - DOT_PRODUCT(h(j, (j + 1):i), y((j + 1):i))) /&
             (h(j, j) + tiny)
      END DO
      !! w = Sum(y_i v(i))
      wP(:) = y(1) * v(:, 1)
      DO j = 2, i
        wP(:) = wP(:) + y(j) * v(:, j)
      END DO
      w(:) = SGS(n, a, diag(:), wP(:))
      !! x = x + w
      x(:) = x(:) + w(:)
      !! convergence test
      IF (ABS(resi - resiOld) < stopVal) THEN
        IF (PRESENT(err)) err = curIter
        EXIT
      END IF
    END DO

    IF (PRESENT(resOut)) resOut = resi
    IF (PRESENT(err)) err = curIter

  END SUBROUTINE GMRES_Solve

  FUNCTION SGS(n, a, diag, vec) RESULT(res)

    !! order of matrix
    INTEGER, INTENT(IN) :: n

    !! non-zero values of the matrix
    REAL(KIND = hp), DIMENSION(n, n), INTENT(IN) :: a

    !! diagonal of matrix
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: diag

    !! right-hand side
    REAL(KIND = hp), DIMENSION(n), INTENT(IN) :: vec

    !! resulting vector
    REAL(KIND = hp), DIMENSION(n) :: res

    INTEGER :: i

    res(:) = vec(:)

    !! Symmetric Gauss Seidel (SGS)
    !! M = (D + L) D^(-1) (D + U)
    !! first part - solve (I + L D^(-1)) x' = y
    DO i = 2, n
      !! calculate resulting rows on current node
      res(i) = res(i) - DOT_PRODUCT(a(i, 1:(i - 1)), res(1:(i - 1)) /&
           diag(1:(i - 1)))
    END DO

    !! second part - solve (D + U) x = x'
    res(n) = res(n) / diag(n)
    DO i = n - 1, 1, -1
      !! calculate resulting rows on current node
      res(i) = (res(i) - DOT_PRODUCT(a(i, (i + 1):n), res((i + 1):n))) /&
           diag(i)
    END DO

  END FUNCTION SGS

END MODULE MOD_VMEC_Mappings
