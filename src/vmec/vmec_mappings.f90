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

!> number of grid points in u direction
!  INTEGER :: nu = 100

!> number of grid points in v direction
!  INTEGER :: nv = 100

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

!> iota on full mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: iotaf

!> pressure on full mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: presf

!> toroidal flux on full mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: phi

!> poloidal flux on full mesh
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: chi

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
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: bmnc_nyq

  !! B^u, B^v (cosine components (read on half mesh, interpolated on full
  !! mesh))
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: bsupumnc_nyq, bsupvmnc_nyq
!
  !! dlambda/du, dlambda/dv (cosine components on full mesh)
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: lUmnc, lVmnc

  !! dkappa/du, dkappa/dv (cosine components on full mesh)
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: kUmnc, kVmnc

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
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradSRmnc, gradSPmns
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradSZmns

!  !! fourier components of grad theta
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradTRmns, gradTPmnc
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradTZmnc

  !! fourier components of |B|
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: absBmnc

  !! fourier components of |B| * |sqrtG|
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: absBsmnc

  !! fourier components of d|B|/ds
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dBdSmnc

  !! fourier components of d(|B| * |sqrtG|)/ds
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dBdSsmnc

  !! fourier components of b
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: bRmns, bZmnc, bPmnc

  !! fourier components of dB/du, dB/dv
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dBdUmns, dBdVmns

  !! fourier components of scaled dB/du, dB/dv
  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: dBdUsmns, dBdVsmns

  !! fourier components of grad B
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradBRmnc, gradBPmns
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: gradBZmns

  !! fourier components of b x grad B / |B|
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: bgradBRmns, bgradBPmnc
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: bgradBZmnc

  !! fourier components of div b
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: divBmns

  !! fourier components of rot B parallel
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: rotBparAbsmnc

  !! fourier components of rot B perpendicular
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: rotBperpRmns, rotBperpPmnc
!  REAL(KIND = hp), DIMENSION(:, :), ALLOCATABLE :: rotBperpZmnc

  !! poloidal and toroidal current
!  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: iPol, jTor

  !! precomputed cosinus and sinus values
  REAL(KIND = hp), DIMENSION(:, :, :), ALLOCATABLE :: vCos, vSin, vCosN, vSinN

  !! precomputed integration points
  REAL(KIND = hp), DIMENSION(:), ALLOCATABLE :: uIv, vJv


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
!    ALLOCATE(kUmnc(mn_mode, radius), kVmnc(mn_mode, radius), stat = aStat)
!    aError = aError + aStat
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
!    ALLOCATE(gradSRmnc(mn_mode, radius), gradTRmns(mn_mode, radius), stat =&
!         aStat)
!    aError = aError + aStat
!    ALLOCATE(gradSPmns(mn_mode, radius), gradTPmnc(mn_mode, radius), stat =&
!         aStat)
!    aError = aError + aStat
!    ALLOCATE(gradSZmns(mn_mode, radius), gradTZmnc(mn_mode, radius), stat =&
!         aStat)
!    aError = aError + aStat
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
!    ALLOCATE(bRmns(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(bZmnc(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(bPmnc(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
    ALLOCATE(dBdUmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dBdVmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dBdUsmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(dBdVsmns(mn_mode_nyq, radius), stat = aStat)
    aError = aError + aStat
!    ALLOCATE(gradBRmnc(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(gradBZmns(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(gradBPmns(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(bgradBRmns(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(bgradBZmnc(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(bgradBPmnc(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(divBmns(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(rotBparAbsmnc(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(rotBperpRmns(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(rotBperpZmnc(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(rotBperpPmnc(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
    ALLOCATE(sVal(stInd:radius), theta(radius), sValVMEC(nFluxVMEC), stat = aStat)
    aError = aError + aStat
    ALLOCATE(axisM(nAxis), axisN(nAxis), stat = aStat)
    aError = aError + aStat
!    ALLOCATE(iPol(radius), jTor(radius), stat = aStat)
!    aError = aError + aStat
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

    !! allocate memory for fourier arrays
    aError = 0
    ALLOCATE(xm(mn_mode), xn(mn_mode), stat = aStat)
    aError = aError + aStat
    ALLOCATE(xm_nyq(mn_mode_nyq), xn_nyq(mn_mode_nyq), stat = aStat)
    aError = aError + aStat
!    ALLOCATE(iotas(stInd:radius), stat = aStat)
!    aError = aError + aStat
    ALLOCATE(iotaf(radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(presf(radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(phi(radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(chi(radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(phipf(radius), stat = aStat)
    aError = aError + aStat
!    ALLOCATE(icur(radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(jcur(radius), stat = aStat)
!    aError = aError + aStat
    ALLOCATE(rmnc(mn_mode, stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(zmns(mn_mode, stInd:radius), stat = aStat)
    aError = aError + aStat
    ALLOCATE(lmns(mn_mode, stInd:radius), stat = aStat)
    aError = aError + aStat
!    ALLOCATE(kmns(mn_mode, radius), stat = aStat)
!    aError = aError + aStat
    ALLOCATE(gmnc(mn_mode_nyq, stInd:radius), stat = aStat)
    aError = aError + aStat
!    ALLOCATE(dgdSmnc(mn_mode_nyq, stInd:radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(bmnc_nyq(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(bsupumnc_nyq(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat
!    ALLOCATE(bsupvmnc_nyq(mn_mode_nyq, radius), stat = aStat)
!    aError = aError + aStat

    IF (aError /= 0) THEN
      PRINT *, "Allocation failure in subroutine ReadVmecOutput!"
      CALL EXIT(4)
    END IF

    !! initialize some arrays
!    iotas(:) = 0
    iotaf(:) = 0
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
!    !! read iotas
!    ioError = NF_INQ_VARID(ncid, "iotas", id)
!    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
!         (/ nFluxVMEC /), iotas(1:))
    !! read iotaf
    ioError = NF_INQ_VARID(ncid, "iotaf", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), iotaf(:))
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
    !! read chi
    ioError = NF_INQ_VARID(ncid, "chi", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), chi(:))
    !! scale poloidal flux to get internal VMEC chi
    chi(:nFluxVMEC) = REAL(signgs, hp) * chi(:nFluxVMEC) / TwoPi
    !! read phipf
    ioError = NF_INQ_VARID(ncid, "phipf", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
         (/ nFluxVMEC /), phipf(:))
    !! scale toroidal flux to get internal VMEC Phi
    phipf(:nFluxVMEC) = REAL(signgs, hp) * phipf(:nFluxVMEC) / TwoPi
!    !! read toroidal current
!    ioError = NF_INQ_VARID(ncid, "jcuru", id)
!    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
!         (/ nFluxVMEC /), icur(:))
!    !! read poloidal current
!    ioError = NF_INQ_VARID(ncid, "jcurv", id)
!    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1 /),&
!         (/ nFluxVMEC /), jcur(:))
    !! read R_mn
    ioError = NF_INQ_VARID(ncid, "rmnc", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), rmnc(:, 1:))
    !! read Z_mn
    ioError = NF_INQ_VARID(ncid, "zmns", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), zmns(:, 1:))
    !! read lambda_mn on HALF MESH
    ioError = NF_INQ_VARID(ncid, "lmns", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/ mn_mode,&
         nFluxVMEC /), lmns(:, 1:))
    !! read jacobian_mn on HALF MESH!!
    ioError = NF_INQ_VARID(ncid, "gmnc", id)
    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
         mn_mode_nyq, nFluxVMEC /), gmnc(:, 1:))
!    !! read b_mn
!    ioError = NF_INQ_VARID(ncid, "bmnc", id)
!    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
!         mn_mode_nyq, nFluxVMEC /), bmnc_nyq(:, :))
!    !! read B^u_mn
!    ioError = NF_INQ_VARID(ncid, "bsupumnc", id)
!    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
!         mn_mode_nyq, nFluxVMEC /), bsupumnc_nyq(:, :))
!    !! read B^v_mn
!    ioError = NF_INQ_VARID(ncid, "bsupvmnc", id)
!    ioError = ioError + NF_GET_VARA_DOUBLE(ncid, id, (/ 1, 1 /), (/&
!         mn_mode_nyq, nFluxVMEC /), bsupvmnc_nyq(:, :))

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
    CALL HalfMeshToFullMesh2D(mn_mode_nyq, nFluxVMEC, bmnc_nyq(:, :nFluxVMEC))
    CALL HalfMeshToFullMesh2D(mn_mode_nyq, nFluxVMEC, bsupumnc_nyq(:, :nFluxVMEC))
    CALL HalfMeshToFullMesh2D(mn_mode_nyq, nFluxVMEC, bsupvmnc_nyq(:, :nFluxVMEC))

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
      dBdUmns(:, curS) = -xm_nyq(:) * bmnc_nyq(:, curS)
      !! dB/dv
      dBdVmns(:, curS) = xn_nyq(:) * bmnc_nyq(:, curS)
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
    sVal(1) = 1E-9

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
    bmnc_nyq(:, :) = absBmnc(:, :) / REAL(intPointsU * intPointsV, hp)
    !! get fourier coefficients of the radial derivative of |B|
    DO curS = 1, radius
      DO i = 1, mn_mode_nyq
        !! scaled |B|
        dBdSsmnc(i, curS) = GetRadialDerivative(radius, curS, absBsmnc(i, :))
        !! |B|
        dBdSmnc(i, curS) = GetRadialDerivative(radius, curS, bmnc_nyq(i, :))
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
