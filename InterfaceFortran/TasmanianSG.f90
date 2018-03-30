!==================================================================================================================================================================================
! Copyright (c) 2017, Miroslav Stoyanov
!
! This file is part of
! Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
!    and the following disclaimer in the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
!    or promote products derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
! OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
! OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
! THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
! COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
! THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
! IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
!==================================================================================================================================================================================
Module TasmanianSG
IMPLICIT NONE
PUBLIC :: tsgInitialize,        &
          tsgFinalize,          &
          tsgNewGridID,         &
          tsgFreeGridID,        &
!===== ABOVE ARE FORTRAN SPECIFIC FUNCTIONS ======!
!===== BELOW ARE WRAPPERS FOR STANDARD C++ CLASS =!
          tsgGetVersionMajor,   &
          tsgGetVersionMinor,   &
          tsgGetLicense,        &
          tsgMakeGlobalGrid,    &
          tsgMakeSequenceGrid,  &
          tsgMakeLocalPolynomialGrid, &
          tsgMakeWaveletGrid,   &
          tsgCopyGrid,          &
          tsgUpdateGlobalGrid,  &
          tsgUpdateSequenceGrid,&
          tsgRead,              &
          tsgWrite,             &
          tsgGetAlpha,          &
          tsgGetBeta,           &
          tsgGetOrder,          &
          tsgGetNumDimensions,  &
          tsgGetNumOutputs,     &
          tsgGetRule,           &
          tsgGetNumLoaded,      &
          tsgGetNumNeeded,      &
          tsgGetNumPoints,      &
          tsgGetLoadedPoints,   &
          tsgGetNeededPoints,   &
          tsgGetPoints,         &
          tsgGetLoadedPointsStatic,   &
          tsgGetNeededPointsStatic,   &
          tsgGetPointsStatic,         &
          tsgGetQuadratureWeights,          &
          tsgGetQuadratureWeightsStatic,    &
          tsgGetInterpolationWeights,       &
          tsgGetInterpolationWeightsStatic, &
          tsgLoadNeededPoints, &
          tsgEvaluate,         &
          tsgEvaluateFast,     &
          tsgEvaluateBatch,    &
          tsgIntegrate,        &
          tsgSetDomainTransform,   &
          tsgIsSetDomainTransfrom, &
          tsgClearDomainTransform, &
          tsgGetDomainTransform,   &
          tsgSetAnisotropicRefinement,        &
          tsgEstimateAnisotropicCoefficients, &
          tsgSetGlobalSurplusRefinement,      &
          tsgSetLocalSurplusRefinement,       &
          tsgClearRefinement,                 &
          tsgSetConformalTransformASIN,   &
          tsgIsSetConformalTransformASIN, &
          tsgClearConformalTransform,     &
          tsgGetConformalTransformASIN,   &
          tsgPrintStats, &
!===== DO NOT USE THESE FUNCTIONS DIRECTLY ======!
!===== THESE ARE NEEDED TO LINK TO C++ ==========!
          !tsgReceiveInt,    &
          !tsgReceiveScalar, &
          tsgReceiveString, &
          tsgReceiveVector, &
          tsgReceiveMatrix, &
          tsg_clenshaw_curtis     ,  tsg_rule_clenshaw_curtis_zero, &
          tsg_chebyshev           ,  tsg_chebyshev_odd        , &
          tsg_gauss_legendre      ,  tsg_gauss_legendreodd    , & 
          tsg_gauss_patterson     ,  tsg_leja                 , &
          tsg_lejaodd             ,  tsg_rleja                , &
          tsg_rleja_odd           ,  tsg_rleja_double_2       , &
          tsg_rleja_double_4      ,  tsg_rleja_shifted        , &
          tsg_rleja_shifted_even  ,  tsg_rleja_shifted_double , &
          tsg_max_lebesgue        ,  tsg_max_lebesgue_odd     , &
          tsg_min_lebesgue        ,  tsg_min_lebesgue_odd     , &
          tsg_min_delta           ,  tsg_min_delta_odd        , &
          tsg_gauss_chebyshev_1   ,  tsg_gauss_chebyshev_1_odd, &
          tsg_gauss_chebyshev_2   ,  tsg_gauss_chebyshev_2_odd, &
          tsg_fejer2              ,  tsg_gauss_gegenbauer     , &
          tsg_gauss_gegenbauer_odd,  tsg_gauss_jacobi         , &
          tsg_gauss_jacobi_odd    ,  tsg_gauss_laguerre       , &
          tsg_gauss_laguerre_odd  ,  tsg_gauss_hermite        , &
          tsg_gauss_hermite_odd   ,  tsg_custom_tabulated     , &
          tsg_localp              ,  tsg_localp_zero          , &
          tsg_semi_localp         ,  tsg_wavelet              , &
          tsg_level       ,  tsg_curved      , &
          tsg_iptotal     ,  tsg_ipcurved    , &
          tsg_qptotal     ,  tsg_qpcurved    , &
          tsg_hyperbolic  ,  tsg_iphyperbolic, &
          tsg_qphyperbolic,  tsg_tensor      , &
          tsg_iptensor    ,  tsg_qptensor    , &
          tsg_classic    ,  tsg_parents_first, &
          tsg_directional,  tsg_fds          , &
          tsg_acc_none    , tsg_acc_cpu_blas , tsg_acc_gpu_cublas, &
          tsg_acc_gpu_cuda, tsg_acc_gpu_magma
                        
  INTEGER, PARAMETER :: tsg_clenshaw_curtis      =  1,  tsg_rule_clenshaw_curtis_zero = 2, &
                        tsg_chebyshev            =  3,  tsg_chebyshev_odd         =  4, &
                        tsg_gauss_legendre       =  5,  tsg_gauss_legendreodd     =  6, & 
                        tsg_gauss_patterson      =  7,  tsg_leja                  =  8, &
                        tsg_lejaodd              =  9,  tsg_rleja                 = 10, &
                        tsg_rleja_odd            = 11,  tsg_rleja_double_2        = 12, &
                        tsg_rleja_double_4       = 13,  tsg_rleja_shifted         = 14, &
                        tsg_rleja_shifted_even   = 15,  tsg_rleja_shifted_double  = 16, &
                        tsg_max_lebesgue         = 17,  tsg_max_lebesgue_odd      = 18, &
                        tsg_min_lebesgue         = 19,  tsg_min_lebesgue_odd      = 20, &
                        tsg_min_delta            = 21,  tsg_min_delta_odd         = 22, &
                        tsg_gauss_chebyshev_1    = 23,  tsg_gauss_chebyshev_1_odd = 24, &
                        tsg_gauss_chebyshev_2    = 25,  tsg_gauss_chebyshev_2_odd = 26, &
                        tsg_fejer2               = 27,  tsg_gauss_gegenbauer      = 28, &
                        tsg_gauss_gegenbauer_odd = 29,  tsg_gauss_jacobi          = 30, &
                        tsg_gauss_jacobi_odd     = 31,  tsg_gauss_laguerre        = 32, &
                        tsg_gauss_laguerre_odd   = 33,  tsg_gauss_hermite         = 34, &
                        tsg_gauss_hermite_odd    = 35,  tsg_custom_tabulated      = 36, &
                        tsg_localp               = 37,  tsg_localp_zero           = 38, &
                        tsg_semi_localp          = 39,  tsg_wavelet               = 40, &
                        tsg_level        =   1,  tsg_curved       =  2, &
                        tsg_iptotal      =   3,  tsg_ipcurved     =  4, &
                        tsg_qptotal      =   5,  tsg_qpcurved     =  6, &
                        tsg_hyperbolic   =   7,  tsg_iphyperbolic =  8, &
                        tsg_qphyperbolic =   9,  tsg_tensor       = 10, &
                        tsg_iptensor     =  11,  tsg_qptensor     = 12, &
                        tsg_classic     = 1,  tsg_parents_first = 2, &
                        tsg_directional = 3,  tsg_fds           = 4, &
                        tsg_acc_none     = 0, tsg_acc_cpu_blas  = 1, tsg_acc_gpu_cublas = 2, &
                        tsg_acc_gpu_cuda = 3, tsg_acc_gpu_magma = 4
PRIVATE
  INTEGER :: rows, cols, length
  DOUBLE PRECISION, pointer :: matrix(:,:), vector(:)
  CHARACTER, pointer :: string(:)
CONTAINS
!=======================================================================
SUBROUTINE tsgInitialize()
  CALL tsgbeg()
END SUBROUTINE tsgInitialize
!=======================================================================
SUBROUTINE tsgFinalize()
  CALL tsgend()
END SUBROUTINE tsgFinalize
!=======================================================================
FUNCTION tsgNewGridID() result(newid)
  INTEGER :: newid
  CALL tsgnew(newid)
END FUNCTION tsgNewGridID
!=======================================================================
SUBROUTINE tsgFreeGridID(gridID)
  INTEGER :: gridID
  CALL tsgfre(gridID)
END SUBROUTINE tsgFreeGridID
!=======================================================================
!   class TasmanianSparseGrid wrapper functions
!=======================================================================
FUNCTION tsgGetVersionMajor() result(ver)
  INTEGER :: ver
  CALL tsggvm(ver)
END FUNCTION tsgGetVersionMajor
!=======================================================================
FUNCTION tsgGetVersionMinor() result(ver)
  INTEGER :: ver
  CALL tsggvn(ver)
END FUNCTION tsgGetVersionMinor
!=======================================================================
FUNCTION tsgGetLicense() result(lic)
  CHARACTER, pointer :: lic(:)
  CALL tsggli()
  lic => string
END FUNCTION tsgGetLicense
!=======================================================================
SUBROUTINE tsgMakeGlobalGrid(gridID, dims, outs, depth, gtype, rule, &
                             aweights, alpha, beta, levelLimits)
  INTEGER, intent(in) :: gridID, dims, outs, depth, gtype, rule
  INTEGER :: i
  INTEGER, optional :: aweights(*), levelLimits(dims)
  DOUBLE PRECISION, optional :: alpha, beta
  DOUBLE PRECISION :: al, be
  INTEGER, allocatable :: aw(:)
  INTEGER, allocatable :: ll(:)
  IF(PRESENT(alpha))then
    al = alpha
  ELSE
    al = 0.0
  ENDIF
  IF(PRESENT(beta))then
    be = beta
  ELSE
    be = 0.0
  ENDIF
  ALLOCATE(ll(dims))
  IF(PRESENT(levelLimits))then
    DO i = 1, dims
      ll(i) = levelLimits(i)
    END DO
  ELSE
    DO i = 1, dims
      ll(i) = -1
    END DO
  ENDIF
  IF(PRESENT(aweights))then
    CALL tsgmg(gridID, dims, outs, depth, gtype, rule, aweights, al, be, ll)
  ELSE
    ALLOCATE(aw(2*dims))
    DO i = 1, dims
      aw(i) = 1
    END DO
    DO i = dims+1, 2*dims
      aw(i) = 0
    END DO
    CALL tsgmg(gridID, dims, outs, depth, gtype, rule, aw, al, be, ll)
    DEALLOCATE(aw)
  ENDIF
  DEALLOCATE(ll)
END SUBROUTINE tsgMakeGlobalGrid
!=======================================================================
SUBROUTINE tsgMakeSequenceGrid(gridID, dims, outs, depth, gtype, rule, &
                               aweights)
  INTEGER :: gridID, dims, outs, depth, gtype, rule, i
  INTEGER, optional :: aweights(*)
  INTEGER, allocatable :: aw(:)
  IF(PRESENT(aweights))then
    CALL tsgms(gridID, dims, outs, depth, gtype, rule, aweights)
  ELSE
    ALLOCATE(aw(2*dims))
    DO i = 1, dims
      aw(i) = 1
    END DO
    DO i = dims+1, 2*dims
      aw(i) = 0
    END DO
    CALL tsgms(gridID, dims, outs, depth, gtype, rule, aw)
    DEALLOCATE(aw)
  ENDIF
END SUBROUTINE tsgMakeSequenceGrid
!=======================================================================
SUBROUTINE tsgMakeLocalPolynomialGrid(gridID, dims, outs, depth, order,&
                                      rule)
  INTEGER :: gridID, dims, outs, depth
  INTEGER, optional :: order, rule
  INTEGER :: or, ru
  IF(PRESENT(order))then
    or = order
  ELSE
    or = 1
  ENDIF
  IF(PRESENT(rule))then
    ru = rule
  ELSE
    ru = 1
  ENDIF
  CALL tsgml(gridID, dims, outs, depth, or, ru)
END SUBROUTINE tsgMakeLocalPolynomialGrid
!=======================================================================
SUBROUTINE tsgMakeWaveletGrid(gridID, dims, outs, depth, order)
  INTEGER :: gridID, dims, outs, depth, or
  INTEGER, optional :: order
  IF(PRESENT(order))THEN
    or = order
  ELSE
    or = 1
  ENDIF
  CALL tsgmw(gridID, dims, outs, depth, or)
END SUBROUTINE tsgMakeWaveletGrid
!=======================================================================
SUBROUTINE tsgCopyGrid(gridID, sourceID)
  INTEGER :: gridID, sourceID
  CALL tsgcp(gridID, sourceID)
END SUBROUTINE tsgCopyGrid
!=======================================================================
SUBROUTINE tsgUpdateGlobalGrid(gridID, depth, gtype, aweights)
  INTEGER, intent(in) :: gridID, depth, gtype
  INTEGER, optional :: aweights(*)
  INTEGER, allocatable :: aw(:)
  INTEGER :: dims, i
  IF(PRESENT(aweights))THEN
    CALL tsgug(gridID, depth, gtype, aweights)
  ELSE
    dims = tsgGetNumDimensions(gridID)
    ALLOCATE(aw(2*dims))
    DO i = 1, dims
      aw(i) = 1
    END DO
    DO i = dims+1, 2*dims
      aw(i) = 0
    END DO
    CALL tsgug(gridID, depth, gtype, aw)
    DEALLOCATE(aw)
  ENDIF
END SUBROUTINE tsgUpdateGlobalGrid
!=======================================================================
SUBROUTINE tsgUpdateSequenceGrid(gridID, depth, gtype, aweights)
  INTEGER, intent(in) :: gridID, depth, gtype
  INTEGER :: dims, i
  INTEGER, optional :: aweights(*)
  INTEGER, allocatable :: aw(:)
  IF(PRESENT(aweights))THEN
    CALL tsgus(gridID, depth, gtype, aweights)
  ELSE
    dims = tsgGetNumDimensions(gridID)
    ALLOCATE(aw(2*dims))
    DO i = 1, dims
      aw(i) = 1
    END DO
    DO i = dims+1, 2*dims
      aw(i) = 0
    END DO
    CALL tsgus(gridID, depth, gtype, aw)
    DEALLOCATE(aw)
  ENDIF
END SUBROUTINE tsgUpdateSequenceGrid
!=======================================================================
SUBROUTINE tsgRead(gridID, filename)
  INTEGER, intent(in) :: gridID
  CHARACTER(len=*), intent(in) :: filename
  INTEGER :: N, i
  CHARACTER, allocatable :: cfilename(:)
  N = len(trim(filename))
  ALLOCATE(cfilename(N+1))
  DO i = 1, N
    cfilename(i) = filename(i:i)
  END DO
  cfilename(N+1) = CHAR(0)
  CALL tsgrea(gridID, cfilename)
  DEALLOCATE(cfilename)
END SUBROUTINE tsgRead
!=======================================================================
SUBROUTINE tsgWrite(gridID, filename, useBinary)
  INTEGER, intent(in) :: gridID
  LOGICAL, intent(in), optional :: useBinary
  INTEGER :: ubin
  CHARACTER(len=*), intent(in) :: filename
  INTEGER :: N, i
  CHARACTER, allocatable :: cfilename(:)
  N = len(trim(filename))
  ALLOCATE(cfilename(N+1))
  DO i = 1, N
    cfilename(i) = filename(i:i)
  END DO
  cfilename(N+1) = CHAR(0)
  IF(PRESENT(useBinary))then
    IF(useBinary)then
      ubin = 1
    ELSE
      ubin = 0
    ENDIF
  ELSE
    ubin = 0
  ENDIF
  CALL tsgwri(gridID, cfilename, ubin)
  DEALLOCATE(cfilename)
END SUBROUTINE tsgWrite
!=======================================================================
FUNCTION tsgGetAlpha(gridID) result(alpha)
  INTEGER :: gridID
  DOUBLE PRECISION :: alpha
  CALL tsggal(gridID, alpha)
END FUNCTION tsgGetAlpha
!=======================================================================
FUNCTION tsgGetBeta(gridID) result(beta)
  INTEGER :: gridID
  DOUBLE PRECISION :: beta
  CALL tsggbe(gridID, beta)
END FUNCTION tsgGetBeta
!=======================================================================
FUNCTION tsgGetOrder(gridID) result(order)
  INTEGER :: gridID, order
  CALL tsggor(gridID, order)
END FUNCTION tsgGetOrder
!=======================================================================
FUNCTION tsgGetNumDimensions(gridID) result(dims)
  INTEGER :: gridID, dims
  CALL tsggnd(gridID, dims)
END FUNCTION tsgGetNumDimensions
!=======================================================================
FUNCTION tsgGetNumOutputs(gridID) result(outs)
  INTEGER :: gridID, outs
  CALL tsggno(gridID, outs)
END FUNCTION tsgGetNumOutputs
!=======================================================================
FUNCTION tsgGetRule(gridID) result(rule)
  INTEGER :: gridID, rule
  CALL tsggor(gridID, rule)
END FUNCTION tsgGetRule
!=======================================================================
FUNCTION tsgGetNumLoaded(gridID) result(num)
  INTEGER :: gridID, num
  CALL tsggnl(gridID, num)
END FUNCTION tsgGetNumLoaded
!=======================================================================
FUNCTION tsgGetNumNeeded(gridID) result(num)
  INTEGER :: gridID, num
  CALL tsggnn(gridID, num)
END FUNCTION tsgGetNumNeeded
!=======================================================================
FUNCTION tsgGetNumPoints(gridID) result(num)
  INTEGER :: gridID, num
  CALL tsggnp(gridID, num)
END FUNCTION tsgGetNumPoints
!=======================================================================
FUNCTION tsgGetLoadedPoints(gridID) result(p)
  DOUBLE PRECISION, pointer :: p(:,:)
  INTEGER :: gridID, rows, cols
  rows = tsgGetNumDimensions(gridID)
  cols = tsgGetNumLoaded(gridID)
  ALLOCATE(p(rows,cols))
  CALL tsgglp(gridID, p)
END FUNCTION tsgGetLoadedPoints
!=======================================================================
FUNCTION tsgGetNeededPoints(gridID) result(p)
  DOUBLE PRECISION, pointer :: p(:,:)
  INTEGER :: gridID, rows, cols
  rows = tsgGetNumDimensions(gridID)
  cols = tsgGetNumNeeded(gridID)
  ALLOCATE(p(rows,cols))
  CALL tsggdp(gridID, p)
END FUNCTION tsgGetNeededPoints
!=======================================================================
FUNCTION tsgGetPoints(gridID) result(p)
  DOUBLE PRECISION, pointer :: p(:,:)
  INTEGER :: gridID, rows, cols
  rows = tsgGetNumDimensions(gridID)
  cols = tsgGetNumPoints(gridID)
  ALLOCATE(p(rows,cols))
  CALL tsggpp(gridID, p)
END FUNCTION tsgGetPoints
!=======================================================================
SUBROUTINE tsgGetLoadedPointsStatic(gridID, dims, points)
  INTEGER :: gridID, dims
  DOUBLE PRECISION :: points(dims,*)
  CALL tsgglp(gridID, points)
END SUBROUTINE tsgGetLoadedPointsStatic
!=======================================================================
SUBROUTINE tsgGetNeededPointsStatic(gridID, dims, points)
  INTEGER :: gridID, dims
  DOUBLE PRECISION :: points(dims,*)
  CALL tsggdp(gridID, points)
END SUBROUTINE tsgGetNeededPointsStatic
!=======================================================================
SUBROUTINE tsgGetPointsStatic(gridID, dims, points)
  INTEGER :: gridID, dims
  DOUBLE PRECISION :: points(dims,*)
  CALL tsggpp(gridID, points)
END SUBROUTINE tsgGetPointsStatic
!=======================================================================
FUNCTION tsgGetQuadratureWeights(gridID) result(w)
  INTEGER :: gridID, length
  DOUBLE PRECISION, pointer :: w(:)
  length = tsgGetNumPoints(gridID)
  ALLOCATE(w(length))
  CALL tsggqw(gridID, w)
END FUNCTION tsgGetQuadratureWeights
!=======================================================================
SUBROUTINE tsgGetQuadratureWeightsStatic(gridID, weights)
  INTEGER :: gridID
  DOUBLE PRECISION :: weights(*)
  CALL tsggqw(gridID, weights)
END SUBROUTINE tsgGetQuadratureWeightsStatic
!=======================================================================
FUNCTION tsgGetInterpolationWeights(gridID, x) result(w)
  INTEGER :: gridID, length
  DOUBLE PRECISION :: x(*)
  DOUBLE PRECISION, pointer :: w(:)
  length = tsgGetNumPoints(gridID)
  ALLOCATE(w(length))
  CALL tsggiw(gridID, x, w)
END FUNCTION tsgGetInterpolationWeights
!=======================================================================
SUBROUTINE tsgGetInterpolationWeightsStatic(gridID, x, weights)
  INTEGER :: gridID
  DOUBLE PRECISION :: x(:)
  DOUBLE PRECISION :: weights(:)
  CALL tsggiw(gridID, x, weights)
END SUBROUTINE tsgGetInterpolationWeightsStatic
!=======================================================================
SUBROUTINE tsgLoadNeededPoints(gridID, values)
  INTEGER :: gridID
  DOUBLE PRECISION :: values(:,:)
  CALL tsglnp(gridID, values)
END SUBROUTINE tsgLoadNeededPoints
!=======================================================================
SUBROUTINE tsgEvaluate(gridID, x, y)
  INTEGER :: gridID
  DOUBLE PRECISION, intent(in) :: x(:)
  DOUBLE PRECISION :: y(:)
  CALL tsgeva(gridID, x, y)
END SUBROUTINE tsgEvaluate
!=======================================================================
SUBROUTINE tsgEvaluateFast(gridID, x, y)
  INTEGER :: gridID
  DOUBLE PRECISION, intent(in) :: x(:)
  DOUBLE PRECISION :: y(:)
  CALL tsgevf(gridID, x, y)
END SUBROUTINE tsgEvaluateFast
!=======================================================================
SUBROUTINE tsgEvaluateBatch(gridID, x, numX, y)
  INTEGER :: gridID, numX
  DOUBLE PRECISION :: x(:,:), y(:,:)
  CALL tsgevb(gridID, x, numX, y)
END SUBROUTINE tsgEvaluateBatch
!=======================================================================
SUBROUTINE tsgIntegrate(gridID, q)
  INTEGER :: gridID
  DOUBLE PRECISION :: q(:)
  CALL tsgint(gridID, q)
END SUBROUTINE tsgIntegrate
!=======================================================================
SUBROUTINE tsgSetDomainTransform(gridID, transformA, transformB)
  INTEGER :: gridID
  DOUBLE PRECISION :: transformA(*), transformB(*)
  CALL tsgsdt(gridID, transformA, transformB)
END SUBROUTINE tsgSetDomainTransform
!=======================================================================
FUNCTION tsgIsSetDomainTransfrom(gridID) result(res)
  INTEGER :: gridID, res
  CALL tsgidt(gridID, res)
END FUNCTION tsgIsSetDomainTransfrom
!=======================================================================
SUBROUTINE tsgClearDomainTransform(gridID)
  INTEGER :: gridID
  CALL tsgcdt(gridID)
END SUBROUTINE tsgClearDomainTransform
!=======================================================================
SUBROUTINE tsgGetDomainTransform(gridID, transformA, transformB)
  INTEGER :: gridID
  DOUBLE PRECISION :: transformA(*), transformB(*)
  CALL tsggdt(gridID, transformA, transformB)
END SUBROUTINE tsgGetDomainTransform
!=======================================================================
SUBROUTINE tsgSetAnisotropicRefinement(gridID, gtype, minGrowth, output)
  INTEGER :: gridID, gtype, minGrowth, output
  CALL tsgsar(gridID, gtype, minGrowth, output)
END SUBROUTINE tsgSetAnisotropicRefinement
!=======================================================================
FUNCTION tsgEstimateAnisotropicCoefficients(gridID, gtype, output) &
                                              result(coeff)
  INTEGER :: gridID, gtype, output, N
  INTEGER, allocatable :: coeff(:)
  N = tsgGetNumDimensions(gridID)
  IF ((gtype .EQ. 2) .OR. (gtype .EQ. 4) .OR. (gtype .EQ. 6))then
    N = N * 2
  ENDIF
  ALLOCATE(coeff(N))
  CALL tsgeac(gridID, gtype, output, coeff)
END FUNCTION tsgEstimateAnisotropicCoefficients
!=======================================================================
SUBROUTINE tsgSetGlobalSurplusRefinement(gridID, tolerance, output)
  INTEGER :: gridID, output
  DOUBLE PRECISION :: tolerance
  CALL tsgssr(gridID, tolerance, output)
END SUBROUTINE tsgSetGlobalSurplusRefinement
!=======================================================================
SUBROUTINE tsgSetLocalSurplusRefinement(gridID, tolerance, rtype, output)
  INTEGER :: gridID, rtype, theout
  INTEGER, optional :: output
  DOUBLE PRECISION :: tolerance
  IF(PRESENT(output))then
    theout = output
  ELSE
    theout = -1
  ENDIF
  CALL tsgshr(gridID, tolerance, rtype, theout)
END SUBROUTINE tsgSetLocalSurplusRefinement
!=======================================================================
SUBROUTINE tsgClearRefinement(gridID)
  INTEGER :: gridID
  CALL tsgcre(gridID)
END SUBROUTINE tsgClearRefinement
!=======================================================================
SUBROUTINE tsgSetConformalTransformASIN(gridID, truncate)
  INTEGER :: gridID, truncate(:)
  CALL tsgsca(gridID, truncate)
END SUBROUTINE tsgSetConformalTransformASIN
!=======================================================================
FUNCTION tsgIsSetConformalTransformASIN(gridID) result(isset)
  INTEGER :: gridID, res
  LOGICAL :: isset
  CALL tsgica(gridID, res)
  IF (res .EQ. 0)then
    isset = .FALSE.
  ELSE
    isset = .TRUE.
  ENDIF
END FUNCTION tsgIsSetConformalTransformASIN
!=======================================================================
SUBROUTINE tsgClearConformalTransform(gridID)
  INTEGER :: gridID
  CALL tsgcct(gridID)
END SUBROUTINE tsgClearConformalTransform
!=======================================================================
FUNCTION tsgGetConformalTransformASIN(gridID) result(truncate)
  INTEGER :: gridID, N
  INTEGER, allocatable :: truncate(:)
  N = tsgGetNumDimensions(gridID)
  ALLOCATE(truncate(N))
  CALL tsgcre(gridID, truncate)
END FUNCTION tsgGetConformalTransformASIN
!=======================================================================
SUBROUTINE tsgPrintStats(gridID)
  INTEGER :: gridID
  CALL tsgpri(gridID)
END SUBROUTINE tsgPrintStats
!=======================================================================
! DO NOT CALL THOSE FUNCTIONS DIRECTLY !
!=======================================================================
SUBROUTINE tsgReceiveString(l, S)
  INTEGER :: l
  CHARACTER, target :: S(l)
  length = l
  string => S(1:length)
END SUBROUTINE tsgReceiveString
!=======================================================================
SUBROUTINE tsgReceiveVector(s, V)
  INTEGER :: s
  DOUBLE PRECISION, target :: V(s)
  length = s
  vector => V(1:length)
END SUBROUTINE tsgReceiveVector
!=======================================================================
SUBROUTINE tsgReceiveMatrix(r, c, M)
  INTEGER :: r, c
  DOUBLE PRECISION, target :: M(r,c)
  rows = r
  cols = c
  matrix => M(1:rows,1:cols)
END SUBROUTINE tsgReceiveMatrix
!=======================================================================
END MODULE
