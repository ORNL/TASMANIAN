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
implicit none
PUBLIC :: tsgClearAll,          &
          tsgNewGridID,         &
          tsgFreeGridID,        &
!===== ABOVE ARE FORTRAN SPECIFIC FUNCTIONS ======!
!===== BELOW ARE WRAPPERS FOR STANDARD C++ CLASS =!
          tsgIsGlobal,          &
          tsgIsSequence,        &
          tsgIsLocalPolynomial, &
          tsgIsWavelet,         &
          tsgIsFourier,         &
          tsgGetVersionMajor,   &
          tsgGetVersionMinor,   &
          tsgGetLicense,        &
          tsgMakeGlobalGrid,    &
          tsgMakeSequenceGrid,  &
          tsgMakeLocalPolynomialGrid, &
          tsgMakeWaveletGrid,   &
          tsgMakeFourierGrid,   &
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
          tsgEvaluateHierarchicalFunctions,        &
          tsgEvaluateComplexHierarchicalFunctions, &
          tsgEvaluateSparseHierarchicalFunctions,  &
          tsgGetHierarchicalCoefficients,          &
          tsgGetHierarchicalCoefficientsStatic,    &
          tsgGetComplexHierarchicalCoefficients,   &
          tsgGetComplexHierarchicalCoefficientsStatic, &
          tsgSetHierarchicalCoefficients, &
          tsgIntegrate,        &
          tsgSetDomainTransform,   &
          tsgIsSetDomainTransform, &
          tsgClearDomainTransform, &
          tsgGetDomainTransform,   &
          tsgSetAnisotropicRefinement,        &
          tsgEstimateAnisotropicCoefficients, &
          tsgSetGlobalSurplusRefinement,      &
          tsgSetLocalSurplusRefinement,       &
          tsgClearRefinement,                 &
          tsgMergeRefinement,                 &
          tsgSetConformalTransformASIN,   &
          tsgIsSetConformalTransformASIN, &
          tsgClearConformalTransform,     &
          tsgGetConformalTransformASIN,   &
          tsgPrintStats, &
          tsgTestInternals,     &
!===== do NOT USE THESE FUNCTIONS DIRECTLY ======!
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

  integer, parameter :: tsg_clenshaw_curtis      =  1,  tsg_rule_clenshaw_curtis_zero = 2, &
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
  integer :: rows, cols, length
  double precision, pointer :: matrix(:,:), vector(:)
  character, pointer :: string(:)
contains
!=======================================================================
function tsgIsGlobal(id) result(res)
  integer :: id, flag
  logical :: res
  call tsgisg(id, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
end function
!=======================================================================
function tsgIsSequence(id) result(res)
  integer :: id, flag
  logical :: res
  call tsgiss(id, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
end function
!=======================================================================
function tsgIsLocalPolynomial(id) result(res)
  integer :: id, flag
  logical :: res
  call tsgisl(id, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
end function
!=======================================================================
function tsgIsWavelet(id) result(res)
  integer :: id, flag
  logical :: res
  call tsgisw(id, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
end function
!=======================================================================
function tsgIsFourier(id) result(res)
  integer :: id, flag
  logical :: res
  call tsgisf(id, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
end function
!=======================================================================
subroutine tsgClearAll()
  call tsgend()
end subroutine tsgClearAll
!=======================================================================
function tsgNewGridID() result(newid)
  integer :: newid
  call tsgnew(newid)
end function tsgNewGridID
!=======================================================================
subroutine tsgFreeGridID(gridID)
  integer :: gridID
  call tsgfre(gridID)
end subroutine tsgFreeGridID
!=======================================================================
!   class TasmanianSparseGrid wrapper functions
!=======================================================================
function tsgGetVersionMajor() result(ver)
  integer :: ver
  call tsggvm(ver)
end function tsgGetVersionMajor
!=======================================================================
function tsgGetVersionMinor() result(ver)
  integer :: ver
  call tsggvn(ver)
end function tsgGetVersionMinor
!=======================================================================
function tsgGetLicense() result(lic)
  character, pointer :: lic(:)
  call tsggli()
  lic => string
end function tsgGetLicense
!=======================================================================
subroutine tsgMakeGlobalGrid(gridID, dims, outs, depth, gtype, rule, &
                             aweights, alpha, beta, customRuleFilename, levelLimits)
  integer, intent(in) :: gridID, dims, outs, depth, gtype, rule
  integer, optional, target  :: aweights(:), levelLimits(dims)
  double precision, optional :: alpha, beta
  character(len=*), optional :: customRuleFilename
  integer           :: opt_flags(5)
  character(len=80) :: cfn = char(0)
  double precision  :: al, be
  integer, pointer  :: aw(:) => null()
  integer, pointer  :: ll(:) => null()

  opt_flags = 0
  if ( present(aweights) ) then
    opt_flags(1) = 1
    aw => aweights
  endif
  if ( present(alpha) ) then
    opt_flags(2) = 1
    al = alpha
  endif
  if ( present(beta) ) then
    opt_flags(3) = 1
    be = beta
    endif
  if ( present(customRuleFilename) ) then
    opt_flags(4) = 1
    cfn = customRuleFilename//char(0)
  endif
  if ( present(levelLimits) ) then
    opt_flags(5) = 1
    ll => levelLimits
  endif

  call tsgmg(gridID, dims, outs, depth, gtype, rule, opt_flags, aw, al, be, cfn, ll)
end subroutine tsgMakeGlobalGrid
!=======================================================================
subroutine tsgMakeSequenceGrid(gridID, dims, outs, depth, gtype, rule, aweights, levelLimits)
integer :: gridID, dims, outs, depth, gtype, rule
integer, optional, target :: aweights(:), levelLimits(dims)
integer          :: opt_flags(2)
integer, pointer :: aw(:) => null()
integer, pointer :: ll(:) => null()

opt_flags = 0
if ( present(aweights) ) then
  opt_flags(1) = 1
  aw => aweights
endif
if ( present(levelLimits) ) then
  opt_flags(2) = 1
  ll => levelLimits
endif

call tsgms(gridID, dims, outs, depth, gtype, rule, opt_flags, aw, ll)
end subroutine tsgMakeSequenceGrid
!=======================================================================
subroutine tsgMakeLocalPolynomialGrid(gridID, dims, outs, depth, order, rule, levelLimits)
  integer :: gridID, dims, outs, depth
  integer, optional :: order, rule
  integer, optional, target :: levelLimits(dims)
  integer          :: opt_flags(3), or, ru
  integer, pointer :: ll(:) => null()

  opt_flags = 0
  if ( present(order) ) then
    opt_flags(1) = 1
    or = order
  endif
  if ( present(rule) ) then
    opt_flags(2) = 1
    ru = rule
  endif
  if ( present(levelLimits) ) then
    opt_flags(3) = 1
    ll => levelLimits
  endif

  call tsgml(gridID, dims, outs, depth, opt_flags, or, ru, ll)
end subroutine tsgMakeLocalPolynomialGrid
!=======================================================================
subroutine tsgMakeWaveletGrid(gridID, dims, outs, depth, order, levelLimits)
  integer :: gridID, dims, outs, depth
  integer, optional, target :: levelLimits(dims)
  integer, optional :: order
  integer           :: opt_flags(2), or
  integer, pointer  :: ll(:) => null()

  opt_flags = 0
  if ( present(order) ) then
    opt_flags(1) = 1
    or = order
  endif
  if ( present(levelLimits) ) then
    opt_flags(2) = 1
    ll => levelLimits
  endif

  call tsgmw(gridID, dims, outs, depth, opt_flags, or, ll)
end subroutine tsgMakeWaveletGrid
!=======================================================================
subroutine tsgMakeFourierGrid(gridID, dims, outs, depth, gtype, aweights, levelLimits)
  integer :: gridID, dims, outs, depth, gtype
  integer, optional, target :: aweights(:), levelLimits(dims)
  integer          :: opt_flags(2)
  integer, pointer :: aw(:) => null()
  integer, pointer :: ll(:) => null()

  opt_flags = 0
  if ( present(aweights) ) then
    opt_flags(1) = 1
    aw => aweights
  endif
  if ( present(levelLimits) ) then
    opt_flags(2) = 1
    ll => levelLimits
  endif

  call tsgmf(gridID, dims, outs, depth, gtype, opt_flags, aw, ll)
end subroutine tsgMakeFourierGrid
!=======================================================================
subroutine tsgCopyGrid(gridID, sourceID)
  integer :: gridID, sourceID
  call tsgcp(gridID, sourceID)
end subroutine tsgCopyGrid
!=======================================================================
subroutine tsgUpdateGlobalGrid(gridID, depth, gtype, aweights)
  integer, intent(in) :: gridID, depth, gtype
  integer, optional, target :: aweights(:)
  integer          :: opt_flags = 0
  integer, pointer :: aw(:) => null()
  if ( present(aweights) ) then
    opt_flags = 1
    aw => aweights
  endif
  call tsgug(gridID, depth, gtype, opt_flags, aw)
end subroutine tsgUpdateGlobalGrid
!=======================================================================
subroutine tsgUpdateSequenceGrid(gridID, depth, gtype, aweights)
  integer, intent(in) :: gridID, depth, gtype
  integer, optional, target :: aweights(:)
  integer          :: opt_flags = 0
  integer, pointer :: aw(:) => null()
  if ( present(aweights) ) then
    opt_flags = 1
    aw => aweights
  endif
  call tsgus(gridID, depth, gtype, opt_flags, aw)
end subroutine tsgUpdateSequenceGrid
!=======================================================================
function tsgRead(gridID, filename) result(res)
  integer, intent(in) :: gridID
  character(len=*), intent(in) :: filename
  logical :: res
  integer :: N, i, flag
  character, allocatable :: cfilename(:)
  N = len(trim(filename))
  allocate(cfilename(N+1))
  do i = 1, N
    cfilename(i) = filename(i:i)
  end do
  cfilename(N+1) = CHAR(0)
  call tsgrea(gridID, cfilename, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
  deallocate(cfilename)
end function tsgRead
!=======================================================================
subroutine tsgWrite(gridID, filename, useBinary)
  integer, intent(in) :: gridID
  logical, intent(in), optional :: useBinary
  integer :: ubin
  character(len=*), intent(in) :: filename
  integer :: N, i
  character, allocatable :: cfilename(:)
  N = len(trim(filename))
  allocate(cfilename(N+1))
  do i = 1, N
    cfilename(i) = filename(i:i)
  end do
  cfilename(N+1) = CHAR(0)
  if(present(useBinary))then
    if(useBinary)then
      ubin = 1
    else
      ubin = 0
    endif
  else
    ubin = 0
  endif
  call tsgwri(gridID, cfilename, ubin)
  deallocate(cfilename)
end subroutine tsgWrite
!=======================================================================
function tsgGetAlpha(gridID) result(alpha)
  integer :: gridID
  double precision :: alpha
  call tsggal(gridID, alpha)
end function tsgGetAlpha
!=======================================================================
function tsgGetBeta(gridID) result(beta)
  integer :: gridID
  double precision :: beta
  call tsggbe(gridID, beta)
end function tsgGetBeta
!=======================================================================
function tsgGetOrder(gridID) result(order)
  integer :: gridID, order
  call tsggor(gridID, order)
end function tsgGetOrder
!=======================================================================
function tsgGetNumDimensions(gridID) result(dims)
  integer :: gridID, dims
  call tsggnd(gridID, dims)
end function tsgGetNumDimensions
!=======================================================================
function tsgGetNumOutputs(gridID) result(outs)
  integer :: gridID, outs
  call tsggno(gridID, outs)
end function tsgGetNumOutputs
!=======================================================================
function tsgGetRule(gridID) result(rule)
  integer :: gridID, rule
  call tsggor(gridID, rule)
end function tsgGetRule
!=======================================================================
function tsgGetNumLoaded(gridID) result(num)
  integer :: gridID, num
  call tsggnl(gridID, num)
end function tsgGetNumLoaded
!=======================================================================
function tsgGetNumNeeded(gridID) result(num)
  integer :: gridID, num
  call tsggnn(gridID, num)
end function tsgGetNumNeeded
!=======================================================================
function tsgGetNumPoints(gridID) result(num)
  integer :: gridID, num
  call tsggnp(gridID, num)
end function tsgGetNumPoints
!=======================================================================
function tsgGetLoadedPoints(gridID) result(p)
  double precision, pointer :: p(:,:)
  integer :: gridID, rows, cols
  rows = tsgGetNumDimensions(gridID)
  cols = tsgGetNumLoaded(gridID)
  allocate(p(rows,cols))
  call tsgglp(gridID, p)
end function tsgGetLoadedPoints
!=======================================================================
function tsgGetNeededPoints(gridID) result(p)
  double precision, pointer :: p(:,:)
  integer :: gridID, rows, cols
  rows = tsgGetNumDimensions(gridID)
  cols = tsgGetNumNeeded(gridID)
  allocate(p(rows,cols))
  call tsggdp(gridID, p)
end function tsgGetNeededPoints
!=======================================================================
function tsgGetPoints(gridID) result(p)
  double precision, pointer :: p(:,:)
  integer :: gridID, rows, cols
  rows = tsgGetNumDimensions(gridID)
  cols = tsgGetNumPoints(gridID)
  allocate(p(rows,cols))
  call tsggpp(gridID, p)
end function tsgGetPoints
!=======================================================================
subroutine tsgGetLoadedPointsStatic(gridID, points)
  integer :: gridID
  double precision :: points(:,:)
  call tsgglp(gridID, points)
end subroutine tsgGetLoadedPointsStatic
!=======================================================================
subroutine tsgGetNeededPointsStatic(gridID, points)
  integer :: gridID
  double precision :: points(:,:)
  call tsggdp(gridID, points)
end subroutine tsgGetNeededPointsStatic
!=======================================================================
subroutine tsgGetPointsStatic(gridID, points)
  integer :: gridID
  double precision :: points(:,:)
  call tsggpp(gridID, points)
end subroutine tsgGetPointsStatic
!=======================================================================
function tsgGetQuadratureWeights(gridID) result(w)
  integer :: gridID, length
  double precision, pointer :: w(:)
  length = tsgGetNumPoints(gridID)
  allocate(w(length))
  call tsggqw(gridID, w)
end function tsgGetQuadratureWeights
!=======================================================================
subroutine tsgGetQuadratureWeightsStatic(gridID, weights)
  integer :: gridID
  double precision :: weights(*)
  call tsggqw(gridID, weights)
end subroutine tsgGetQuadratureWeightsStatic
!=======================================================================
function tsgGetInterpolationWeights(gridID, x) result(w)
  integer :: gridID, length
  double precision :: x(*)
  double precision, pointer :: w(:)
  length = tsgGetNumPoints(gridID)
  allocate(w(length))
  call tsggiw(gridID, x, w)
end function tsgGetInterpolationWeights
!=======================================================================
subroutine tsgGetInterpolationWeightsStatic(gridID, x, weights)
  integer :: gridID
  double precision :: x(:)
  double precision :: weights(:)
  call tsggiw(gridID, x, weights)
end subroutine tsgGetInterpolationWeightsStatic
!=======================================================================
subroutine tsgLoadNeededPoints(gridID, values)
  integer :: gridID
  double precision :: values(:,:)
  call tsglnp(gridID, values)
end subroutine tsgLoadNeededPoints
!=======================================================================
subroutine tsgEvaluate(gridID, x, y)
  integer :: gridID
  double precision, intent(in) :: x(:)
  double precision :: y(:)
  call tsgeva(gridID, x, y)
end subroutine tsgEvaluate
!=======================================================================
subroutine tsgEvaluateFast(gridID, x, y)
  integer :: gridID
  double precision, intent(in) :: x(:)
  double precision :: y(:)
  call tsgevf(gridID, x, y)
end subroutine tsgEvaluateFast
!=======================================================================
subroutine tsgEvaluateBatch(gridID, x, numX, y)
  integer :: gridID, numX
  double precision :: x(:,:), y(:,:)
  call tsgevb(gridID, x, numX, y)
end subroutine tsgEvaluateBatch
!=======================================================================
subroutine tsgEvaluateHierarchicalFunctions(gridID, x, numX, y)
  integer :: gridID, numX
  double precision :: x(:,:), y(:,:)
  if ( .not. tsgIsFourier(gridID) ) then
    call tsgehf(gridID, x, numX, y)
  else
    write(*,*) "ERROR: called tsgEvaluateHierarchicalFunctions() on a Fourier grid, "
    write(*,*) "       use tsgEvaluateComplexHierarchicalFunctions() instead"
  endif
end subroutine tsgEvaluateHierarchicalFunctions
!=======================================================================
subroutine tsgEvaluateComplexHierarchicalFunctions(gridID, x, numX, y)
  integer :: gridID, numX
  double precision :: x(:,:)
  double complex   :: y(:,:)
  double precision, allocatable :: y_c_style(:)
  integer :: i, j
  if ( tsgIsFourier(gridID) ) then
    allocate(y_c_style(2*size(y,1)*size(y,2)))
    call tsgehf(gridID, x, numX, y_c_style)
    do j = 1,size(y,2)
      do i = 1,size(y,1)
        y(i,j) = cmplx( y_c_style( 2*size(y,1)*(j-1) + 2*i-1 ), y_c_style( 2*size(y,1)*(j-1) + 2*i), kind(y) )
      enddo
    enddo
  else
    write(*,*) "ERROR: called tsgEvaluateComplexHierarchicalFunctions() on a non-Fourier grid, "
    write(*,*) "       use tsgEvaluateHierarchicalFunctions() instead"
  endif
end subroutine tsgEvaluateComplexHierarchicalFunctions
!=======================================================================
subroutine tsgEvaluateSparseHierarchicalFunctions(gridID, x, numX, pntr, indx, y)
  integer :: gridID, numX
  double precision :: x(:,:)
  integer, pointer :: pntr(:), indx(:)
  double precision, pointer :: y(:)
  integer :: numNZ
  if ( .not. tsgIsFourier(gridID) ) then
    call tsgehz(gridID, x, numX, numNZ)
    allocate( pntr(numX+1), indx(numNZ), y(numNZ) )
    call tsgehs(gridID, x, numX, pntr, indx, y)
  else
    write(*,*) "ERROR: called tsgEvaluateSparseHierarchicalFunctions() on a Fourier grid"
  endif
end subroutine tsgEvaluateSparseHierarchicalFunctions
!=======================================================================
function tsgGetHierarchicalCoefficients(gridID) result(c)
  integer :: gridID
  double precision, pointer :: c(:)
  if ( .not. tsgIsFourier(gridID) ) then
    allocate(c(tsgGetNumOutputs(gridID)*tsgGetNumPoints(gridID)))
    call tsgghc(gridID, c)
  else
    write(*,*) "ERROR: called tsgGetHierarchicalCoefficients() on a Fourier grid, "
    write(*,*) "       use tsgGetComplexHierarchicalCoefficients() instead"
  endif
end function tsgGetHierarchicalCoefficients
!=======================================================================
subroutine tsgGetHierarchicalCoefficientsStatic(gridID, c)
  integer :: gridID
  double precision :: c(:)
  if ( .not. tsgIsFourier(gridID) ) then
    call tsgghc(gridID, c)
  else
    write(*,*) "ERROR: called tsgGetHierarchicalCoefficientsStatic() on a Fourier grid, "
    write(*,*) "       use tsgGetComplexHierarchicalCoefficientsStatic() instead"
  endif
end subroutine tsgGetHierarchicalCoefficientsStatic
!=======================================================================
function tsgGetComplexHierarchicalCoefficients(gridID) result(c)
  integer :: gridID
  double complex, pointer       :: c(:)
  double precision, allocatable :: c_real(:)
  integer :: i
  if ( tsgIsFourier(gridID) ) then
    allocate(c(tsgGetNumOutputs(gridID)*tsgGetNumPoints(gridID)))
    allocate(c_real(2*tsgGetNumOutputs(gridID)*tsgGetNumPoints(gridID)))
    call tsgghc(gridID, c_real)
    do i = 1,size(c)
      c(i) = cmplx( c_real(i), c_real(i + size(c)), kind(c) )
    enddo
    deallocate(c_real)
  else
    write(*,*) "ERROR: called tsgGetComplexHierarchicalCoefficients() on a non-Fourier grid, "
    write(*,*) "       use tsgGetHierarchicalCoefficients() instead"
  endif
end function tsgGetComplexHierarchicalCoefficients
!=======================================================================
subroutine tsgGetComplexHierarchicalCoefficientsStatic(gridID, c)
  integer :: gridID
  double complex   :: c(:)
  double precision :: c_real(2*size(c))
  integer :: i
  if ( tsgIsFourier(gridID) ) then
    call tsgghc(gridID, c_real)
    do i = 1,size(c)
      c(i) = cmplx( c_real(i), c_real(i + size(c)), kind(c) )
    enddo
  else
    write(*,*) "ERROR: called tsgGetComplexHierarchicalCoefficientsStatic() on a non-Fourier grid, "
    write(*,*) "       use tsgGetHierarchicalCoefficientsStatic() instead"
  endif
end subroutine tsgGetComplexHierarchicalCoefficientsStatic
!=======================================================================
subroutine tsgSetHierarchicalCoefficients(gridID,c)
  integer :: gridID
  double precision :: c(:)
  call tsgshc(gridID, c)
end subroutine tsgSetHierarchicalCoefficients
!=======================================================================
subroutine tsgIntegrate(gridID, q)
  integer :: gridID
  double precision :: q(:)
  call tsgint(gridID, q)
end subroutine tsgIntegrate
!=======================================================================
subroutine tsgSetDomainTransform(gridID, transformA, transformB)
  integer :: gridID
  double precision :: transformA(*), transformB(*)
  call tsgsdt(gridID, transformA, transformB)
end subroutine tsgSetDomainTransform
!=======================================================================
function tsgIsSetDomainTransform(gridID) result(res)
  integer :: gridID
  logical :: res
  call tsgidt(gridID, res)
end function tsgIsSetDomainTransform
!=======================================================================
subroutine tsgClearDomainTransform(gridID)
  integer :: gridID
  call tsgcdt(gridID)
end subroutine tsgClearDomainTransform
!=======================================================================
subroutine tsgGetDomainTransform(gridID, transformA, transformB)
  integer :: gridID
  double precision :: transformA(*), transformB(*)
  call tsggdt(gridID, transformA, transformB)
end subroutine tsgGetDomainTransform
!=======================================================================
subroutine tsgSetAnisotropicRefinement(gridID, gtype, minGrowth, output, levelLimits)
  integer :: gridID, gtype, minGrowth, output
  integer, optional, target :: levelLimits(:)
  integer          :: opt_flags = 0
  integer, pointer :: ll(:) => null()
  if (present(levelLimits)) then
    opt_flags = 1
    ll => levelLimits
  endif
  call tsgsar(gridID, gtype, minGrowth, output-1, opt_flags, ll)
end subroutine tsgSetAnisotropicRefinement
!=======================================================================
function tsgEstimateAnisotropicCoefficients(gridID, gtype, output) result(coeff)
  integer :: gridID, gtype, output, N
  integer, pointer :: coeff(:)
  N = tsgGetNumDimensions(gridID)
  if ((gtype .EQ. tsg_curved) .OR. (gtype .EQ. tsg_ipcurved) .OR. (gtype .EQ. tsg_qpcurved))then
    N = N * 2
  endif
  allocate(coeff(N))
  call tsgeac(gridID, gtype, output-1, coeff)
end function tsgEstimateAnisotropicCoefficients
!=======================================================================
subroutine tsgSetGlobalSurplusRefinement(gridID, tolerance, output, levelLimits)
  integer :: gridID, output
  integer, optional, target :: levelLimits(:)
  double precision :: tolerance
  integer          :: opt_flags = 0
  integer, pointer :: ll(:) => null()
  if (present(levelLimits)) then
    opt_flags = 1
    ll => levelLimits
  endif
  call tsgssr(gridID, tolerance, output-1, opt_flags, ll)
end subroutine tsgSetGlobalSurplusRefinement
!=======================================================================
subroutine tsgSetLocalSurplusRefinement(gridID, tolerance, rtype, output, levelLimits)
  integer :: gridID, rtype
  integer, optional :: output
  integer, optional, target :: levelLimits(:)
  double precision :: tolerance
  integer          :: opt_flags(2), theout
  integer, pointer :: ll(:) => null()

  opt_flags = 0
  if (present(output)) then
    opt_flags(1) = 1
    theout = output-1
  endif
  if (present(levelLimits)) then
    opt_flags(2) = 1
    ll => levelLimits
  endif

  call tsgshr(gridID, tolerance, rtype, opt_flags, theout, ll)
end subroutine tsgSetLocalSurplusRefinement
!=======================================================================
subroutine tsgClearRefinement(gridID)
  integer :: gridID
  call tsgcre(gridID)
end subroutine tsgClearRefinement
!=======================================================================
subroutine tsgMergeRefinement(gridID)
  integer :: gridID
  call tsgmre(gridID)
end subroutine tsgMergeRefinement
!=======================================================================
subroutine tsgSetConformalTransformASIN(gridID, truncate)
  integer :: gridID, truncate(:)
  call tsgsca(gridID, truncate)
end subroutine tsgSetConformalTransformASIN
!=======================================================================
function tsgIsSetConformalTransformASIN(gridID) result(isset)
  integer :: gridID, res
  logical :: isset
  call tsgica(gridID, res)
  if (res .EQ. 0)then
    isset = .FALSE.
  else
    isset = .TRUE.
  endif
end function tsgIsSetConformalTransformASIN
!=======================================================================
subroutine tsgClearConformalTransform(gridID)
  integer :: gridID
  call tsgcct(gridID)
end subroutine tsgClearConformalTransform
!=======================================================================
function tsgGetConformalTransformASIN(gridID) result(truncate)
  integer :: gridID, N
  integer, allocatable :: truncate(:)
  N = tsgGetNumDimensions(gridID)
  allocate(truncate(N))
  call tsgcre(gridID, truncate)
end function tsgGetConformalTransformASIN
!=======================================================================
subroutine tsgPrintStats(gridID)
  integer :: gridID
  call tsgpri(gridID)
end subroutine tsgPrintStats
!=======================================================================
function tsgTestInternals() result(res)
  integer, parameter :: i_a = 8, i_b = 4
  integer :: i, num_ag, num_ag_in
  integer :: int_1d_a(i_a), int_1d_b(i_a) = (/1,8,6,4,7,3,2,5/)
  logical :: res
  logical :: verb
  ! NOTE: calling tsgTestInternals() will clear all grids and invalidate all grid IDs
  write(*,*) "Perfroming internal consistency test in Tasmanian Fortran 90/95 module, all existing grids will be cleared."

  res = .true.
  call tsggag(num_ag_in)
  do i = 1, i_a
    int_1d_a(i) = tsgNewGridID()
    call tsggag(num_ag)
    if ( num_ag .ne. i+num_ag_in ) then
      write(*,*) "Mismatch in number of active grids, adding grids: num_ag = ", num_ag
      res = .false.
    endif
  enddo
  do i = 1, i_a
    call tsgFreeGridID(int_1d_a(int_1d_b(i)))
    call tsggag(num_ag)
    if ( num_ag .ne. i_a-i+num_ag_in ) then
      write(*,*) "Mismatch in number of active grids, removing grids: num_ag = ", num_ag
      res = .false.
    endif
  enddo
  int_1d_a(1) = tsgNewGridID()
  int_1d_a(2) = tsgNewGridID()
  call tsgClearAll()
  call tsggag(num_ag)
  if ( num_ag .ne. 0 ) then
    write(*,*) "Mismatch in number of active grids after clearall: num_ag = ", num_ag
    res = .false.
  endif
end function tsgTestInternals
!=======================================================================
! do NOT call THOSE FUNCTIONS DIRECTLY !
!=======================================================================
subroutine tsgReceiveString(l, S)
  integer :: l
  character, target :: S(l)
  length = l
  string => S(1:length)
end subroutine tsgReceiveString
!=======================================================================
subroutine tsgReceiveVector(s, V)
  integer :: s
  double precision, target :: V(s)
  length = s
  vector => V(1:length)
end subroutine tsgReceiveVector
!=======================================================================
subroutine tsgReceiveMatrix(r, c, M)
  integer :: r, c
  double precision, target :: M(r,c)
  rows = r
  cols = c
  matrix => M(1:rows,1:cols)
end subroutine tsgReceiveMatrix
!=======================================================================
end MODULE
