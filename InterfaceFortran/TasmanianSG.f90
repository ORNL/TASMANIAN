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
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN if ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
! THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
! COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
! THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
! IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
!==================================================================================================================================================================================
Module TasmanianSG
implicit none
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
subroutine tsgInitialize()
  call tsgbeg()
end subroutine tsgInitialize
!=======================================================================
subroutine tsgFinalize()
  call tsgend()
end subroutine tsgFinalize
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
                             aweights, alpha, beta, levelLimits)
  integer, intent(in) :: gridID, dims, outs, depth, gtype, rule
  integer :: i
  integer, optional :: aweights(*), levelLimits(dims)
  double precision, optional :: alpha, beta
  double precision :: al, be
  integer, allocatable :: aw(:)
  integer, allocatable :: ll(:)
  if(present(alpha))then
    al = alpha
  else
    al = 0.0
  endif
  if(present(beta))then
    be = beta
  else
    be = 0.0
  endif
  allocate(ll(dims))
  if(present(levelLimits))then
    do i = 1, dims
      ll(i) = levelLimits(i)
    end do
  else
    do i = 1, dims
      ll(i) = -1
    end do
  endif
  if(present(aweights))then
    call tsgmg(gridID, dims, outs, depth, gtype, rule, aweights, al, be, ll)
  else
    allocate(aw(2*dims))
    do i = 1, dims
      aw(i) = 1
    end do
    do i = dims+1, 2*dims
      aw(i) = 0
    end do
    call tsgmg(gridID, dims, outs, depth, gtype, rule, aw, al, be, ll)
    deallocate(aw)
  endif
  deallocate(ll)
end subroutine tsgMakeGlobalGrid
!=======================================================================
subroutine tsgMakeSequenceGrid(gridID, dims, outs, depth, gtype, rule, &
                               aweights)
  integer :: gridID, dims, outs, depth, gtype, rule, i
  integer, optional :: aweights(*)
  integer, allocatable :: aw(:)
  if(present(aweights))then
    call tsgms(gridID, dims, outs, depth, gtype, rule, aweights)
  else
    allocate(aw(2*dims))
    do i = 1, dims
      aw(i) = 1
    end do
    do i = dims+1, 2*dims
      aw(i) = 0
    end do
    call tsgms(gridID, dims, outs, depth, gtype, rule, aw)
    deallocate(aw)
  endif
end subroutine tsgMakeSequenceGrid
!=======================================================================
subroutine tsgMakeLocalPolynomialGrid(gridID, dims, outs, depth, order,&
                                      rule)
  integer :: gridID, dims, outs, depth
  integer, optional :: order, rule
  integer :: or, ru
  if(present(order))then
    or = order
  else
    or = 1
  endif
  if(present(rule))then
    ru = rule
  else
    ru = 1
  endif
  call tsgml(gridID, dims, outs, depth, or, ru)
end subroutine tsgMakeLocalPolynomialGrid
!=======================================================================
subroutine tsgMakeWaveletGrid(gridID, dims, outs, depth, order)
  integer :: gridID, dims, outs, depth, or
  integer, optional :: order
  if(present(order))THEN
    or = order
  else
    or = 1
  endif
  call tsgmw(gridID, dims, outs, depth, or)
end subroutine tsgMakeWaveletGrid
!=======================================================================
subroutine tsgMakeFourierGrid(gridID, dims, outs, depth, gtype, aweights, levelLimits)
  integer :: gridID, dims, outs, depth, gtype
  integer, optional :: aweights(*), levelLimits(dims)
  integer, allocatable :: aw(:)
  integer, allocatable :: ll(:)
  allocate(ll(dims))
  if(present(levelLimits)) then
    ll = levelLimits
  else
    ll = -1
  endif
  if(present(aweights)) then
    call tsgmf(gridID, dims, outs, depth, gtype, aweights, ll)
  else
    allocate(aw(2*dims))
    aw(1:dims)  = 1
    aw(dims+1:) = 0
    call tsgmf(gridID, dims, outs, depth, gtype, aw, ll)
    deallocate(aw)
  endif
  deallocate(ll)
end subroutine tsgMakeFourierGrid
!=======================================================================
subroutine tsgCopyGrid(gridID, sourceID)
  integer :: gridID, sourceID
  call tsgcp(gridID, sourceID)
end subroutine tsgCopyGrid
!=======================================================================
subroutine tsgUpdateGlobalGrid(gridID, depth, gtype, aweights)
  integer, intent(in) :: gridID, depth, gtype
  integer, optional :: aweights(*)
  integer, allocatable :: aw(:)
  integer :: dims, i
  if(present(aweights))THEN
    call tsgug(gridID, depth, gtype, aweights)
  else
    dims = tsgGetNumDimensions(gridID)
    allocate(aw(2*dims))
    do i = 1, dims
      aw(i) = 1
    end do
    do i = dims+1, 2*dims
      aw(i) = 0
    end do
    call tsgug(gridID, depth, gtype, aw)
    deallocate(aw)
  endif
end subroutine tsgUpdateGlobalGrid
!=======================================================================
subroutine tsgUpdateSequenceGrid(gridID, depth, gtype, aweights)
  integer, intent(in) :: gridID, depth, gtype
  integer :: dims, i
  integer, optional :: aweights(*)
  integer, allocatable :: aw(:)
  if(present(aweights))THEN
    call tsgus(gridID, depth, gtype, aweights)
  else
    dims = tsgGetNumDimensions(gridID)
    allocate(aw(2*dims))
    do i = 1, dims
      aw(i) = 1
    end do
    do i = dims+1, 2*dims
      aw(i) = 0
    end do
    call tsgus(gridID, depth, gtype, aw)
    deallocate(aw)
  endif
end subroutine tsgUpdateSequenceGrid
!=======================================================================
subroutine tsgRead(gridID, filename)
  integer, intent(in) :: gridID
  character(len=*), intent(in) :: filename
  integer :: N, i
  character, allocatable :: cfilename(:)
  N = len(trim(filename))
  allocate(cfilename(N+1))
  do i = 1, N
    cfilename(i) = filename(i:i)
  end do
  cfilename(N+1) = CHAR(0)
  call tsgrea(gridID, cfilename)
  deallocate(cfilename)
end subroutine tsgRead
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
subroutine tsgGetLoadedPointsStatic(gridID, dims, points)
  integer :: gridID, dims
  double precision :: points(dims,*)
  call tsgglp(gridID, points)
end subroutine tsgGetLoadedPointsStatic
!=======================================================================
subroutine tsgGetNeededPointsStatic(gridID, dims, points)
  integer :: gridID, dims
  double precision :: points(dims,*)
  call tsggdp(gridID, points)
end subroutine tsgGetNeededPointsStatic
!=======================================================================
subroutine tsgGetPointsStatic(gridID, dims, points)
  integer :: gridID, dims
  double precision :: points(dims,*)
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
  double precision :: y_c_style(size(y,2),size(y,1))
  call tsgehf(gridID, x, numX, y_c_style)
  y = transpose(y_c_style)
end subroutine tsgEvaluateHierarchicalFunctions
!=======================================================================
subroutine tsgEvaluateComplexHierarchicalFunctions(gridID, x, numX, y)
  integer :: gridID, numX
  double precision :: x(:,:)
  double complex   :: y(:,:)
  double precision :: y_c_style(2*size(y,2),size(y,1))
  integer :: i, j
  call tsgehf(gridID, x, numX, y_c_style)
  do i = 1,size(y,1)
    do j = 1,size(y,2)
      y(i,j) = complex( y_c_style(2*j-1,i), y_c_style(2*j,i) )
    enddo
  enddo
end subroutine tsgEvaluateComplexHierarchicalFunctions
!=======================================================================
subroutine tsgEvaluateSparseHierarchicalFunctions(gridID, x, numX, pntr, indx, y)
  integer :: gridID, numX
  double precision :: x(:,:)
  integer, pointer :: pntr(:), indx(:)
  double precision, pointer :: y(:)
  integer :: numNZ
  call tsgehz(gridID, x, numX, numNZ)
  allocate( pntr(numX+1), indx(numNZ), y(numNZ) )
  call tsgehs(gridID, x, numX, pntr, indx, y)
end subroutine tsgEvaluateSparseHierarchicalFunctions
!=======================================================================
function tsgGetHierarchicalCoefficients(gridID) result(c)
  integer :: gridID
  double precision, pointer :: c(:)
  allocate(c(tsgGetNumOutputs(gridID)*tsgGetNumPoints(gridID)))
  call tsgghc(gridID, c)
end function tsgGetHierarchicalCoefficients
!=======================================================================
subroutine tsgGetHierarchicalCoefficientsStatic(gridID, c)
  integer :: gridID
  double precision :: c(:)
  call tsgghc(gridID, c)
end subroutine tsgGetHierarchicalCoefficientsStatic
!=======================================================================
function tsgGetComplexHierarchicalCoefficients(gridID) result(c)
  integer :: gridID
  double complex, pointer :: c(:)
  allocate(c(tsgGetNumOutputs(gridID)*tsgGetNumPoints(gridID)))
  call tsggchc(gridID, c)
end function tsgGetComplexHierarchicalCoefficients
!=======================================================================
subroutine tsgGetComplexHierarchicalCoefficientsStatic(gridID, c)
  integer :: gridID
  double complex :: c(:)
  call tsggchc(gridID, c)
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
function tsgIsSetDomainTransfrom(gridID) result(res)
  integer :: gridID, res
  call tsgidt(gridID, res)
end function tsgIsSetDomainTransfrom
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
subroutine tsgSetAnisotropicRefinement(gridID, gtype, minGrowth, output)
  integer :: gridID, gtype, minGrowth, output
  call tsgsar(gridID, gtype, minGrowth, output)
end subroutine tsgSetAnisotropicRefinement
!=======================================================================
function tsgEstimateAnisotropicCoefficients(gridID, gtype, output) &
                                              result(coeff)
  integer :: gridID, gtype, output, N
  integer, allocatable :: coeff(:)
  N = tsgGetNumDimensions(gridID)
  if ((gtype .EQ. 2) .OR. (gtype .EQ. 4) .OR. (gtype .EQ. 6))then
    N = N * 2
  endif
  allocate(coeff(N))
  call tsgeac(gridID, gtype, output, coeff)
end function tsgEstimateAnisotropicCoefficients
!=======================================================================
subroutine tsgSetGlobalSurplusRefinement(gridID, tolerance, output)
  integer :: gridID, output
  double precision :: tolerance
  call tsgssr(gridID, tolerance, output)
end subroutine tsgSetGlobalSurplusRefinement
!=======================================================================
subroutine tsgSetLocalSurplusRefinement(gridID, tolerance, rtype, output)
  integer :: gridID, rtype, theout
  integer, optional :: output
  double precision :: tolerance
  if(present(output))then
    theout = output
  else
    theout = -1
  endif
  call tsgshr(gridID, tolerance, rtype, theout)
end subroutine tsgSetLocalSurplusRefinement
!=======================================================================
subroutine tsgClearRefinement(gridID)
  integer :: gridID
  call tsgcre(gridID)
end subroutine tsgClearRefinement
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
