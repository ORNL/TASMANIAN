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
PUBLIC ! default include all symbols, list only exceptions as PRIVATE ::

  type TasmanianSparseGrid
    double complex :: pntr
  end type TasmanianSparseGrid

  integer, parameter :: tsg_clenshaw_curtis      =  1,  tsg_clenshaw_curtis_zero  =  2, &
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
                        tsg_fourier              = 41,  tsg_localpb               = 42, &
                        tsg_level        =   1,  tsg_curved       =  2, &
                        tsg_iptotal      =   3,  tsg_ipcurved     =  4, &
                        tsg_qptotal      =   5,  tsg_qpcurved     =  6, &
                        tsg_hyperbolic   =   7,  tsg_iphyperbolic =  8, &
                        tsg_qphyperbolic =   9,  tsg_tensor       = 10, &
                        tsg_iptensor     =  11,  tsg_qptensor     = 12, &
                        tsg_classic     = 1,  tsg_parents_first = 2, &
                        tsg_directional = 3,  tsg_fds           = 4, &
                        tsg_stable      = 5,                         &
                        tsg_acc_none     = 0, tsg_acc_cpu_blas  = 1, tsg_acc_gpu_cublas = 2, &
                        tsg_acc_gpu_rocblas = 2, tsg_acc_gpu_hip = 3, &
                        tsg_acc_gpu_cuda = 3, tsg_acc_gpu_magma = 4

  integer :: rows, cols, length
  double precision, pointer :: matrix(:,:), vector(:)
  character, pointer :: string(:)

PRIVATE :: rows, cols, length, matrix, vector, string

contains
!=======================================================================
subroutine tsgAllocateGrid(grid)
  type(TasmanianSparseGrid) :: grid
  call tsgalloc(grid%pntr);
end subroutine tsgAllocateGrid
subroutine tsgDeallocateGrid(grid)
  type(TasmanianSparseGrid) :: grid
  call tsgfree(grid%pntr);
end subroutine tsgDeallocateGrid
!=======================================================================
function tsgIsGlobal(grid) result(res)
  type(TasmanianSparseGrid) :: grid
  integer :: flag
  logical :: res
  call tsgisg(grid%pntr, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
end function
!=======================================================================
function tsgIsSequence(grid) result(res)
  type(TasmanianSparseGrid) :: grid
  integer :: flag
  logical :: res
  call tsgiss(grid%pntr, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
end function
!=======================================================================
function tsgIsLocalPolynomial(grid) result(res)
  type(TasmanianSparseGrid) :: grid
  integer :: flag
  logical :: res
  call tsgisl(grid%pntr, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
end function
!=======================================================================
function tsgIsWavelet(grid) result(res)
  type(TasmanianSparseGrid) :: grid
  integer :: flag
  logical :: res
  call tsgisw(grid%pntr, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
end function
!=======================================================================
function tsgIsFourier(grid) result(res)
  type(TasmanianSparseGrid) :: grid
  integer :: flag
  logical :: res
  call tsgisf(grid%pntr, flag)
  if (flag .ne. 0) then
    res = .true.
  else
    res = .false.
  endif
end function
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
subroutine tsgMakeGlobalGrid(grid, dims, outs, depth, gtype, rule, &
                             aweights, alpha, beta, customRuleFilename, levelLimits)
  type(TasmanianSparseGrid) :: grid
  integer, intent(in) :: dims, outs, depth, gtype, rule
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

  call tsgmg(grid%pntr, dims, outs, depth, gtype, rule, opt_flags, aw, al, be, cfn, ll)
end subroutine tsgMakeGlobalGrid
!=======================================================================
subroutine tsgMakeSequenceGrid(grid, dims, outs, depth, gtype, rule, aweights, levelLimits)
  type(TasmanianSparseGrid) :: grid
  integer :: dims, outs, depth, gtype, rule
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

  call tsgms(grid%pntr, dims, outs, depth, gtype, rule, opt_flags, aw, ll)
end subroutine tsgMakeSequenceGrid
!=======================================================================
subroutine tsgMakeLocalPolynomialGrid(grid, dims, outs, depth, order, rule, levelLimits)
  type(TasmanianSparseGrid) :: grid
  integer :: dims, outs, depth
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

  call tsgml(grid%pntr, dims, outs, depth, opt_flags, or, ru, ll)
end subroutine tsgMakeLocalPolynomialGrid
!=======================================================================
subroutine tsgMakeWaveletGrid(grid, dims, outs, depth, order, levelLimits)
  type(TasmanianSparseGrid) :: grid
  integer :: dims, outs, depth
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

  call tsgmw(grid%pntr, dims, outs, depth, opt_flags, or, ll)
end subroutine tsgMakeWaveletGrid
!=======================================================================
subroutine tsgMakeFourierGrid(grid, dims, outs, depth, gtype, aweights, levelLimits)
  type(TasmanianSparseGrid) :: grid
  integer :: dims, outs, depth, gtype
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

  call tsgmf(grid%pntr, dims, outs, depth, gtype, opt_flags, aw, ll)
end subroutine tsgMakeFourierGrid
!=======================================================================
subroutine tsgCopyGrid(grid, source)
  type(TasmanianSparseGrid) :: grid, source
  call tsgcp(grid%pntr, source)
end subroutine tsgCopyGrid
!=======================================================================
subroutine tsgUpdateGlobalGrid(grid, depth, gtype, aweights)
  type(TasmanianSparseGrid) :: grid
  integer, intent(in) :: depth, gtype
  integer, optional, target :: aweights(:)
  integer          :: opt_flags = 0
  integer, pointer :: aw(:) => null()
  if ( present(aweights) ) then
    opt_flags = 1
    aw => aweights
  endif
  call tsgug(grid%pntr, depth, gtype, opt_flags, aw)
end subroutine tsgUpdateGlobalGrid
!=======================================================================
subroutine tsgUpdateSequenceGrid(grid, depth, gtype, aweights)
  type(TasmanianSparseGrid) :: grid
  integer, intent(in) :: depth, gtype
  integer, optional, target :: aweights(:)
  integer          :: opt_flags = 0
  integer, pointer :: aw(:) => null()
  if ( present(aweights) ) then
    opt_flags = 1
    aw => aweights
  endif
  call tsgus(grid%pntr, depth, gtype, opt_flags, aw)
end subroutine tsgUpdateSequenceGrid
!=======================================================================
function tsgRead(grid, filename) result(res)
  type(TasmanianSparseGrid) :: grid
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
  call tsgrea(grid%pntr, cfilename, flag)
  res = (flag .ne. 0)
  deallocate(cfilename)
end function tsgRead
!=======================================================================
subroutine tsgWrite(grid, filename, useBinary)
  type(TasmanianSparseGrid) :: grid
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
    ubin = 1
  endif
  call tsgwri(grid%pntr, cfilename, ubin)
  deallocate(cfilename)
end subroutine tsgWrite
!=======================================================================
function tsgGetAlpha(grid) result(alpha)
  type(TasmanianSparseGrid) :: grid
  double precision :: alpha
  call tsggal(grid%pntr, alpha)
end function tsgGetAlpha
!=======================================================================
function tsgGetBeta(grid) result(beta)
  type(TasmanianSparseGrid) :: grid
  double precision :: beta
  call tsggbe(grid%pntr, beta)
end function tsgGetBeta
!=======================================================================
function tsgGetOrder(grid) result(order)
  type(TasmanianSparseGrid) :: grid
  integer :: order
  call tsggor(grid%pntr, order)
end function tsgGetOrder
!=======================================================================
function tsgGetNumDimensions(grid) result(dims)
  type(TasmanianSparseGrid) :: grid
  integer :: dims
  call tsggnd(grid%pntr, dims)
end function tsgGetNumDimensions
!=======================================================================
function tsgGetNumOutputs(grid) result(outs)
  type(TasmanianSparseGrid) :: grid
  integer :: outs
  call tsggno(grid%pntr, outs)
end function tsgGetNumOutputs
!=======================================================================
function tsgGetRule(grid) result(rule)
  type(TasmanianSparseGrid) :: grid
  integer :: rule
  call tsggru(grid%pntr, rule)
end function tsgGetRule
!=======================================================================
function tsgGetNumLoaded(grid) result(num)
  type(TasmanianSparseGrid) :: grid
  integer :: num
  call tsggnl(grid%pntr, num)
end function tsgGetNumLoaded
!=======================================================================
function tsgGetNumNeeded(grid) result(num)
  type(TasmanianSparseGrid) :: grid
  integer :: num
  call tsggnn(grid%pntr, num)
end function tsgGetNumNeeded
!=======================================================================
function tsgGetNumPoints(grid) result(num)
  type(TasmanianSparseGrid) :: grid
  integer :: num
  call tsggnp(grid%pntr, num)
end function tsgGetNumPoints
!=======================================================================
function tsgGetLoadedPoints(grid) result(p)
  type(TasmanianSparseGrid) :: grid
  double precision, pointer :: p(:,:)
  integer :: rows, cols
  rows = tsgGetNumDimensions(grid)
  cols = tsgGetNumLoaded(grid)
  allocate(p(rows,cols))
  call tsgglp(grid%pntr, p)
end function tsgGetLoadedPoints
!=======================================================================
function tsgGetNeededPoints(grid) result(p)
  type(TasmanianSparseGrid) :: grid
  double precision, pointer :: p(:,:)
  integer :: rows, cols
  rows = tsgGetNumDimensions(grid)
  cols = tsgGetNumNeeded(grid)
  allocate(p(rows,cols))
  call tsggdp(grid%pntr, p)
end function tsgGetNeededPoints
!=======================================================================
function tsgGetPoints(grid) result(p)
  type(TasmanianSparseGrid) :: grid
  double precision, pointer :: p(:,:)
  integer :: rows, cols
  rows = tsgGetNumDimensions(grid)
  cols = tsgGetNumPoints(grid)
  allocate(p(rows,cols))
  call tsggpp(grid%pntr, p)
end function tsgGetPoints
!=======================================================================
subroutine tsgGetLoadedPointsStatic(grid, points)
  type(TasmanianSparseGrid) :: grid
  double precision :: points(:,:)
  call tsgglp(grid%pntr, points)
end subroutine tsgGetLoadedPointsStatic
!=======================================================================
subroutine tsgGetNeededPointsStatic(grid, points)
  type(TasmanianSparseGrid) :: grid
  double precision :: points(:,:)
  call tsggdp(grid%pntr, points)
end subroutine tsgGetNeededPointsStatic
!=======================================================================
subroutine tsgGetPointsStatic(grid, points)
  type(TasmanianSparseGrid) :: grid
  double precision :: points(:,:)
  call tsggpp(grid%pntr, points)
end subroutine tsgGetPointsStatic
!=======================================================================
function tsgGetQuadratureWeights(grid) result(w)
  type(TasmanianSparseGrid) :: grid
  integer :: length
  double precision, pointer :: w(:)
  length = tsgGetNumPoints(grid)
  allocate(w(length))
  call tsggqw(grid%pntr, w)
end function tsgGetQuadratureWeights
!=======================================================================
subroutine tsgGetQuadratureWeightsStatic(grid, weights)
  type(TasmanianSparseGrid) :: grid
  double precision :: weights(*)
  call tsggqw(grid%pntr, weights)
end subroutine tsgGetQuadratureWeightsStatic
!=======================================================================
function tsgGetInterpolationWeights(grid, x) result(w)
  type(TasmanianSparseGrid) :: grid
  integer :: length
  double precision :: x(*)
  double precision, pointer :: w(:)
  length = tsgGetNumPoints(grid)
  allocate(w(length))
  call tsggiw(grid%pntr, x, w)
end function tsgGetInterpolationWeights
!=======================================================================
subroutine tsgGetInterpolationWeightsStatic(grid, x, weights)
  type(TasmanianSparseGrid) :: grid
  double precision :: x(:)
  double precision :: weights(:)
  call tsggiw(grid%pntr, x, weights)
end subroutine tsgGetInterpolationWeightsStatic
!=======================================================================
subroutine tsgLoadNeededPoints(grid, values)
  type(TasmanianSparseGrid) :: grid
  double precision :: values(:,:)
  call tsglnp(grid%pntr, values)
end subroutine tsgLoadNeededPoints
!=======================================================================
subroutine tsgEvaluate(grid, x, y)
  type(TasmanianSparseGrid) :: grid
  double precision, intent(in) :: x(:)
  double precision :: y(:)
  call tsgeva(grid%pntr, x, y)
end subroutine tsgEvaluate
!=======================================================================
subroutine tsgEvaluateFast(grid, x, y)
  type(TasmanianSparseGrid) :: grid
  double precision, intent(in) :: x(:)
  double precision :: y(:)
  call tsgevf(grid%pntr, x, y)
end subroutine tsgEvaluateFast
!=======================================================================
subroutine tsgEvaluateBatch(grid, x, numX, y)
  type(TasmanianSparseGrid) :: grid
  integer :: numX
  double precision :: x(:,:), y(:,:)
  call tsgevb(grid%pntr, x, numX, y)
end subroutine tsgEvaluateBatch
!=======================================================================
subroutine tsgEvaluateHierarchicalFunctions(grid, x, numX, y)
  type(TasmanianSparseGrid) :: grid
  integer :: numX
  double precision :: x(:,:), y(:,:)
  if ( .not. tsgIsFourier(grid) ) then
    call tsgehf(grid%pntr, x, numX, y)
  else
    write(*,*) "ERROR: called tsgEvaluateHierarchicalFunctions() on a Fourier grid, "
    write(*,*) "       use tsgEvaluateComplexHierarchicalFunctions() instead"
  endif
end subroutine tsgEvaluateHierarchicalFunctions
!=======================================================================
subroutine tsgEvaluateComplexHierarchicalFunctions(grid, x, numX, y)
  type(TasmanianSparseGrid) :: grid
  integer :: numX
  double precision :: x(:,:)
  double complex   :: y(:,:)
  double precision, allocatable :: y_c_style(:)
  integer :: i, j
  if ( tsgIsFourier(grid) ) then
    allocate(y_c_style(2*size(y,1)*size(y,2)))
    call tsgehf(grid%pntr, x, numX, y_c_style)
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
subroutine tsgEvaluateSparseHierarchicalFunctions(grid, x, numX, pntr, indx, y)
  type(TasmanianSparseGrid) :: grid
  integer :: numX
  double precision :: x(:,:)
  integer, pointer :: pntr(:), indx(:)
  double precision, pointer :: y(:)
  integer :: numNZ
  if ( .not. tsgIsFourier(grid) ) then
    call tsgehz(grid%pntr, x, numX, numNZ)
    allocate( pntr(numX+1), indx(numNZ), y(numNZ) )
    call tsgehs(grid%pntr, x, numX, pntr, indx, y)
  else
    write(*,*) "ERROR: called tsgEvaluateSparseHierarchicalFunctions() on a Fourier grid"
  endif
end subroutine tsgEvaluateSparseHierarchicalFunctions
!=======================================================================
function tsgGetHierarchicalCoefficients(grid) result(c)
  type(TasmanianSparseGrid) :: grid
  double precision, pointer :: c(:)
  if ( .not. tsgIsFourier(grid) ) then
    allocate(c(tsgGetNumOutputs(grid)*tsgGetNumPoints(grid)))
    call tsgghc(grid%pntr, c)
  else
    write(*,*) "ERROR: called tsgGetHierarchicalCoefficients() on a Fourier grid, "
    write(*,*) "       use tsgGetComplexHierarchicalCoefficients() instead"
  endif
end function tsgGetHierarchicalCoefficients
!=======================================================================
subroutine tsgGetHierarchicalCoefficientsStatic(grid, c)
  type(TasmanianSparseGrid) :: grid
  double precision :: c(:)
  if ( .not. tsgIsFourier(grid) ) then
    call tsgghc(grid%pntr, c)
  else
    write(*,*) "ERROR: called tsgGetHierarchicalCoefficientsStatic() on a Fourier grid, "
    write(*,*) "       use tsgGetComplexHierarchicalCoefficientsStatic() instead"
  endif
end subroutine tsgGetHierarchicalCoefficientsStatic
!=======================================================================
function tsgGetComplexHierarchicalCoefficients(grid) result(c)
  type(TasmanianSparseGrid) :: grid
  double complex, pointer       :: c(:)
  double precision, allocatable :: c_real(:)
  integer :: i
  if ( tsgIsFourier(grid) ) then
    allocate(c(tsgGetNumOutputs(grid)*tsgGetNumPoints(grid)))
    allocate(c_real(2*tsgGetNumOutputs(grid)*tsgGetNumPoints(grid)))
    call tsgghc(grid%pntr, c_real)
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
subroutine tsgGetComplexHierarchicalCoefficientsStatic(grid, c)
  type(TasmanianSparseGrid) :: grid
  double complex   :: c(:)
  double precision :: c_real(2*size(c))
  integer :: i
  if ( tsgIsFourier(grid) ) then
    call tsgghc(grid%pntr, c_real)
    do i = 1,size(c)
      c(i) = cmplx( c_real(i), c_real(i + size(c)), kind(c) )
    enddo
  else
    write(*,*) "ERROR: called tsgGetComplexHierarchicalCoefficientsStatic() on a non-Fourier grid, "
    write(*,*) "       use tsgGetHierarchicalCoefficientsStatic() instead"
  endif
end subroutine tsgGetComplexHierarchicalCoefficientsStatic
!=======================================================================
subroutine tsgSetHierarchicalCoefficients(grid,c)
  type(TasmanianSparseGrid) :: grid
  double precision :: c(:)
  call tsgshc(grid%pntr, c)
end subroutine tsgSetHierarchicalCoefficients
!=======================================================================
function tsgGetHierarchicalSupport(grid) result(c)
  type(TasmanianSparseGrid) :: grid
  double precision, pointer :: c(:,:)
  allocate(c(tsgGetNumDimensions(grid), tsgGetNumPoints(grid)))
  call tsghsu(grid%pntr, c)
end function tsgGetHierarchicalSupport
!=======================================================================
subroutine tsgIntegrate(grid, q)
  type(TasmanianSparseGrid) :: grid
  double precision :: q(:)
  call tsgint(grid%pntr, q)
end subroutine tsgIntegrate
!=======================================================================
subroutine tsgSetDomainTransform(grid, transformA, transformB)
  type(TasmanianSparseGrid) :: grid
  double precision :: transformA(*), transformB(*)
  call tsgsdt(grid%pntr, transformA, transformB)
end subroutine tsgSetDomainTransform
!=======================================================================
function tsgIsSetDomainTransform(grid) result(res)
  type(TasmanianSparseGrid) :: grid
  logical :: res
  call tsgidt(grid%pntr, res)
end function tsgIsSetDomainTransform
!=======================================================================
subroutine tsgClearDomainTransform(grid)
  type(TasmanianSparseGrid) :: grid
  call tsgcdt(grid%pntr)
end subroutine tsgClearDomainTransform
!=======================================================================
subroutine tsgGetDomainTransform(grid, transformA, transformB)
  type(TasmanianSparseGrid) :: grid
  double precision :: transformA(*), transformB(*)
  call tsggdt(grid%pntr, transformA, transformB)
end subroutine tsgGetDomainTransform
!=======================================================================
subroutine tsgSetAnisotropicRefinement(grid, gtype, minGrowth, output, levelLimits)
  type(TasmanianSparseGrid) :: grid
  integer :: gtype, minGrowth, output
  integer, optional, target :: levelLimits(:)
  integer          :: opt_flags = 0
  integer, pointer :: ll(:) => null()
  if (present(levelLimits)) then
    opt_flags = 1
    ll => levelLimits
  endif
  call tsgsar(grid%pntr, gtype, minGrowth, output-1, opt_flags, ll)
end subroutine tsgSetAnisotropicRefinement
!=======================================================================
function tsgEstimateAnisotropicCoefficients(grid, gtype, output) result(coeff)
  type(TasmanianSparseGrid) :: grid
  integer :: gtype, output, N
  integer, pointer :: coeff(:)
  N = tsgGetNumDimensions(grid)
  if ((gtype .EQ. tsg_curved) .OR. (gtype .EQ. tsg_ipcurved) .OR. (gtype .EQ. tsg_qpcurved))then
    N = N * 2
  endif
  allocate(coeff(N))
  call tsgeac(grid%pntr, gtype, output-1, coeff)
end function tsgEstimateAnisotropicCoefficients
!=======================================================================
subroutine tsgSetGlobalSurplusRefinement(grid, tolerance, output, levelLimits)
  type(TasmanianSparseGrid) :: grid
  integer :: output
  integer, optional, target :: levelLimits(:)
  double precision :: tolerance
  integer          :: opt_flags = 0
  integer, pointer :: ll(:) => null()
  if (present(levelLimits)) then
    opt_flags = 1
    ll => levelLimits
  endif
  call tsgssr(grid%pntr, tolerance, output-1, opt_flags, ll)
end subroutine tsgSetGlobalSurplusRefinement
!=======================================================================
subroutine tsgSetLocalSurplusRefinement(grid, tolerance, rtype, output, levelLimits)
  type(TasmanianSparseGrid) :: grid
  integer :: rtype
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

  call tsgshr(grid%pntr, tolerance, rtype, opt_flags, theout, ll)
end subroutine tsgSetLocalSurplusRefinement
!=======================================================================
subroutine tsgClearRefinement(grid)
  type(TasmanianSparseGrid) :: grid
  call tsgcre(grid%pntr)
end subroutine tsgClearRefinement
!=======================================================================
subroutine tsgMergeRefinement(grid)
  type(TasmanianSparseGrid) :: grid
  call tsgmre(grid%pntr)
end subroutine tsgMergeRefinement
!=======================================================================
subroutine tsgSetConformalTransformASIN(grid, truncate)
  type(TasmanianSparseGrid) :: grid
  integer :: truncate(:)
  call tsgsca(grid%pntr, truncate)
end subroutine tsgSetConformalTransformASIN
!=======================================================================
function tsgIsSetConformalTransformASIN(grid) result(isset)
  type(TasmanianSparseGrid) :: grid
  integer :: res
  logical :: isset
  call tsgica(grid%pntr, res)
  if (res .EQ. 0)then
    isset = .FALSE.
  else
    isset = .TRUE.
  endif
end function tsgIsSetConformalTransformASIN
!=======================================================================
subroutine tsgClearConformalTransform(grid)
  type(TasmanianSparseGrid) :: grid
  call tsgcct(grid%pntr)
end subroutine tsgClearConformalTransform
!=======================================================================
function tsgGetConformalTransformASIN(grid) result(truncate)
  type(TasmanianSparseGrid) :: grid
  integer :: N
  integer, allocatable :: truncate(:)
  N = tsgGetNumDimensions(grid)
  allocate(truncate(N))
  call tsggca(grid%pntr, truncate)
end function tsgGetConformalTransformASIN
!=======================================================================
subroutine tsgPrintStats(grid)
  type(TasmanianSparseGrid) :: grid
  call tsgpri(grid%pntr)
end subroutine tsgPrintStats
!=======================================================================
! Addon/MPI methods !
!=======================================================================
subroutine tsgMPIGridSend(grid, destination, tag, comm, ierr)
  type(TasmanianSparseGrid) :: grid
  integer :: destination, tag, comm, ierr
  call tsgmpigsend(grid, destination, tag, comm, ierr)
end subroutine tsgMPIGridSend
subroutine tsgMPIGridRecv(grid, source, tag, comm, stats, ierr)
  type(TasmanianSparseGrid) :: grid
  integer :: source, tag, comm, stats(:), ierr
  call tsgmpigrecv(grid, source, tag, comm, ierr)
end subroutine tsgMPIGridRecv
subroutine tsgMPIGridBcast(grid, root, comm, ierr)
  type(TasmanianSparseGrid) :: grid
  integer :: root, comm, ierr
  call tsgmpigcast(grid, root, comm, ierr)
end subroutine tsgMPIGridBcast
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
