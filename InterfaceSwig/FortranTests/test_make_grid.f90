!==================================================================================================================================================================================
! Copyright (c) 2018, Miroslav Stoyanov
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

subroutine test_make_global_grid()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE), dimension(:,:), allocatable :: points, ref_points

    grid = TasmanianSparseGrid()
    call grid%makeGlobalGrid(2, 1, 1, tsg_type_level, tsg_rule_clenshawcurtis)

    allocate(points(grid%getNumDimensions(), grid%getNumPoints()))
    call grid%getPoints(points(:, 1))

    ! correct set of points, the generated set above must match the ref_points
    ref_points = reshape([ 0.0D-0, 0.0D-0, 0.0D-0, -1.0D-0, 0.0D-0, 1.0D-0, -1.0D-0, 0.0D-0, 1.0D-0, 0.0D-0 ], [ 2, 5 ])

    call approx2d(2, 5, points, ref_points)

    call tassert(grid%isGlobal())
    call tassert(.NOT. grid%isLocalPolynomial())
    call tassert(grid%getAlpha() .EQ. 0.0D-0)
    call tassert(grid%getBeta()  .EQ. 0.0D-0)
    call tassert(grid%getOrder() .EQ. -1)

    deallocate(points)
    call grid%release()
end subroutine

subroutine test_make_localp_grid()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE), dimension(:,:), allocatable :: points, ref_points

    grid = TasmanianSparseGrid()
    call grid%makeLocalPolynomialGrid(3, 0, 1, 3, tsg_rule_localp0)

    allocate(points(grid%getNumDimensions(), grid%getNumPoints()))
    call grid%getLoadedPoints(points(:, 1))

    ! correct set of points, the generated set above must match the ref_points
    ref_points = reshape([ 0.0D-0,  0.0D-0,  0.0D-0, &
                           0.0D-0,  0.0D-0, -0.5D-0, &
                           0.0D-0,  0.0D-0,  0.5D-0, &
                           0.0D-0, -0.5D-0,  0.0D-0, &
                           0.0D-0,  0.5D-0,  0.0D-0, &
                          -0.5D-0,  0.0D-0,  0.0D-0, &
                           0.5D-0,  0.0D-0,  0.0D-0  &
                          ], [ 3, 7 ])

    call approx2d(2, 7, points, ref_points)

    call tassert(grid%getOrder() .EQ. 3)

    deallocate(points)
    call grid%release()
end subroutine

subroutine test_make_sequence_grid()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid, grid_read
    real(C_DOUBLE), dimension(:,:), allocatable :: points, ref_points

    grid = TasmanianSparseGrid()
    grid_read = TasmanianSparseGrid()

    call grid%makeSequenceGrid(1, 1, 3, tsg_type_iptotal, tsg_rule_leja)
    call grid%write("f03_test_file", .false.)
    call grid_read%read("f03_test_file")

    ref_points = reshape([0.0D+0, 1.0D+0, -1.0D+0, sqrt(1.0D+0/3.0D+0)], [1, 4])
    call tassert(grid_read%getNumPoints() .EQ. 4)

    allocate(points(grid%getNumDimensions(), grid%getNumPoints()))
    call grid_read%getNeededPoints(points(:, 1))

    call approx2d(1, 4, points, ref_points)
    deallocate(points)
    call grid%release()
    call grid_read%release()
end subroutine

subroutine test_make_fourier_grid()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    integer :: i
    type(TasmanianSparseGrid) :: grid, grid_copy
    real(C_DOUBLE), dimension(:), allocatable :: weights, ref_weights

    grid = TasmanianSparseGrid()
    grid_copy = TasmanianSparseGrid()

    call grid%makeFourierGrid(1, 1, 2, tsg_type_level)
    call grid_copy%copyGrid(grid)

    call tassert(grid_copy%getNumPoints() .EQ. 9)

    allocate(weights(grid_copy%getNumPoints()))
    allocate(ref_weights(grid_copy%getNumPoints()))
    call grid_copy%getQuadratureWeights(weights)

    do i = 1, 9
        ref_weights(i) = 1.0D-0 / 9.0D-0
    enddo

    call approx1d(9, weights, ref_weights)

    deallocate(weights, ref_weights)
    call grid%release()
    call grid_copy%release()
end subroutine

subroutine test_make_wavelet_grid()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    integer :: i
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE) :: weights(3), ref_weights(3), x(1)

    grid = TasmanianSparseGrid()
    call grid%makeWaveletGrid(1, 1, 0)

    call tassert(grid%getNumPoints() .EQ. 3)
    call tassert(grid%getNumDimensions() .EQ. 1)
    call tassert(grid%getRule() .EQ. tsg_rule_wavelet)
    call tassert(.NOT. grid%isSetDomainTransfrom())

    x = [-0.5]
    ref_weights = [0.5, 0.5, 0.0]
    call grid%getInterpolationWeights(x, weights)

    call approx1d(3, weights, ref_weights)

    call grid%release()
end subroutine

subroutine test_make_grid()
    call test_make_global_grid()
    call test_make_localp_grid()
    call test_make_sequence_grid()
    call test_make_fourier_grid()
    call test_make_wavelet_grid()
    write(*,*) "  Performing tests on make***Grid() methods:       PASS"
end subroutine
