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

subroutine test_exact_eval_linear()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE), dimension(:,:), pointer :: points
    real(C_DOUBLE), dimension(:,:), pointer :: values
    real(C_DOUBLE) :: x(2, 3), y(2, 3), y_ref(2, 3)
    integer :: i, num_points

    grid = TasmanianSparseGrid()
    call grid%makeLocalPolynomialGrid(2, 2, 3)

    points => grid%returnNeededPoints()
    num_points = grid%getNumNeeded()
    allocate(values(2, num_points))

    do i = 1, num_points
    ! constant and linear models are captures exactly by the linear grid
        values(1, i) = 3.0D-0
        values(2, i) = 2.0D-0 + points(1, i) + 2.0D+0 * points(2, i)
    enddo
    call grid%loadNeededPoints(values(:, 1))

    x = reshape([ 0.3D-0, 0.3D-0, 0.7D-0, -0.3D-0, 0.44D-0, -0.11D-0 ], [2, 3])
    do i = 1, 3
        y_ref(1, i) = 3.0D-0
        y_ref(2, i) = 2.0D-0 + x(1, i) + 2.0D+0 * x(2, i)
    enddo

    call grid%evaluate(x(:,1), y(:,1))
    call approx2d(2, 1, y, y_ref)

    y = 0.0D-0
    call grid%evaluateFast(x(:,1), y(:,1))
    call approx2d(2, 1, y, y_ref)

    y(1, 1) = 0.0D-0
    y(2, 1) = 0.0D-0
    call grid%evaluateBatch(x(:,1), 3, y(:,1))
    call approx2d(2, 3, y, y_ref)

    deallocate(points, values)
    call grid%release()
end subroutine

subroutine test_exact_eval_cubic()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE), dimension(:,:), pointer :: points, lpoints
    real(C_DOUBLE), dimension(:,:), pointer :: values
    real(C_DOUBLE) :: x(2, 4), y(3, 4), y_ref(3, 4)
    integer :: i, num_points

    grid = TasmanianSparseGrid()
    call grid%makeLocalPolynomialGrid(2, 3, 3, 3)

    points => grid%returnNeededPoints()
    num_points = grid%getNumNeeded()
    allocate(values(3, num_points))

    do i = 1, num_points
    ! constant and linear models are captures exactly by the linear grid
        values(1, i) = 3.0D-0
        values(2, i) = points(1, i) + 2.0D+0 * points(2, i)
        values(3, i) = 5.0D-0 + points(1, i)**3 + 2.0D+0 * points(2, i)**2
    enddo
    call grid%loadNeededValues(values(:, 1))
    lpoints => grid%returnLoadedPoints()
    call approx2d(2, num_points, points, lpoints)

    x = reshape([ 0.3D-0, 0.3D-0, 0.7D-0, -0.3D-0, 0.44D-0, -0.11D-0, -0.17D-0, -0.83D-0 ], [2, 4])
    do i = 1, 4
        y_ref(1, i) = 3.0D-0
        y_ref(2, i) = x(1, i) + 2.0D+0 * x(2, i)
        y_ref(3, i) = 5.0D-0 + x(1, i)**3 + 2.0D+0 * x(2, i)**2
    enddo

    call grid%evaluate(x(:,1), y(:,1))
    call approx2d(3, 1, y, y_ref)

    y(1:3, 1) = 0.0D-0
    call grid%evaluateBatch(x(:,1), 4, y(:,1))
    call approx2d(3, 4, y, y_ref)

    y(1:3, 1) = 0.0D-0
    call grid%integrate(y(:,1))
    call tassert(abs(y(1,1) - 12.0D-0) < 1.0D-11)

    deallocate(points, lpoints, values)
    call grid%release()
end subroutine

subroutine test_eval_sequence()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid, grid1, grid2
    real(C_DOUBLE), dimension(:,:), pointer :: points
    real(C_DOUBLE), dimension(:,:), allocatable :: values
    real(C_DOUBLE) :: x(2, 4), y1(1, 4), y2(2, 4), y_ref1(1, 4), y_ref2(2, 4)
    real(C_FLOAT) :: xf(2, 4), yf(1, 4), y_ref1f(1, 4)
    integer :: i, num_points

    grid  = TasmanianSparseGrid()
    grid1 = TasmanianSparseGrid()
    grid2 = TasmanianSparseGrid()

    call grid%makeSequenceGrid(2, 3, 6, tsg_type_level, tsg_rule_mindelta)
    points => grid%returnNeededPoints()
    num_points = grid%getNumNeeded()

    allocate(values(3, num_points))
    do i = 1, num_points
        values(1, i) = 1.0D-0 + points(1, i)
        values(2, i) = 2.0D-0 * points(2, i)
        values(3, i) = points(1, i)**2 + points(2, i)**4 + points(1, i) * points(2, i)
    enddo
    call grid%loadNeededPoints(values(:,1))

    call grid1%copyGrid(grid, 0, 1)
    call grid2%copyGrid(grid, 1)

    x = reshape([ 0.3D-0, 0.4D-0, 0.7D-0, -0.7D-0, 0.44D-0, -0.11D-0, -0.13D-0, -0.83D-0 ], [2, 4])
    do i = 1, 4
        y_ref1(1, i) = 1.0D-0 + x(1, i)
        y_ref2(1, i) = 2.0D-0 * x(2, i)
        y_ref2(2, i) = x(1, i)**2 + x(2, i)**4 + x(1, i) * x(2, i)
        y_ref1f(1, i) = y_ref1(1, i)
    enddo

    call grid1%evaluate(x(:,1), y1(:,1))
    call approx2d(1, 1, y1, y_ref1)

    call grid1%enableAcceleration(tsg_accel_none)
    call grid1%evaluateBatch(x(:,1), 4, y1(:,1))
    call approx2d(1, 4, y1, y_ref1)

    call grid2%enableAcceleration(tsg_accel_cpu_blas)
    call grid2%evaluateBatch(x(:,1), 4, y2(:,1))
    call approx2d(2, 4, y2, y_ref2)

    if (grid2%isAccelerationAvailable(tsg_accel_cpu_blas)) then
        call tassert(grid2%getAccelerationType() == tsg_accel_cpu_blas)
    else
        call tassert(grid2%getAccelerationType() == tsg_accel_none)
    endif

    if (grid1%isAccelerationAvailable(tsg_accel_gpu_cuda)) then
        call grid1%enableAcceleration(tsg_accel_gpu_cuda)
        do i = 1, 4
            xf(1,i) = x(1,i)
            xf(2,i) = x(2,i)
        enddo
        call grid1%evaluateBatch(xf(:,1), 4, yf(:,1))
        call approx2df(1, 4, yf, y_ref1f)

        ! reset yf to call eval fast in single precision
        yf(1, 1) = -333.33
        call grid1%evaluateFast(xf(:,1), yf(:,1))
        call approx2df(1, 1, yf, y_ref1f)

        call tassert(grid1%getGPUID() == 0)
        call grid1%setGPUID(0) ! just for coverage
    endif

    deallocate(points, values)
    call grid%release()
    call grid1%release()
    call grid2%release()
end subroutine

subroutine test_eval_surrogate()
    call test_exact_eval_linear()
    call test_exact_eval_cubic()
    call test_eval_sequence()
    write(*,*) "  Performing tests on evaluate methods:            PASS"
end subroutine
