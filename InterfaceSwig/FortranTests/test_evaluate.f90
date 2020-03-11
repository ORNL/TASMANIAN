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

subroutine test_exact_eval()
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

    points => tsgGetNeededPoints(grid)
    num_points = grid%getNumNeeded()
    allocate(values(2, num_points))

    do i = 1, num_points
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

    y(1, 1) = 0.0D-0
    y(2, 1) = 0.0D-0
    call grid%evaluateBatch(x(:,1), 3, y(:,1))
    call approx2d(2, 3, y, y_ref)

    deallocate(points, values)
    call grid%release()
end subroutine


subroutine test_eval_surrogate()
    call test_exact_eval()
    write(*,*) "  Performing tests on evaluate methods:            PASS"
end subroutine
