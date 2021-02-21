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

! Tests if two matrices with dimension n by m are approximately the same
! Calls "error stop" if the entries don't match to 1.D-12
subroutine approx2d(n, m, A, B)
    integer :: n, m
    integer :: i, j
    double precision, dimension(n,m) :: A, B

    do j = 1, m
        do i = 1, n
            if ( abs(A(i, j) - B(i, j)) > 1.D-12 ) then
                error stop
            endif
        enddo
    enddo
end subroutine

! Same as approx2d() but uses a 1D array
subroutine approx1d(n, x, y)
    integer :: n, i
    double precision, dimension(n) :: x, y

    do i = 1, n
        if ( abs(x(i) - y(i)) > 1.D-12 ) then
            write(*,*) x(i), y(i)
            error stop
        endif
    enddo
end subroutine

! Same as approx2d() but works in single precision
subroutine approx2df(n, m, A, B)
    integer :: n, m
    integer :: i, j
    real, dimension(n,m) :: A, B

    do j = 1, m
        do i = 1, n
            if ( abs(A(i, j) - B(i, j)) > 1.E-4 ) then
                write(*,*) "Error exceeds tolerance"
                error stop
            endif
        enddo
    enddo
end subroutine

! similar to cassert, exits with "error stop" if the variable is false
subroutine tassert(x)
    logical :: x

    if (.NOT. x) then
        error stop
    endif
end subroutine

! compare grid points and weights
subroutine approx_grid_pw(grid, grid_ref)
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid), intent(in) :: grid, grid_ref
    real(C_DOUBLE), dimension(:), pointer :: weights, weights_ref
    real(C_DOUBLE), dimension(:,:), pointer :: points, points_ref

    weights => tsgGetQuadratureWeights(grid)
    points => grid%returnPoints()

    weights_ref => tsgGetQuadratureWeights(grid_ref)
    points_ref => tsgGetPoints(grid_ref)

    call approx1d(grid_ref%getNumPoints(), weights, weights_ref)
    call approx2d(grid_ref%getNumDimensions(), grid_ref%getNumPoints(), points, points_ref)

    deallocate(weights, weights_ref, points, points_ref)
end subroutine

! print the points of the grid
subroutine print_points(grid)
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid), intent(in) :: grid
    integer :: i, j
    real(C_DOUBLE), dimension(:,:), pointer :: points

    points => tsgGetPoints(grid)
    do i = 1, grid%getNumPoints()
        do j = 1, grid%getNumDimensions()
            write(*, "(ES15.4)", advance="no") points(j, i)
        enddo
        write(*,*)
    enddo
    deallocate(points)

end subroutine
