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
program FORSWIGTESTER
    use Tasmanian
    implicit none

write(*,*)
write(*,'(a)') 'Testing Tasmanian Fortran 2003-SWIG interface'

call test_make_global_grid()

write(*,*)


contains
subroutine test_make_global_grid()
    type(TasmanianSparseGrid) :: grid
    double precision, dimension(:,:), allocatable :: points, ref_points
    integer :: i, j

    grid = TasmanianSparseGrid()
    call grid%makeGlobalGrid(2, 1, 1, tsg_type_level, tsg_rule_clenshawcurtis)
    allocate(points(grid%getNumDimensions(), grid%getNumPoints()))

    call grid%getPoints(reshape(points, [size(points)]))

    ! correct set of points, the generated set above must match the ref_points
    ref_points = reshape((/ 0.0D-0, 0.0D-0, 0.0D-0, -1.0D-0, 0.0D-0, 1.0D-0, -1.0D-0, 0.0D-0, 1.0D-0, 0.0D-0 /), (/ 2, 5 /) )

    ! Points don't seem to be initialized here, some of them contain random garbage numbers
    write(*,*) "Dimensions should be 2, 5 = ", grid%getNumDimensions(), grid%getNumPoints()
    do j = 1, 5
        do i = 1, 2
            write(*,*) points(i, j)
        enddo
    enddo

    call approx2d(2, 5, points, ref_points)

    deallocate(points)
    call grid%release()
end subroutine


! Tests if two matrices with dimension n by m are approximately the same
! Calls "error stop" if the entries don't match to 1.D-12
subroutine approx2d(n, m, A, B)
    integer :: n, m
    integer :: i, j
    double precision, dimension(:,:) :: A, B

    do j = 1, m
        do i = 1, n
            if ( abs(A(i, j) - B(i, j)) > 1.D-12 ) then
                error stop
            endif
        enddo
    enddo
end subroutine

end program FORSWIGTESTER
