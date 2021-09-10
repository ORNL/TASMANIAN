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

subroutine test_set_get_rcoeff()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: gridA, gridB
    real(C_DOUBLE), dimension(:,:), pointer :: points, vals, coeffs
    integer :: i, nump

    gridA = TasmanianSequenceGrid(2, 1, 11, tsg_type_level, tsg_rule_rleja)
    gridB = TasmanianSequenceGrid(2, 1, 11, tsg_type_level, tsg_rule_rleja)

    points => gridA%returnPoints()
    nump = gridA%getNumPoints()
    allocate(vals(1, nump))
    do i = 1, nump
        vals(1,i) = exp(points(1,i)**2 - points(2,i))
    enddo
    call gridA%loadNeededPoints(vals(:,1))

    coeffs => gridA%returnHierarchicalCoefficients()

    call gridB%setHierarchicalCoefficients(coeffs(:,1))
    call approx_grid_pv(gridA, gridB)

    deallocate(points, vals, coeffs)

    call gridA%release()
    call gridB%release()
end subroutine

subroutine test_set_get_ccoeff()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: gridA, gridB
    real(C_DOUBLE), dimension(:,:), pointer :: points, vals
    complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: coeffs
    integer :: i, j, nump

    gridA = TasmanianFourierGrid(2, 2, 1, tsg_type_level)
    gridB = TasmanianFourierGrid(2, 2, 1, tsg_type_level)

    call tassert(gridA%isFourier())

    points => gridA%returnPoints()
    nump = gridA%getNumPoints()
    allocate(vals(2, nump))
    do i = 1, nump
        vals(1,i) = points(1,i) ** 3 - points(1,i)
        vals(2,i) = 0.25 * points(1,i) ** 4 - 0.5 * points(1,i) ** 2
    enddo
    call gridA%loadNeededPoints(vals(:,1))

    coeffs => gridA%returnComplexHierarchicalCoefficients()

    call gridB%setComplexHierarchicalCoefficients(coeffs)
    call approx_grid_pv(gridA, gridB)

    deallocate(points, vals, coeffs)

    call gridA%release()
    call gridB%release()
end subroutine

subroutine test_hierarchy_transforms()
    call test_set_get_rcoeff()
    call test_set_get_ccoeff()
    write(*,*) "  Performing tests on hierarchy coefficients:      PASS"
end subroutine
