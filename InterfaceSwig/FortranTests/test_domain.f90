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

subroutine test_domain_range()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE), dimension(:), pointer :: weights
    real(C_DOUBLE), dimension(:,:), pointer :: points
    real(C_DOUBLE), dimension(3) :: lower, upper
    real(C_DOUBLE), dimension(3) :: tlower, tupper
    integer :: i

    lower = [3.0D-0, -7.0D-0, -12.0D-0]
    upper = [5.0D-0, -6.0D-0,  17.0D-0]

    grid = TasmanianGlobalGrid(3, 0, 2, tsg_type_level, tsg_rule_clenshawcurtis, [1, 1, 1])
    call grid%setDomainTransform(lower, upper)
    call tassert(grid%isSetDomainTransfrom())

    weights => tsgGetQuadratureWeights(grid)
    points => tsgGetPoints(grid)

    call tassert(abs(sum(weights) - 58.0D-0) < 1.D-11)
    do i = 1, 3
        call tassert(abs(minval(points(i,:))-lower(i))  < 1.D-11)
        call tassert(abs(maxval(points(i,:))-upper(i))  < 1.D-11)
    enddo

    call grid%getDomainTransform(tlower, tupper)
    call approx1d(3, tlower, lower)
    call approx1d(3, tupper, upper)

    call grid%clearDomainTransform()
    call tassert(.not. grid%isSetDomainTransfrom())

    deallocate(points, weights)
    call grid%release()
end subroutine

subroutine test_domain_gauss_hermite()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE), dimension(:), pointer :: weights
    double precision, parameter :: pi = 4.0D+0 * atan(1.0D+0)

    grid = TasmanianGlobalGrid(1, 0, 4, tsg_type_level, tsg_rule_gausshermite, alpha=2.0D+0)
    weights => tsgGetQuadratureWeights(grid)

    call tassert(abs(sum(weights) - 0.5D+0 * sqrt(pi)) < 1.D-11)
    call tassert(grid%getAlpha() == 2.0D-0)
    call tassert(grid%getBeta() == 0.0D-0)

    deallocate(weights)
    call grid%release()
end subroutine

subroutine test_domain_gauss_jacobi()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE), dimension(:), pointer :: weights

    grid = TasmanianGlobalGrid(1, 0, 4, tsg_type_level, tsg_rule_gaussjacobi, beta=1.0D+0)
    weights => tsgGetQuadratureWeights(grid)

    call tassert(abs(sum(weights) - 2.0D-0) < 1.D-11)
    call tassert(grid%getAlpha() == 0.0D-0)
    call tassert(grid%getBeta() == 1.0D-0)

    deallocate(weights)
    call grid%release()
end subroutine

subroutine test_domain_aniso()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE), dimension(:), pointer :: weights
    real(C_DOUBLE), dimension(:,:), pointer :: points
    real(C_DOUBLE) :: points_ref(2, 4)

    points_ref = reshape([ 0.0D+0, 0.0D+0, 0.0D+0, 1.0D+0, 0.0D+0, -1.0D+0, 1.0D+0, 0.0D+0 ], [2, 4])
    grid = TasmanianGlobalGrid(2, 1, 2, tsg_type_level, tsg_rule_leja, [2, 1])

    weights => tsgGetQuadratureWeights(grid)
    points =>tsgGetNeededPoints(grid)

    call tassert(abs(sum(weights) - 4.0D-0) < 1.E-11)
    call approx2d(points, points_ref)

    deallocate(weights, points)
    call grid%release()
end subroutine

subroutine test_domain_transforms()
    call test_domain_range()
    call test_domain_gauss_hermite()
    call test_domain_gauss_jacobi()
    call test_domain_aniso()
    write(*,*) "  Performing tests on domain transforms:           PASS"
end subroutine
