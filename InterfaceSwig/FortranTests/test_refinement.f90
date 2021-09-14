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

subroutine test_anisotropic_refinement()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE), dimension(:,:), pointer :: points
    real(C_DOUBLE), dimension(:), pointer :: values

    grid = TasmanianSequenceGrid(3, 1, 3, tsg_type_level, tsg_rule_leja, level_limits=(/3,2,1/))
    points => grid%returnPoints()

    call approx2d(3, grid%getNumPoints(), points, &
             reshape( (/ 0.d0, 0.d0, 0.d0, &
                         0.d0, 0.d0, 1.d0, &
                         0.d0, 1.d0, 0.d0, &
                         0.d0, 1.d0, 1.d0, &
                         0.d0,-1.d0, 0.d0, &
                         0.d0,-1.d0, 1.d0, &
                         1.d0, 0.d0, 0.d0, &
                         1.d0, 0.d0, 1.d0, &
                         1.d0, 1.d0, 0.d0, &
                         1.d0, 1.d0, 1.d0, &
                         1.d0,-1.d0, 0.d0, &
                        -1.d0, 0.d0, 0.d0, &
                        -1.d0, 0.d0, 1.d0, &
                        -1.d0, 1.d0, 0.d0, &
                         1.d0/sqrt(3.d0), 0.d0, 0.d0 /), &
                      (/3, grid%getNumPoints()/) ))

    allocate( values(grid%getNumPoints()) )
    values(:) = exp( -points(1,:)**2 - points(2,:)**2 )
    call grid%loadNeededPoints(values)

    call grid%setAnisotropicRefinement(tsg_type_iptotal,5,0)
    if ( grid%getNumNeeded() .eq. 0 ) then
        write(*,*) "Failed to generate new points in refinement"
        error stop
    endif
    deallocate(points)
    points => grid%returnNeededPoints()

    call approx2d(3, grid%getNumNeeded(), points, &
             reshape( (/ 1.d0, -1.d0, 1.d0, &
                        -1.d0,  1.d0, 1.d0, &
                        -1.d0, -1.d0, 0.d0, &
                         1.d0/sqrt(3.d0), 0.d0, 1.d0, &
                         1.d0/sqrt(3.d0), 1.d0, 0.d0 /), &
                      (/3, grid%getNumNeeded()/) ))

    call grid%clearRefinement()
    if ( .not. grid%getNumNeeded() .eq. 0 ) then
        write(*,*) "Failed to clear refinement"
        error stop
    endif

    call grid%release()
    deallocate(points, values)
end subroutine

subroutine test_surplus_refinement()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    real(C_DOUBLE), dimension(:,:), pointer :: points
    real(C_DOUBLE), dimension(:), pointer :: values

    grid = TasmanianLocalPolynomialGrid(3, 1, 1, 1, tsg_rule_localp, level_limits=(/1,2,3/))
    points => grid%returnPoints();

    call approx2d(3, grid%getNumPoints(), points, &
             reshape( (/ 0.d0, 0.d0, 0.d0, &
                         0.d0, 0.d0,-1.d0, &
                         0.d0, 0.d0, 1.d0, &
                         0.d0,-1.d0, 0.d0, &
                         0.d0, 1.d0, 0.d0, &
                        -1.d0, 0.d0, 0.d0, &
                         1.d0, 0.d0, 0.d0 /), &
                      (/3, grid%getNumPoints()/) ))

    allocate( values(grid%getNumPoints()) )
    values(:) = exp( -points(1,:)**2 - points(2,:)**2 )
    call grid%loadNeededPoints(values)

    call grid%setSurplusRefinement(1.D-3, tsg_refine_classic)
    if ( grid%getNumNeeded() .eq. 0 ) then
        write(*,*) "Refinement failed to set new points"
        error stop
    endif

    deallocate(points)
    points => grid%returnNeededPoints()

    call approx2d(3, grid%getNumNeeded(), points, &
             reshape( (/ 0.d0, -1.d0, -1.d0, &
                         0.d0, -1.d0,  1.d0, &
                         0.d0,  1.d0, -1.d0, &
                         0.d0,  1.d0,  1.d0, &
                         0.d0, -5.d-1,  0.d0, &
                         0.d0,  5.d-1,  0.d0, &
                        -1.d0,  0.d0, -1.d0, &
                        -1.d0,  0.d0,  1.d0, &
                        -1.d0, -1.d0,  0.d0, &
                        -1.d0,  1.d0,  0.d0, &
                         1.d0,  0.d0, -1.d0, &
                         1.d0,  0.d0,  1.d0, &
                         1.d0, -1.d0,  0.d0, &
                         1.d0,  1.d0,  0.d0 /), &
                      (/3, grid%getNumNeeded()/) ))

    call grid%clearLevelLimits()
    call grid%setSurplusRefinement(1.D-3, tsg_refine_classic, 0)

    call tassert(grid%getNumNeeded() == 16)

    call grid%setSurplusRefinement(1.D-3, tsg_refine_classic, -1, level_limits=(/3,2,0/))

    deallocate(points)
    points => grid%returnNeededPoints()

    call approx2d(3, grid%getNumNeeded(), points, &
             reshape( (/ 0.d0, -5.d-1, 0.d0, &
                         0.d0,  5.d-1, 0.d0, &
                        -1.d0, -1.d0,  0.d0, &
                        -1.d0,  1.d0,  0.d0, &
                         1.d0, -1.d0,  0.d0, &
                         1.d0,  1.d0,  0.d0, &
                        -5.d-1, 0.d0,  0.d0, &
                         5.d-1, 0.d0,  0.d0 /), &
                      (/3, grid%getNumNeeded()/) ))

    call grid%setSurplusRefinement(1.D-3, tsg_refine_classic, -1, level_limits=(/3,2,0/), &
                                   scale_correction=(/0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0/))

    deallocate(points)
    points => grid%returnNeededPoints()

    call approx2d(3, grid%getNumNeeded(), points, &
             reshape( (/ 0.d0, -5.d-1, 0.d0, &
                        -1.d0, -1.d0,  0.d0, &
                         1.d0, -1.d0,  0.d0 /), &
                      (/3, grid%getNumNeeded()/) ))

    call grid%release()
    deallocate(points, values)
end subroutine

subroutine test_refinement()
    call test_anisotropic_refinement()
    call test_surplus_refinement()
    write(*,*) "  Performing tests on refinement:                  PASS"
end subroutine
