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

subroutine test_update_global()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid, grid_ref
    integer(C_INT) :: aweights(2), level_limits(2)

    aweights = [ 2, 1 ]
    level_limits = [ 2, 1 ]

    grid_ref = TasmanianGlobalGrid(2, 1, 3, tsg_type_level, tsg_rule_clenshawcurtis)
    grid = TasmanianGlobalGrid(2, 1, 1, tsg_type_level, tsg_rule_clenshawcurtis)

    call grid%updateGlobalGrid(3, tsg_type_level)
    call approx_grid_pw(grid, grid_ref)

    call grid_ref%release()
    call grid%release()

    grid_ref = TasmanianGlobalGrid(2, 1, 4, tsg_type_level, tsg_rule_fejer2, aweights)
    grid = TasmanianGlobalGrid(2, 1, 1, tsg_type_level, tsg_rule_fejer2)

    call grid%updateGlobalGrid(4, tsg_type_level, aweights)
    call approx_grid_pw(grid, grid_ref)

    call grid_ref%release()
    call grid%release()

    grid_ref = TasmanianGlobalGrid(2, 1, 6, tsg_type_level, tsg_rule_rleja, aweights, level_limits=level_limits)
    grid = TasmanianGlobalGrid(2, 1, 1, tsg_type_level, tsg_rule_rleja)

    call grid%updateGlobalGrid(6, tsg_type_level, aweights, level_limits)
    call approx_grid_pw(grid, grid_ref)

    call grid_ref%release()
    call grid%release()
end subroutine

subroutine test_update_sequence()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid, grid_ref
    integer(C_INT) :: aweights(4), level_limits(2)

    aweights = [ 2, 2, 1, 1 ]
    level_limits = [ 2, 1 ]

    grid_ref = TasmanianSequenceGrid(2, 1, 3, tsg_type_level, tsg_rule_minlebesgue)
    grid = TasmanianSequenceGrid(2, 1, 1, tsg_type_level, tsg_rule_minlebesgue)

    call grid%updateSequenceGrid(3, tsg_type_level)
    call approx_grid_pw(grid, grid_ref)

    call grid_ref%release()
    call grid%release()

    grid_ref = TasmanianSequenceGrid(2, 1, 4, tsg_type_level, tsg_rule_leja, aweights)
    grid = TasmanianSequenceGrid(2, 1, 1, tsg_type_level, tsg_rule_leja)

    call grid%updateSequenceGrid(4, tsg_type_level, aweights)
    call approx_grid_pw(grid, grid_ref)

    call grid_ref%release()
    call grid%release()

    grid_ref = TasmanianSequenceGrid(2, 1, 8, tsg_type_ipcurved, tsg_rule_rleja, aweights, level_limits=level_limits)
    grid = TasmanianSequenceGrid(2, 1, 1, tsg_type_ipcurved, tsg_rule_rleja)

    call grid%updateSequenceGrid(8, tsg_type_ipcurved, aweights, level_limits)
    call approx_grid_pw(grid, grid_ref)

    call grid_ref%release()
    call grid%release()
end subroutine

subroutine test_update_fourier()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid, grid_ref
    integer(C_INT) :: aweights(2), level_limits(2)

    aweights = [ 2, 1 ]
    level_limits = [ 2, 1 ]

    grid_ref = TasmanianFourierGrid(2, 1, 2, tsg_type_iphyperbolic)
    grid = TasmanianFourierGrid(2, 1, 1, tsg_type_iphyperbolic)

    call grid%updateFourierGrid(2, tsg_type_iphyperbolic)
    call approx_grid_pw(grid, grid_ref)

    call grid_ref%release()
    call grid%release()

    grid_ref = TasmanianFourierGrid(2, 1, 4, tsg_type_level, aweights)
    grid = TasmanianFourierGrid(2, 1, 1, tsg_type_level)

    call grid%updateFourierGrid(4, tsg_type_level, aweights)
    call approx_grid_pw(grid, grid_ref)

    call grid_ref%release()
    call grid%release()

    grid_ref = TasmanianFourierGrid(2, 1, 8, tsg_type_iphyperbolic, aweights, level_limits=level_limits)
    grid = TasmanianFourierGrid(2, 1, 1, tsg_type_iphyperbolic)

    call grid%updateFourierGrid(8, tsg_type_iphyperbolic, aweights, level_limits)
    call approx_grid_pw(grid, grid_ref)

    call grid_ref%release()
    call grid%release()
end subroutine

subroutine test_update_grids()
    call test_update_global()
    call test_update_sequence()
    call test_update_fourier()
    write(*,*) "  Performing tests on update grid methods:         PASS"
end subroutine
