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
program MPITESTER
    use Tasmanian
    use Tasmanian_mpi
    use mpi
    use, intrinsic :: iso_c_binding
    implicit none

    real(C_DOUBLE), dimension(:,:), pointer :: points, values
    integer :: mpi_ierr, me, i

    type(TasmanianSparseGrid) :: grid, scattered_grid, reference_grid

! ======= begin code ======= !
call MPI_Init(mpi_ierr)

call MPI_Comm_rank(MPI_COMM_WORLD, me, mpi_ierr)

if (ierr /= MPI_SUCCESS) then
  write(*,*) "ERROR: failed to initialize MPI"
  stop 1
endif

reference_grid = TasmanianSequenceGrid(3, 1, 3, tsg_type_level, tsg_rule_rleja)
grid = TasmanianSparseGrid()

if (me == 0) then
    mpi_ierr = tsgMPIGridSend(reference_grid, 1, 11, MPI_COMM_WORLD)
endif
if (me == 1) then
    mpi_ierr = tsgMPIGridRecv(grid, 0, 11, MPI_COMM_WORLD)
    if (grid%getNumPoints() /= reference_grid%getNumPoints() .or. &
        grid%getNumDimensions() /= reference_grid%getNumDimensions()) then
        write(*,*) "ERROR: failed to receive the grid"
        stop 1
    endif
endif

if (mpi_ierr /= MPI_SUCCESS) then
    write(*,*) "ERROR: failed at send/recv stage with mpi code: ", mpi_ierr
    error stop
endif

call grid%release()

call MPI_Barrier(MPI_COMM_WORLD, mpi_ierr)

if (me == 0) then
    write(*,*) ""
    write(*,*) " Tasmanian Fortran 2003 MPI Tests:"
    write(*,*) ""
    write(*,*) "    MPI Send/Recv         PASS"
endif

!
! Testing tsgMPIGridBcast()
!

! make a silly grid to overwrite later
grid = TasmanianSequenceGrid(2, 1, 0, tsg_type_level, tsg_rule_rleja)

if (me == 2) then
    mpi_ierr = tsgMPIGridBcast(reference_grid, 2, MPI_COMM_WORLD)
    call grid%copyGrid(reference_grid)
endif
if (me /= 2) then
    mpi_ierr = tsgMPIGridBcast(grid, 2, MPI_COMM_WORLD)
endif

if (grid%getNumPoints() /= reference_grid%getNumPoints() .or. &
    grid%getNumDimensions() /= reference_grid%getNumDimensions()) then
    write(*,*) "ERROR: failed to broadcast the grid"
    error stop
endif

if (mpi_ierr /= MPI_SUCCESS) then
    write(*,*) "ERROR: failed at broadcast stage with mpi code: ", mpi_ierr
    error stop
endif

if (me == 0) then
    write(*,*) "    MPI Bcast             PASS"
endif

call MPI_Barrier(MPI_COMM_WORLD, mpi_ierr)

! Cleanup before scatter
call grid%release()
call reference_grid%release()

! scatter_grid is the destination, grid is the source
scattered_grid = TasmanianSparseGrid()

if (me == 1) then
! make the source grid
    grid = TasmanianLocalPolynomialGrid(3, 2, 4, 1, tsg_rule_localp)
    points => grid%returnPoints()
    allocate(values(2, grid%getNumPoints()))
    do i = 1, grid%getNumPoints()
        values(1, i) = exp(points(1,i) + points(2, i))
        values(2, i) = exp(points(2,i) - points(3, i))
    enddo
    call grid%loadNeededPoints(values(:, 1))
    deallocate(points, values)
else
    grid = TasmanianSparseGrid() ! non-null grid is needed outside of root
endif

reference_grid = TasmanianLocalPolynomialGrid(3, 1, 4, 1, tsg_rule_localp)
if (me == 0) then
! rank 0 gets the first output
    points => reference_grid%returnPoints()
    allocate(values(1, reference_grid%getNumPoints()))
    do i = 1, reference_grid%getNumPoints()
        values(1, i) = exp(points(1,i) + points(2, i))
    enddo
    call reference_grid%loadNeededPoints(values(:, 1))
    deallocate(points, values)
elseif (me == 1) then
! rank 1 gets the second output
    points => reference_grid%returnPoints()
    allocate(values(1, reference_grid%getNumPoints()))
    do i = 1, reference_grid%getNumPoints()
        values(1, i) = exp(points(2,i) - points(3, i))
    enddo
    call reference_grid%loadNeededPoints(values(:, 1))
    deallocate(points, values)
endif
! rank 2 gets and empty grid

mpi_ierr = tsgMPIGridScatterOutputs(grid, scattered_grid, 1, 47, MPI_COMM_WORLD)
if (mpi_ierr /= 0) then
    write(*,*) "ERROR: calling tsgMPIGridScatterOutputs() with code: ", mpi_ierr
    error stop
endif

if (me == 2) then
    if (.not. scattered_grid%isEmpty()) then
        write(*,*) "rank 2 grid is not empty()"
        error stop
    endif
else
    call approx_grid_pv(scattered_grid, reference_grid)
endif

! clean up the grids
call reference_grid%release()
call scattered_grid%release()
call grid%release()

if (me == 0) then
    write(*,*) "    MPI Scatter Outputs   PASS"
    write(*,*) ""
endif

call MPI_Finalize(mpi_ierr)

end program MPITESTER
