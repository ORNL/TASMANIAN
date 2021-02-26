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
    type(TasmanianSparseGrid) :: grid
    integer :: i

write(*,*)
write(*,'(a,i1,a,i1)') 'Testing Tasmanian Fortran 2003-SWIG interface: version ', &
                       grid%getVersionMajor(), ".", grid%getVersionMinor()
write(*,*)

write(*,'(a,a)') '    module version: ', grid%getVersion()
write(*,'(a,a)') '    license:        ', grid%getLicense()
write(*,'(a,a)') '    git hash:       ', grid%getGitCommitHash()
write(*,'(a,a)') '    cxx flags:      ', grid%getCmakeCxxFlags()
write(*,'(a,l)') '    openmp enabled: ', grid%isOpenMPEnabled()
write(*,'(a,l)') '    blas   enabled: ', grid%isAccelerationAvailable(tsg_accel_cpu_blas)
write(*,'(a,l)') '    magma  enabled: ', grid%isAccelerationAvailable(tsg_accel_gpu_magma)
write(*,'(a,l)') '    gpu    enabled: ', grid%isAccelerationAvailable(tsg_accel_gpu_cuda)

if (grid%isAccelerationAvailable(tsg_accel_gpu_cuda)) then
    do i = 1, grid%getNumGPUs()
        write(*,"(a,i1,a,a20,a,i6,a)") "      device ", i-1, ": ", &
            grid%getGPUName(i-1), " with ", grid%getGPUMemory(i-1),"MB of RAM"
    enddo
endif

write(*,*)

call test_make_grid()
call test_domain_transforms()
call test_eval_surrogate()
call test_update_grids()
call test_refinement()

write(*,*)

end program FORSWIGTESTER
