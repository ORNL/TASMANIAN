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
program FORTESTER
  USE TasmanianSG
  implicit none

integer,          parameter :: seed = 86456
double precision, parameter :: pi = 4.0D+0 * atan(1.0D+0)

character, pointer :: licence(:)
character(14), parameter :: filename = 'fortranio.test'

type(TasmanianSparseGrid) :: grid, grid_II

double precision, pointer :: rnd(:,:)
double precision, pointer :: points(:,:), weights(:), pointsb(:,:), weightsb(:)


integer          :: i_a, i_b, i_c, i_d, i, j
double precision :: d_a, d_b, d_c, d_d, int, int2
logical          :: bool

integer,          allocatable :: int_1d_a(:),    int_1d_b(:),    int_1d_c(:),    int_1d_d(:)
double precision, allocatable :: double_1d_a(:), double_1d_b(:), double_1d_c(:), double_1d_d(:)
double complex,   allocatable :: dcmplx_1d_a(:), dcmplx_1d_b(:), dcmplx_1d_c(:), dcmplx_1d_d(:)

integer,          allocatable :: int_2d_a(:,:),    int_2d_b(:,:),    int_2d_c(:,:),    int_2d_d(:,:)
double precision, allocatable :: double_2d_a(:,:), double_2d_b(:,:), double_2d_c(:,:), double_2d_d(:,:)
double complex,   allocatable :: dcmplx_2d_a(:,:), dcmplx_2d_b(:,:), dcmplx_2d_c(:,:), dcmplx_2d_d(:,:)

integer,          pointer :: int_pnt_1d_a(:),    int_pnt_1d_b(:),    int_pnt_1d_c(:)
double precision, pointer :: double_pnt_1d_a(:), double_pnt_1d_b(:), double_pnt_1d_c(:)
double complex,   pointer :: dcmplx_pnt_1d_a(:), dcmplx_pnt_1d_b(:), dcmplx_pnt_1d_c(:)



write(*,'(a)') 'Testing TASMANIAN FORTRAN interface'
write(*,*)

! read library meta-data
licence => tsgGetLicense()
write(*,"(A,I2,A,I1)") "Tasmanian Sparse Grid module: ", tsgGetVersionMajor(), ".", tsgGetVersionMinor()
write(*,"(A,40A)") "Licence: ", licence


! initialize random number generator
call srand(seed)


call tsgAllocateGrid(grid)
call tsgAllocateGrid(grid_II)


!=======================================================================
!       tsgIs***()
!=======================================================================
! test type of grid
call tsgMakeGlobalGrid(grid, 2, 0, 1, tsg_level, tsg_clenshaw_curtis)
if ( (.not. tsgIsGlobal(grid)) .or. tsgIsSequence(grid) .or. tsgIsLocalPolynomial(grid) &
     .or. tsgIsFourier(grid) .or. tsgIsWavelet(grid) ) then
  write(*,*) "Mismatch in tsgIsGlobal"
  stop 1
endif
call tsgMakeSequenceGrid(grid, 2, 1, 3, tsg_level, tsg_min_lebesgue)
if ( .not. tsgIsSequence(grid) ) then
  write(*,*) "Mismatch in tsgIsSequence"
  stop 1
endif
call tsgMakeLocalPolynomialGrid(grid, 3, 1, 2, 2, tsg_localp)
if ( (.not. tsgIsLocalPolynomial(grid)) .or. tsgIsSequence(grid) .or. tsgIsGlobal(grid) &
     .or. tsgIsFourier(grid) .or. tsgIsWavelet(grid) ) then
  write(*,*) "Mismatch in tsgIsLocalPolynomial"
  stop 1
endif
call tsgMakeWaveletGrid(grid, 3, 1, 2, 1)
if ( .not. tsgIsWavelet(grid) ) then
  write(*,*) "Mismatch in tsgIsWavelet"
  stop 1
endif
call tsgMakeFourierGrid(grid, 3, 1, 1, tsg_level)
if ( .not. tsgIsFourier(grid) ) then
  write(*,*) "Mismatch in tsgIsFourier"
  stop 1
endif
! test transform
if ( tsgIsSetDomainTransform(grid) ) then
  write(*,*) "Mismatch in tsgIsSetDomainTransform"
  stop 1
endif
call tsgMakeGlobalGrid(grid, 3, 0, 2, tsg_level, tsg_clenshaw_curtis)
call tsgSetDomainTransform(grid, transformA=(/ 3.0d0, -7.0d0, -12.0d0 /), transformB=(/ 5.0d0, -6.0d0,  17.0d0 /))
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( .not. tsgIsSetDomainTransform(grid) ) then
  write(*,*) "Mismatch in tsgIsSetDomainTransform"
  stop 1
endif
deallocate(points, weights)







!=======================================================================
!       tsgMakeGlobalGrid()
!=======================================================================
! Check for correct points and weights
allocate( pointsb(2,5), weightsb(5) )
pointsb  = reshape((/0.0D+0, 0.0D+0, 0.0D+0, -1.0D+0, 0.0D+0, 1.0D+0, -1.0D+0, 0.0D+0, 1.0D+0, 0.0D+0/), shape(pointsb))
weightsb = reshape((/4.0D+0/3.0D+0, 2.0D+0/3.0D+0, 2.0D+0/3.0D+0, 2.0D+0/3.0D+0, 2.0D+0/3.0D+0/), shape(weightsb))
call tsgMakeGlobalGrid(grid, 2, 0, 1, tsg_level, tsg_clenshaw_curtis)
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( norm(pointsb-points,size(points)) > 1.d-11 .OR. norm1d(weightsb-weights) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgMakeGlobal: core case 1", sum(abs(pointsb-points)), sum(abs(weightsb-weights))
  stop 1
end if
deallocate( pointsb, weightsb, points, weights )


call tsgMakeGlobalGrid(grid, 2, 0, 2, tsg_level, tsg_clenshaw_curtis)
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ((abs(points(2, 4) + sqrt(0.5D+0)) > 1.0D-11) .OR. &
    (abs(points(2, 5) - sqrt(0.5D+0)) > 1.0D-11) .OR. &
    (abs(weights(7) - 1.0D+0/9.0D+0)  > 1.0D-11))then
  write(*,*) "Mismatch in tsgMakeGlobal: core case 2", &
             abs(points(2, 4) - sqrt(0.5D+0)), abs(points(2, 5) - sqrt(0.5D+0)), abs(weights(7) - 1.0D+0/9.0D+0)
  stop 1
end if
deallocate(points, weights)


call tsgMakeGlobalGrid(grid, 3, 0, 4, tsg_level, tsg_fejer2)
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( (abs(sum(points)) > 1.0D-11) .OR. (abs(sum(weights)-8.0D+0) > 1.0D-11) )then
  write(*,*) "Mismatch in tsgMakeGlobal: core case 3"
  write(*,*) abs(sum(points)), abs(sum(weights) - 8.0D+0)
  stop 1
end if
deallocate(points, weights)


allocate( pointsb(1,4), weightsb(4) )
pointsb(1,:) = (/0.0D+0, 1.0D+0, -1.0D+0, sqrt(1.0D+0/3.0D+0)/)
weightsb     = (/4.0D+0/3.0D+0, 1.0D+0/3.0D+0, 1.0D+0/3.0D+0, 0.0D+0/)
call tsgMakeGlobalGrid(grid, 1, 0, 3, tsg_level, tsg_leja)
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( (norm2d(pointsb-points) > 1.0D-11) .OR. (norm1d(weightsb-weights) > 1.0D-11) )then
  write(*,*) "Mismatch in tsgMakeGlobal: core case 4", norm2d(pointsb-points), norm1d(weightsb-weights)
  stop 1
end if
deallocate(pointsb, weightsb, points, weights)


! test transform
call tsgMakeGlobalGrid(grid, 3, 0, 2, tsg_level, tsg_clenshaw_curtis)
call tsgSetDomainTransform(grid, transformA=(/ 3.0d0, -7.0d0, -12.0d0 /), transformB=(/ 5.0d0, -6.0d0,  17.0d0 /))
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( (abs(sum(weights)-58.d0)        > 1.d-11) .or.  &
     (abs(maxval(points(1,:))-5.d0)  > 1.d-11) .or. (abs(minval(points(1,:))-3.d0)  > 1.d-11) .or. &
     (abs(maxval(points(2,:))+6.d0)  > 1.d-11) .or. (abs(minval(points(2,:))+7.d0)  > 1.d-11) .or. &
     (abs(maxval(points(3,:))-17.d0) > 1.d-11) .or. (abs(minval(points(3,:))+12.d0) > 1.d-11) ) then
  write(*,*) "Mismatch in tsgMakeGlobal: transform "
  stop 1
endif
deallocate(points, weights)


! test alpha/beta
call tsgMakeGlobalGrid(grid, 1, 0, 4, tsg_level, tsg_gauss_hermite, alpha=2.0D+0)
weights => tsgGetQuadratureWeights(grid)
if ( (abs(sum(weights)-0.5d0*sqrt(pi)) > 1.D-11) .OR. (abs(tsgGetBeta(grid)) > 1.D-11) &
      .OR. (abs(tsgGetAlpha(grid) - 2.0D+0) > 1.D-11) .OR. (tsgGetOrder(grid) /= -1)) then
  write(*,*) "Mismatch in tsgMakeGlobal: alpha/beta", sum(weights), tsgGetBeta(grid)
  stop 1
endif
deallocate(weights)


! test anisotropy
allocate( pointsb(2,4) )
pointsb = reshape((/0.0D+0, 0.0D+0, 0.0D+0, 1.0D+0, 0.0D+0, -1.0D+0, 1.0D+0, 0.0D+0/), shape(pointsb))
call tsgMakeGlobalGrid(grid, 2, 0, 2, tsg_level, tsg_leja, aweights=(/2,1/))
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( (norm2d(pointsb-points) > 1.D-11) .OR. (abs(sum(weights)-4.0D+0) > 1.D-11))then
  write(*,*) "Mismatch in tsgMakeGlobal: anisotropy", norm2d(pointsb-points), sum(weights)
  stop 1
endif
deallocate(pointsb, points, weights)


! make custom rule file
open(unit=17,iostat=i,file='tasmanianFortranTestCustomRule.table',status='replace')
if (i .ne. 0) then
  print*, "Cannot open file with custom rule"
  stop 1
endif
write(17,*) "description: test Gauss-Legendre"
write(17,*) "levels: 5"
write(17,*) "1 1 2 3 3 5 4 7 5 9"
do i = 0,4
  call tsgMakeGlobalGrid(grid, 1, 0, i, tsg_level, tsg_gauss_legendre)
  points  => tsgGetPoints(grid)
  weights => tsgGetQuadratureWeights(grid)
  do j = 1,tsgGetNumPoints(grid)
    write(17,*) weights(j), points(1,j)
  enddo
  deallocate(points,weights)
enddo
close(unit=17)
! test custom rule
call tsgMakeGlobalGrid(grid, 2, 0, 4, tsg_level, tsg_custom_tabulated, &
                       customRuleFilename="tasmanianFortranTestCustomRule.table")
call tsgMakeGlobalGrid(grid_II, 2, 0, 4, tsg_level, tsg_gauss_legendre)
points   => tsgGetPoints(grid)
pointsb  => tsgGetPoints(grid_II)
weights  => tsgGetQuadratureWeights(grid)
weightsb => tsgGetQuadratureWeights(grid_II)
if ( (norm1d(weights-weightsb) > 1.d-11) .or. (norm2d(points-pointsb) > 1.d-11) ) then
  write(*,*) "Mismatch in tsgMakeGlobal: custom rule"
endif
deallocate(points,pointsb,weights,weightsb)
! delete custom rule file
open(unit=17,file='tasmanianFortranTestCustomRule.table',status='replace')
close(unit=17,status='delete')


! test conformal
do i = 7,8
  call tsgMakeGlobalGrid(grid,    2, 0, i, tsg_qptotal, tsg_gauss_patterson)
  call tsgMakeGlobalGrid(grid_II, 2, 0, i, tsg_qptotal, tsg_gauss_patterson)
  call tsgSetConformalTransformASIN(grid_II, (/4,4/))
  points   => tsgGetPoints(grid)
  pointsb  => tsgGetPoints(grid_II)
  weights  => tsgGetQuadratureWeights(grid)
  weightsb => tsgGetQuadratureWeights(grid_II)

  int  = sum( weights  * ( 1.d0/((1.d0+5.d0*points(1,:)**2)*(1.d0+5.d0*points(2,:)**2  )) ) )
  int2 = sum( weightsb * ( 1.d0/((1.d0+5.d0*pointsb(1,:)**2)*(1.d0+5.d0*pointsb(2,:)**2)) ) )

  deallocate(points, pointsb, weights, weightsb)
  if ( abs(int-1.028825601981092d0**2) < abs(int2 - 1.028825601981092d0**2) ) then
    write(*,*) "Mismatch in tsgMakeGlobal: conformal map"
    stop 1
  endif
enddo


! test level limits
call tsgMakeGlobalGrid(grid, 2, 0, 20, tsg_qptotal, tsg_clenshaw_curtis,  (/1,1/), 0.0D+0, 0.0D+0, levelLimits=(/1,3/))
call tsgMakeGlobalGrid(grid_II, 2, 0, 1, tsg_tensor, tsg_clenshaw_curtis, (/1,3/))
points   => tsgGetPoints(grid)
pointsb  => tsgGetPoints(grid_II)
weights  => tsgGetQuadratureWeights(grid)
weightsb => tsgGetQuadratureWeights(grid_II)
if ( (tsgGetNumPoints(grid_II) .ne. tsgGetNumPoints(grid)) .or. &
     (norm2d(points-pointsb) > 1.0D-11) .or. (norm1d(weights-weightsb) > 1.0D-11) ) then
  write(*,*) "Mismatch in tsgMakeGlobal: level limits", tsgGetNumPoints(grid_II), norm2d(points-pointsb)
  stop 1
end if
deallocate(points, pointsb, weights, weightsb)








!=======================================================================
!       tsgCopyGrid()
!=======================================================================
call tsgMakeGlobalGrid(grid, 2, 1, 3, tsg_iptotal, tsg_clenshaw_curtis)
call tsgSetConformalTransformASIN(grid, (/4,4/))
points  => tsgGetPoints(grid)
call tsgCopyGrid(grid_II,grid)
call tsgDeallocateGrid(grid)
pointsb => tsgGetPoints(grid_II)
if ( norm2d(points-pointsb) > 1.0D-11 ) then
  write(*,*) "Mismatch in tsgMakeGlobal: copy grid"
  stop 1
end if
call tsgAllocateGrid(grid)
deallocate(points, pointsb)




!=======================================================================
!       tsgRead/Write()
!=======================================================================
call tsgMakeFourierGrid(grid, 3, 1, 3, tsg_level, (/1, 2, 3/))
points  => tsgGetPoints(grid)
call tsgWrite(grid, filename)
call tsgDeallocateGrid(grid)
call tsgAllocateGrid(grid) ! reset the grid
if (.not. tsgRead(grid, filename)) then
  write(*,*) "File read failed: 1"
  stop 1
end if
pointsb => tsgGetPoints(grid)
if ( norm2d(points-pointsb) > 1.0D-11 ) then
  write(*,*) "Mismatch in tsgMakeGlobal: read/write 1"
  stop 1
end if
deallocate(pointsb)

call tsgWrite(grid, filename, .false.)
call tsgDeallocateGrid(grid)
call tsgAllocateGrid(grid) ! reset the grid
if (.not. tsgRead(grid, filename)) then
  write(*,*) "File read failed: 2"
  stop 1
end if
pointsb => tsgGetPoints(grid)
if ( norm2d(points-pointsb) > 1.0D-11 ) then
  write(*,*) "Mismatch in tsgMakeGlobal: read/write 2"
  stop 1
end if
deallocate(points, pointsb)

call tsgWrite(grid, filename, .true.)
call tsgDeallocateGrid(grid)
call tsgAllocateGrid(grid) ! reset the grid
if (.not. tsgRead(grid, filename)) then
  write(*,*) "File read failed: 3"
  stop 1
end if
if ( tsgGetRule(grid) /= tsg_fourier ) then
  write(*,*) "Mismatch in tsgMakeGlobal: read/write 3"
  stop 1
end if

if (tsgRead(grid, 'nonafile')) then
  write(*,*) "File read failed: 4"
  stop 1
end if




!=======================================================================
!       tsgMakeSequenceGrid()
!=======================================================================
call tsgMakeSequenceGrid(grid, 2, 1, 3, tsg_level, tsg_min_lebesgue)
points => tsgGetPoints(grid)
allocate(pointsb(2,10))
pointsb = reshape((/ 0.D0, 0.D0, 0.D0,  1.D0,  0.D0, -1.D0,  0.D0, sqrt(1.D0/3.D0), 1.D0, 0.D0, &
                     1.D0, 1.D0, 1.D0, -1.D0, -1.D0,  0.D0, -1.D0, 1.D0, sqrt(1.D0/3.D0), 0.D0 /), shape(pointsb))
if ( norm2d(points-pointsb) > 1.0D-11 ) then
  write(*,*) "Mismatch in tsgMakeSequenceGrid: core case 1"
  stop 1
endif
deallocate(points, pointsb)


! test transform
call tsgMakeSequenceGrid(grid, 3, 1, 2, tsg_level, tsg_rleja)
call tsgSetDomainTransform(grid, transformA=(/3.0d0,-7.0d0,-12.0d0/), transformB=(/5.0d0,-6.0d0,17.0d0/))
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( (abs(sum(weights)-58.d0) > 1.d-11) .or.  &
     (abs(maxval(points(1,:))-5.d0)  > 1.d-11) .or. (abs(minval(points(1,:))-3.d0)  > 1.d-11) .or. &
     (abs(maxval(points(2,:))+6.d0)  > 1.d-11) .or. (abs(minval(points(2,:))+7.d0)  > 1.d-11) .or. &
     (abs(maxval(points(3,:))-17.d0) > 1.d-11) .or. (abs(minval(points(3,:))+12.d0) > 1.d-11) ) then
  write(*,*) "Mismatch in tsgMakeSequenceGrid: transform "
  stop 1
endif

allocate( double_1d_a(tsgGetNumPoints(grid)) )
call tsgGetQuadratureWeightsStatic(grid, double_1d_a)
if ( (abs(sum(weights)-58.d0) > 1.d-11) ) then
    write(*,*) "Mismatch in tsgGetQuadratureWeightsStatic"
  stop 1
endif
deallocate( double_1d_a )

allocate(double_1d_a(3), double_1d_b(3))
call tsgGetDomainTransform(grid, double_1d_a, double_1d_b)
if ( (norm1d(double_1d_a - (/ 3.0d0, -7.0d0, -12.0d0 /)) > 1.d-11) .or. &
     (norm1d(double_1d_b - (/ 5.0d0, -6.0d0,  17.0d0 /)) > 1.d-11) ) then
    write(*,*) "Mismatch in tsgGetDomainTransform "
  stop 1
endif
deallocate(points, weights, double_1d_a, double_1d_b)

call tsgClearDomainTransform(grid)
points  => tsgGetPoints(grid)
if ( tsgIsSetDomainTransform(grid) .or.  &
     (abs(maxval(points(1,:))-1.d0) > 1.d-11) .or. (abs(minval(points(1,:))+1.d0) > 1.d-11) .or. &
     (abs(maxval(points(2,:))-1.d0) > 1.d-11) .or. (abs(minval(points(2,:))+1.d0) > 1.d-11) .or. &
     (abs(maxval(points(3,:))-1.d0) > 1.d-11) .or. (abs(minval(points(3,:))+1.d0) > 1.d-11) ) then
  write(*,*) "Mismatch in tsgMakeSequenceGrid: cleared transform "
  stop 1
endif
deallocate(points)


! test anisotropy
allocate(pointsb(2,4))
pointsb = reshape((/0.0D+0, 0.0D+0, 0.0D+0, 1.0D+0, 0.0D+0, -1.0D+0, 1.0D+0, 0.0D+0/), shape(pointsb))
call tsgMakeSequenceGrid(grid, 2, 1, 2, tsg_level, tsg_leja, aweights=(/2,1/))
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( (norm2d(pointsb-points) > 1.D-11) .or. (abs(sum(weights)-4.0D+0) > 1.D-11) ) then
  write(*,*) "Mismatch in tsgMakeSequenceGrid: anisotropy", norm2d(pointsb-points), sum(weights)
  stop 1
endif
deallocate(pointsb, points, weights)


! test conformal
allocate( pointsb(2,100), double_2d_a(1,100), double_2d_b(1,100), double_2d_c(1,100) )
rnd => random(2,100)
pointsb          = -1.d0 + 2.d0 * rnd
double_2d_a(1,:) = 1.d0/((1.d0+5.d0*pointsb(1,:)**2)*(1.d0+5.d0*pointsb(2,:)**2))
do i = 7,8
  call tsgMakeSequenceGrid(grid, 2, 1, i, tsg_iptotal, tsg_rleja)
  points => tsgGetPoints(grid)
  allocate(double_2d_d(1, tsgGetNumPoints(grid)))
  double_2d_d(1,:) = 1.d0/((1.d0+5.d0*points(1,:)**2)*(1.d0+5.d0*points(2,:)**2))
  call tsgLoadNeededPoints(grid, double_2d_d)
  call tsgEvaluateBatch(grid,pointsb,100,double_2d_b)
  deallocate(points, double_2d_d)

  call tsgMakeSequenceGrid(grid_II, 2, 1, i, tsg_iptotal, tsg_rleja)
  call tsgSetConformalTransformASIN(grid_II, (/4,4/))
  points => tsgGetPoints(grid_II)
  allocate(double_2d_d(1, tsgGetNumPoints(grid_II)))
  double_2d_d(1,:) = 1.d0/((1.d0+5.d0*points(1,:)**2)*(1.d0+5.d0*points(2,:)**2))
  call tsgLoadNeededPoints(grid_II, double_2d_d)
  call tsgEvaluateBatch(grid_II,pointsb,100,double_2d_c)
  deallocate(points, double_2d_d)

  if ( norm2d(double_2d_a-double_2d_b) < norm2d(double_2d_a-double_2d_c) ) then
    write(*,*) "Mismatch in tsgMakeSequenceGrid: conformal map"
    stop 1
  endif
enddo
deallocate(rnd, pointsb, double_2d_a, double_2d_b, double_2d_c)


! test level limits
call tsgMakeSequenceGrid(grid, 2, 1, 20, tsg_iptotal, tsg_min_delta, levelLimits=(/1,3/))
call tsgMakeSequenceGrid(grid_II, 2, 1, 1, tsg_tensor, tsg_min_delta, aweights=(/1,3/))
points   => tsgGetPoints(grid)
pointsb  => tsgGetPoints(grid_II)
weights  => tsgGetQuadratureWeights(grid)
weightsb => tsgGetQuadratureWeights(grid_II)
if ( (tsgGetNumPoints(grid_II) .ne. tsgGetNumPoints(grid)) .or. &
     (norm2d(points-pointsb) > 1.0D-11)  .or. (norm1d(weights-weightsb) > 1.0D-11) ) then
  write(*,*) "Mismatch in tsgMakeSequenceGrid: level limits", tsgGetNumPoints(grid_II), norm2d(points-pointsb)
  stop 1
end if
deallocate(points, pointsb, weights, weightsb)








!=======================================================================
!       tsgMakeLocalPolynomialGrid()
!=======================================================================
! test transform
call tsgMakeLocalPolynomialGrid(grid, 3, 1, 2, 2, tsg_localp)
call tsgSetDomainTransform(grid, transformA=(/3.0d0,-7.0d0,-12.0d0/), transformB=(/5.0d0,-6.0d0,17.0d0/))
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( (abs(sum(weights)-58.d0)        > 1.d-11) .or.  &
     (abs(maxval(points(1,:))-5.d0)  > 1.d-11) .or. (abs(minval(points(1,:))-3.d0)  > 1.d-11) .or. &
     (abs(maxval(points(2,:))+6.d0)  > 1.d-11) .or. (abs(minval(points(2,:))+7.d0)  > 1.d-11) .or. &
     (abs(maxval(points(3,:))-17.d0) > 1.d-11) .or. (abs(minval(points(3,:))+12.d0) > 1.d-11) ) then
  write(*,*) "Mismatch in tsgMakeLocalPolynomialGrid: transform "
  stop 1
endif
deallocate(points, weights)


! test conformal
allocate( pointsb(2,100), double_2d_a(1,100), double_2d_b(1,100), double_2d_c(1,100) )
rnd => random(2,100)
pointsb          = -1.d0 + 2.d0 * rnd
double_2d_a(1,:) = 1.d0/((1.d0+5.d0*pointsb(1,:)**2)*(1.d0+5.d0*pointsb(2,:)**2))
do i = 3,4
  call tsgMakeLocalPolynomialGrid(grid, 2, 1, i, 2, tsg_semi_localp)
  points => tsgGetPoints(grid)
  allocate(double_2d_d(1, tsgGetNumPoints(grid)))
  double_2d_d(1,:) = 1.d0/((1.d0+5.d0*points(1,:)**2)*(1.d0+5.d0*points(2,:)**2))
  call tsgLoadNeededPoints(grid, double_2d_d)
  call tsgEvaluateBatch(grid,pointsb,100,double_2d_b)
  deallocate(points, double_2d_d)

  call tsgMakeLocalPolynomialGrid(grid_II, 2, 1, i, 2, tsg_semi_localp)
  call tsgSetConformalTransformASIN(grid_II, (/4,4/))
  points => tsgGetPoints(grid_II)
  allocate(double_2d_d(1, tsgGetNumPoints(grid_II)))
  double_2d_d(1,:) = 1.d0/((1.d0+5.d0*points(1,:)**2)*(1.d0+5.d0*points(2,:)**2))
  call tsgLoadNeededPoints(grid_II, double_2d_d)
  call tsgEvaluateBatch(grid_II,pointsb,100,double_2d_c)
  deallocate(points, double_2d_d)

  if ( norm2d(double_2d_a-double_2d_b) < norm2d(double_2d_a-double_2d_c) ) then
    write(*,*) "Mismatch in tsgMakeSequenceGrid: conformal map"
    stop 1
  endif
enddo
deallocate(rnd, pointsb, double_2d_a, double_2d_b, double_2d_c)


! test level limits
call tsgMakeLocalPolynomialGrid(grid, 3, 1, 3, 2, tsg_semi_localp, levelLimits=(/1,2,3/))
points => tsgGetPoints(grid)
if ( minval(abs(points(1,:)-0.5)) < 1.e-8 ) then
  write(*,*) "Mismatch in tsgMakeLocalPolynomialGrid: level limits, dim 1"
  stop 1
endif
if ( minval(abs(points(2,:)-0.75)) < 1.e-8 ) then
  write(*,*) "Mismatch in tsgMakeLocalPolynomialGrid: level limits, dim 2"
  stop 1
endif
if ( minval(abs(points(3,:)-0.125)) < 1.e-8 ) then
  write(*,*) "Mismatch in tsgMakeLocalPolynomialGrid: level limits, dim 3"
  stop 1
endif
deallocate(points)








!=======================================================================
!       tsgMakeWaveletGrid()
!=======================================================================
! test transform
call tsgMakeWaveletGrid(grid, 3, 1, 2, 1)
call tsgSetDomainTransform(grid, transformA=(/3.0d0,-7.0d0,-12.0d0/), transformB=(/5.0d0,-6.0d0,17.0d0/))
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( (abs(sum(weights)-58.d0)        > 1.d-11) .or.  &
     (abs(maxval(points(1,:))-5.d0)  > 1.d-11) .or. (abs(minval(points(1,:))-3.d0)  > 1.d-11) .or. &
     (abs(maxval(points(2,:))+6.d0)  > 1.d-11) .or. (abs(minval(points(2,:))+7.d0)  > 1.d-11) .or. &
     (abs(maxval(points(3,:))-17.d0) > 1.d-11) .or. (abs(minval(points(3,:))+12.d0) > 1.d-11) ) then
  write(*,*) "Mismatch in tsgMakeWaveletGrid: transform "
  stop 1
endif
deallocate(points, weights)


! correctness of 1-D
call tsgMakeWaveletGrid(grid, 1, 1, 2, 1)
call tsgMakeLocalPolynomialGrid(grid_II, 1, 1, 3, 1, tsg_localp)
points  => tsgGetPoints(grid)
pointsb => tsgGetPoints(grid_II)
if ( norm2d(points-pointsb) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgMakeWaveletGrid: points "
  stop 1
endif
deallocate(points, pointsb)


! test conformal
do i = 3,4
  call tsgMakeWaveletGrid(grid,    2, 1, i, 1)
  call tsgMakeWaveletGrid(grid_II, 2, 1, i, 1)
  call tsgSetConformalTransformASIN(grid_II, (/4,4/))
  points   => tsgGetPoints(grid)
  pointsb  => tsgGetPoints(grid_II)
  weights  => tsgGetQuadratureWeights(grid)
  weightsb => tsgGetQuadratureWeights(grid_II)

  int  = sum( weights  * ( 1.d0/((1.d0+5.d0*points(1,:)**2)*(1.d0+5.d0*points(2,:)**2)) ) )
  int2 = sum( weightsb * ( 1.d0/((1.d0+5.d0*pointsb(1,:)**2)*(1.d0+5.d0*pointsb(2,:)**2)) ) )

  if ( abs(int-1.028825601981092d0**2) < abs(int2 - 1.028825601981092d0**2) ) then
    write(*,*) "Mismatch in tsgMakeWaveletGrid: conformal map"
    stop 1
  endif
  deallocate(points, pointsb, weights, weightsb)
enddo


! test level limits
call tsgMakeWaveletGrid(grid, 3, 1, 2, 1, levelLimits=(/0,1,2/))
points => tsgGetPoints(grid)
if ( minval(abs(points(1,:)-0.5)) < 1.e-8 ) then
  write(*,*) "Mismatch in tsgMakeWaveletGrid: level limits, dim 1"
  stop 1
endif
if ( minval(abs(points(2,:)-0.75)) < 1.e-8 ) then
  write(*,*) "Mismatch in tsgMakeWaveletGrid: level limits, dim 2"
  stop 1
endif
if ( minval(abs(points(3,:)-0.125)) < 1.e-8 ) then
  write(*,*) "Mismatch in tsgMakeWaveletGrid: level limits, dim 3"
  stop 1
endif
deallocate(points)








!=======================================================================
!       tsgMakeFourierGrid()
!=======================================================================
! test transform
call tsgMakeFourierGrid(grid, 3, 1, 1, tsg_level)
call tsgSetDomainTransform(grid, transformA=(/3.0d0,-7.0d0,-12.0d0/), transformB=(/5.0d0,-6.0d0,17.0d0/))
points  => tsgGetPoints(grid)
weights => tsgGetQuadratureWeights(grid)
if ( (abs(sum(weights)-58.d0)           > 1.d-11) .or.  &
     (abs(maxval(points(1,:))-13./3.d0) > 1.d-11) .or. (abs(minval(points(1,:))-3.d0)  > 1.d-11) .or. &
     (abs(maxval(points(2,:))+19./3.d0) > 1.d-11) .or. (abs(minval(points(2,:))+7.d0)  > 1.d-11) .or. &
     (abs(maxval(points(3,:))-22./3.d0) > 1.d-11) .or. (abs(minval(points(3,:))+12.d0) > 1.d-11) ) then
  write(*,*) "Mismatch in tsgMakeFourierGrid: transform "
  stop 1
endif
deallocate(points, weights)


! correctness of 1-D
call tsgMakeFourierGrid(grid, 1, 1, 2, tsg_level)
allocate(pointsb(1,9))
pointsb = reshape( (/ 0.d0, 1.0/3.d0, 2.0/3.d0, 1.0/9.d0, 2.0/9.d0, &
                            4.0/9.d0, 5.0/9.d0, 7.0/9.d0, 8.0/9.d0 /), (/1, 9/) )
points => tsgGetPoints(grid)
if ( norm2d(points-pointsb) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgMakeFourierGrid: points "
  stop 1
endif
deallocate(points, pointsb)


! test level limits
call tsgMakeFourierGrid(grid, 3, 1, 3, tsg_level, levelLimits=(/0,1,2/))
points => tsgGetPoints(grid)
if ( maxval(abs(points(1,:))) > 1.e-8 ) then
  write(*,*) "Mismatch in tsgMakeFourierGrid: level limits, dim 1"
  stop 1
endif
if ( minval(abs(points(2,:)-1./9.d0)) < 1.e-8 ) then
  write(*,*) "Mismatch in tsgMakeFourierGrid: level limits, dim 2"
  stop 1
endif
if ( minval(abs(points(3,:)-1./27.d0)) < 1.e-8 ) then
  write(*,*) "Mismatch in tsgMakeFourierGrid: level limits, dim 3"
  stop 1
endif
deallocate(points)



write(*,*) "tsgMake* functions:       PASS"




!=======================================================================
!       tsgGet***Points()
!=======================================================================
allocate(pointsb(2,5))
pointsb = reshape( (/ 0.d0, 0.d0, 0.d0, -1.d0, 0.d0, &
                      1.d0, -1.d0, 0.d0, 1.d0, 0.d0 /), (/ 2, 5 /) )
call tsgMakeGlobalGrid(grid, 2, 1, 1, tsg_level, tsg_clenshaw_curtis)

points => tsgGetPoints(grid)
if ( norm2d(points-pointsb) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetPoints: core case 1"
  stop 1
endif
deallocate(points)

points => tsgGetNeededPoints(grid)
if ( norm2d(points-pointsb) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetNeededPoints: core case 1"
  stop 1
endif
deallocate(points)

allocate(points(tsgGetNumDimensions(grid),tsgGetNumPoints(grid)))
call tsgGetPointsStatic(grid, points )
if ( norm2d(points-pointsb) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetPointsStatic: core case 1"
  stop 1
endif
deallocate(points)

allocate(points(tsgGetNumDimensions(grid),tsgGetNumPoints(grid)))
call tsgGetNeededPointsStatic(grid, points )
if ( norm2d(points-pointsb) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetNeededPointsStatic: core case 1"
  stop 1
endif
deallocate(points, pointsb)








!=======================================================================
!       tsgLoadNeededPoints(), tsgGetLoadedPoints()
!=======================================================================
call tsgMakeGlobalGrid(grid, 2, 2, 4, tsg_level, tsg_min_delta)
points  => tsgGetPoints(grid)
pointsb => tsgGetNeededPoints(grid)
if ( norm2d(points-pointsb) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgLoadNeededPoints: tsgGetNeededPoints"
  stop 1
endif
deallocate(pointsb)
allocate(double_2d_a(2,tsgGetNumPoints(grid)))
double_2d_a(1,:) = exp( -points(1,:)**2 - points(2,:)**2 )
double_2d_a(2,:) = cos( -points(1,:) - 2.d0 * points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
pointsb => tsgGetLoadedPoints(grid)
if ( norm2d(points-pointsb) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgLoadNeededPoints: tsgGetLoadedPoints"
  stop 1
endif
deallocate(pointsb)
allocate(pointsb(tsgGetNumDimensions(grid),tsgGetNumPoints(grid)))
call tsgGetLoadedPointsStatic(grid, pointsb)
if ( norm2d(points-pointsb) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgLoadNeededPoints: tsgGetLoadedPointsStatic"
  stop 1
endif
if ( tsgGetNumNeeded(grid) .ne. 0 ) then
  write(*,*) "Mismatch in tsgLoadNeededPoints: tsgGetNumNeeded"
  stop 1
endif
deallocate(points, pointsb, double_2d_a)








!=======================================================================
!       tsgEvaluate(), tsgEvaluateBatch()
!=======================================================================
allocate(double_2d_a(4,2), double_2d_c(4,2), double_1d_a(4), double_1d_b(4), pointsb(2,2))
pointsb = reshape( (/ 1./3.d0, 1./3.d0, pi/6.d0, -sqrt(2.d0)/2.d0 /), shape(pointsb) )
double_2d_a(1,:) = 0.3d0
double_2d_a(2,:) = pointsb(1,:)    + pointsb(2,:)
double_2d_a(3,:) = pointsb(1,:)**2 + pointsb(2,:)**2 + pointsb(1,:)*pointsb(2,:)
double_2d_a(4,:) = pointsb(1,:)**3 + pointsb(2,:)**3 + pointsb(1,:)*(pointsb(2,:)**2)
call tsgMakeLocalPolynomialGrid(grid, 2, 4, 1, 1, tsg_localp)
points => tsgGetPoints(grid)
allocate(double_2d_b(4,tsgGetNumPoints(grid)))
double_2d_b(1,:) = 0.3d0
double_2d_b(2,:) = points(1,:)    + points(2,:)
double_2d_b(3,:) = points(1,:)**2 + points(2,:)**2 + points(1,:)*points(2,:)
double_2d_b(4,:) = points(1,:)**3 + points(2,:)**3 + points(1,:)*(points(2,:)**2)
call tsgLoadNeededPoints(grid, double_2d_b)
call tsgEvaluate(grid,reshape(pointsb(:,1),(/2/)),double_1d_a)
call tsgEvaluateFast(grid,reshape(pointsb(:,2),(/2/)),double_1d_b)
call tsgEvaluateBatch(grid,pointsb,2,double_2d_c)
do i=1,2
  if ( ( norm1d(double_2d_a(i,:)-double_2d_c(i,:)) > 1.d-11 ) .or. &
       ( abs(double_2d_a(i,1)-double_1d_a(i)) + abs(double_2d_a(i,2)-double_1d_b(i)) ) > 1.d-11 ) then
    write(*,*) "Mismatch in tsgEvaluateBatch: case 1, output ", i
    stop 1
  endif
enddo
do i=3,4
  if ( ( norm1d(double_2d_a(i,:)-double_2d_c(i,:)) < 1.d-8 ) .or. &
       ( abs(double_2d_a(i,1)-double_1d_a(i)) + abs(double_2d_a(i,2)-double_1d_b(i)) ) < 1.d-8 ) then
    write(*,*) "Mismatch in tsgEvaluateBatch: case 1, output ", i
    stop 1
  endif
enddo
deallocate(double_2d_a, double_2d_b, double_2d_c, double_1d_a, double_1d_b, pointsb, points)


allocate(double_2d_a(4,2), double_2d_c(4,2), double_1d_a(4), double_1d_b(4), pointsb(2,2))
pointsb = reshape( (/ 1./3.d0, 1./3.d0, pi/6.d0, -sqrt(2.d0)/2.d0 /), shape(pointsb) )
double_2d_a(1,:) = 0.3d0
double_2d_a(2,:) = pointsb(1,:)    + pointsb(2,:)
double_2d_a(3,:) = pointsb(1,:)**2 + pointsb(2,:)**2 + pointsb(1,:)*pointsb(2,:)
double_2d_a(4,:) = pointsb(1,:)**3 + pointsb(2,:)**3 + pointsb(1,:)*(pointsb(2,:)**2)
call tsgMakeLocalPolynomialGrid(grid, 2, 4, 1, 2, tsg_localp)
points => tsgGetPoints(grid)
allocate(double_2d_b(4,tsgGetNumPoints(grid)))
double_2d_b(1,:) = 0.3d0
double_2d_b(2,:) = points(1,:)    + points(2,:)
double_2d_b(3,:) = points(1,:)**2 + points(2,:)**2 + points(1,:)*points(2,:)
double_2d_b(4,:) = points(1,:)**3 + points(2,:)**3 + points(1,:)*(points(2,:)**2)
call tsgLoadNeededPoints(grid, double_2d_b)
call tsgEvaluate(grid,reshape(pointsb(:,1),(/2/)),double_1d_a)
call tsgEvaluate(grid,reshape(pointsb(:,2),(/2/)),double_1d_b)
call tsgEvaluateBatch(grid,pointsb,2,double_2d_c)
do i=1,2
  if ( norm1d(double_2d_a(i,:)-double_2d_c(i,:)) > 1.d-11 .or. &
     ( abs(double_2d_a(i,1)-double_1d_a(i)) + abs(double_2d_a(i,2)-double_1d_b(i)) ) > 1.d-11 ) then
    write(*,*) "Mismatch in tsgEvaluateBatch: case 2, output ", i
    stop 1
  endif
enddo
do i=3,4
  if ( ( norm1d(double_2d_a(i,:)-double_2d_c(i,:)) < 1.d-8 ) .or. &
       ( abs(double_2d_a(i,1)-double_1d_a(i)) + abs(double_2d_a(i,2)-double_1d_b(i)) ) < 1.d-8 ) then
    write(*,*) "Mismatch in tsgEvaluateBatch: case 2, output ", i
    stop 1
  endif
enddo
deallocate(double_2d_a, double_2d_b, double_2d_c, double_1d_a, double_1d_b, pointsb, points)


allocate(double_2d_a(4,2), double_2d_c(4,2), double_1d_a(4), double_1d_b(4), pointsb(2,2))
pointsb = reshape( (/ 1./3.d0, 1./3.d0, pi/6.d0, -sqrt(2.d0)/2.d0 /), shape(pointsb) )
double_2d_a(1,:) = 0.3d0
double_2d_a(2,:) = pointsb(1,:)    + pointsb(2,:)
double_2d_a(3,:) = pointsb(1,:)**2 + pointsb(2,:)**2
double_2d_a(4,:) = pointsb(1,:)**3 + pointsb(2,:)**3 + pointsb(1,:)*(pointsb(2,:)**2)
call tsgMakeLocalPolynomialGrid(grid, 2, 4, 1, 2, tsg_semi_localp)
points => tsgGetPoints(grid)
allocate(double_2d_b(4,tsgGetNumPoints(grid)))
double_2d_b(1,:) = 0.3d0
double_2d_b(2,:) = points(1,:)    + points(2,:)
double_2d_b(3,:) = points(1,:)**2 + points(2,:)**2
double_2d_b(4,:) = points(1,:)**3 + points(2,:)**3 + points(1,:)*(points(2,:)**2)
call tsgLoadNeededPoints(grid, double_2d_b)
call tsgEvaluate(grid,reshape(pointsb(:,1),(/2/)),double_1d_a)
call tsgEvaluate(grid,reshape(pointsb(:,2),(/2/)),double_1d_b)
call tsgEvaluateBatch(grid,pointsb,2,double_2d_c)
do i=1,3
  if ( norm1d(double_2d_a(i,:)-double_2d_c(i,:)) > 1.d-11 .or. &
      ( abs(double_2d_a(i,1)-double_1d_a(i)) + abs(double_2d_a(i,2)-double_1d_b(i)) ) > 1.d-11 ) then
    write(*,*) "Mismatch in tsgEvaluateBatch: case 3, output ", i
    stop 1
  endif
enddo
do i=4,4
  if ( ( norm1d(double_2d_a(i,:)-double_2d_c(i,:)) < 1.d-8 ) .or. &
       ( abs(double_2d_a(i,1)-double_1d_a(i)) + abs(double_2d_a(i,2)-double_1d_b(i)) ) < 1.d-8 ) then
    write(*,*) "Mismatch in tsgEvaluateBatch: case 3, output ", i
    stop 1
  endif
enddo
deallocate(double_2d_a, double_2d_b, double_2d_c, double_1d_a, double_1d_b, pointsb, points)


allocate(double_2d_a(4,2), double_2d_c(4,2), double_1d_a(4), double_1d_b(4), pointsb(2,2))
pointsb = reshape( (/ 1./3.d0, 1./3.d0, pi/6.d0, -sqrt(2.d0)/2.d0 /), shape(pointsb) )
double_2d_a(1,:) = 0.3d0
double_2d_a(2,:) = pointsb(1,:)    + pointsb(2,:)
double_2d_a(3,:) = pointsb(1,:)**2 + pointsb(2,:)**2 + pointsb(1,:)*pointsb(2,:)
double_2d_a(4,:) = pointsb(1,:)**3 + pointsb(2,:)**3 + pointsb(1,:)*(pointsb(2,:)**2)
call tsgMakeLocalPolynomialGrid(grid, 2, 4, 2, 2, tsg_localp)
points => tsgGetPoints(grid)
allocate(double_2d_b(4,tsgGetNumPoints(grid)))
double_2d_b(1,:) = 0.3d0
double_2d_b(2,:) = points(1,:)    + points(2,:)
double_2d_b(3,:) = points(1,:)**2 + points(2,:)**2 + points(1,:)*points(2,:)
double_2d_b(4,:) = points(1,:)**3 + points(2,:)**3 + points(1,:)*(points(2,:)**2)
call tsgLoadNeededPoints(grid, double_2d_b)
call tsgEvaluate(grid,reshape(pointsb(:,1),(/2/)),double_1d_a)
call tsgEvaluate(grid,reshape(pointsb(:,2),(/2/)),double_1d_b)
call tsgEvaluateBatch(grid,pointsb,2,double_2d_c)
do i=1,3
  if ( norm1d(double_2d_a(i,:)-double_2d_c(i,:)) > 1.d-11 .or. &
      ( abs(double_2d_a(i,1)-double_1d_a(i)) + abs(double_2d_a(i,2)-double_1d_b(i)) ) > 1.d-11 ) then
    write(*,*) "Mismatch in tsgEvaluateBatch: case 4, output ", i
    stop 1
  endif
enddo
do i=4,4
  if ( ( norm1d(double_2d_a(i,:)-double_2d_c(i,:)) < 1.d-8 ) .or. &
       ( abs(double_2d_a(i,1)-double_1d_a(i)) + abs(double_2d_a(i,2)-double_1d_b(i)) ) < 1.d-8 ) then
    write(*,*) "Mismatch in tsgEvaluateBatch: case 4, output ", i
    stop 1
  endif
enddo
deallocate(double_2d_a, double_2d_b, double_2d_c, double_1d_a, double_1d_b, pointsb, points)


allocate(double_2d_a(4,2), double_2d_c(4,2), double_1d_a(4), double_1d_b(4), pointsb(2,2))
pointsb = reshape( (/ 1./3.d0, 1./3.d0, pi/6.d0, -sqrt(2.d0)/2.d0 /), shape(pointsb) )
double_2d_a(1,:) = 0.3d0
double_2d_a(2,:) = pointsb(1,:)    + pointsb(2,:)
double_2d_a(3,:) = pointsb(1,:)**2 + pointsb(2,:)**2 + pointsb(1,:)*pointsb(2,:)
double_2d_a(4,:) = pointsb(1,:)**3 + pointsb(2,:)**3 + pointsb(1,:)*(pointsb(2,:)**2)
call tsgMakeLocalPolynomialGrid(grid, 2, 4, 3, 3, tsg_localp)
points => tsgGetPoints(grid)
allocate(double_2d_b(4,tsgGetNumPoints(grid)))
double_2d_b(1,:) = 0.3d0
double_2d_b(2,:) = points(1,:)    + points(2,:)
double_2d_b(3,:) = points(1,:)**2 + points(2,:)**2 + points(1,:)*points(2,:)
double_2d_b(4,:) = points(1,:)**3 + points(2,:)**3 + points(1,:)*(points(2,:)**2)
call tsgLoadNeededPoints(grid, double_2d_b)
call tsgEvaluate(grid,reshape(pointsb(:,1),(/2/)),double_1d_a)
call tsgEvaluate(grid,reshape(pointsb(:,2),(/2/)),double_1d_b)
call tsgEvaluateBatch(grid,pointsb,2,double_2d_c)
do i=1,4
  if ( ( norm1d(double_2d_a(i,:)-double_2d_c(i,:)) > 1.d-11 ) .or. &
       ( abs(double_2d_a(i,1)-double_1d_a(i)) + abs(double_2d_a(i,2)-double_1d_b(i)) ) > 1.d-11 ) then
    write(*,*) "Mismatch in tsgEvaluateBatch: case 5, output ", i
    stop 1
  endif
enddo
deallocate(double_2d_a, double_2d_b, double_2d_c, double_1d_a, double_1d_b, pointsb, points)


call tsgMakeGlobalGrid(grid, 2, 1, 22, tsg_iptotal, tsg_chebyshev)
points => tsgGetPoints(grid)
allocate(double_2d_a(1,tsgGetNumPoints(grid)))
double_2d_a(1,:) = exp( -points(1,:)**2 - points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
allocate(pointsb(2,1000),double_2d_b(1,1000),double_2d_c(1,1000))
rnd     => random(2,1000)
pointsb = -1.d0 + 2.d0 * rnd
double_2d_b(1,:) = exp( -pointsb(1,:)**2 - pointsb(2,:)**2 )
call tsgEvaluateBatch(grid,pointsb,1000,double_2d_c)
if ( norm2d(double_2d_b-double_2d_c) > 1.d-8 ) then
  write(*,*) "Mismatch in tsgEvaluateBatch: global grid with chebyshev points, output ", norm2d(double_2d_b-double_2d_c)
  stop 1
endif
deallocate(points,pointsb,double_2d_a,double_2d_b,double_2d_c,rnd)








!=======================================================================
!       tsgEvaluateHierarchicalFunctions()
!=======================================================================
call tsgMakeGlobalGrid(grid, 3, 1, 4, tsg_level, tsg_fejer2)
points => tsgGetPoints(grid); i_a = tsgGetNumPoints(grid)
allocate(double_2d_a(i_a,i_a))
call tsgEvaluateHierarchicalFunctions(grid,points,i_a,double_2d_a)
forall(i = 1:i_a) double_2d_a(i,i) = double_2d_a(i,i) - 1
if ( norm2d(double_2d_a) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgEvaluateHierarchy: lagrange polynomials do not form identity"
  stop 1
endif
deallocate(points, double_2d_a)


call tsgMakeSequenceGrid(grid, 2, 1, 2, tsg_level, tsg_leja)
allocate(points(2,6), double_2d_a(6,6), double_2d_b(6,6))
points = reshape( (/ 0.33d0, 0.25d0, -0.270d0, 0.39d0,  0.970d0, -0.760d0, &
                    -0.44d0, 0.21d0, -0.813d0, 0.03d0, -0.666d0,  0.666d0 /), shape(points) )
double_2d_a(1,:) = 1
double_2d_a(2,:) = points(2,:)
double_2d_a(3,:) = 0.5d0 * points(2,:) * ( points(2,:) - 1.d0 )
double_2d_a(4,:) = points(1,:)
double_2d_a(5,:) = points(1,:) * points(2,:)
double_2d_a(6,:) = 0.5d0 * points(1,:) * ( points(1,:) - 1.d0 )
call tsgEvaluateHierarchicalFunctions(grid,points,6,double_2d_b)
if ( norm2d(double_2d_a-double_2d_b) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgEvaluateHierarchicalFunctions: sequence grid test"
  stop 1
endif
deallocate(points, double_2d_a, double_2d_b)


call tsgMakeLocalPolynomialGrid(grid, 2, 1, 4, 1, tsg_localp)
points => tsgGetPoints(grid)
i_a = tsgGetNumPoints(grid)
i_b = 13
allocate(double_2d_a(1,i_a), double_2d_b(1,i_b), double_2d_c(i_a,i_b), double_1d_b(i_b), pointsb(2,i_b))
double_2d_a(1,:) = exp( -points(1,:)**2 - 2.d0 * points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
rnd     => random(2,i_b)
pointsb = -1.d0 + 2.d0 * rnd
call tsgEvaluateBatch(grid,pointsb,i_b,double_2d_b)
call tsgEvaluateHierarchicalFunctions(grid,pointsb,i_b,double_2d_c)
weights => tsgGetHierarchicalCoefficients(grid)
double_1d_b = matmul( weights, double_2d_c )
if ( norm1d(double_1d_b-double_2d_b(1,:)) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgEvaluateHierarchicalFunctions: local grid test"
  stop 1
endif
allocate( double_1d_a(tsgGetNumPoints(grid)) )
call tsgGetHierarchicalCoefficientsStatic(grid, double_1d_a)
double_1d_b = matmul( double_1d_a, double_2d_c )
if ( norm1d(double_1d_b-double_2d_b(1,:)) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetHierarchicalCoefficientsStatic: local grid test"
  stop 1
endif
deallocate(points,double_1d_a,double_2d_a,double_2d_b,double_2d_c,weights,double_1d_b,rnd,pointsb)


call tsgMakeLocalPolynomialGrid(grid, 2, 1, 4, 1, tsg_localp)
points => tsgGetPoints(grid)
i_a = tsgGetNumPoints(grid)
i_b = 13
allocate(double_2d_a(1,i_a), double_2d_b(1,i_b), double_1d_b(i_b), pointsb(2,i_b))
double_2d_a(1,:) = exp( -points(1,:)**2 - 2.d0 * points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
rnd     => random(2,i_b)
pointsb = -1.d0 + 2.d0 * rnd
call tsgEvaluateBatch(grid,pointsb,i_b,double_2d_b)
call tsgEvaluateSparseHierarchicalFunctions(grid,pointsb, i_b, int_pnt_1d_a, int_pnt_1d_b, double_pnt_1d_a)
weights => tsgGetHierarchicalCoefficients(grid)
double_1d_b = sparse_matmul( int_pnt_1d_a, int_pnt_1d_b, double_pnt_1d_a, weights )
if ( norm1d(double_1d_b-double_2d_b(1,:)) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgEvaluateSparseHierarchicalFunctions: local grid test"
  stop 1
endif
deallocate(points,double_2d_a,double_2d_b,weights,double_1d_b,rnd,pointsb,int_pnt_1d_a,int_pnt_1d_b,double_pnt_1d_a)


! test reading a complex matrix
call tsgMakeFourierGrid(grid, 2, 1, 4, tsg_level)
points => tsgGetPoints(grid)
i_a = tsgGetNumPoints(grid)
i_b = 13
allocate(double_2d_a(1,i_a), double_2d_b(1,i_b), double_2d_c(i_b,2*i_a), &
         double_1d_a(i_b),   double_1d_b(i_b),   double_1d_c(i_b), pointsb(2,i_b), double_1d_d(2*i_a))
double_2d_a(1,:) = exp( -points(1,:)**2 - 2.d0 * points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
rnd     => random(2,i_b)
pointsb = -1.d0 + 2.d0 * rnd
call tsgEvaluateBatch(grid,pointsb,i_b,double_2d_b)

allocate(dcmplx_1d_a(i_a),dcmplx_2d_c(i_a,i_b))
call tsgEvaluateComplexHierarchicalFunctions(grid,pointsb,i_b,dcmplx_2d_c)
call tsgGetComplexHierarchicalCoefficientsStatic(grid, dcmplx_1d_a)
double_1d_b = real(matmul(dcmplx_1d_a,dcmplx_2d_c))

dcmplx_pnt_1d_a => tsgGetComplexHierarchicalCoefficients(grid)
double_1d_c = real(matmul(dcmplx_pnt_1d_a,dcmplx_2d_c))

if ( norm1d(double_1d_b - double_2d_b(1,:)) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetComplexHierarchicalCoefficientsStatic: fourier grid test 1"
  stop 1
endif
if ( norm1d(double_1d_c - double_2d_b(1,:)) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetComplexHierarchicalCoefficientsStatic: fourier grid test 2"
  stop 1
endif
deallocate(points,pointsb,double_2d_a,double_2d_b,double_2d_c,double_1d_a,double_1d_b,double_1d_c,double_1d_d,rnd)
deallocate(dcmplx_pnt_1d_a,dcmplx_1d_a,dcmplx_2d_c)

! load coefficients
i_a = 13
allocate( double_2d_b(1,i_a), double_2d_c(1,i_a), pointsb(2,i_a) )
call tsgMakeLocalPolynomialGrid(grid,    2, 1, 5, 1, tsg_semi_localp)
call tsgMakeLocalPolynomialGrid(grid_II, 2, 1, 5, 1, tsg_semi_localp)
points  => tsgGetPoints(grid);  i_b = tsgGetNumPoints(grid)
rnd     => random(2,i_a)
pointsb =  -1.d0 + 2.d0 * rnd
allocate( double_2d_a(1,i_b), double_2d_d(i_b,i_b) )
double_2d_a(1,:) = exp( -points(1,:)**2 - 2.d0 * points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
double_pnt_1d_b => tsgGetHierarchicalCoefficients(grid)
call tsgSetHierarchicalCoefficients(grid_II,double_pnt_1d_b)
call tsgEvaluateBatch(grid,pointsb,i_a,double_2d_b)
call tsgEvaluateBatch(grid_II,pointsb,i_a,double_2d_c)
if ( norm2d(double_2d_b - double_2d_c) > 1.d-11 ) then
  write(*,*) "Mismatch in setHierarchicalCoefficients: localp grid solve ", norm2d(double_2d_b - double_2d_c)
  stop 1
endif
deallocate(points,pointsb,double_2d_a,double_2d_b,double_2d_c,double_2d_d,double_pnt_1d_b,rnd)


call tsgMakeLocalPolynomialGrid(grid, 2, 1, 5, 1, tsg_semi_localp)
call tsgMakeLocalPolynomialGrid(grid_II, 2, 1, 5, 1, tsg_semi_localp)
call tsgSetDomainTransform(grid,    (/-1.d0, 7.d0/), (/2.d0, 9.d0/))
call tsgSetDomainTransform(grid_II, (/-1.d0, 7.d0/), (/2.d0, 9.d0/))
points => tsgGetPoints(grid);  i_b = tsgGetNumPoints(grid)
i_a = 13
allocate(double_2d_a(1,i_b), double_2d_b(1,i_a), double_2d_c(1,i_a), double_2d_d(i_b,i_b), pointsb(2,i_a))
double_2d_a(1,:) = exp( -points(1,:)**2 - 2.d0 * points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
double_pnt_1d_b => tsgGetHierarchicalCoefficients(grid)
call tsgSetHierarchicalCoefficients(grid_II,double_pnt_1d_b)
rnd     => random(2,i_a)
pointsb = -1.d0 + 2.d0 * rnd
call tsgEvaluateBatch(grid,    pointsb, i_a, double_2d_b)
call tsgEvaluateBatch(grid_II, pointsb, i_a, double_2d_c)
if ( norm2d(double_2d_b - double_2d_c) > 1.d-11 ) then
  write(*,*) "Mismatch in setHierarchicalCoefficients: local grid solve, transform"
  stop 1
endif
deallocate(points,pointsb,double_2d_a,double_2d_b,double_2d_c,double_2d_d,double_pnt_1d_b,rnd)


call tsgMakeSequenceGrid(grid,    2, 1, 5, tsg_level, tsg_rleja)
call tsgMakeSequenceGrid(grid_II, 2, 1, 5, tsg_level, tsg_rleja)
call tsgSetDomainTransform(grid, (/1.d0, 1.d0/), (/2.d0, 2.d0/))
call tsgSetDomainTransform(grid_II, (/1.d0, 1.d0/), (/2.d0, 2.d0/))
points => tsgGetPoints(grid);  i_b = tsgGetNumPoints(grid)
i_a = 13
allocate(double_2d_a(1,i_b), double_2d_b(1,i_a), double_2d_c(1,i_a), double_2d_d(i_b,i_b), pointsb(2,i_a))
double_2d_a(1,:) = exp( -points(1,:)**2 - 2.d0 * points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
double_pnt_1d_b => tsgGetHierarchicalCoefficients(grid)
call tsgSetHierarchicalCoefficients(grid_II,double_pnt_1d_b)
rnd     => random(2,i_a)
pointsb = -1.d0 + 2.d0 * rnd
call tsgEvaluateBatch(grid,    pointsb, i_a, double_2d_b)
call tsgEvaluateBatch(grid_II, pointsb, i_a, double_2d_c)
if ( norm2d(double_2d_b - double_2d_c) > 1.d-11 ) then
  write(*,*) "Mismatch in setHierarchicalCoefficients: sequence grid solve, transform"
  stop 1
endif
deallocate(points,pointsb,double_2d_a,double_2d_b,double_2d_c,double_2d_d,double_pnt_1d_b,rnd)


call tsgMakeWaveletGrid(grid, 1, 1, 1, 3)
points => tsgGetHierarchicalSupport(grid)
double_2d_b = reshape( (/ 2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0, 2.d0, &
                          1.d4, 1.d4, 1.d4, 1.d4, 1.d4, 1.d4, 1.d4, 1.d4 /), shape(points) )
if ( norm2d(points - double_2d_b) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetHierarchicalSupport: wavelet grid"
  stop 1
endif
deallocate(points, double_2d_b)





!=======================================================================
!       tsgIntegrate()
!=======================================================================
call tsgMakeGlobalGrid(grid, 1, 1, 2, tsg_level, tsg_gauss_hermite, alpha=0.d0, beta=0.d0)
points => tsgGetPoints(grid)
i_a = tsgGetNumPoints(grid)
allocate(double_2d_a(1,i_a),double_1d_a(tsgGetNumOutputs(grid)))
double_2d_a(1,:) = points(1,:)**2
call tsgLoadNeededPoints(grid, double_2d_a)
call tsgIntegrate(grid, double_1d_a)
if ( abs(double_1d_a(1) - sqrt(pi)/2.d0) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgIntegrate: case 1"
  stop 1
endif
deallocate(points, double_2d_a, double_1d_a)


call tsgMakeGlobalGrid(grid, 1, 1, 2, tsg_level, tsg_gauss_hermite, alpha=2.d0, beta=0.d0)
points => tsgGetPoints(grid)
i_a = tsgGetNumPoints(grid)
allocate(double_2d_a(1,i_a),double_1d_a(tsgGetNumOutputs(grid)))
double_2d_a(1,:) = sqrt(2.d0)
call tsgLoadNeededPoints(grid, double_2d_a)
call tsgIntegrate(grid, double_1d_a)
if ( abs(double_1d_a(1) - sqrt(2.d0*pi)/2.d0) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgIntegrate: case 2"
  stop 1
endif
deallocate(points, double_2d_a, double_1d_a)




write(*,*) "Core I/O and evaluate:    PASS"



!=======================================================================
!       tsgGetInterpolationWeights()
!=======================================================================
i_a = 32
call tsgMakeGlobalGrid(grid,    2, 1, 4, tsg_level, tsg_fejer2)
call tsgMakeGlobalGrid(grid_II, 2, 1, 4, tsg_level, tsg_fejer2)
points  => tsgGetPoints(grid);  i_b = tsgGetNumPoints(grid)
allocate( pointsb(2,i_a), double_2d_a(1,i_b), double_2d_b(1,i_a), double_2d_c(1,i_a) )
rnd     => random(2,i_a)
pointsb =  -1.d0 + 2.d0 * rnd
double_2d_a(1,:) = exp( -points(1,:)**2 - 2.d0 * points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
call tsgEvaluateBatch(grid, pointsb, i_a, double_2d_b)
do i = 1,i_a
  double_pnt_1d_a  => tsgGetInterpolationWeights(grid_II, pointsb(:,i))
  double_2d_c(:,i) =  matmul(double_2d_a,double_pnt_1d_a)
  deallocate(double_pnt_1d_a)
enddo
if ( norm2d(double_2d_b - double_2d_c) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetInterpolationWeights: global case"
  stop 1
endif
do i = 1,i_a
  allocate( double_1d_a(tsgGetNumPoints(grid)) )
  call tsgGetInterpolationWeightsStatic(grid, pointsb(:,i), double_1d_a)
  double_2d_c(:,i) =  matmul(double_2d_a,double_1d_a)
  deallocate(double_1d_a)
enddo
if ( norm2d(double_2d_b - double_2d_c) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetInterpolationWeights: global case"
  stop 1
endif
deallocate(points, pointsb, double_2d_a, double_2d_b, double_2d_c, rnd)


i_a = 32
call tsgMakeSequenceGrid(grid,    2, 1, 7, tsg_level, tsg_min_delta)
call tsgMakeSequenceGrid(grid_II, 2, 1, 7, tsg_level, tsg_min_delta)
points  => tsgGetPoints(grid);  i_b = tsgGetNumPoints(grid)
allocate( pointsb(2,i_a), double_2d_a(1,i_b), double_2d_b(1,i_a), double_2d_c(1,i_a) )
rnd     => random(2,i_a)
pointsb =  -1.d0 + 2.d0 * rnd
double_2d_a(1,:) = exp( -points(1,:)**2 - 2.d0 * points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
call tsgEvaluateBatch(grid, pointsb, i_a, double_2d_b)
do i = 1,i_a
  double_pnt_1d_a  => tsgGetInterpolationWeights(grid_II, pointsb(:,i))
  double_2d_c(:,i) =  matmul(double_2d_a,double_pnt_1d_a)
  deallocate(double_pnt_1d_a)
enddo
if ( norm2d(double_2d_b - double_2d_c) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgGetInterpolationWeights: global case"
  stop 1
endif
deallocate(points, pointsb, double_2d_a, double_2d_b, double_2d_c, rnd)








!=======================================================================
!       tsgEstimateAnisotropicCoefficients()
!=======================================================================
call tsgMakeGlobalGrid(grid, 2, 1, 9, tsg_level, tsg_rleja)
points => tsgGetPoints(grid);  i_b = tsgGetNumPoints(grid)
allocate( double_2d_a(1,i_b) )
double_2d_a(1,:) = exp( points(1,:) + points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
int_pnt_1d_a => tsgEstimateAnisotropicCoefficients(grid,tsg_iptotal,1)
if ( abs( int_pnt_1d_a(1)/dble(int_pnt_1d_a(2)) - 2.0 ) > 0.2 ) then
  write(*,*) "Mismatch in tsgEstimateAnisotropicCoefficients: total degree"
  stop 1
endif
deallocate(int_pnt_1d_a)
int_pnt_1d_a => tsgEstimateAnisotropicCoefficients(grid,tsg_ipcurved,1)
if ( size(int_pnt_1d_a) .ne. 4 ) then
  write(*,*) "Mismatch in tsgEstimateAnisotropicCoefficients: curved dimensions"
  stop 1
endif
if ( (abs( int_pnt_1d_a(1)/dble(int_pnt_1d_a(2)) - 2.0 ) > 0.2) &
     .or. (int_pnt_1d_a(3)>0) .or. (int_pnt_1d_a(4)>0) ) then
  write(*,*) "Mismatch in tsgEstimateAnisotropicCoefficients: curved"
  stop 1
endif
deallocate(int_pnt_1d_a)
deallocate(points, double_2d_a)








!=======================================================================
!       tsgSetAnisotropicRefinement()
!=======================================================================
call tsgMakeSequenceGrid(grid, 3, 1, 3, tsg_level, tsg_leja, levelLimits=(/3,2,1/))
points => tsgGetPoints(grid);  i_b = tsgGetNumPoints(grid)
if ( check_points( points(1,:), (/0.d0,-1.d0,1.d0,1.d0/sqrt(3.d0)/), (/12345.d0/) ) .or. &
     check_points( points(2,:), (/0.d0,-1.d0,1.d0/), (/1.d0/sqrt(3.d0)/) ) .or. &
     check_points( points(3,:), (/0.d0,1.d0/), (/-1.d0,1.d0/sqrt(3.d0)/) ) ) then
  write(*,*) "Mismatch in tsgSetAnisotropicRefinement: limits in make"
  stop 1
endif
allocate( double_2d_a(1,i_b) )
double_2d_a(1,:) = exp( -points(1,:)**2 - points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
call tsgSetAnisotropicRefinement( grid, tsg_iptotal, 5, 1 )
pointsb => tsgGetNeededPoints(grid)
if ( size(pointsb,2) .eq. 0 ) then
  write(*,*) "Mismatch in tsgSetAnisotropicRefinement: did not refine"
  stop 1
endif
if ( check_points( points(2,:), (/12345.d0/), (/1.d0/sqrt(3.d0)/) ) .or. &
     check_points( points(3,:), (/12345.d0/), (/-1.d0,1.d0/sqrt(3.d0)/) ) ) then
  write(*,*) "Mismatch in tsgSetAnisotropicRefinement: limits refine using existing limits"
  stop 1
endif
deallocate( pointsb )
call tsgSetAnisotropicRefinement(grid,tsg_iptotal,10,1,levelLimits=(/3,2,2/))
pointsb => tsgGetNeededPoints(grid)
if ( size(pointsb) .eq. 0 ) then
  write(*,*) "Mismatch in tsgSetAnisotropicRefinement: did not refine"
  stop 1
endif
if ( check_points( points(2,:), (/12345.d0/), (/1.d0/sqrt(3.d0)/) ) .or. &
     check_points( points(3,:), (/-1.d0/), (/1.d0/sqrt(3.d0)/) ) ) then
  write(*,*) "Mismatch in tsgSetAnisotropicRefinement: limits refine using new limits"
  stop 1
endif
deallocate( points, pointsb, double_2d_a )








!=======================================================================
!       tsgSetLocalSurplusRefinement()
!=======================================================================
call tsgMakeLocalPolynomialGrid(grid, 3, 1, 3, 1, tsg_localp, levelLimits=(/1,2,3/))
points => tsgGetPoints(grid);  i_b = tsgGetNumPoints(grid)
if ( check_points( points(1,:), (/0.d0,-1.d0,1.d0/), (/0.5d0,-0.5d0,-0.75d0,-0.25d0,0.25d0,0.75d0/) ) .or. &
     check_points( points(2,:), (/0.d0,-1.d0,1.d0,0.5d0,-0.5d0/), (/-0.75d0,-0.25d0,0.25d0,0.75d0/) ) .or. &
     check_points( points(3,:), (/0.d0,-1.d0,1.d0,0.5d0,-0.5d0,-0.75d0,-0.25d0,0.25d0,0.75d0/), (/12345.d0/) ) ) then
  write(*,*) "Mismatch in tsgSetLocalSurplusRefinement: limits in make"
  stop 1
endif
allocate( double_2d_a(1,i_b) )
double_2d_a(1,:) = exp( -points(1,:)**2 - points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
call tsgSetLocalSurplusRefinement(grid,1.d-8,tsg_classic,1)
if ( tsgGetNumNeeded(grid) .eq. 0 ) then
  write(*,*) "Mismatch in tsgSetLocalSurplusRefinement: did not refine local polynomial"
  stop 1
endif
if ( check_points( points(1,:), (/12345.d0/), (/0.5d0,-0.5d0,-0.75d0,-0.25d0,0.25d0,0.75d0/) ) .or. &
     check_points( points(2,:), (/12345.d0/), (/-0.75d0,-0.25d0,0.25d0,0.75d0/) ) .or. &
     check_points( points(3,:), (/12345.d0/), (/12345.d0/) ) ) then
  write(*,*) "Mismatch in tsgSetSurplusRefinement: limits refine using existing limits"
  stop 1
endif
call tsgSetLocalSurplusRefinement(grid,1.d-8,tsg_classic,1,(/2,2,3/))
if ( tsgGetNumNeeded(grid) .eq. 0 ) then
  write(*,*) "Mismatch in tsgSetLocalSurplusRefinement: did not refine on second pass"
  stop 1
endif
if ( check_points( points(1,:), (/0.5d0,-0.5d0/), (/-0.75d0,-0.25d0,0.25d0,0.75d0/) ) .or. &
     check_points( points(2,:), (/12345.d0/), (/-0.75d0,-0.25d0,0.25d0,0.75d0/) ) ) then
  write(*,*) "Mismatch in tsgSetSurplusRefinement: limits refine using new limits"
  stop 1
endif
deallocate(points, double_2d_a)








!=======================================================================
!       tsgClearRefinement()
!=======================================================================
call tsgMakeLocalPolynomialGrid(grid, 3, 1, 4, 2, tsg_localp, levelLimits=(/300,300,300/))
points => tsgGetPoints(grid);  i_b = tsgGetNumPoints(grid)
allocate( double_2d_a(1,i_b) )
double_2d_a(1,:) = exp( -points(1,:)**2 - 0.5*points(2,:)**2 - 2.d0*points(3,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
if ( tsgGetNumNeeded(grid) > 0 ) then
  write(*,*) "Mismatch in cancel refine: did not load values"
  stop 1
endif
call tsgSetLocalSurplusRefinement(grid,1.d-4,tsg_directional)
if ( tsgGetNumNeeded(grid) .eq. 0 ) then
  write(*,*) "Mismatch in cancel refine: did not set refinement at output -1"
  stop 1
endif
call tsgClearRefinement(grid)
if ( tsgGetNumNeeded(grid) > 0 ) then
  write(*,*) "Mismatch in cancel refine: did not cancel the refinement"
  stop 1
endif
deallocate(points, double_2d_a)








!=======================================================================
!       tsgMergeRefinement()
!=======================================================================
call tsgMakeGlobalGrid(grid,    2, 1, 4, tsg_level, tsg_fejer2)
call tsgMakeGlobalGrid(grid_II, 2, 1, 4, tsg_level, tsg_fejer2)

points => tsgGetPoints(grid);  i_b = tsgGetNumPoints(grid)
allocate( double_2d_a(1,i_b) )
double_2d_a(1,:) = exp( -points(1,:)**2 - points(2,:)**2 )
call tsgLoadNeededPoints(grid,    double_2d_a)
call tsgLoadNeededPoints(grid_II, double_2d_a)
call tsgSetAnisotropicRefinement( grid,    tsg_iptotal, 10, 1 )
call tsgSetAnisotropicRefinement( grid_II, tsg_iptotal, 10, 1 )
deallocate(points, double_2d_a)

points => tsgGetNeededPoints(grid);  i_b = tsgGetNumNeeded(grid)
allocate( double_2d_a(1,i_b) )
double_2d_a(1,:) = exp( -points(1,:)**2 - points(2,:)**2 )
call tsgLoadNeededPoints(grid, double_2d_a)
call tsgMergeRefinement(grid_II)
deallocate(points, double_2d_a)

i_a = 20
points  => tsgGetPoints(grid);     i_b = tsgGetNumPoints(grid)
pointsb => tsgGetPoints(grid_II);  i_c = tsgGetNumPoints(grid_II)
if ( norm2d(points-pointsb) > 1.d-11 )then
  write(*,*) "Mismatch in tsgMergeRefine(): case 2, tsgGetPoints()"
  stop 1
endif
deallocate(points, pointsb)

allocate(points(2,i_a), double_2d_a(1,i_a))
rnd    => random(2,i_a)
points = -1.d0 + 2.d0 * rnd
call tsgEvaluateBatch(grid_II,points,i_a,double_2d_a)
if ( norm2d(double_2d_a) > 1.d-11 ) then
  write(*,*) "Mismatch in tsgMergeRefine(): case 3, tsgEvaluate() not zero"
  stop 1
endif
deallocate(rnd, points, double_2d_a)

points  => tsgGetPoints(grid_II);     i_b = tsgGetNumPoints(grid_II)
allocate( double_2d_a(1,i_b) )
double_2d_a(1,:) = exp( -points(1,:)**2 - points(2,:)**2 )
call tsgLoadNeededPoints(grid_II, double_2d_a)
deallocate(points, double_2d_a)

i_a = 30
allocate( points(2,i_a), double_2d_a(1,i_a), double_2d_b(1,i_a) )
rnd     => random(2,i_a)
points = -1.d0 + 2.d0 * rnd
call tsgEvaluateBatch(grid,    points, i_a, double_2d_a)
call tsgEvaluateBatch(grid_II, points, i_a, double_2d_b)
if ( norm2d(double_2d_a-double_2d_b) > 1.d-11 )then
  write(*,*) "Mismatch in tsgMergeRefine(): case 3, tsgEvaluate() not equal"
  stop 1
endif
deallocate(points, double_2d_a, double_2d_b, rnd)


write(*,*) "Refinement functions:     PASS"




call tsgDeallocateGrid(grid)
call tsgDeallocateGrid(grid_II)





contains

function random(n,m) result(res)
  integer :: n, m
  double precision, pointer :: res(:,:)

  allocate(res(n,m))
  call random_number(res)
end function

function norm(x,n) result(res)
  integer :: n
  double precision :: x(n)
  double precision :: res

  res = sqrt(sum(x**2))
end function

function norm1d(x) result(res)
  double precision :: x(:)
  double precision :: res

  res = sqrt(sum(x**2))
end function

function norm2d(x) result(res)
  double precision :: x(:,:)
  double precision :: res

  res = sqrt(sum(x**2))
end function

function check_points(points,mustHave,mustNotHave) result(res)
  double precision :: points(:), mustHave(:), mustNotHave(:)
  logical :: res
  integer :: i, j

  res = .false.
  do i = 1,size(points)
    if ( any(points(i).ne.mustHave) .and. any(points(i).eq.mustNotHave) ) then
      res = .true.
      exit
    endif
  enddo
end function

subroutine print_vector(vec)
  double precision :: vec(:)
  integer :: i

  do i = 1,size(vec)
    write(*,'(f6.3 a)',advance='no') vec(i), " "
  enddo
  write(*,*)
end subroutine

subroutine print_array(arr)
  double precision :: arr(:,:)
  integer :: shp(2), i, j

  shp = shape(arr)

  do i = 1,shp(1)
    do j = 1,shp(2)
      write(*,'(f6.3 a)',advance='no') arr(i,j), " "
    enddo
    write(*,*)
  enddo
end subroutine

function sparse_matmul(pntr,indx,val,vec) result(res)
  integer :: pntr(:), indx(:)
  double precision :: val(:), vec(:)
  integer :: i, j
  double precision :: res(size(pntr)-1)

  res = 0.d0
  do i = 1,size(pntr)-1
    do j = pntr(i)+1,pntr(i+1)
      res(i) = res(i) + val(j) * vec(indx(j)+1)
    enddo
  enddo
end function

function Re_complex_matmul(mat,vec) result(res)
  double precision :: mat(:,:), vec(:)
  integer :: i, j
  integer :: shp(2)
  double precision :: res(size(mat,1))

  shp = shape(mat)

  res = 0.d0
  do i = 1,shp(1)
    do j = 1,shp(2)/2
      res(i) = res(i) + mat(i,2*j-1) * vec(2*j-1) - mat(i,2*j) * vec(2*j)
    enddo
  enddo
end function

end program FORTESTER
