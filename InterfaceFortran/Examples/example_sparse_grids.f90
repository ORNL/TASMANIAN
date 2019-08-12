!==================================================================================================================================================================================
! Copyright (c) 2017, Miroslav Stoyanov
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
PROGRAM TasmanianSGExample
  USE TasmanianSG
IMPLICIT NONE
  type(TasmanianSparseGrid) :: grid, grid1, grid2, grid3
  INTEGER :: dims, outs, level
  INTEGER :: N1, N2, N3
  INTEGER :: N, i, j, verm, vern
  DOUBLE PRECISION :: err1, err2, err3, err4, exact
  REAL :: cpuStart, cpuEnd, stages(2,3)
  INTEGER :: conformal(3)
  CHARACTER, pointer :: string(:)
  DOUBLE PRECISION, pointer :: points(:,:), weights(:)
  DOUBLE PRECISION :: x, y, integ, E
  DOUBLE PRECISION, allocatable :: transformA(:), transformB(:), values(:,:), tvalues(:,:)
  DOUBLE PRECISION, allocatable :: res(:), res2(:,:)
  DOUBLE PRECISION :: desired_x(2)
  DOUBLE PRECISION :: randPoints(4,1000), randPoints2(2,1000), randPoints3(3,1000)
  DOUBLE PRECISION :: PI = 4.D0*DATAN(1.D0)


! This is the sound "Glaucodon Ballaratensis" makes :)
!  WRITE(*,*) "Ghurrrrrphurrr"

! ============ reference table of rules ============ !
!  1: clenshaw-curtis             2: clenshaw-curtis-zero
!  3: chebyshev                   4: chebyshev-odd
!  5: gauss-legendre              6: gauss-legendreodd
!  7: gauss-patterson             8: leja
!  9: lejaodd                    10: rleja
! 11: rleja-odd                  12: rleja-double-2
! 13: rleja-double-4             14: rleja-shifted
! 15: rleja-shifted-even         16: rleja-shifted-double
! 17: max-lebesgue               18: max-lebesgue-odd
! 19: min-lebesgue               20: min-lebesgue-odd
! 21: min-delta                  22: min-delta-odd
! 23: gauss-chebyshev-1          24: gauss-chebyshev-1-odd
! 25: gauss-chebyshev-2          26: gauss-chebyshev-2-odd
! 27: fejer2                     28: gauss-gegenbauer
! 29: gauss-gegenbauer-odd       30: gauss-jacobi
! 31: gauss-jacobi-odd           32: gauss-laguerre
! 33: gauss-laguerre-odd         34: gauss-hermite
! 35: gauss-hermite-odd          36: custom-tabulated
! 37: localp                     38: localp-zero
! 39: semi-localp                40: wavelet

! ============ reference table of grid types ============ !
!  1: level                       2: curved
!  3: iptotal                     4: ipcurved
!  5: qptotal                     6: qpcurved
!  7: hyperbolic                  8: iphyperbolic
!  9: qphyperbolic               10: tensor
! 11: iptensor                   12: qptensor

! ============ reference table of refinement types ============ !
!  1: classic            2: parents first
!  3: directional        4: FDS (both parents and directions)

! ============ reference table of acceleration types ============ !
!  0: none       1: CPU BLAS    2: GPU cuBLAS
!  3: GPU CUDA   4: GPU MAGMA

! ==================================================================== !
! EXAMPLE 1: integrate: f(x,y) = exp(-x^2) * cos(y) over [-1,1] x [-1,1]
! using classical Smolyak grid with Clenshaw-Curtis points and weights


  verm = tsgGetVersionMajor()
  vern = tsgGetVersionMinor()

  ! WARNING: do not DEALLOCATE the string pointer, it is const char*
  string => tsgGetLicense()
  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Tasmanian Sparse Grids Fortran Module (TasmanianSG)"
  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,"(A,I2,A,I1)") " Tasmanian Sparse Grid module, version: ", verm, ".", vern
  WRITE(*,"(A,40A)")       "                               license: ", string
  WRITE(*,*)

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 1:  integrate f(x,y) = exp(-x^2) * cos(y), using clenshaw-curtis level nodes"

  dims = 2
  level = 6

! before you use a grid, you must allocate a new grid
  call tsgAllocateGrid(grid)

  CALL tsgMakeGlobalGrid(grid, dims, 0, level, tsg_level, tsg_clenshaw_curtis)
  points => tsgGetPoints(grid)
  weights => tsgGetQuadratureWeights(grid)

  N = tsgGetNumPoints(grid)
  integ = 0.0

  DO i = 1, N
    x = points(1, i)
    y = points(2, i)
    integ = integ + weights(i) * exp(-x*x) * cos(y)
  END DO

  E = abs(integ - 2.513723354063905D+00)

  WRITE(*,"(A,I4)")   "      at level:     ", level
  WRITE(*,"(A,I4,A)") "      the grid has: ", N, " points"
  WRITE(*,"(A,E25.16)") "      integral:  ", integ
  WRITE(*,"(A,E25.16)") "      error:     ", E
  WRITE(*,*)

  level = 7
! no need to ask for a new ID when remaking an existing grid
  CALL tsgMakeGlobalGrid(grid, dims, 0, level, tsg_level, tsg_clenshaw_curtis)

! do not forget to release the memory associated with points and weights
  DEALLOCATE(points)
  DEALLOCATE(weights)
  points => tsgGetPoints(grid)
  weights => tsgGetQuadratureWeights(grid)

  N = tsgGetNumPoints(grid)
  integ = 0.0

  DO i = 1, N
    x = points(1, i)
    y = points(2, i)
    integ = integ + weights(i) * exp(-x*x) * cos(y)
  END DO

  E = abs(integ - 2.513723354063905D+00)
  WRITE(*,"(A,I4)")   "      at level:     ", level
  WRITE(*,"(A,I4,A)") "      the grid has: ", N, " points"
  WRITE(*,"(A,E25.16)") "      integral:  ", integ
  WRITE(*,"(A,E25.16)") "      error:     ", E
  WRITE(*,*)

  DEALLOCATE(points)
  DEALLOCATE(weights)

! after calling tsgDeallocateGrid(), we can no longer use this grid
! until we call tsgAllocateGrid() again
  CALL tsgDeallocateGrid(grid)

! ==================================================================== !
! EXAMPLE 2: integrate: f(x,y) = exp(-x^2) * cos(y)
!                       over (x,y) in [-5,5] x [-2,3]
! using Gauss-Patterson rules chosen to integrate exactly polynomials of
! total degree up to degree specified by prec

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 2: integrate f(x,y) = exp(-x^2) * cos(y) over [-5,5] x [-2,3] using  Gauss-Patterson nodes"

  dims = 2
  level = 20

  ALLOCATE(transformA(dims))
  ALLOCATE(transformB(dims))
  transformA(1) = -5.0
  transformA(2) = -2.0
  transformB(1) =  5.0
  transformB(2) =  3.0

! need new grid, since we freed this earlier
  call tsgAllocateGrid(grid)

! gauss-patterson = 7, type_qptotal = 5
  CALL tsgMakeGlobalGrid(grid, dims, 0, level, tsg_qptotal, tsg_gauss_patterson)
  CALL tsgSetDomainTransform(grid, transformA, transformB)

  points => tsgGetPoints(grid)
  weights => tsgGetQuadratureWeights(grid)

  N = tsgGetNumPoints(grid)
  integ = 0.0

  DO i = 1, N
    x = points(1, i)
    y = points(2, i)
    integ = integ + weights(i) * exp(-x*x) * cos(y)
  END DO

  E = abs(integ - 1.861816427518323D+00)
  WRITE(*,"(A,I4)")   "      at precision:   ", level
  WRITE(*,"(A,I4,A)") "      the grid has:   ", N, " points"
  WRITE(*,"(A,E25.16)") "      integral:    ", integ
  WRITE(*,"(A,E25.16)") "      error:       ", E
  WRITE(*,*)

  level = 40
! no need to ask for a new ID when remaking an existing grid
  CALL tsgMakeGlobalGrid(grid, dims, 0, level, tsg_qptotal, tsg_gauss_patterson)
  CALL tsgSetDomainTransform(grid, transformA, transformB)

! do not forget to release the memory associated with points and weights
  DEALLOCATE(points)
  DEALLOCATE(weights)
  points => tsgGetPoints(grid)
  weights => tsgGetQuadratureWeights(grid)

  N = tsgGetNumPoints(grid)
  integ = 0.0

  DO i = 1, N
    x = points(1, i)
    y = points(2, i)
    integ = integ + weights(i) * exp(-x*x) * cos(y)
  END DO

  E = abs(integ - 1.861816427518323D+00)
  WRITE(*,"(A,I4)")   "      at precision:   ", level
  WRITE(*,"(A,I4,A)") "      the grid has:   ", N, " points"
  WRITE(*,"(A,E25.16)") "      integral:    ", integ
  WRITE(*,"(A,E25.16)") "      error:       ", E
  WRITE(*,*)

  DEALLOCATE(points)
  DEALLOCATE(weights)
  ! keep transformA and transformB for the next example
  CALL tsgDeallocateGrid(grid)

! ==================================================================== !
! EXAMPLE 3: integrate: f(x,y) = exp(-x^2) * cos(y)
!                       over (x,y) in [-5,5] x [-2,3]
! using different rules

  call tsgAllocateGrid(grid1)
  call tsgAllocateGrid(grid2)
  call tsgAllocateGrid(grid3)

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 3: integrate f(x,y) = exp(-x^2) * cos(y) over [-5,5] x [-2,3] using different rules"
  WRITE(*,*)

  WRITE(*,*) "               Clenshaw-Curtis      Gauss-Legendre     Gauss-Patterson"
  WRITE(*,*) " precision    points     error    points     error    points     error"

  DO level = 9, 30, 4
    CALL tsgMakeGlobalGrid(grid1, dims, 0, level, tsg_qptotal, tsg_clenshaw_curtis)
    CALL tsgSetDomainTransform(grid1, transformA, transformB)
    CALL tsgMakeGlobalGrid(grid2, dims, 0, level, tsg_qptotal, tsg_gauss_legendre)
    CALL tsgSetDomainTransform(grid2, transformA, transformB)
    CALL tsgMakeGlobalGrid(grid3, dims, 0, level, tsg_qptotal, tsg_gauss_patterson)
    CALL tsgSetDomainTransform(grid3, transformA, transformB)

    points => tsgGetPoints(grid1)
    weights => tsgGetQuadratureWeights(grid1)
    N1 = tsgGetNumPoints(grid1)
    integ = 0.0
    DO i = 1, N1
      x = points(1, i)
      y = points(2, i)
      integ = integ + weights(i) * exp(-x*x) * cos(y)
    END DO
    err1 = abs(integ - 1.861816427518323D+00)
    DEALLOCATE(points)
    DEALLOCATE(weights)

    points => tsgGetPoints(grid2)
    weights => tsgGetQuadratureWeights(grid2)
    N2 = tsgGetNumPoints(grid2)
    integ = 0.0
    DO i = 1, N2
      x = points(1, i)
      y = points(2, i)
      integ = integ + weights(i) * exp(-x*x) * cos(y)
    END DO
    err2 = abs(integ - 1.861816427518323D+00)
    DEALLOCATE(points)
    DEALLOCATE(weights)

    points => tsgGetPoints(grid3)
    weights => tsgGetQuadratureWeights(grid3)
    N3 = tsgGetNumPoints(grid3)
    integ = 0.0
    DO i = 1, N3
      x = points(1, i)
      y = points(2, i)
      integ = integ + weights(i) * exp(-x*x) * cos(y)
    END DO
    err3 = abs(integ - 1.861816427518323D+00)
    DEALLOCATE(points)
    DEALLOCATE(weights)

    WRITE(*,"(I10,I10,E11.3,I9,E11.3,I9,E11.3)") level, N1, err1, N2, err2, N3, err3

  END DO
  WRITE(*,*)

  CALL tsgDeallocateGrid(grid1)
  CALL tsgDeallocateGrid(grid2)
  CALL tsgDeallocateGrid(grid3)

  DEALLOCATE(transformA)
  DEALLOCATE(transformB)


! ==================================================================== !
! EXAMPLE 4: interpolate: f(x,y) = exp(-x^2) * cos(y)
! with a rule that exactly interpolates polynomials of total degree

  call tsgAllocateGrid(grid)

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 4: interpolate f(x,y) = exp(-x^2) * cos(y), using clenshaw-curtis iptotal rule"
  WRITE(*,*)

  dims = 2
  outs = 1
  level = 10

  desired_x(1) = 0.3
  desired_x(2) = 0.7

  ! desired value
  exact = exp(-desired_x(1)**2) * cos(desired_x(2))

  CALL tsgMakeGlobalGrid(grid, dims, outs, level, tsg_iptotal, tsg_clenshaw_curtis)

  N = tsgGetNumNeeded(grid)
  points => tsgGetNeededPoints(grid)
  ALLOCATE(values(outs,N))

  DO i = 1, N
    x = points(1, i)
    y = points(2, i)
    values(1,i) = exp(-x**2) * cos(y)
  END DO

  CALL tsgLoadNeededPoints(grid, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  ALLOCATE(res(outs)) ! will DEALLOCATE later
  CALL tsgEvaluate(grid, desired_x, res)
  E = abs(res(1) - exact)

  WRITE(*,"(A,I4)")   "  using polynomials of total degree:  ", level
  WRITE(*,"(A,I4,A)") "      the grid has:                   ", N, " points"
  WRITE(*,"(A,E25.16)") "      interpolant at (0.3,0.7):    ", res(1)
  WRITE(*,"(A,E25.16)") "      error:                       ", E
  WRITE(*,*)

  ! do the same with level = 12
  level = 12

  CALL tsgMakeGlobalGrid(grid, dims, outs, level, tsg_iptotal, tsg_clenshaw_curtis)

  N = tsgGetNumNeeded(grid)
  points => tsgGetNeededPoints(grid)
  ALLOCATE(values(outs,N))

  DO i = 1, N
    x = points(1, i)
    y = points(2, i)
    values(1,i) = exp(-x**2) * cos(y)
  END DO

  CALL tsgLoadNeededPoints(grid, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  CALL tsgEvaluate(grid, desired_x, res)
  E = abs(res(1) - exact)

  WRITE(*,"(A,I4)")   "  using polynomials of total degree:  ", level
  WRITE(*,"(A,I4,A)") "      the grid has:                   ", N, " points"
  WRITE(*,"(A,E25.16)") "      interpolant at (0.3,0.7):    ", res(1)
  WRITE(*,"(A,E25.16)") "      error:                       ", E
  WRITE(*,*)
  DEALLOCATE(res)

  CALL tsgDeallocateGrid(grid)

! prepare random smaples for future tests
  call random_seed()

  call random_number(randPoints)

! ==================================================================== !
! EXAMPLE 5:
! interpolate: f(x1,x2,x3,x4) = exp(-x1^2) * cos(x2) * exp(-x3^2) * cos(x4)
! with Global and Sequence Leja rules
! Skipping example 5: it takes time and slows down testing, but it dosn't show/test any API not covered in other cases

!  dims = 4
!  outs = 1
!  level = 15
!
!  CALL cpu_time(cpuStart)
!  CALL tsgMakeGlobalGrid(grid, dims, outs, level, tsg_level, tsg_leja)
!  CALL cpu_time(cpuEnd)
!  stages(1,1) = cpuEnd - cpuStart
!
!  N = tsgGetNumPoints(grid)
!
!  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
!  WRITE(*,*) "Example 5: interpolate f(x1,x2,x3,x4) = exp(-x1^2) * cos(x2) * exp(-x3^2) * cos(x4)"
!  WRITE(*,*) "           comparign the performance of Global and Sequence grids with leja nodes"
!  WRITE(*,"(A,I4)") "            using polynomials of total degree up to: ", level
!  WRITE(*,"(A,I4,A)") "            the grids have:                          ", N, " points"
!  WRITE(*,*) "           both grids are evaluated at 1000 random points "
!  WRITE(*,*)
!
!  points => tsgGetNeededPoints(grid)
!  ALLOCATE(values(outs,N))
!  DO i = 1, N
!    values(1,i) = exp(-points(1,i)**2) * cos(points(2,i)) * exp(-points(3,i)**2) * cos(points(4,i))
!  END DO
!  DEALLOCATE(points)
!
!  CALL cpu_time(cpuStart)
!  CALL tsgLoadNeededPoints(grid, values)
!  CALL cpu_time(cpuEnd)
!  stages(1,2) = cpuEnd - cpuStart
!
!  ALLOCATE(res2(outs,1000)) ! 2-D result
!
!  CALL cpu_time(cpuStart)
!  CALL tsgEvaluateBatch(grid, randPoints, 1000, res2)
!  CALL cpu_time(cpuEnd)
!  stages(1,3) = cpuEnd - cpuStart
!
!  CALL cpu_time(cpuStart)
!  CALL tsgMakeSequenceGrid(grid, dims, outs, level, tsg_level, tsg_leja)
!  CALL cpu_time(cpuEnd)
!  stages(2,1) = cpuEnd - cpuStart
!
!  ! points are the same, no need to recompue values
!  CALL cpu_time(cpuStart)
!  CALL tsgLoadNeededPoints(grid, values)
!  CALL cpu_time(cpuEnd)
!  stages(2,2) = cpuEnd - cpuStart
!
!  DEALLOCATE(values)
!
!  CALL cpu_time(cpuStart)
!  CALL tsgEvaluateBatch(grid, randPoints, 1000, res2)
!  CALL cpu_time(cpuEnd)
!  stages(2,3) = cpuEnd - cpuStart
!
!  WRITE(*,*) "Stage            Global Grid         Sequence Grid"
!  WRITE(*,"(A,E20.8,E20.8)") " make grid  ", stages(1,1), stages(2,1)
!  WRITE(*,"(A,E20.8,E20.8)") " load needed", stages(1,2), stages(2,2)
!  WRITE(*,"(A,E20.8,E20.8)") " evaluate   ", stages(1,3), stages(2,3)
!
!  CALL tsgDeallocateGrid(grid)
!  DEALLOCATE(res2)

! ==================================================================== !
! EXAMPLE 6:
! interpolate: f(x,y) = exp(-x^2) * cos(y)
! using different refinement schemes

  ALLOCATE(tvalues(1,1000)) ! true values of f(x,y)
  DO i = 1, 1000
    randPoints2(1,i) = randPoints(1,i)
    randPoints2(2,i) = randPoints(2,i)
    tvalues(1,i) = exp(-randPoints(1,i)**2) * cos(randPoints(2,i))
  END DO

  call tsgAllocateGrid(grid1)
  call tsgAllocateGrid(grid2)
  call tsgAllocateGrid(grid3)

  dims = 2
  outs = 1

  CALL tsgMakeGlobalGrid(grid1, dims, outs, 3, tsg_iptotal, tsg_leja)
  CALL tsgMakeGlobalGrid(grid2, dims, outs, 3, tsg_iptotal, tsg_leja)
  CALL tsgMakeGlobalGrid(grid3, dims, outs, 3, tsg_iptotal, tsg_leja)

  N = tsgGetNumNeeded(grid1)
  points => tsgGetNeededPoints(grid1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
  END DO
  CALL tsgLoadNeededPoints(grid1, values)
  CALL tsgLoadNeededPoints(grid2, values)
  CALL tsgLoadNeededPoints(grid3, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 6: interpolate: f(x,y) = exp(-x^2) * cos(y)"
  WRITE(*,*) "   using leja nodes and different refinement schemes "
  WRITE(*,*) "   the error is estimated as the maximum from 1000 random points"

  WRITE(*,*) "           Total Degree            Curved               Surplus"
  WRITE(*,*) "iteration  points     error     points     error     points     error"

  ALLOCATE(res2(outs,1000)) ! 2-D result

  ! iptotal: 3, ipcurved: 4
  DO j = 1, 10
    CALL tsgSetAnisotropicRefinement(grid1, tsg_iptotal, 10, 1)

    N = tsgGetNumNeeded(grid1)
    points => tsgGetNeededPoints(grid1)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
    END DO
    CALL tsgLoadNeededPoints(grid1, values)
    DEALLOCATE(values)
    DEALLOCATE(points)

    CALL tsgEvaluateBatch(grid1, randPoints2, 1000, res2)
    err1 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
        err1 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO

    CALL tsgSetAnisotropicRefinement(grid2, tsg_ipcurved, 10, 1)

    N = tsgGetNumNeeded(grid2)
    points => tsgGetNeededPoints(grid2)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
    END DO
    CALL tsgLoadNeededPoints(grid2, values)
    DEALLOCATE(values)
    DEALLOCATE(points)

    CALL tsgEvaluateBatch(grid2, randPoints2, 1000, res2)
    err2 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err2)then
        err2 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO

    CALL tsgSetGlobalSurplusRefinement(grid3, 1.D-10, 1)

    N = tsgGetNumNeeded(grid3)
    points => tsgGetNeededPoints(grid3)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
    END DO
    CALL tsgLoadNeededPoints(grid3, values)
    DEALLOCATE(values)
    DEALLOCATE(points)

    CALL tsgEvaluateBatch(grid3, randPoints2, 1000, res2)
    err3 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err3)then
        err3 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO

    N1 = tsgGetNumPoints(grid1)
    N2 = tsgGetNumPoints(grid2)
    N3 = tsgGetNumPoints(grid3)

    WRITE(*,"(I9,I9,E12.4,I9,E12.4,I9,E12.4)") j, N1, err1, N2, err2, N3, err3

  END DO
  WRITE(*,*)
  CALL tsgDeallocateGrid(grid1)
  CALL tsgDeallocateGrid(grid2)
  CALL tsgDeallocateGrid(grid3)

! ==================================================================== !
! EXAMPLE 7:
! interpolate: f(x,y) = exp(-x^2) * cos(y)
! using localp and semilocalp grids

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 7: interpolate: f(x,y) = exp(-x^2) * cos(y)"
  WRITE(*,*) "       using localp and semi-localp rules with depth 7"
  WRITE(*,*) "       the error is estimated as the maximum from 1000 random points"
  WRITE(*,*)

  call tsgAllocateGrid(grid1)
  call tsgAllocateGrid(grid2)

  dims = 2
  outs = 1

  CALL tsgMakeLocalPolynomialGrid(grid1, dims, outs, 7, 2, tsg_localp)
  CALL tsgMakeLocalPolynomialGrid(grid2, dims, outs, 7, 2, tsg_semi_localp)

  N = tsgGetNumNeeded(grid1)
  points => tsgGetNeededPoints(grid1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
  END DO
  CALL tsgLoadNeededPoints(grid1, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  CALL tsgEvaluateBatch(grid1, randPoints2, 1000, res2)
  err1 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
      err1 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO

  points => tsgGetNeededPoints(grid2)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
  END DO
  CALL tsgLoadNeededPoints(grid2, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  CALL tsgEvaluateBatch(grid2, randPoints2, 1000, res2)
  err2 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err2)then
      err2 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO

  WRITE(*,"(A,I5)") "     Number of points: ", N
  WRITE(*,"(A,E12.4)") "     Error for      rule_localp: ", err1
  WRITE(*,"(A,E12.4)") "     Error for  rule_semilocalp: ", err2
  WRITE(*,*) " Note: semi-localp wins this competition because the function is very smooth"
  WRITE(*,*)

! ==================================================================== !
! EXAMPLE 8:
! interpolate: f(x,y) = cos(0.5 * pi * x) * cos(0.5 * pi * y)
! using localp and semilocalp grids

! remake the true values for the next example
  DO i = 1, 1000
    tvalues(1,i) = cos(0.5 * PI * randPoints2(1,i)) &
                 * cos(0.5 * PI * randPoints2(2,i))
  END DO

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 8: interpolate f(x,y) = cos(0.5 * pi * x) * cos(0.5 * pi * y)"
  WRITE(*,*) "       using localp and localp-zero rules with depths 7 and 6"
  WRITE(*,*) "       the error is estimated as the maximum from 1000 random points"
  WRITE(*,*)

  dims = 2
  outs = 1

  CALL tsgMakeLocalPolynomialGrid(grid1, dims, outs, 7, 2, tsg_localp)
  CALL tsgMakeLocalPolynomialGrid(grid2, dims, outs, 6, 2, tsg_localp_zero)

  N = tsgGetNumNeeded(grid1)
  points => tsgGetNeededPoints(grid1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = cos(0.5 * PI * points(1,i)) &
                * cos(0.5 * PI * points(2,i))
  END DO
  CALL tsgLoadNeededPoints(grid1, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  CALL tsgEvaluateBatch(grid1, randPoints2, 1000, res2)
  err1 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
      err1 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO

  N = tsgGetNumNeeded(grid2)
  points => tsgGetNeededPoints(grid2)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = cos(0.5 * PI * points(1,i)) * cos(0.5 * PI * points(2,i))
  END DO
  CALL tsgLoadNeededPoints(grid2, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  CALL tsgEvaluateBatch(grid2, randPoints2, 1000, res2)
  err2 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err2)then
      err2 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO

  N1 = tsgGetNumPoints(grid1)
  N2 = tsgGetNumPoints(grid2)

  WRITE(*,"(A,I5,E12.4)") " For rule_localp   Number of points: ", N1, err1
  WRITE(*,"(A,I5,E12.4)") " For rule_localp0  Number of points: ", N2, err2
  WRITE(*,*) " Note: localp-zero wins this competition because the function is zero at the boundary"
  WRITE(*,*)

! ==================================================================== !
! EXAMPLE 9:
! interpolate: f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))
! using different refinement schemes

! remake the true values for the next example
  DO i = 1, 1000
    tvalues(1,i) = exp(-randPoints2(1,i)) / (1.0 + 100.0 * exp(-10.0 * randPoints2(2,i)))
  END DO

  dims = 2
  outs = 1

  CALL tsgMakeLocalPolynomialGrid(grid1, dims, outs, 2, -1, tsg_localp)
  CALL tsgMakeLocalPolynomialGrid(grid2, dims, outs, 2, -1, tsg_localp)

  N = tsgGetNumNeeded(grid1)
  points => tsgGetNeededPoints(grid1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
  END DO
  CALL tsgLoadNeededPoints(grid1, values)
  CALL tsgLoadNeededPoints(grid2, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 9: interpolate f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))"
  WRITE(*,*) "   the error is estimated as the maximum from 1000 random points"
  WRITE(*,*) "   tolerance is set at 1.E-5 and maximal order polynomials are used"
  WRITE(*,*)

  WRITE(*,*) "             Classic              FDS"
  WRITE(*,*) "iteration  points     error     points     error"

  DO j = 1, 7
    ! 1 below corresponds to classic refinement
    CALL tsgSetLocalSurplusRefinement(grid1, 1.D-5, tsg_classic)

    N = tsgGetNumNeeded(grid1)
    points => tsgGetNeededPoints(grid1)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
    END DO
    CALL tsgLoadNeededPoints(grid1, values)
    DEALLOCATE(values)
    DEALLOCATE(points)

    CALL tsgEvaluateBatch(grid1, randPoints2, 1000, res2)
    err1 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
        err1 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO

    CALL tsgSetLocalSurplusRefinement(grid2, 1.D-5, tsg_fds)

    N = tsgGetNumNeeded(grid2)
    points => tsgGetNeededPoints(grid2)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
    END DO
    CALL tsgLoadNeededPoints(grid2, values)
    DEALLOCATE(values)
    DEALLOCATE(points)

    CALL tsgEvaluateBatch(grid2, randPoints2, 1000, res2)
    err2 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err2)then
        err2 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO

    N1 = tsgGetNumPoints(grid1)
    N2 = tsgGetNumPoints(grid2)

    WRITE(*,"(I9,I9,E12.4,I9,E12.4)") j, N1, err1, N2, err2
  END DO
  WRITE(*,*)

! ==================================================================== !
! EXAMPLE 10:
! interpolate: f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))
! using local polynomails and wavelets

  dims = 2
  outs = 1

  CALL tsgMakeLocalPolynomialGrid(grid1, dims, outs, 3, 1, tsg_localp)
  CALL tsgMakeWaveletGrid(grid2, dims, outs, 1, 1)

  N = tsgGetNumNeeded(grid1)
  points => tsgGetNeededPoints(grid1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
  END DO
  CALL tsgLoadNeededPoints(grid1, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  N = tsgGetNumNeeded(grid2)
  points => tsgGetNeededPoints(grid2)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
  END DO
  CALL tsgLoadNeededPoints(grid2, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 10: interpolate f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))"
  WRITE(*,*) "   the error is estimated as the maximum from 1000 random points"
  WRITE(*,*) "   using local polynomials and wavelets"
  WRITE(*,*)

  WRITE(*,*) "                Polynomials             Wavelets"
  WRITE(*,*) "iteration  points     error     points     error"

  DO j = 1, 8
    CALL tsgSetLocalSurplusRefinement(grid1, 1.D-5, tsg_fds)

    N = tsgGetNumNeeded(grid1)
    points => tsgGetNeededPoints(grid1)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
    END DO
    CALL tsgLoadNeededPoints(grid1, values)
    DEALLOCATE(values)
    DEALLOCATE(points)

    CALL tsgEvaluateBatch(grid1, randPoints2, 1000, res2)
    err1 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
        err1 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO

    CALL tsgSetLocalSurplusRefinement(grid2, 1.D-5, tsg_fds)

    N = tsgGetNumNeeded(grid2)
    points => tsgGetNeededPoints(grid2)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
    END DO
    CALL tsgLoadNeededPoints(grid2, values)
    DEALLOCATE(values)
    DEALLOCATE(points)

    CALL tsgEvaluateBatch(grid2, randPoints2, 1000, res2)
    err2 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err2)then
        err2 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO

    N1 = tsgGetNumPoints(grid1)
    N2 = tsgGetNumPoints(grid2)

    WRITE(*,"(I9,I9,E12.4,I9,E12.4)") j, N1, err1, N2, err2
  END DO
  WRITE(*,*)

  CALL tsgDeallocateGrid(grid1)
  CALL tsgDeallocateGrid(grid2)

! ==================================================================== !
! EXAMPLE 11: interpolate: f(x,y,z) = 1/((1+4x^2)*(1+5y^2)*(1+6z^2))
! using classical and conformal transformation

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 11: interpolate f(x,y,z) = 1/((1+4x^2)*(1+5y^2)*(1+6z^2))"
  WRITE(*,*) "            using conformal transformation"
  WRITE(*,*) "            the error is estimated as the maximum from 1000 random points"

  dims = 3
  outs = 1
  level = 12
  conformal(1) = 4
  conformal(2) = 4
  conformal(3) = 4

  DO i = 1, 1000
    randPoints3(1,i) = randPoints(1,i)
    randPoints3(2,i) = randPoints(2,i)
    randPoints3(3,i) = randPoints(3,i)
    tvalues(1,i) = 1.0 / ((1.0 + 4.0 * randPoints3(1,i)**2) * &
                          (1.0 + 5.0 * randPoints3(2,i)**2) * &
                          (1.0 + 6.0 * randPoints3(3,i)**2))
  END DO

  call tsgAllocateGrid(grid1)

  CALL tsgMakeGlobalGrid(grid1, dims, outs, level, tsg_iptotal, tsg_clenshaw_curtis)

  N = tsgGetNumNeeded(grid1)
  points => tsgGetNeededPoints(grid1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = 1.0 / ((1.0 + 4.0 * points(1,i)**2) * &
                          (1.0 + 5.0 * points(2,i)**2) * &
                          (1.0 + 6.0 * points(3,i)**2))
  END DO
  CALL tsgLoadNeededPoints(grid1, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  CALL tsgEvaluateBatch(grid1, randPoints3, 1000, res2)
  err1 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
      err1 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO
  N1 = tsgGetNumPoints(grid1)

  CALL tsgMakeGlobalGrid(grid1, dims, outs, level, tsg_iptotal, tsg_clenshaw_curtis)
  CALL tsgSetConformalTransformASIN(grid1, conformal)

  N = tsgGetNumNeeded(grid1)
  points => tsgGetNeededPoints(grid1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = 1.0 / ((1.0 + 4.0 * points(1,i)**2) * &
                          (1.0 + 5.0 * points(2,i)**2) * &
                          (1.0 + 6.0 * points(3,i)**2))
  END DO
  CALL tsgLoadNeededPoints(grid1, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  CALL tsgEvaluateBatch(grid1, randPoints3, 1000, res2)
  err2 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err2)then
      err2 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO

  CALL tsgMakeLocalPolynomialGrid(grid1, dims, outs, level-4, 2, tsg_localp)

  N = tsgGetNumNeeded(grid1)
  points => tsgGetNeededPoints(grid1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = 1.0 / ((1.0 + 4.0 * points(1,i)**2) * &
                          (1.0 + 5.0 * points(2,i)**2) * &
                          (1.0 + 6.0 * points(3,i)**2))
  END DO
  CALL tsgLoadNeededPoints(grid1, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  CALL tsgEvaluateBatch(grid1, randPoints3, 1000, res2)
  err3 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err3)then
      err3 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO
  N2 = tsgGetNumPoints(grid1)

  CALL tsgMakeLocalPolynomialGrid(grid1, dims, outs, level-4, 2, tsg_localp)
  CALL tsgSetConformalTransformASIN(grid1, conformal)

  N = tsgGetNumNeeded(grid1)
  points => tsgGetNeededPoints(grid1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = 1.0 / ((1.0 + 4.0 * points(1,i)**2) * &
                          (1.0 + 5.0 * points(2,i)**2) * &
                          (1.0 + 6.0 * points(3,i)**2))
  END DO
  CALL tsgLoadNeededPoints(grid1, values)
  DEALLOCATE(values)
  DEALLOCATE(points)

  CALL tsgEvaluateBatch(grid1, randPoints3, 1000, res2)
  err4 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err4)then
      err4 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO

  WRITE(*,*) "Grid Type    nodes     error regular   error conformal"
  WRITE(*,"(A,I8,E18.4,E18.4)") " Global    ", N1, err1, err2
  WRITE(*,"(A,I8,E18.4,E18.4)") " Localp    ", N2, err3, err4
  WRITE(*,*)
  WRITE(*,*) "Note: conformal maps address specific problems with the region of analyticity of a function"
  WRITE(*,*) "      the map can accelerate or slow down convergence depending on the problem"
  WRITE(*,*)

  CALL tsgDeallocateGrid(grid1)
! ==================================================================== !

! cleanup
  DEALLOCATE(res2)
  DEALLOCATE(tvalues)

  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*)

END PROGRAM TasmanianSGExample

