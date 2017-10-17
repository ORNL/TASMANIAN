PROGRAM TasmanianSGExample
  USE TasmanianSG, ONLY: tsgInitialize, tsgFinalize, tsgNewGridID, tsgFreeGridID, &
       tsgGetVersionMajor, tsgGetVersionMinor, tsgGetLicense, &
       tsgMakeGlobalGrid, tsgMakeSequenceGrid, tsgMakeLocalPolynomialGrid, tsgMakeWaveletGrid, &
       tsgUpdateGlobalGrid, tsgUpdateSequenceGrid, tsgRead, tsgWrite, &
       tsgGetAlpha, tsgGetBeta, tsgGetOrder, tsgGetNumDimensions, tsgGetNumOutputs, tsgGetRule, &
       tsgGetNumLoaded, tsgGetNumNeeded, tsgGetNumPoints, &
       tsgGetLoadedPoints, tsgGetNeededPoints, tsgGetPoints, &
       tsgGetLoadedPointsStatic, tsgGetNeededPointsStatic, tsgGetPointsStatic, &
       tsgLoadNeededPoints, tsgEvaluate, tsgEvaluateFast, tsgEvaluateBatch, tsgIntegrate, &
       tsgGetQuadratureWeights, tsgGetQuadratureWeightsStatic, &
       tsgGetInterpolationWeights, tsgGetInterpolationWeightsStatic, &
       tsgSetDomainTransform, tsgIsSetDomainTransfrom, tsgClearDomainTransform, tsgGetDomainTransform, &
       tsgSetAnisotropicRefinement, tsgSetGlobalSurplusRefinement, tsgSetLocalSurplusRefinement, tsgClearRefinement
IMPLICIT NONE
  INTEGER :: gridID, dims, outs, level
  INTEGER :: gridID1, gridID2, gridID3, N1, N2, N3
  INTEGER :: N, i, j, verm, vern
  DOUBLE PRECISION :: err1, err2, err3, exact
  REAL :: cpuStart, cpuEnd, stages(2,3)
  !INTEGER :: aweights(3)
  CHARACTER, pointer :: string(:)
  DOUBLE PRECISION, pointer :: points(:,:), weights(:)
  DOUBLE PRECISION :: x, y, integ, E
  DOUBLE PRECISION, allocatable :: transformA(:), transformB(:), values(:,:), tvalues(:,:)
  DOUBLE PRECISION, allocatable :: res(:), res2(:,:)
  DOUBLE PRECISION :: desired_x(2)
  DOUBLE PRECISION :: randPoints(4,1000), randPoints2(2,1000)
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

! ==================================================================== !  
! EXAMPLE 1: integrate: f(x,y) = exp(-x^2) * cos(y) over [-1,1] x [-1,1]
! using classical Smolyak grid with Clenshaw-Curtis points and weights

! must call tsgInitialize() once per program
  CALL tsgInitialize()
  
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
  
! before you use a grid, you must ask for a new valid grid ID
  gridID = tsgNewGridID()
  
! clenshaw-curtis = 1, type_level = 1
  CALL tsgMakeGlobalGrid(gridID, dims, 0, level, 1, 1)
  
  points => tsgGetPoints(gridID)
  weights => tsgGetQuadratureWeights(gridID)
  
  N = tsgGetNumPoints(gridID)
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
  CALL tsgMakeGlobalGrid(gridID, dims, 0, level, 1, 1)
  
! do not forget to release the memory associated with points and weights
  DEALLOCATE(points)
  DEALLOCATE(weights)
  points => tsgGetPoints(gridID)
  weights => tsgGetQuadratureWeights(gridID)
  
  N = tsgGetNumPoints(gridID)
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
  
! after calling tsgFreeGridID(), we can no longer use this gridID
! until we call tsgNewGridID()
  CALL tsgFreeGridID(gridID)

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
  
! need new gridID, since we freed this earlier
  gridID = tsgNewGridID()
  
! gauss-patterson = 7, type_qptotal = 5
  CALL tsgMakeGlobalGrid(gridID, dims, 0, level, 5, 7)
  CALL tsgSetDomainTransform(gridID, transformA, transformB)
  
  points => tsgGetPoints(gridID)
  weights => tsgGetQuadratureWeights(gridID)
  
  N = tsgGetNumPoints(gridID)
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
  CALL tsgMakeGlobalGrid(gridID, dims, 0, level, 5, 7)
  CALL tsgSetDomainTransform(gridID, transformA, transformB)
  
! do not forget to release the memory associated with points and weights
  DEALLOCATE(points)
  DEALLOCATE(weights)
  points => tsgGetPoints(gridID)
  weights => tsgGetQuadratureWeights(gridID)
  
  N = tsgGetNumPoints(gridID)
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
  CALL tsgFreeGridID(gridID)

! ==================================================================== !
! EXAMPLE 3: integrate: f(x,y) = exp(-x^2) * cos(y) 
!                       over (x,y) in [-5,5] x [-2,3]
! using different rules

  gridID1 = tsgNewGridID()
  gridID2 = tsgNewGridID()
  gridID3 = tsgNewGridID()
  
  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 3: integrate f(x,y) = exp(-x^2) * cos(y) over [-5,5] x [-2,3] using different rules"
  WRITE(*,*)
  
  WRITE(*,*) "               Clenshaw-Curtis      Gauss-Legendre     Gauss-Patterson"
  WRITE(*,*) " precision    points     error    points     error    points     error"
  
  DO level = 9, 30, 4
    ! clenshaw-curtis = 1, gauss-legendre = 5, gauss-patterson = 7
    ! type_qptotal = 5
    CALL tsgMakeGlobalGrid(gridID1, dims, 0, level, 5, 1)
    CALL tsgSetDomainTransform(gridID1, transformA, transformB)
    CALL tsgMakeGlobalGrid(gridID2, dims, 0, level, 5, 5)
    CALL tsgSetDomainTransform(gridID2, transformA, transformB)
    CALL tsgMakeGlobalGrid(gridID3, dims, 0, level, 5, 7)
    CALL tsgSetDomainTransform(gridID3, transformA, transformB)
    
    points => tsgGetPoints(gridID1)
    weights => tsgGetQuadratureWeights(gridID1)
    N1 = tsgGetNumPoints(gridID1)
    integ = 0.0
    DO i = 1, N1
      x = points(1, i)
      y = points(2, i)
      integ = integ + weights(i) * exp(-x*x) * cos(y)
    END DO
    err1 = abs(integ - 1.861816427518323D+00)
    DEALLOCATE(points)
    DEALLOCATE(weights)
    
    points => tsgGetPoints(gridID2)
    weights => tsgGetQuadratureWeights(gridID2)
    N2 = tsgGetNumPoints(gridID2)
    integ = 0.0
    DO i = 1, N2
      x = points(1, i)
      y = points(2, i)
      integ = integ + weights(i) * exp(-x*x) * cos(y)
    END DO
    err2 = abs(integ - 1.861816427518323D+00)
    DEALLOCATE(points)
    DEALLOCATE(weights)
    
    points => tsgGetPoints(gridID3)
    weights => tsgGetQuadratureWeights(gridID3)
    N3 = tsgGetNumPoints(gridID3)
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
  
  CALL tsgFreeGridID(gridID1)
  CALL tsgFreeGridID(gridID2)
  CALL tsgFreeGridID(gridID3)
  
  DEALLOCATE(transformA)
  DEALLOCATE(transformB)


! ==================================================================== !
! EXAMPLE 4: interpolate: f(x,y) = exp(-x^2) * cos(y)
! with a rule that exactly interpolates polynomials of total degree

  gridID = tsgNewGridID()
  
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
  
  ! iptotal = 3, clenshaw-curtis = 1
  CALL tsgMakeGlobalGrid(gridID, dims, outs, level, 3, 1)
  
  N = tsgGetNumNeeded(gridID)
  points => tsgGetNeededPoints(gridID)
  ALLOCATE(values(outs,N))
  
  DO i = 1, N
    x = points(1, i)
    y = points(2, i)
    values(1,i) = exp(-x**2) * cos(y)
  END DO
  
  CALL tsgLoadNeededPoints(gridID, values)
  DEALLOCATE(values)
  DEALLOCATE(points)
  
  ALLOCATE(res(outs)) ! will DEALLOCATE later
  CALL tsgEvaluate(gridID, desired_x, res)
  E = abs(res(1) - exact)
  
  WRITE(*,"(A,I4)")   "  using polynomials of total degree:  ", level
  WRITE(*,"(A,I4,A)") "      the grid has:                   ", N, " points"
  WRITE(*,"(A,E25.16)") "      interpolant at (0.3,0.7):    ", integ
  WRITE(*,"(A,E25.16)") "      error:                       ", E
  WRITE(*,*)
  
  ! do the same with level = 12
  level = 12
  
  ! iptotal = 3, clenshaw-curtis = 1
  CALL tsgMakeGlobalGrid(gridID, dims, outs, level, 3, 1)
  
  N = tsgGetNumNeeded(gridID)
  points => tsgGetNeededPoints(gridID)
  ALLOCATE(values(outs,N))
  
  DO i = 1, N
    x = points(1, i)
    y = points(2, i)
    values(1,i) = exp(-x**2) * cos(y)
  END DO
  
  CALL tsgLoadNeededPoints(gridID, values)
  DEALLOCATE(values)
  DEALLOCATE(points)
  
  CALL tsgEvaluate(gridID, desired_x, res)
  E = abs(res(1) - exact)
  DEALLOCATE(res)
  
  WRITE(*,"(A,I4)")   "  using polynomials of total degree:  ", level
  WRITE(*,"(A,I4,A)") "      the grid has:                   ", N, " points"
  WRITE(*,"(A,E25.16)") "      interpolant at (0.3,0.7):    ", integ
  WRITE(*,"(A,E25.16)") "      error:                       ", E
  WRITE(*,*)

! ==================================================================== !
! prepare random smaples for future tests
  call srand(TIME())

  DO i = 1, 1000
    DO j = 1, 4
      randPoints(j,i) = rand()
    END DO
  END DO

!! ==================================================================== !
!! EXAMPLE 5:
!! interpolate: f(x1,x2,x3,x4) = exp(-x1^2) * cos(x2) * exp(-x3^2) * cos(x4)
!! with Global and Sequence Leja rules
!  
!  dims = 4
!  outs = 1
!  level = 15
!  
!  ! 8: rule leja,    1: type level
!  CALL cpu_time(cpuStart)
!  CALL tsgMakeGlobalGrid(gridID, dims, outs, level, 1, 8)
!  CALL cpu_time(cpuEnd)
!  stages(1,1) = cpuEnd - cpuStart
!  
!  N = tsgGetNumPoints(gridID)
!  
!  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
!  WRITE(*,*) "Example 5: interpolate f(x1,x2,x3,x4) = exp(-x1^2) * cos(x2) * exp(-x3^2) * cos(x4)"
!  WRITE(*,*) "       comparign the performance of Global and Sequence grids with leja nodes"
!  WRITE(*,"(A,I4)") "       using polynomials of total degree up to: ", level
!  WRITE(*,"(A,I4,A)") "      the grids have:                   ", N, " points"
!  WRITE(*,*) "       both grids are evaluated at 1000 random points "
!  WRITE(*,*)
!  
!  points => tsgGetNeededPoints(gridID)
!  ALLOCATE(values(outs,N))
!  DO i = 1, N
!    values(1,i) = exp(-points(1,i)**2) * cos(points(2,i)) * exp(-points(3,i)**2) * cos(points(4,i))
!  END DO
!  DEALLOCATE(points)
!  
!  CALL cpu_time(cpuStart)
!  CALL tsgLoadNeededPoints(gridID, values)
!  CALL cpu_time(cpuEnd)
!  stages(1,2) = cpuEnd - cpuStart
!  
!  ALLOCATE(res2(outs,1000)) ! 2-D result
!  
!  CALL cpu_time(cpuStart)
!  CALL tsgEvaluateBatch(gridID, randPoints, 1000, res2)
!  CALL cpu_time(cpuEnd)
!  stages(1,3) = cpuEnd - cpuStart
!  
!  ! 8: rule leja,    1: type level
!  CALL cpu_time(cpuStart)
!  CALL tsgMakeSequenceGrid(gridID, dims, outs, level, 1, 8)
!  CALL cpu_time(cpuEnd)
!  stages(2,1) = cpuEnd - cpuStart
!  
!  ! points are the same, no need to recompue values
!  CALL cpu_time(cpuStart)
!  CALL tsgLoadNeededPoints(gridID, values)
!  CALL cpu_time(cpuEnd)
!  stages(2,2) = cpuEnd - cpuStart
!  
!  DEALLOCATE(values)
!  
!  CALL cpu_time(cpuStart)
!  CALL tsgEvaluateBatch(gridID, randPoints, 1000, res2)
!  CALL cpu_time(cpuEnd)
!  stages(2,3) = cpuEnd - cpuStart
!  
!  WRITE(*,*) "Stage        Global Grid      Sequence Grid"
!  WRITE(*,"(A,E20.8,E20.8)") " make grid  ", stages(1,1), stages(2,1)
!  WRITE(*,"(A,E20.8,E20.8)") " load needed  ", stages(1,2), stages(2,2)
!  WRITE(*,"(A,E20.8,E20.8)") " evaluate  ", stages(1,3), stages(2,3)
!  WRITE(*,*) "WARNING: I have not figured out how to time execution under Fortran"
!  
!  CALL tsgFreeGridID(gridID)
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
  
  gridID1 = tsgNewGridID()
  gridID2 = tsgNewGridID()
  gridID3 = tsgNewGridID()
  
  dims = 2
  outs = 1
  ! 8: rule leja,    1: type level
  CALL tsgMakeGlobalGrid(gridID1, dims, outs, 3, 1, 8)
  CALL tsgMakeGlobalGrid(gridID2, dims, outs, 3, 1, 8)
  CALL tsgMakeGlobalGrid(gridID3, dims, outs, 3, 1, 8)

  N = tsgGetNumNeeded(gridID1)
  points => tsgGetNeededPoints(gridID1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
  END DO
  CALL tsgLoadNeededPoints(gridID1, values)
  CALL tsgLoadNeededPoints(gridID2, values)
  CALL tsgLoadNeededPoints(gridID3, values)
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
    CALL tsgSetAnisotropicRefinement(gridID1, 3, 10, 0)
  
    N = tsgGetNumNeeded(gridID1)
    points => tsgGetNeededPoints(gridID1)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
    END DO
    CALL tsgLoadNeededPoints(gridID1, values)
    DEALLOCATE(values)
    DEALLOCATE(points)
    
    CALL tsgEvaluateBatch(gridID1, randPoints2, 1000, res2)
    err1 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
        err1 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO
    
    CALL tsgSetAnisotropicRefinement(gridID2, 4, 10, 0)
  
    N = tsgGetNumNeeded(gridID2)
    points => tsgGetNeededPoints(gridID2)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
    END DO
    CALL tsgLoadNeededPoints(gridID2, values)
    DEALLOCATE(values)
    DEALLOCATE(points)
    
    CALL tsgEvaluateBatch(gridID2, randPoints2, 1000, res2)
    err2 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err2)then
        err2 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO
    
    CALL tsgSetGlobalSurplusRefinement(gridID3, 1.D-10, 0)
  
    N = tsgGetNumNeeded(gridID3)
    points => tsgGetNeededPoints(gridID3)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
    END DO
    CALL tsgLoadNeededPoints(gridID3, values)
    DEALLOCATE(values)
    DEALLOCATE(points)
    
    CALL tsgEvaluateBatch(gridID3, randPoints2, 1000, res2)
    err3 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err3)then
        err3 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO
    
    N1 = tsgGetNumPoints(gridID1)
    N2 = tsgGetNumPoints(gridID2)
    N3 = tsgGetNumPoints(gridID3)
    
    WRITE(*,"(I9,I9,E12.4,I9,E12.4,I9,E12.4)") j, N1, err1, N2, err2, N3, err3
    
  END DO
  WRITE(*,*)
  CALL tsgFreeGridID(gridID1)
  CALL tsgFreeGridID(gridID2)
  CALL tsgFreeGridID(gridID3)

! ==================================================================== !
! EXAMPLE 7:
! interpolate: f(x,y) = exp(-x^2) * cos(y)
! using localp and semilocalp grids
  
  WRITE(*,*) "-------------------------------------------------------------------------------------------------"
  WRITE(*,*) "Example 7: interpolate: f(x,y) = exp(-x^2) * cos(y)"
  WRITE(*,*) "       using localp and semi-localp rules with depth 7"
  WRITE(*,*) "       the error is estimated as the maximum from 1000 random points"
  WRITE(*,*)
  
  ! localp: 37, semi-localp: 39
  gridID1 = tsgNewGridID()
  gridID2 = tsgNewGridID()
  
  dims = 2
  outs = 1
  
  CALL tsgMakeLocalPolynomialGrid(gridID1, dims, outs, 7, 2, 37)
  CALL tsgMakeLocalPolynomialGrid(gridID2, dims, outs, 7, 2, 39)
  
  N = tsgGetNumNeeded(gridID1)
  points => tsgGetNeededPoints(gridID1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
  END DO
  CALL tsgLoadNeededPoints(gridID1, values)
  DEALLOCATE(values)
  DEALLOCATE(points)
  
  CALL tsgEvaluateBatch(gridID1, randPoints2, 1000, res2)
  err1 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
      err1 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO
  
  points => tsgGetNeededPoints(gridID2)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)**2) * cos(points(2,i))
  END DO
  CALL tsgLoadNeededPoints(gridID2, values)
  DEALLOCATE(values)
  DEALLOCATE(points)
  
  CALL tsgEvaluateBatch(gridID2, randPoints2, 1000, res2)
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
  
  ! localp: 37, localp-zero: 38
  dims = 2
  outs = 1
  
  CALL tsgMakeLocalPolynomialGrid(gridID1, dims, outs, 7, 2, 37)
  CALL tsgMakeLocalPolynomialGrid(gridID2, dims, outs, 6, 2, 38)
  
  N = tsgGetNumNeeded(gridID1)
  points => tsgGetNeededPoints(gridID1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = cos(0.5 * PI * points(1,i)) &
                * cos(0.5 * PI * points(2,i))
  END DO
  CALL tsgLoadNeededPoints(gridID1, values)
  DEALLOCATE(values)
  DEALLOCATE(points)
  
  CALL tsgEvaluateBatch(gridID1, randPoints2, 1000, res2)
  err1 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
      err1 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO
  
  N = tsgGetNumNeeded(gridID2)
  points => tsgGetNeededPoints(gridID2)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = cos(0.5 * PI * points(1,i)) * cos(0.5 * PI * points(2,i))
  END DO
  CALL tsgLoadNeededPoints(gridID2, values)
  DEALLOCATE(values)
  DEALLOCATE(points)
  
  CALL tsgEvaluateBatch(gridID2, randPoints2, 1000, res2)
  err2 = 0.0
  DO i = 1, 1000
    IF(abs(res2(1,i) - tvalues(1,i)) .GT. err2)then
      err2 = abs(res2(1,i) - tvalues(1,i))
    ENDIF
  END DO
  
  N1 = tsgGetNumPoints(gridID1)
  N2 = tsgGetNumPoints(gridID2)
  
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
  
  ! localp: 37,  classic: 1, FDS: 4
  CALL tsgMakeLocalPolynomialGrid(gridID1, dims, outs, 2, -1, 37)
  CALL tsgMakeLocalPolynomialGrid(gridID2, dims, outs, 2, -1, 37)
  
  N = tsgGetNumNeeded(gridID1)
  points => tsgGetNeededPoints(gridID1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
  END DO
  CALL tsgLoadNeededPoints(gridID1, values)
  CALL tsgLoadNeededPoints(gridID2, values)
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
    CALL tsgSetLocalSurplusRefinement(gridID1, 1.D-5, 1)
    
    N = tsgGetNumNeeded(gridID1)
    points => tsgGetNeededPoints(gridID1)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
    END DO
    CALL tsgLoadNeededPoints(gridID1, values)
    DEALLOCATE(values)
    DEALLOCATE(points)
    
    CALL tsgEvaluateBatch(gridID1, randPoints2, 1000, res2)
    err1 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
        err1 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO
    
    CALL tsgSetLocalSurplusRefinement(gridID2, 1.D-5, 4)
  
    N = tsgGetNumNeeded(gridID2)
    points => tsgGetNeededPoints(gridID2)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
    END DO
    CALL tsgLoadNeededPoints(gridID2, values)
    DEALLOCATE(values)
    DEALLOCATE(points)
    
    CALL tsgEvaluateBatch(gridID2, randPoints2, 1000, res2)
    err2 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err2)then
        err2 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO
    
    N1 = tsgGetNumPoints(gridID1)
    N2 = tsgGetNumPoints(gridID2)
    
    WRITE(*,"(I9,I9,E12.4,I9,E12.4)") j, N1, err1, N2, err2
  END DO
  WRITE(*,*)
  
! ==================================================================== !
! EXAMPLE 10:
! interpolate: f(x,y) = exp(-x) / (1 + 100 * exp(-10 * y))
! using local polynomails and wavelets

  dims = 2
  outs = 1
  
  ! localp: 37
  CALL tsgMakeLocalPolynomialGrid(gridID1, dims, outs, 3, 1, 37)
  CALL tsgMakeWaveletGrid(gridID2, dims, outs, 1, 1)
  
  N = tsgGetNumNeeded(gridID1)
  points => tsgGetNeededPoints(gridID1)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
  END DO
  CALL tsgLoadNeededPoints(gridID1, values)
  DEALLOCATE(values)
  DEALLOCATE(points)
  
  N = tsgGetNumNeeded(gridID2)
  points => tsgGetNeededPoints(gridID2)
  ALLOCATE(values(outs,N))
  DO i = 1, N
    values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
  END DO
  CALL tsgLoadNeededPoints(gridID2, values)
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
    CALL tsgSetLocalSurplusRefinement(gridID1, 1.D-5, 4)
  
    N = tsgGetNumNeeded(gridID1)
    points => tsgGetNeededPoints(gridID1)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
    END DO
    CALL tsgLoadNeededPoints(gridID1, values)
    DEALLOCATE(values)
    DEALLOCATE(points)
    
    CALL tsgEvaluateBatch(gridID1, randPoints2, 1000, res2)
    err1 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err1)then
        err1 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO
    
    CALL tsgSetLocalSurplusRefinement(gridID2, 1.D-5, 4)
  
    N = tsgGetNumNeeded(gridID2)
    points => tsgGetNeededPoints(gridID2)
    ALLOCATE(values(outs,N))
    DO i = 1, N
      values(1,i) = exp(-points(1,i)) / (1.0 + 100.0 * exp(-10.0 * points(2,i)))
    END DO
    CALL tsgLoadNeededPoints(gridID2, values)
    DEALLOCATE(values)
    DEALLOCATE(points)
    
    CALL tsgEvaluateBatch(gridID2, randPoints2, 1000, res2)
    err2 = 0.0
    DO i = 1, 1000
      IF(abs(res2(1,i) - tvalues(1,i)) .GT. err2)then
        err2 = abs(res2(1,i) - tvalues(1,i))
      ENDIF
    END DO
    
    N1 = tsgGetNumPoints(gridID1)
    N2 = tsgGetNumPoints(gridID2)
    
    WRITE(*,"(I9,I9,E12.4,I9,E12.4)") j, N1, err1, N2, err2
  END DO
  WRITE(*,*)
  
  
  DEALLOCATE(res2)
  DEALLOCATE(tvalues)

! Add Example 11




! no need to free grid IDs before tsgFinalize(),
! all memory will be freed regardless
! must deallocate points, weights, etc.
  CALL tsgFinalize()

END PROGRAM TasmanianSGExample

