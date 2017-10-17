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
PROGRAM FORTESTER
  USE TasmanianSG, ONLY: tsgInitialize, tsgFinalize, tsgNewGridID, tsgFreeGridID, &
       tsgGetVersionMajor, tsgGetVersionMinor, &
       tsgMakeGlobalGrid, tsgMakeSequenceGrid, tsgMakeLocalPolynomialGrid, tsgMakeWaveletGrid, &
       tsgGetAlpha, tsgGetBeta, tsgGetOrder, tsgGetNumDimensions, tsgGetNumOutputs, tsgGetRule, &
       tsgGetNumLoaded, tsgGetNumNeeded, tsgGetNumPoints, &
       tsgGetLoadedPoints, tsgGetNeededPoints, tsgGetPoints, &
       tsgGetLoadedPointsStatic, tsgGetNeededPointsStatic, tsgGetPointsStatic, &
       tsgGetQuadratureWeights, tsgGetQuadratureWeightsStatic, &
       tsgGetInterpolationWeights, tsgGetInterpolationWeightsStatic
IMPLICIT NONE
  INTEGER :: gridID, dims, level
  INTEGER :: N, i, verm, vern
  INTEGER :: aweights(3)
  DOUBLE PRECISION :: alpha, beta
  DOUBLE PRECISION, pointer :: points(:,:), weights(:)
  DOUBLE PRECISION :: x, y, integ, E
  DOUBLE PRECISION, allocatable :: transformA(:), transformB(:)

! This is the sound "Glaucodon Ballaratensis" makes :)
!  WRITE(*,*) "Ghurrrrrphurrr"

  
! EXAMPLE 1: integrate: f(x,y) = exp(-x^2) * cos(y) over [-1,1] x [-1,1]
! using classical Smolyak grid with Clenshaw-Curtis points and weights

! must call tsgInitialize() once per program
  CALL tsgInitialize()
  
  verm = tsgGetVersionMajor()
  vern = tsgGetVersionMinor()
  WRITE(*,*) "Tasmanian Sparse Grid module: ", verm, vern
  
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
  WRITE(*,*) "      at level:       ", level
  WRITE(*,*) "      the grid has:   ", N, " points"
  WRITE(*,*) "      integral:     ", integ
  WRITE(*,*) "      error:        ", E
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
  WRITE(*,*) "      at level:       ", level
  WRITE(*,*) "      the grid has:   ", N, " points"
  WRITE(*,*) "      integral:     ", integ
  WRITE(*,*) "      error:        ", E
  WRITE(*,*)
  
  DEALLOCATE(points)
  DEALLOCATE(weights)
  
! after calling tsgFreeGridID(), we can no longer use this gridID
! until we call tsgNewGridID()
  CALL tsgFreeGridID(gridID)

! ==================================================================== !
! EXAMPLE 2: integrate: f(x,y) = exp(-x^2) * cos(y) over (x,y) in [-5,5] x [-2,3]
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
  
! clenshaw-curtis = 1, type_level = 1
  CALL tsgMakeGlobalGrid(gridID, dims, 0, level, 1, 1)
  

  CALL tsgFinalize()

END PROGRAM FORTESTER

