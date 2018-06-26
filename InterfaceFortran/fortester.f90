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
  USE TasmanianSG
  IMPLICIT NONE
  INTEGER :: verm, vern!, i
  CHARACTER, pointer :: licence(:)
  INTEGER :: gridID, gridID_II
  DOUBLE PRECISION :: tp1(2, 5), tw1(5), tp2(1, 4), tw2(4), tp3(2, 4)
  INTEGER :: n
  DOUBLE PRECISION, pointer :: points(:,:), weights(:), pointsb(:,:)
  DOUBLE PRECISION :: getErrorMat, getErrorVec, sumMat, sumVec
  INTEGER :: anisoWeights(4), levelLimits(4)
!  INTEGER :: N, i, 
!  DOUBLE PRECISION, pointer :: points(:,:), weights(:)
!  DOUBLE PRECISION :: x, y, integ, E
!  DOUBLE PRECISION, allocatable :: transformA(:), transformB(:), testTransA(:), testTransB(:)

! This is the sound "Glaucodon Ballaratensis" makes :)
!  WRITE(*,*) "Ghurrrrrphurrr"

! must call tsgInitialize() once per program
  CALL tsgInitialize()

! read library meta-data  
  verm = tsgGetVersionMajor()
  vern = tsgGetVersionMinor()
  licence => tsgGetLicense()
  WRITE(*,"(A,I2,A,I1)") "Tasmanian Sparse Grid module: ", verm, ".", vern
  WRITE(*,"(A,40A)") "Licence: ", licence

!=======================================================================
!       tsgMakeGlobalGrid()
!=======================================================================
  gridID = tsgNewGridID()
  CALL tsgMakeGlobalGrid(gridID, 2, 0, 1, tsg_level, tsg_clenshaw_curtis)
  tp1 = reshape((/0.0D+0, 0.0D+0, 0.0D+0, -1.0D+0, 0.0D+0, 1.0D+0, -1.0D+0, 0.0D+0, 1.0D+0, 0.0D+0/), shape(tp1))
  tw1 = reshape((/4.0D+0/3.0D+0, 2.0D+0/3.0D+0, 2.0D+0/3.0D+0, 2.0D+0/3.0D+0, 2.0D+0/3.0D+0/), shape(tw1))
  points  => tsgGetPoints(gridID)
  weights => tsgGetQuadratureWeights(gridID)
  IF ((getErrorMat(2, 5, tp1, points) > 1.0D-11) .OR. (getErrorVec(5, tw1, weights) > 1.0D-11)) THEN
    WRITE(*,*) "Mismatch in tsgMakeGlobal: core case 1", getErrorMat(2, 5, tp1, points), getErrorVec(5, tw1, weights)
    STOP 1
  END IF
  DEALLOCATE(points)
  DEALLOCATE(weights)
  
  CALL tsgMakeGlobalGrid(gridID, 2, 0, 2, tsg_level, tsg_clenshaw_curtis)
  points  => tsgGetPoints(gridID)
  weights => tsgGetQuadratureWeights(gridID)
  IF ((abs(points(2, 4) + sqrt(0.5D+0)) > 1.0D-11) .OR. &
      (abs(points(2, 5) - sqrt(0.5D+0)) > 1.0D-11) .OR. &
      (abs(weights(7) - 1.0D+0/9.0D+0)  > 1.0D-11))then
    WRITE(*,*) "Mismatch in tsgMakeGlobal: core case 2"
    WRITE(*,*) abs(points(2, 4) - sqrt(0.5D+0)), abs(points(2, 5) - sqrt(0.5D+0)), abs(weights(7) - 1.0D+0/9.0D+0)
    STOP 1
  END IF
  DEALLOCATE(points)
  DEALLOCATE(weights)
  
  CALL tsgMakeGlobalGrid(gridID, 3, 0, 4, tsg_level, tsg_fejer2)
  points  => tsgGetPoints(gridID)
  weights => tsgGetQuadratureWeights(gridID)
  n = tsgGetNumPoints(gridID)
  IF ((abs(sumMat(3, n, points)) > 1.0D-11) .OR. (abs(sumVec(n, weights) - 8.0D+0) > 1.0D-11))then
    WRITE(*,*) "Mismatch in tsgMakeGlobal: core case 3"
    WRITE(*,*) abs(sumMat(3, n, points)), abs(sumVec(n, weights) - 8.0D+0)
    STOP 1
  END IF
  DEALLOCATE(points)
  DEALLOCATE(weights)
  
  CALL tsgMakeGlobalGrid(gridID, 1, 0, 3, tsg_level, tsg_leja)
  tp2 = reshape((/0.0D+0, 1.0D+0, -1.0D+0, sqrt(1.0D+0/3.0D+0)/), shape(tp2))
  tw2 = reshape((/4.0D+0/3.0D+0, 1.0D+0/3.0D+0, 1.0D+0/3.0D+0, 0.0D+0/), shape(tw2))
  points  => tsgGetPoints(gridID)
  weights => tsgGetQuadratureWeights(gridID)
  IF ((getErrorMat(1, 4, tp2, points) > 1.0D-11) .OR. (getErrorVec(1, tw2, weights) > 1.0D-11))then
    WRITE(*,*) "Mismatch in tsgMakeGlobal: core case 4", getErrorMat(1, 4, tp2, points), getErrorVec(1, tw2, weights)
    STOP 1
  END IF
  DEALLOCATE(points)
  DEALLOCATE(weights)

  anisoWeights(1) = 1
  CALL tsgMakeGlobalGrid(gridID, 1, 0, 4, tsg_level, tsg_gauss_hermite, anisoWeights, 2.0D+0)
  weights => tsgGetQuadratureWeights(gridID)
  n = tsgGetNumPoints(gridID)
  IF ((abs(sumVec(n, weights) - 0.5D+0 * sqrt(4.0D+0 * atan(1.0D+0))) > 1.D-11) .OR. (abs(tsgGetBeta(gridID)) > 1.D-11))then
    WRITE(*,*) "Mismatch in tsgMakeGlobal: core case 5", sumVec(n, weights), tsgGetBeta(gridID)
    STOP 1
  END IF
  DEALLOCATE(weights)
  
  anisoWeights(1) = 2
  anisoWeights(2) = 1
  CALL tsgMakeGlobalGrid(gridID, 2, 0, 2, tsg_level, tsg_leja, anisoWeights)
  tp3 = reshape((/0.0D+0, 0.0D+0, 0.0D+0, 1.0D+0, 0.0D+0, -1.0D+0, 1.0D+0, 0.0D+0/), shape(tp3))
  points  => tsgGetPoints(gridID)
  weights => tsgGetQuadratureWeights(gridID)
  IF ((getErrorMat(2, 4, tp3, points) > 1.D-11) .OR. (abs(sumVec(4, weights) - 4.0D+0) > 1.D-11))then
    WRITE(*,*) "Mismatch in tsgMakeGlobal: core case 6", getErrorMat(2, 4, tp3, points), sumVec(4, weights)
    STOP 1
  END IF
  DEALLOCATE(points)
  DEALLOCATE(weights)
  
  gridID_II = tsgNewGridID()
  levelLimits(1)  = 1
  levelLimits(2)  = 3
  anisoWeights(1) = 1
  anisoWeights(2) = 1
  CALL tsgMakeGlobalGrid(gridID, 2, 0, 20, tsg_qptotal, tsg_clenshaw_curtis, anisoWeights, 0.0D+0, 0.0D+0, levelLimits)
  anisoWeights(1) = 1
  anisoWeights(2) = 3
  CALL tsgMakeGlobalGrid(gridID_II, 2, 0, 1, tsg_tensor, tsg_clenshaw_curtis, anisoWeights)
  points  => tsgGetPoints(gridID)
  pointsb => tsgGetPoints(gridID_II)
  n = tsgGetNumPoints(gridID)
  IF ((tsgGetNumPoints(gridID_II) .NE. n) .OR. (getErrorMat(2, n, points, pointsb) > 1.0D-11))then
    WRITE(*,*) "Mismatch in tsgMakeGlobal: core case 7", tsgGetNumPoints(gridID_II), getErrorMat(2, n, points, pointsb)
    STOP 1
  END IF
  DEALLOCATE(points)
  DEALLOCATE(pointsb)
  
  CALL tsgFreeGridID(gridID)
  CALL tsgFreeGridID(gridID_II)
  


! Tasmanian holds to some RAM until tsgFinalize() is called
  CALL tsgFinalize()

END PROGRAM FORTESTER

FUNCTION getErrorMat(m, n, x, y) result(error)
  INTEGER, intent(in) :: m, n
  INTEGER :: i, j
  DOUBLE PRECISION :: x(m, n), y(m, n)
  DOUBLE PRECISION :: error
  error = 0.0
  DO i = 1, m
    DO j = 1, n
      error = error + abs(x(i,j) - y(i,j))
    END DO
  END DO
END FUNCTION getErrorMat

FUNCTION getErrorVec(n, x, y) result(error)
  INTEGER, intent(in) :: n
  INTEGER :: i
  DOUBLE PRECISION :: x(n), y(n)
  DOUBLE PRECISION :: error
  error = 0.0
  DO i = 1, n
    error = error + abs(x(i) - y(i))
  END DO
END FUNCTION getErrorVec

FUNCTION sumMat(m, n, x) result(s)
  INTEGER, intent(in) :: m, n
  INTEGER :: i, j
  DOUBLE PRECISION :: x(m, n)
  DOUBLE PRECISION :: s
  s = 0.0
  DO i = 1, m
    DO j = 1, n
      s = s + x(i,j)
    END DO
  END DO
END FUNCTION sumMat

FUNCTION sumVec(n, x) result(s)
  INTEGER, intent(in) :: n
  INTEGER :: i
  DOUBLE PRECISION :: x(n)
  DOUBLE PRECISION :: s
  s = 0.0
  DO i = 1, n
    s = s + x(i)
  END DO
END FUNCTION sumVec
