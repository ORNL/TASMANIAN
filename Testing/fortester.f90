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
  INTEGER :: verm, vern
  CHARACTER, pointer :: licence(:)
  INTEGER :: gridID
  DOUBLE PRECISION :: tp1(2, 5)
  DOUBLE PRECISION, pointer :: points(:,:)
  DOUBLE PRECISION :: getError
!  INTEGER :: N, i, 
!  DOUBLE PRECISION, pointer :: points(:,:), weights(:)
!  DOUBLE PRECISION :: x, y, integ, E
!  DOUBLE PRECISION, allocatable :: transformA(:), transformB(:)

! This is the sound "Glaucodon Ballaratensis" makes :)
!  WRITE(*,*) "Ghurrrrrphurrr"


! must call tsgInitialize() once per program
  CALL tsgInitialize()
  
  verm = tsgGetVersionMajor()
  vern = tsgGetVersionMinor()
  licence => tsgGetLicense()
  WRITE(*,"(A,I2,A,I1)") "Tasmanian Sparse Grid module: ", verm, ".", vern
  WRITE(*,"(A,40A)") "Licence: ", licence
  
  gridID = tsgNewGridID()
  CALL tsgMakeGlobalGrid(gridID, 2, 1, 1, 1, tsg_clenshaw_curtis)
  tp1 = reshape((/ 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, -1.0, 0.0, 1.0, 0.0 /), shape(tp1))
  points => tsgGetPoints(gridID)
  IF (getError(2, 5, tp1, points) > 1.0E-11) THEN
    WRITE(*,*) "Mismatch in tsgMakeGlobal: core case 1", getError(2, 5, tp1, points)
    STOP 1
  END IF


! Tasmanian holds to some RAN until tsgFinalize() is called
  CALL tsgFinalize()

END PROGRAM FORTESTER

FUNCTION getError(m, n, x, y) result(error)
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
END FUNCTION getError
