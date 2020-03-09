
subroutine example_sparse_grid_02()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    integer :: dimension = 2, exactness = 20
    integer :: i, num_points
    ! using the static array API requires only allocatable variables
    ! but the allocation has to be made before the call
    real(C_DOUBLE), dimension(:), allocatable :: weights
    real(C_DOUBLE), dimension(:,:), allocatable :: points
    ! pointer variables can be returned through the functional API
    ! allowing for more expressive calls
    real(C_DOUBLE), dimension(:), pointer :: pweights
    real(C_DOUBLE), dimension(:,:), pointer :: ppoints
    real(C_DOUBLE), dimension(2) :: lower, upper
    double precision :: exact = 1.861816427518323D+00;
    double precision :: integral, err

    write(*,*)
    write(*,*) "-------------------------------------------------------------------------------------------------"
    write(*,*) "Example 2: integrate f(x,y) = exp(-x^2) * cos(y) over [-5,5] x [-2,3]"
    write(*,*) "           using  Gauss-Patterson nodes and total degree polynomial space"

    grid = TasmanianGlobalGrid(dimension, 0, exactness, tsg_type_qptotal, tsg_rule_gausspatterson)
    num_points = grid%getNumPoints()

    lower = [ -5.0D-0, -2.0D-0 ]
    upper = [  5.0D-0,  3.0D-0 ]
    call grid%setDomainTransform(lower, upper)

    allocate(weights(num_points))
    call grid%getQuadratureWeights(weights)
    allocate(points(dimension, num_points))
    call grid%getPoints(points(:,1))

    integral = 0.0D-0
    do i = 1, num_points
        integral = integral + weights(i) * exp( -points(1,i)**2 ) * cos( points(2,i) )
    enddo
    err = abs(integral - exact)

    write(*,*) "     at polynomial exactness: ", exactness
    write(*,*) "     the grid has:            ", num_points
    write(*,*) "     integral: ", integral
    write(*,*) "     error:    ", err
    write(*,*)
    deallocate(weights, points)
    call grid%release()

    exactness = 40
    grid = TasmanianGlobalGrid(dimension, 0, exactness, tsg_type_qptotal, tsg_rule_gausspatterson)
    num_points = grid%getNumPoints()

    call grid%setDomainTransform(lower, upper)

    ! using the functional API, saves the need to manually allocate the sizes
    pweights => tsgGetQuadratureWeights(grid)
    ppoints => tsgGetPoints(grid)

    integral = 0.0D-0
    do i = 1, num_points
        integral = integral + pweights(i) * exp( -ppoints(1,i)**2 ) * cos( ppoints(2,i) )
    enddo
    err = abs(integral - exact)

    write(*,*) "     at polynomial exactness: ", exactness
    write(*,*) "     the grid has:            ", num_points
    write(*,*) "     integral: ", integral
    write(*,*) "     error:    ", err
    write(*,*)

    deallocate(pweights, ppoints)
    call grid%release()
end subroutine
