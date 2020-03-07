
subroutine example_sparse_grid_01()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    type(TasmanianSparseGrid) :: grid
    integer :: dimension = 2, level = 6
    integer :: i, num_points
    real(C_DOUBLE), dimension(:), allocatable :: weights
    real(C_DOUBLE), dimension(:,:), allocatable :: points
    double precision :: exact = 2.513723354063905D+0;
    double precision :: integral, err

    write(*,*)
    write(*,*) "-------------------------------------------------------------------------------------------------"
    write(*,*) "Example 1:  integrate f(x,y) = exp(-x^2) * cos(y)"
    write(*,*) "            using clenshaw-curtis nodes and grid of type level"

    grid = TasmanianSparseGrid()
    call grid%makeGlobalGrid(dimension, 0, level, tsg_type_level, tsg_rule_clenshawcurtis)
    num_points = grid%getNumPoints()

    allocate(weights(num_points))
    call grid%getQuadratureWeights(weights)
    allocate(points(dimension, num_points))
    call grid%getPoints(points(:,1))

    integral = 0.0D-0
    do i = 1, num_points
        integral = integral + weights(i) * exp( -points(1,i)**2 ) * cos( points(2,i) )
    enddo
    err = abs(integral - exact)

    write(*,*) "     at level:     ", level
    write(*,*) "     the grid has: ", num_points
    write(*,*) "     integral: ", integral
    write(*,*) "     error:    ", err
    write(*,*)
    deallocate(weights, points)

    level = 7
    call grid%makeGlobalGrid(dimension, 0, level, tsg_type_level, tsg_rule_clenshawcurtis)
    num_points = grid%getNumPoints()

    allocate(weights(num_points))
    call grid%getQuadratureWeights(weights)
    allocate(points(dimension, num_points))
    call grid%getPoints(points(:,1))

    integral = 0.0D-0
    do i = 1, num_points
        integral = integral + weights(i) * exp( -points(1,i)**2 ) * cos( points(2,i) )
    enddo
    err = abs(integral - exact)

    write(*,*) "     at level:     ", level
    write(*,*) "     the grid has: ", num_points
    write(*,*) "     integral: ", integral
    write(*,*) "     error:    ", err
    write(*,*)

    deallocate(weights, points)
    call grid%release()
end subroutine
