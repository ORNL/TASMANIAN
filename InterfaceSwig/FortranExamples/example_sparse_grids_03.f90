subroutine example_sparse_grid_03()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid_cc, grid_gl, grid_gp
    integer :: prec
    real(C_DOUBLE) :: r
    real(C_DOUBLE), external :: ex3_get_error

    write(*,*)
    write(*,*) "-------------------------------------------------------------------------------------------------"
    write(*,*) "Example 3: integrate exp(-x1^2 - x2^2) * cos(x3) * cos(x4)"
    write(*,*) "           for x1, x2 in [-5,5]; x3, x4 in [-2,3];"
    write(*,*) "           using different rules and total degree polynomial space"
    write(*,*)

    write(*,*) "               Clenshaw-Curtis      Gauss-Legendre    Gauss-Patterson"
    write(*,*) " precision    points     error    points     error    points    error"

    do prec = 5, 40, 5
        grid_cc = TasmanianSparseGrid()
        grid_gl = TasmanianSparseGrid()
        grid_gp = TasmanianSparseGrid()

        call ex3_make_grid(grid_cc, prec, tsg_rule_clenshawcurtis)
        call ex3_make_grid(grid_gl, prec, tsg_rule_gausslegendreodd)
        call ex3_make_grid(grid_gp, prec, tsg_rule_gausspatterson)

        write(*,"(i10,i10,ES10.2,i10,ES10.2,i10,ES10.2)")  prec, &
            grid_cc%getNumPoints(), ex3_get_error(grid_cc), &
            grid_gl%getNumPoints(), ex3_get_error(grid_gl), &
            grid_gp%getNumPoints(), ex3_get_error(grid_gp)

        call grid_cc%release()
        call grid_gl%release()
        call grid_gp%release()
    enddo

    write(*,*)
    write(*,*) "see the example comments in the on-line manual (or the equivalent C++ example)"
    write(*,*)

end subroutine

function ex3_get_error(grid) result(err)
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    real(C_DOUBLE) :: err
    type(TasmanianSparseGrid), intent(in) :: grid
    real(C_DOUBLE), dimension(:), pointer :: weights
    real(C_DOUBLE), dimension(:,:), pointer :: points
    real(C_DOUBLE) :: integ
    integer :: i

    weights => tsgGetQuadratureWeights(grid)
    points  => tsgGetPoints(grid)

    integ = 0.0D-0
    do i = 1, grid%getNumPoints()
        integ = integ + weights(i) * exp(-points(1,i)**2 -points(2,i)**2) &
                                   * cos(points(3,i)) * cos(points(4,i))
    enddo
    err = abs(1.861816427518323D+00 * 1.861816427518323D+00 - integ)

    deallocate(weights, points)
end function

subroutine ex3_make_grid(grid, prec, rule)
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: grid
    integer(C_INT), intent(in) :: prec, rule
    real(C_DOUBLE) :: lower(4), upper(4)

    grid = TasmanianGlobalGrid(4, 0, prec, tsg_type_qptotal, rule)

    lower = [ -5.0D-0, -5.0D-0, -2.0D-0, -2.0D-0 ]
    upper = [  5.0D-0,  5.0D-0,  3.0D-0,  3.0D-0 ]
    call grid%setDomainTransform(lower, upper)
end subroutine
