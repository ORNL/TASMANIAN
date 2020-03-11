subroutine example_sparse_grid_04()
    use Tasmanian
    use, intrinsic :: iso_c_binding
    implicit none
    type(TasmanianSparseGrid) :: surrogate
    integer :: prec, i
    integer(C_INT), parameter :: num_inputs = 2
    integer(C_INT), parameter :: num_outputs = 1
    real(C_DOUBLE), dimension(:,:), allocatable :: points
    real(C_DOUBLE), dimension(:,:), allocatable :: values
    real(C_DOUBLE) :: x(2), y(1)

    write(*,*) "-------------------------------------------------------------------------------------------------"
    write(*,*) "Example 4: interpolate f(x,y) = exp(-x^2) * cos(y), using clenshaw-curtis iptotal rule"
    write(*,*)

    x = [ 0.3D-0, 0.7D-0 ]

    do prec = 6, 12, 6
        surrogate = TasmanianGlobalGrid(num_inputs, num_outputs, prec, &
                                        tsg_type_iptotal, tsg_rule_clenshawcurtis)

        allocate(points(num_inputs, surrogate%getNumNeeded()))
        call surrogate%getNeededPoints(points(:,1))
        allocate(values(num_outputs, surrogate%getNumNeeded()))

        do i = 1, surrogate%getNumNeeded()
            ! fill the values with the model values at the points
            values(1, i) = exp(-points(1,i)**2) * cos(points(2,i))
        enddo
        call surrogate%loadNeededPoints(values(:,1))
        ! after the load call, the surrogate model is ready

        ! let y = surrogate(x)
        call surrogate%evaluate(x, y)

        write(*,"(A,i2)") "  using polynomials of total degree up to: ", prec
        write(*,"(A,i3,A)") "                             the grid has: ", &
                                surrogate%getNumPoints(), " points"
        write(*,"(A,ES10.4)") "                 interpolant at (0.3,0.7): ", y(1)
        write(*,"(A,ES10.4)") "                                    error: ", &
                                abs(y(1) - exp(-x(1)**2) * cos(x(2)))
        write(*,*)

        call surrogate%release()
        deallocate(points, values)
    enddo
end subroutine
