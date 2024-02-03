program P1

    implicit none
   
    integer :: i, j, fu
    ! character(len=*), parameter :: OUT_FILE = 'data-p1.txt' ! Output file.
    character(len=*), parameter :: PLT_SAVE = 'plot-p1.plt'   ! Save
    
    double precision :: x_vec(100) = 0
    double precision :: y_vec(100) = 0
    double precision :: z_vec(100, 100) = 0
    

    ! Computing the sin(2x) fucntion
    do i = 1, 100
        x_vec(i) = dble(i) / 20.0
        y_vec(i) = sin(2*x_vec(i))
    end do

    ! print *, x_vec
    
    
    ! Writing results in ".txt" files
    open (action='write', file="data-p1-sin.txt", newunit=fu, status='replace')
    do i = 1, 100
        write (fu, '(F10.5 F10.5)') x_vec(i), y_vec(i)
    end do
    close(fu)

    ! Now we create a 3D plot of tanh(x*y) with nestded loops
    open (action='write', file="data-p1-surface.txt", newunit=fu, status='replace')
    do i = 1, 100
        x_vec(i) = dble(i) / 20.0 - 2.5
        do j = 1, 100
            y_vec(j) = dble(j) / 20.0 - 2.5
            z_vec(i,j) = tanh(x_vec(i)*y_vec(j))

            write (fu, '(F10.5 F10.5 F10.5)') x_vec(i), y_vec(j), z_vec(i,j)
        end do
        write (fu, *) ""
    end do
    close(fu)


    print *, "Ploting started"
    call execute_command_line('gnuplot -p ' // PLT_SAVE)    ! Plot and save
    
    
    print *, "The program is finished"
end program P1
