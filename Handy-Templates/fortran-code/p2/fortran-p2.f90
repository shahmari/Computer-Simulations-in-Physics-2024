program P2

    implicit none
   
    integer :: i, fu
    character(len=*), parameter :: PLT_SAVE = 'plot-p2.plt'     ! Save
    integer, parameter :: dp = selected_real_kind(15, 307)      ! Double-precision type
    
    integer, parameter :: N = 1000                              ! Dimension of the matrix
    real(dp) :: A_real(N,N), A_img(N,N)
    complex(dp) :: A(N,N), A_t(N,N)
    complex(dp) :: H_sum(N,N), H_mul(N,N)
    real(dp) :: H_sum_eigen(N), H_mul_eigen(N)
    


    call random_seed()                    ! Setting the random seed
    call random_number(A_real)            ! Random real part
    call random_number(A_img)             ! Random imagenary part

    ! Creating A and adjoint of it
    A = A_real + A_img*complex(0.0, 1.0)
    A_t = conjg(transpose(A))

    ! Making Hermitian matrices 
    H_sum = A + A_t                         
    H_mul = matmul(A, A_t)
    
    ! Computing the eigenvalues 
    H_sum_eigen = calc_eigen(H_sum, N)
    H_mul_eigen = calc_eigen(H_mul, N)
    
    ! Writing results in ".txt" files
    open (action='write', file="data-p2-sum.txt", newunit=fu, status='replace')
    do i = 1, N-1
        write (fu, '(F10.5)') H_sum_eigen(i)
    end do
    close(fu)

    open (action='write', file="data-p2-mul.txt", newunit=fu, status='replace')
    do i = 1, N-1
        write (fu, '(F10.5)') H_mul_eigen(i)
    end do
    close(fu)



    print *, "Ploting started"
    call execute_command_line('gnuplot -p ' // PLT_SAVE)    ! Plot and save
    
    
    print *, "The program is finished"


    contains

    function  calc_eigen(mat, N_mat)     ! Calculates eigenvalues of hermitian matrices

        ! This functions uses the ZHEEV subroutine from LAPACK lib.

        implicit none

        integer :: N_mat
        integer :: LWORK
        complex(dp) :: mat(N_mat, N_mat)
        complex(dp) :: mat_copy(N_mat, N_mat)
        real(dp) :: eigen_out(N_mat), calc_eigen(N_mat)
        complex(dp) :: WORK(6*N_mat)
        real(dp) :: RWORK(3*N_mat - 2)
        integer :: info

        LWORK = 6*N_mat

        mat_copy = mat
        call ZHEEV('N', 'U', N_mat, mat_copy, N_mat, eigen_out, WORK, LWORK, RWORK, info)

        print *, "ZHEEV info:", info, "\n"
        calc_eigen = eigen_out

    end function calc_eigen


end program P2
