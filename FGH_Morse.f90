! Aditya Barman, Project Associate, MNIT Jaipur, India; April 10, 2025
! Email: atomicadi2023@gmail.com

program Fourier_Grid_Hamiltonia_Morse_pot
       implicit none
       integer :: i, j, p, l, N,  Xmax, Xmin, lda, lwork, info
          real(kind=8) :: val_0, val_1, val_2, val_3, val_4, val_5, val_6, pow, pi, pi_val, delx, delk, beta, D
          real(kind=8), allocatable :: H(:,:), T(:,:), V(:,:), x(:), v_pot(:), eigval(:), work(:)

        pi = 3.141592653589793 
        pi_val = 2.0d0 * pi
        Xmin = 0
        Xmax = 2
           N = 10 
        delx = real(Xmax - Xmin)/real(N-1)
        delk = pi_val/(N * delx)
         lda = N
       lwork = (3 * N) - 1
           D = 0.1744d0   ! Dissociation energy for H2 in a. u. unit
         beta = 1.02764d0 ! beta for H2 in a. u. unit

        call find_odd(N,l)

        allocate(H(N,N))
        allocate(T(N,N))
        allocate(V(N,N))
        allocate(work(lwork))
        allocate(eigval(N))
        allocate(x(N), v_pot(N))
        

!=========================================================
!              Kinetic energy matrix                    ||
!=========================================================
              val_0 = 0.0d0                
              val_5 = 0.0d0

              do i = 1, N
               do j = 1, N
                if (i == j) then
                 do p = 1, l
                    val_1 = (0.5d0/918.0764d0) * (p*delk)**2
                    val_0 = val_1 + val_0
                  end do
                    val_2 = val_0 * (2.0/real(N))
                    T(i,j) = val_2
                    val_0 = 0.0d0
                 else 
                    do p = 1, l
                        val_1 = (0.5d0/918.0764d0) * (p * delk) ** 2 !here 918.0764 is reduced mass of H2 in a. u. unit
                        pow = (p * pi_val * (i-j))/N
                        val_3 = cos(pow) 
                        val_4 = val_3 * val_1
                        val_5 = val_4 + val_5
                     end do
                      val_6 = (val_5 * (2.0d0/(real(N))))
                      T(i,j) = val_6
                      val_5 = 0.0d0
                end if
               end do
             end do
         print*, "Kinetic energy matrix (first 10 rows and 3 column):"
          do i= 1, 10
           print*, (T(i,j), j= 1, 3)
         end do

!=========================================================
!              potential energy matrix                  ||
!=========================================================
        do i = 1, N
           x(i) = Xmin + (i-1) * delx
           v_pot(i) = (D * (1.0d0 - exp(- beta * (x(i) - 1.40201d0))) ** 2)
         end do

         do i = 1, N
            do j = 1,N
              if (i == j) then
                 V(i,j) = v_pot(i)
              else
                  V(i,j) = 0.0
              end if
            end do
         end do

          print*, "Potential energy matrix (first 10 rows and 3 column):"
         do i= 1, 10
           print*, (V(i,j), j= 1, 3)
         end do 

!=========================================================
!             Hamiltonian matrix                        ||
!========================================================= 
         H = T + V 
         print*, "Hamiltonian matrix (first 10 rows and 3 column):"
         do i= 1, 10
           print*, (H(i,j), j= 1, 3)
         end do 

!=========================================================
!             Diagonalize Hamiltonian matrix            ||
!========================================================= 
        call dsyev('V', 'U', N, H, lda, eigval, work, lwork, info)

        print*, "Eigenvalues (first 10):"
        do i = 1, 10
          write(*,*) eigval(i)
        end do
          
         deallocate(H)
         deallocate(T)
         deallocate(x,v_pot)
         deallocate(V)
         deallocate(work)
         deallocate(eigval)

  contains

     subroutine find_odd(input_val,output_val)
          implicit none
          integer :: input_val, output_val

          if (mod(input_val,2) == 1) then
               output_val = (input_val - 1)/2
          else 
              output_val = input_val/2
          end if
     end subroutine 

end program Fourier_Grid_Hamiltonia_Morse_pot

!command line [gfortran -o file file.f90 -llapack -lblas] file = whatever_you_give_the_name_of_your_file 
              
 