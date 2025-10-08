! Aditya Barman, Project Associate, MNIT Jaipur India; April 12, 2025
! Email: atomicadi2023@gmail.com


program Fourier_Grid_Hamiltonia_SHO
       implicit none
       integer :: i, j, p, l, N, Xmax, Xmin, lda, lwork, info
          real(kind=8) :: val_0, val_1, val_2, val_3, val_4, val_5, val_6, pow, pi, pi_val, delx, delk, b, Vmax, mu, re 
          real(kind=8), allocatable :: H(:,:), T(:,:), V(:,:), x(:), v_pot(:), eigval(:), work(:)

        pi = 3.141592653589793 
        pi_val = 2.0d0 * pi
        Xmin = 50  
        Xmax = 130 
           N = 7 ! No. of grid
        delx = real(Xmax - Xmin)/real(N-1)
        delk = pi_val/(N * delx)
         lda = N
       lwork = 3 * N -1
        Vmax = 0.0080929597d0 !kcal/mol to a.u. 
           b = 23.2749d0
          re = 90.0002d0    
          mu = 155.5d0  !Reduced mass for NH3 (in atomic units) at the umbrella inversion mode (2506 cm-1 (in B3LYP/pVDZ level of theory) and k(force constant) = 0.7589 a.u.
            

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
              val_0 = 0.0                
              val_5 = 0.0

              do i = 1, N
               do j = 1, N
                 do p = 1, l
                        val_1 =  (0.5d0/mu) * (p * delk) ** 2
                        pow = (p * pi_val * (i-j))/N
                        val_3 = cos(pow) 
                        val_4 = val_3 * val_1
                        val_5 = val_4 + val_5
                     end do
                      val_6 = val_5 * (2.0d0/(real(N)))
                      T(i,j) = val_6
                      val_5 = 0.0d0
                end do
             end do
         print*, "Kinetic energy matrix (first 10 rows):"
          do i= 1, 10
           print*, (T(i,j), j= 1, 3)
         end do

!=========================================================
!              potential energy matrix                  ||
!=========================================================
        do i = 1, N
           x(i) = Xmin + (i-1) * delx
           v_pot(i) = (Vmax/b**4) * (((x(i) - re)**2) - b**2)**2 - Vmax 
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

          print*, "Potential energy matrix (first 10 rows):"
         do i= 1, 10
           print*, (V(i,j), j= 1, 3)
         end do 

!=========================================================
!             Hamiltonian matrix                        ||
!========================================================= 
         H = T + V
         print*, "Hamiltonian matrix (first 10 rows):"
         do i= 1, 7
           print*, (H(i,j), j= 1, 3)
         end do 

!=========================================================
!             Diagonalize Hamiltonian matrix            ||
!========================================================= 
        call dsyev('V', 'U', N, H, lda, eigval, work, lwork, info)

        print*, "Eigenvalues (first 10):"
        do i = 1, 11
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

end program Fourier_Grid_Hamiltonia_SHO

!command line [gfortran -o file file.f90 -llapack -lblas] file = whatever_you_give_the_name_of_your_file    
 