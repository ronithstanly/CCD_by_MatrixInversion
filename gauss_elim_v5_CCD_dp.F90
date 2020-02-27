MODULE matrices
  INTEGER, PARAMETER :: DP = KIND(1.0d0) 
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a
  REAL(DP), DIMENSION(:), ALLOCATABLE :: b,sol
  REAL(DP), DIMENSION(:), ALLOCATABLE :: f,x,exact_first,exact_second,sol_first,sol_second
  REAL(DP) :: h, error_first, error_second
  INTEGER :: n, n_grid !Size of matrix n x n
    
    CONTAINS

    SUBROUTINE A_matrix
        IMPLICIT NONE
        INTEGER :: i,j
        !ALLOCATE(a(1:n,1:n),b(1:n))
        ALLOCATE(a(1:n,1:n))

        a(1:n,1:n) = 0.0d0  

        a(1,1)=1.0d0; a(1,2)=0.0d0; a(1,3)=2.0d0; a(1,4)=-h/2.0d0;
        a(2,1)=0.0d0; a(2,2)=h;     a(2,3)=-6.0d0;a(2,4)=5.0d0*h; 

        j=1;
        DO i=3,n-2
            IF (MOD(i,2).NE.0) THEN
                a(i,j)=7.0d0; a(i,j+1)=h; a(i,j+2)=16.0d0; a(i,j+3)=0.0d0; a(i,j+4)=7.0d0; a(i,j+5)=-h;
            ELSE IF (MOD(i,2)==0) THEN
                a(i,j)=-9.0d0; a(i,j+1)=-h; a(i,j+2)=0.0d0; a(i,j+3)=8.0d0*h; a(i,j+4)=9.0d0; a(i,j+5)=-h; a(i,j+6)=0.0d0;
                j = j+2
            END IF    
        END DO

        a(n-1,n-3)=2.0d0; a(n-1,n-2)=h/2.0d0; a(n-1,n-1)=1.0d0; a(n-1,n)=0.0d0;
        a(n,n-3)=-6.0d0; a(n,n-2)=-5.0d0*h; a(n,n-1)=0.0d0; a(n,n)=-h;

    OPEN(20,file='test_matrix_a.dat',status='replace')
            DO i=1, n
                WRITE(20,'(50f5.1)') a(i,:) !Printing out array assuming max of 50 elements in each row
            END DO
    CLOSE(20)
    END SUBROUTINE A_Matrix
    
    
    SUBROUTINE B_Matrix
        IMPLICIT NONE
        INTEGER :: i,j

        ALLOCATE(b(1:n))
        ALLOCATE(x(1:n_grid),f(1:n_grid),exact_first(1:n_grid),exact_second(1:n_grid))

        b(1:n) = 0.0d0
    
        x(1:n_grid) = 0.0d0
        f(1:n_grid) = 0.0d0
        exact_first(1:n_grid)   = 0.0d0
        exact_second(1:n_grid)  = 0.0d0
       
        DO i=1, n_grid
            x(i)            = 2.0d0 + DFLOAT(i)*h    ! Grid
            !f(i)            = x(i)**3!DSIN(x(i))/(x(i)**3)  ! Function
            !exact_first(i)  = 3.0d0*x(i)**2!(x(i)*DCOS(x(i)) - 3.0d0*DSIN(x(i)))/(x(i)**4) ! Exact first derivative
            !exact_second(i) = 6.0d0*x(i)!(-6.0d0*x(i)*DCOS(x(i)) - (x(i)**2-12.0d0)*DSIN(x(i)) ) / (x(i)**5) 
            f(i)            = DSIN(x(i))/(x(i)**3)  ! Function
            exact_first(i)  = (x(i)*DCOS(x(i)) - 3.0d0*DSIN(x(i)))/(x(i)**4) ! Exact first derivative
            exact_second(i) = (-6.0d0*x(i)*DCOS(x(i)) - (x(i)**2-12.0d0)*DSIN(x(i)) ) / (x(i)**5) 
        END DO
!    OPEN(25,file='test_matrix_x.dat',status='replace')
!            DO i=1, n
!                IF (MOD(i,2).NE.0) THEN ! Odd
!                    WRITE(25,*) exact_first(i) !Printing out array assuming max of 50 elements in each row
!                ELSE IF (MOD(i,2)==0) THEN ! Even
!                    WRITE(25,*) exact_second(i) !Printing out array assuming max of 50 elements in each row
!                END IF 
!            END DO
!    CLOSE(25)
              

        !B Matrix
        b(1) = 3.0d0*(f(2)-f(1))
        b(2) = 9.0d0*f(1) - 12.0d0*f(2) + 3.0d0*f(3)

        j = 2
        DO i=3, n-2
            IF (MOD(i,2).NE.0) THEN ! Odd
                !b(i) = 15.0d0 * (f(i+1)-f(i-1))
                b(i) = 15.0d0 * (f(j+1)-f(j-1))
            ELSE IF (MOD(i,2)==0) THEN ! Even
                !b(i) = 24.0d0 * (f(i-1)-2.0d0*f(i)+f(i+1))
                b(i) = 24.0d0 * (f(j-1)-2.0d0*f(j)+f(j+1))
                j = j+1
            END IF 
        END DO

        b(n-1)  = 3.0d0* (f(j)-f(j-1))
        b(n)    = -9.0d0*f(j) + 12.0d0*f(j-1) - 3.0d0*f(j-2)   
        
        b(:)    = b(:)/h
    OPEN(24,file='test_matrix_b.dat',status='replace')
            DO i=1, n
                WRITE(24,*) b(i) !Printing out array assuming max of 50 elements in each row
            END DO
    CLOSE(24)
              
    OPEN(21,file='test_out.dat',status='replace')
            WRITE(21,*) 'x, f, exact_first, exact_second'
            DO i=1, n_grid
                WRITE(21,'(100g15.5)') x(i),f(i),exact_first(i),exact_second(i)
            END DO
    CLOSE(21)

    END SUBROUTINE B_Matrix


    SUBROUTINE print_sol
        IMPLICIT NONE
        INTEGER  :: i,j
        REAL(DP) :: sum_first,sum_second,diff_first,diff_second 

        ALLOCATE(sol_first(1:n_grid),sol_second(1:n_grid))
        sol_first(:)    = 0.0d0
        sol_second(:)   = 0.0d0 
    
    OPEN(22,file='test_sol_first_der.dat',status='replace')
            j = 1
            DO i=1, n
                IF (MOD(i,2).NE.0) THEN ! Odd
                    WRITE(22,*) sol(i)
                    sol_first(j) = sol(i)
                    j = j+1
                END IF 
            END DO
    CLOSE(22)

    OPEN(23,file='test_sol_sec_der.dat',status='replace')
            j = 1
            DO i=1, n
                IF (MOD(i,2)==0) THEN ! Even
                    WRITE(23,*) sol(i) ! change this
                    sol_second(j) = sol(i)
                    j = j+1
                END IF 
            END DO
    CLOSE(23)

    error_first = 0.0d0
    error_second = 0.0d0

    sum_first = 0.0d0
    diff_first = 0.0d0
    
    DO i=1, n_grid
    !DO i=3, n_grid-3 ! To check accuracy in the interior points
        diff_first = DABS(sol_first(i) - exact_first(i))
        sum_first = sum_first + diff_first
    END DO
    error_first = sum_first/DFLOAT(n_grid)
    !error_first = sum_first/DFLOAT(n_grid-6) ! To check accuracy in the interior points
    
    sum_second = 0.0d0
    diff_second = 0.0d0

    DO i=1, n_grid
    !DO i=3, n_grid-3 ! To check accuracy in the interior points
        diff_second = DABS(sol_second(i) - exact_second(i))
        sum_second = sum_first + diff_second
    END DO
    error_second = sum_second/DFLOAT(n_grid)
    !error_second = sum_second/DFLOAT(n_grid-6) ! To check accuracy in the interior points

    END SUBROUTINE print_sol

END MODULE matrices

!###########################
PROGRAM gauss_elim
  USE matrices

    n = 200 ! Double of grid points (Matrices are n/2 X n/2)
    n_grid = n/2 ! Grid points
    h = 1/(n/2.0d0)!4.0d0/FLOAT(n_grid)!0.1d0 ! Grid spacing

    CALL A_Matrix
    CALL B_Matrix

        ALLOCATE(sol(1:n))
        sol(1:n) = 0.0d0

        sol = solve_wbs(ge_wpp(a,b))
    CALL print_sol
   
        !WRITE(*,*) '*** Solution Vector ***' 
        !PRINT'(f15.7)',sol
        WRITE(*,*) '*** First Derivatives ***'
        PRINT'(f15.7)',sol_first
        WRITE(*,*) '*** Second Derivatives ***'
        PRINT'(f15.7)',sol_second
        WRITE(*,*) '*** Error-First Derivative ***' 
        WRITE(*,*) error_first!ABS(SUM(sol_first)-SUM(exact_first))/FLOAT(n_grid)
        WRITE(*,*) '*** Error-Second Derivative ***' 
        WRITE(*,*) error_second!ABS(SUM(sol_second)-SUM(exact_second))/FLOAT(n_grid)
        WRITE(*,*) '*** Grid Spacing ***'
        WRITE(*,*) 'h=',h           


    !PRINT'(f15.7)',solve_wbs(ge_wpp(a,b))
 
	CONTAINS
 
	FUNCTION solve_wbs(u) result(x) ! solve with backward substitution
    IMPLICIT NONE
      REAL(DP)                 :: u(:,:)
      INTEGER              :: i,n
      REAL(DP)   , ALLOCATABLE :: x(:)
      n = SIZE(u,1)
      ALLOCATE(x(n))
        FORALL (i=n:1:-1) x(i) = ( u(i,n+1) - SUM(u(i,i+1:n)*x(i+1:n)) ) / u(i,i)
    END FUNCTION
 
    FUNCTION  ge_wpp(a,b) result(u) ! gaussian eliminate with partial pivoting
     IMPLICIT NONE
      real(DP)                 :: a(:,:),b(:),upi
      integer              :: i,j,n,p
      real(DP)   , allocatable :: u(:,:)
      n = SIZE(a,1)
      u = reshape( [a,b], [n,n+1] )
      DO j=1,n
        p = MAXLOC(DABS(u(j:n,j)),1) + j-1 ! maxloc returns indices between (1,n-j+1)
        IF (p /= j) u([p,j],j) = u([j,p],j)
        u(j+1:,j) = u(j+1:,j)/u(j,j)
        DO i=j+1,n+1
          upi = u(p,i)
          IF (p /= j) u([p,j],i) = u([j,p],i)
          u(j+1:n,i) = u(j+1:n,i) - upi*u(j+1:n,j)
        END DO
      END DO
    END FUNCTION
 
END PROGRAM gauss_elim
 
