PROGRAM fd1d_heat_explicit_prb
  USE :: types_mod, ONLY: dp

  IMPLICIT NONE

  INTEGER :: t_num
  PARAMETER (t_num=201)
  INTEGER :: x_num
  PARAMETER (x_num=21)

  REAL (KIND=dp) :: cfl
  REAL (KIND=dp) :: dt
  REAL (KIND=dp) :: h(x_num)
  REAL (KIND=dp) :: h_new(x_num)
! the "matrix" stores all x-values for all t-values
! remember Fortran is column major, meaning that rows are contiguous
  REAL (KIND=dp) :: hmat(x_num, t_num)
  INTEGER :: i
  INTEGER :: j
  REAL (KIND=dp) :: k

  REAL (KIND=dp) :: t(t_num)
  REAL (KIND=dp) :: t_max
  REAL (KIND=dp) :: t_min
  REAL (KIND=dp) :: x(x_num)
  REAL (KIND=dp) :: x_max
  REAL (KIND=dp) :: x_min

  WRITE (*, '(a)') ' '
  WRITE (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
  WRITE (*, '(a)') '  FORTRAN77 version.'
  WRITE (*, '(a)') '  Test the FD1D_HEAT_EXPLICIT library.'

  WRITE (*, '(a)') ' '
  WRITE (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
  WRITE (*, '(a)') '  Normal end of execution.'
  WRITE (*, '(a)') ' '

  WRITE (*, '(a)') ' '
  WRITE (*, '(a)') 'FD1D_HEAT_EXPLICIT_TEST01:'
  WRITE (*, '(a)') '  Compute an approximate solution to the time-dependent'
  WRITE (*, '(a)') '  one dimensional heat equation:'
  WRITE (*, '(a)') ' '
  WRITE (*, '(a)') '    dH/dt - K * d2H/dx2 = f(x,t)'
  WRITE (*, '(a)') ' '
  WRITE (*, '(a)') '  Run a simple test case.'

! heat coefficient
  k = 0.002E+00_dp

! the x-range values
  x_min = 0.0E+00_dp
  x_max = 1.0E+00_dp
! x_num is the number of intervals in the x-direction
  CALL r8vec_linspace(x_num, x_min, x_max, x)

! the t-range values. integrate from t_min to t_max
  t_min = 0.0E+00_dp
  t_max = 80.0E+00_dp

! t_num is the number of intervals in the t-direction
  dt = (t_max-t_min)/real(t_num-1, kind=dp)
  CALL r8vec_linspace(t_num, t_min, t_max, t)

! get the CFL coefficient
  CALL fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, x_max, &
    cfl)

  IF (0.5E+00_dp<=cfl) THEN
    WRITE (*, '(a)') ' '
    WRITE (*, '(a)') 'FD1D_HEAT_EXPLICIT_CFL - Fatal error!'
    WRITE (*, '(a)') '  CFL condition failed.'
    WRITE (*, '(a)') '  0.5 <= K * dT / dX / dX = CFL.'
    STOP
  END IF

! set the initial condition
  DO j = 1, x_num
    h(j) = 50.0E+00_dp
  END DO

! set the bounday condition
  h(1) = 90.0E+00_dp
  h(x_num) = 70.0E+00_dp

! initialise the matrix to the initial condition
  DO i = 1, x_num
    hmat(i, 1) = h(i)
  END DO

! the main time integration loop 
  DO j = 2, t_num
    CALL fd1d_heat_explicit(x_num, x, t(j-1), dt, cfl, h, h_new)

    DO i = 1, x_num
      hmat(i, j) = h_new(i)
      h(i) = h_new(i)
    END DO
  END DO

! write data to files
  CALL r8mat_write('h_test01.txt', x_num, t_num, hmat)
  CALL r8vec_write('t_test01.txt', t_num, t)
  CALL r8vec_write('x_test01.txt', x_num, x)

CONTAINS

  FUNCTION func(j, x_num, x) RESULT (d)
    IMPLICIT NONE

    INTEGER :: j, x_num
    REAL (KIND=dp) :: d
    REAL (KIND=dp) :: x(x_num)

    d = 0.0E+00_dp
  END FUNCTION

  SUBROUTINE fd1d_heat_explicit(x_num, x, t, dt, cfl, h, h_new)
    IMPLICIT NONE

    INTEGER :: x_num

    REAL (KIND=dp) :: cfl
    REAL (KIND=dp) :: dt
    REAL (KIND=dp) :: h(x_num)
    REAL (KIND=dp) :: h_new(x_num)
    INTEGER :: j
    REAL (KIND=dp) :: t
    REAL (KIND=dp) :: x(x_num)
    REAL (KIND=dp) :: f(x_num)

    DO j = 1, x_num
      f(j) = func(j, x_num, x)
    END DO

    h_new(1) = 0.0E+00_dp

    DO j = 2, x_num - 1
      h_new(j) = h(j) + dt*f(j) + cfl*(h(j-1)-2.0E+00_dp*h(j)+h(j+1))
    END DO

! set the boundary conditions again
    h_new(1) = 90.0E+00_dp
    h_new(x_num) = 70.0E+00_dp
  END SUBROUTINE

  SUBROUTINE fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, &
    x_max, cfl)

    IMPLICIT NONE

    REAL (KIND=dp) :: cfl
    REAL (KIND=dp) :: dx
    REAL (KIND=dp) :: dt
    REAL (KIND=dp) :: k
    REAL (KIND=dp) :: t_max
    REAL (KIND=dp) :: t_min
    INTEGER :: t_num
    REAL (KIND=dp) :: x_max
    REAL (KIND=dp) :: x_min
    INTEGER :: x_num

    dx = (x_max-x_min)/real(x_num-1, kind=dp)
    dt = (t_max-t_min)/real(t_num-1, kind=dp)

    cfl = k*dt/dx/dx

    WRITE (*, '(a)') ' '
    WRITE (*, '(a,g14.6)') '  CFL stability criterion value = ', cfl

  END SUBROUTINE

  SUBROUTINE r8mat_write(output_filename, m, n, table)
    IMPLICIT NONE

    INTEGER :: m
    INTEGER :: n

    INTEGER :: j
    CHARACTER (LEN=*) :: output_filename
    INTEGER :: output_unit_id
    CHARACTER (LEN=30) :: string
    REAL (KIND=dp) :: table(m, n)

    output_unit_id = 10
    OPEN (UNIT=output_unit_id, FILE=output_filename, STATUS='replace')

    WRITE (string, '(a1,i8,a1,i8,a1,i8,a1)') '(', m, 'g', 24, '.', 16, ')'

    DO j = 1, n
      WRITE (output_unit_id, string) table(1:m, j)
    END DO

    CLOSE (UNIT=output_unit_id)
  END SUBROUTINE

  SUBROUTINE r8vec_linspace(n, a_first, a_last, a)

    IMPLICIT NONE

    INTEGER :: n
    REAL (KIND=dp) :: a(n)
    REAL (KIND=dp) :: a_first
    REAL (KIND=dp) :: a_last
    INTEGER :: i

    DO i = 1, n
      a(i) = (real(n-i,kind=dp)*a_first+real(i-1,kind=dp)*a_last)/ &
        real(n-1, kind=dp)
    END DO

  END SUBROUTINE

  SUBROUTINE r8vec_write(output_filename, n, x)

    IMPLICIT NONE

    INTEGER :: m
    INTEGER :: n

    INTEGER :: j
    CHARACTER (LEN=*) :: output_filename
    INTEGER :: output_unit_id
    REAL (KIND=dp) :: x(n)

    output_unit_id = 11
    OPEN (UNIT=output_unit_id, FILE=output_filename, STATUS='replace')

    DO j = 1, n
      WRITE (output_unit_id, '(2x,g24.16)') x(j)
    END DO

    CLOSE (UNIT=output_unit_id)
  END SUBROUTINE

END PROGRAM
