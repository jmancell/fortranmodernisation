    Program fd1d_heat_explicit_prb
      Use types_mod, Only: dp

      Implicit None

      Integer :: t_num
      Parameter (t_num=201)
      Integer :: x_num
      Parameter (x_num=21)

      Real (Kind=dp) :: cfl
      Real (Kind=dp) :: dt
      Real (Kind=dp) :: h(x_num)
      Real (Kind=dp) :: h_new(x_num)
! the "matrix" stores all x-values for all t-values
! remember Fortran is column major, meaning that rows are contiguous
      Real (Kind=dp) :: hmat(x_num, t_num)
      Integer :: i
      Integer :: j
      Real (Kind=dp) :: k

      Real (Kind=dp) :: t(t_num)
      Real (Kind=dp) :: t_max
      Real (Kind=dp) :: t_min
      Real (Kind=dp) :: x(x_num)
      Real (Kind=dp) :: x_max
      Real (Kind=dp) :: x_min

      Write (*, '(a)') ' '
      Write (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
      Write (*, '(a)') '  FORTRAN77 version.'
      Write (*, '(a)') '  Test the FD1D_HEAT_EXPLICIT library.'

      Write (*, '(a)') ' '
      Write (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
      Write (*, '(a)') '  Normal end of execution.'
      Write (*, '(a)') ' '

      Write (*, '(a)') ' '
      Write (*, '(a)') 'FD1D_HEAT_EXPLICIT_TEST01:'
      Write (*, '(a)') &
        '  Compute an approximate solution to the time-dependent'
      Write (*, '(a)') '  one dimensional heat equation:'
      Write (*, '(a)') ' '
      Write (*, '(a)') '    dH/dt - K * d2H/dx2 = f(x,t)'
      Write (*, '(a)') ' '
      Write (*, '(a)') '  Run a simple test case.'

! heat coefficient
      k = 0.002E+00_dp

! the x-range values
      x_min = 0.0E+00_dp
      x_max = 1.0E+00_dp
! x_num is the number of intervals in the x-direction
      Call r8vec_linspace(x_num, x_min, x_max, x)

! the t-range values. integrate from t_min to t_max
      t_min = 0.0E+00_dp
      t_max = 80.0E+00_dp

! t_num is the number of intervals in the t-direction
      dt = (t_max-t_min)/real(t_num-1, kind=dp)
      Call r8vec_linspace(t_num, t_min, t_max, t)

! get the CFL coefficient
      Call fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, x_max, &
        cfl)

      If (0.5E+00_dp<=cfl) Then
        Write (*, '(a)') ' '
        Write (*, '(a)') 'FD1D_HEAT_EXPLICIT_CFL - Fatal error!'
        Write (*, '(a)') '  CFL condition failed.'
        Write (*, '(a)') '  0.5 <= K * dT / dX / dX = CFL.'
        Stop
      End If

! set the initial condition
      Do j = 1, x_num
        h(j) = 50.0E+00_dp
      End Do

! set the bounday condition
      h(1) = 90.0E+00_dp
      h(x_num) = 70.0E+00_dp

! initialise the matrix to the initial condition
      Do i = 1, x_num
        hmat(i, 1) = h(i)
      End Do

! the main time integration loop 
      Do j = 2, t_num
        Call fd1d_heat_explicit(x_num, x, t(j-1), dt, cfl, h, h_new)

        Do i = 1, x_num
          hmat(i, j) = h_new(i)
          h(i) = h_new(i)
        End Do
      End Do

! write data to files
      Call r8mat_write('h_test01.txt', x_num, t_num, hmat)
      Call r8vec_write('t_test01.txt', t_num, t)
      Call r8vec_write('x_test01.txt', x_num, x)

    Contains

      Function func(j, x_num, x) Result (d)
        Implicit None

        Integer :: j, x_num
        Real (Kind=dp) :: d
        Real (Kind=dp) :: x(x_num)

        d = 0.0E+00_dp
      End Function

      Subroutine fd1d_heat_explicit(x_num, x, t, dt, cfl, h, h_new)
        Implicit None

        Integer :: x_num

        Real (Kind=dp) :: cfl
        Real (Kind=dp) :: dt
        Real (Kind=dp) :: h(x_num)
        Real (Kind=dp) :: h_new(x_num)
        Integer :: j
        Real (Kind=dp) :: t
        Real (Kind=dp) :: x(x_num)
        Real (Kind=dp) :: f(x_num)

        Do j = 1, x_num
          f(j) = func(j, x_num, x)
        End Do

        h_new(1) = 0.0E+00_dp

        Do j = 2, x_num - 1
          h_new(j) = h(j) + dt*f(j) + cfl*(h(j-1)-2.0E+00_dp*h(j)+h(j+1))
        End Do

! set the boundary conditions again
        h_new(1) = 90.0E+00_dp
        h_new(x_num) = 70.0E+00_dp
      End Subroutine

      Subroutine fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, &
        x_max, cfl)

        Implicit None

        Real (Kind=dp) :: cfl
        Real (Kind=dp) :: dx
        Real (Kind=dp) :: dt
        Real (Kind=dp) :: k
        Real (Kind=dp) :: t_max
        Real (Kind=dp) :: t_min
        Integer :: t_num
        Real (Kind=dp) :: x_max
        Real (Kind=dp) :: x_min
        Integer :: x_num

        dx = (x_max-x_min)/real(x_num-1, kind=dp)
        dt = (t_max-t_min)/real(t_num-1, kind=dp)

        cfl = k*dt/dx/dx

        Write (*, '(a)') ' '
        Write (*, '(a,g14.6)') '  CFL stability criterion value = ', cfl

      End Subroutine

      Subroutine r8mat_write(output_filename, m, n, table)
        Implicit None

        Integer :: m
        Integer :: n

        Integer :: j
        Character (*) :: output_filename
        Integer :: output_unit_id
        Character (30) :: string
        Real (Kind=dp) :: table(m, n)

        output_unit_id = 10
        Open (Unit=output_unit_id, File=output_filename, Status='replace')

        Write (string, '(a1,i8,a1,i8,a1,i8,a1)') '(', m, 'g', 24, '.', 16, ')'

        Do j = 1, n
          Write (output_unit_id, string) table(1:m, j)
        End Do

        Close (Unit=output_unit_id)
      End Subroutine

      Subroutine r8vec_linspace(n, a_first, a_last, a)

        Implicit None

        Integer :: n
        Real (Kind=dp) :: a(n)
        Real (Kind=dp) :: a_first
        Real (Kind=dp) :: a_last
        Integer :: i

        Do i = 1, n
          a(i) = (real(n-i,kind=dp)*a_first+real(i-1,kind=dp)*a_last)/ &
            real(n-1, kind=dp)
        End Do

      End Subroutine

      Subroutine r8vec_write(output_filename, n, x)

        Implicit None

        Integer :: m
        Integer :: n

        Integer :: j
        Character (*) :: output_filename
        Integer :: output_unit_id
        Real (Kind=dp) :: x(n)

        output_unit_id = 11
        Open (Unit=output_unit_id, File=output_filename, Status='replace')

        Do j = 1, n
          Write (output_unit_id, '(2x,g24.16)') x(j)
        End Do

        Close (Unit=output_unit_id)
      End Subroutine

    End Program
