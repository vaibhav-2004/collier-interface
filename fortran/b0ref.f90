program b0ref
  implicit none
  integer, parameter :: dp = kind(1d0)
  real(dp) :: p2, m02, m12, mu2
  complex(dp) :: val
  character(len=64) :: a1, a2, a3, a4

  if (command_argument_count() /= 4) then
     write(0,*) 'usage: b0ref p2 m02 m12 mu2'
     stop 1
  end if
  call get_command_argument(1,a1); read(a1,*) p2
  call get_command_argument(2,a2); read(a2,*) m02
  call get_command_argument(3,a3); read(a3,*) m12
  call get_command_argument(4,a4); read(a4,*) mu2

  val = b0fin(p2, m02, m12, mu2)
  write(*,'(2ES26.17)') real(val), aimag(val)

contains

  complex(dp) function b0fin(p2, m02, m12, mu2)
    real(dp), intent(in) :: p2, m02, m12, mu2
    real(dp), parameter  :: eps = 1d-20
    complex(dp) :: M0, M1, b, sq, xp, xm, r

    M0 = cmplx(m02, -eps, dp)
    M1 = cmplx(m12, -eps, dp)

    if (abs(p2) < 1d-10) then
       if (abs(m02 - m12) < 1d-10) then
          b0fin = -log(M0/mu2)
       else
          b0fin = 1d0 - (M0*log(M0/mu2) - M1*log(M1/mu2))/(M0 - M1)
       end if
       return
    end if

    b  = p2 + M0 - M1
    sq = sqrt(b*b - 4d0*p2*M0)
    xp = (b + sq)/(2d0*p2)
    xm = (b - sq)/(2d0*p2)

    r = 2d0 - log(cmplx(p2, -eps, dp)/mu2)
    r = r - ((1d0 - xp)*log(1d0 - xp) + xp*log(-xp))
    r = r - ((1d0 - xm)*log(1d0 - xm) + xm*log(-xm))
    b0fin = r
  end function b0fin

end program b0ref
