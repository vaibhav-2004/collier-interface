program loopdrv
  use collier
  implicit none
  integer, parameter :: dp = kind(1d0)
  integer :: n, i, nexp
  real(dp) :: a(16), mu2
  complex(dp) :: v0, v1, v2, c1, c2
  character(len=8)  :: what
  character(len=64) :: buf

  n = command_argument_count()
  if (n < 3) then
     write(0,*) 'usage: loopdrv {A0|B0|C0|D0} <args...> mu2'
     stop 1
  end if
  call get_command_argument(1, buf); what = trim(buf)
  do i = 2, n
     call get_command_argument(i, buf); read(buf,*) a(i-1)
  end do

  select case (trim(what))
  case ('A0'); nexp = 3
  case ('B0'); nexp = 5
  case ('C0'); nexp = 8
  case ('D0'); nexp = 12
  case default
     write(0,*) 'unknown integral: ', trim(what); stop 1
  end select
  if (n /= nexp) then
     write(0,*) trim(what), ' expects ', nexp-1, ' numbers, got ', n-1
     stop 1
  end if
  mu2 = a(n-1)

  call Init_cll(4, 0, 'collier_output', .true.)
  call SetMode_cll(1)
  call SetMuUV2_cll(mu2)
  call SetMuIR2_cll(mu2)

  call SetDeltaUV_cll(0d0)
  call SetDeltaIR_cll(0d0, 0d0)
  v0 = evaluate()

  call SetDeltaUV_cll(1d0)
  call SetDeltaIR_cll(1d0, 0d0)
  v1 = evaluate()

  call SetDeltaUV_cll(0d0)
  call SetDeltaIR_cll(0d0, 1d0)
  v2 = evaluate()

  c1 = v1 - v0
  c2 = v2 - v0

  write(*,'(6ES26.17)') real(v0), aimag(v0), &
                        real(c1), aimag(c1), &
                        real(c2), aimag(c2)

contains

  complex(dp) function evaluate()
    complex(dp) :: r
    select case (trim(what))
    case ('A0')
       call A0_cll(r, cx(a(1)))
    case ('B0')
       call B0_cll(r, cx(a(1)), cx(a(2)), cx(a(3)))
    case ('C0')
       call C0_cll(r, cx(a(1)), cx(a(2)), cx(a(3)), &
                      cx(a(4)), cx(a(5)), cx(a(6)))
    case ('D0')
       call D0_cll(r, cx(a(1)), cx(a(2)), cx(a(3)), cx(a(4)), &
                      cx(a(5)), cx(a(6)), &
                      cx(a(7)), cx(a(8)), cx(a(9)), cx(a(10)))
    end select
    evaluate = r
  end function evaluate

  complex(dp) function cx(x)
    real(dp), intent(in) :: x
    cx = cmplx(x, 0d0, dp)
  end function cx

end program loopdrv
