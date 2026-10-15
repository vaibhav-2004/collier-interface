! loopdrv  --  A0 / B0 / C0 / D0 through Collier
!   ./loopdrv A0 m02 mu2
!   ./loopdrv B0 p10 m02 m12 mu2
!   ./loopdrv C0 p10 p21 p20 m02 m12 m22 mu2
!   ./loopdrv D0 p10 p21 p32 p30 p20 p31 m02 m12 m22 m32 mu2
! pij means (p_i - p_j)^2 ; for D0 the last two momenta are the diagonals.
program loopdrv
  use collier
  implicit none
  integer, parameter :: dp = kind(1d0)
  integer :: n, i
  real(dp) :: a(16), mu2
  complex(dp) :: val
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
  ! explicit argument counts, including mu2 as the final number
  select case (trim(what))
  case ('A0'); if (n /= 3)  then; write(0,*) 'A0 needs: m02 mu2';                stop 1; end if
  case ('B0'); if (n /= 5)  then; write(0,*) 'B0 needs: p10 m02 m12 mu2';        stop 1; end if
  case ('C0'); if (n /= 8)  then; write(0,*) 'C0 needs: p10 p21 p20 m02 m12 m22 mu2'; stop 1; end if
  case ('D0'); if (n /= 12) then; write(0,*) 'D0 needs: p10 p21 p32 p30 p20 p31 m02 m12 m22 m32 mu2'; stop 1; end if
  end select
  mu2 = a(n-1)

  call Init_cll(4, 0, 'collier_output', .true.)
  call SetMode_cll(1)
  call SetMuUV2_cll(mu2)
  call SetMuIR2_cll(mu2)
  call SetDeltaUV_cll(0d0)
  call SetDeltaIR_cll(0d0, 0d0)

  select case (what)
  case ('A0')
     call A0_cll(val, cx(a(1)))
  case ('B0')
     call B0_cll(val, cx(a(1)), cx(a(2)), cx(a(3)))
  case ('C0')
     call C0_cll(val, cx(a(1)), cx(a(2)), cx(a(3)), &
                      cx(a(4)), cx(a(5)), cx(a(6)))
  case ('D0')
     call D0_cll(val, cx(a(1)), cx(a(2)), cx(a(3)), cx(a(4)), &
                      cx(a(5)), cx(a(6)), &
                      cx(a(7)), cx(a(8)), cx(a(9)), cx(a(10)))
  case default
     write(0,*) 'unknown: ', what; stop 1
  end select

  write(*,'(2ES26.17)') real(val), aimag(val)

contains
  complex(dp) function cx(x)
    real(dp), intent(in) :: x
    cx = cmplx(x, 0d0, dp)
  end function cx
end program loopdrv
