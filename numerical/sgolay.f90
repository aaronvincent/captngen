!! Savitzky-Golay smoothing routine
subroutine sgolay(xin,xlen,n,order,xout)
  implicit none
  !xin: data in
  !n: left and right length: your filter width is 2n+1
  !order: 0 = smooth the function, 1 = first derivative
  !xout: output
  integer, intent(in) :: n, order,xlen
  double precision :: xin(xlen), xout(xlen), xbig(xlen+2*n)
  double precision :: coeffs(-n:n), norm
  integer ::  ii, jj

  call get_sgolay_coeffs(n,order,coeffs,norm)

  ! print*,n, order, coeffs, norm

  xbig(1:n) = xin(1)
  xbig(n+xlen+1:2*n+xlen) = xin(xlen)
  xbig(n+1:n+xlen)= xin(:)

  ! print*,xbig

  xout = 0.
  do ii = 1,xlen
    do jj = -n,n
      xout(ii) = xout(ii)+coeffs(jj)*xbig(ii+n+jj+1)
    end do
  end do
  xout = xout/norm

return
end subroutine sgolay


subroutine get_sgolay_coeffs(n,order,coeffs,norm)
  implicit none
  integer, intent(in) :: n, order
  double precision :: coeffs(2*n+1), norm

  select case (order)
  case (0)
    select case (n)
      case (2)
        coeffs = (/ -3., 12., 17., 12., -3./)
        norm = 35.
      case (3)
        coeffs = (/ -2., 3., 6., 7., 6., 3., -2./)
        norm = 21.
      case (4)
        coeffs = (/-21., 14., 39., 54., 59., 54., 39., 14., -21./)
        norm = 231
      case default
        stop "Problem with sgolay"
    end select
  case (1)
    select case (n)
    case (2)
      coeffs = (/ -2., -1., 0., 1., 2./)
      norm = 10.
    case (3)
      coeffs = (/ -3., -2., -1., 0., 2., 1., 3./)
      norm = 28
    case (4)
      coeffs = (/ -4., -3., -2., -1., 0., 1., 2., 3., 4./)
      norm = 60.
    case default
      stop "Problem with sgolay"
    end select
  end select
return
end subroutine get_sgolay_coeffs
