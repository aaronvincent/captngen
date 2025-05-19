!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This is for getting & storing the alpha and kappa tables

module alpha_kappa_mod
  implicit none

  double precision :: dm_target_ratio(100), alpha_loaded(100), kappavect(100)
  integer, parameter :: alpha_kappa_length =100 !because lazy
contains





  subroutine interpolate_1d(x, y, length, interp_point, interp_value)
    !1d interpolation assuming monotonically increasing vector
    integer, intent(in) :: length
    double precision, intent(in) :: x(length), y(length), interp_point
    double precision :: interp_value
    integer :: i

    if (interp_point .lt. x(1)) then
      print*,"xout ", interp_point, "xin min ", x(1)
      stop "Error in interpolation: xout < min(xin)"
    end if
    if (interp_point .gt. x(length)) then
    print*,"xout ", interp_point, "xin max ", x(length), length
    stop "Error in interpolation: xout > max(xin)"
  end if

    i = 1
    do while (interp_point .gt. x(i))
      i = i+1
    end do

    interp_value = y(i-1)+(y(i)-y(i-1))/(x(i)-x(i-1))*(interp_point-x(i-1))
    return
  end subroutine interpolate_1d

    


end module alpha_kappa_mod


subroutine read_alpha_kappa(q_pow, v_pow)
  use alpha_kappa_mod
  implicit none
  !v^{2*v_pow}, q^{2*q_pow}
  integer, intent(in) :: q_pow, v_pow
  integer i
  character*300 afilename, kfilename

  open(99,file="ak_files/mVect.dat")
  do i=1,alpha_kappa_length
  read(99,*) dm_target_ratio(i)
  end do
  close(99)

  if ((q_pow .eq. 0) .and. (v_pow .eq. 0)) then
     afilename = "ak_files/aVect_0.dat"
     kfilename = "ak_files/kVect_0.dat"
   else if (q_pow == 1) then
     afilename = "ak_files/aVect_q2.dat"
     kfilename = "ak_files/kVect_q2.dat"
   else if (q_pow == 2) then
     afilename = "ak_files/aVect_q4.dat"
     kfilename = "ak_files/kVect_q4.dat"
   else if (q_pow == -1) then
     afilename = "ak_files/aVect_qm2.dat"
     kfilename = "ak_files/kVect_qm2.dat"
   else if (v_pow == 1) then
     afilename = "ak_files/aVect_v2.dat"
     kfilename = "ak_files/kVect_v2.dat"
   else if (v_pow == 2) then
     afilename = "ak_files/aVect_v4.dat"
     kfilename = "ak_files/kVect_v4.dat"
   else if (v_pow == -1) then
     afilename = "ak_files/aVect_vm2.dat"
     kfilename = "ak_files/kVect_vm2.dat"
   end if


  open(99,file=afilename)
  open(95,file=kfilename)
  do i=1,alpha_kappa_length
  read(99,*) alpha_loaded(i)
  read(95,*) kappavect(i)
  end do
  close(99)
  close(95)

  return
end subroutine read_alpha_kappa
