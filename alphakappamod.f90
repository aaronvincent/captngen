!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This is for getting & storing the alpha and kappa tables

module akmod
  implicit none

  double precision :: muvect(100), alphavect(100), kappavect(100)
  integer, parameter :: nlinesinaktable =100 !because lazy
contains





  subroutine interp1(xin,yin,lenxin,xout,yout)
    !1d interpolation assuming monotonically increasing vector
    integer, intent(in) :: lenxin
    double precision, intent(in) :: xin(lenxin), yin(lenxin), xout
    double precision :: yout
    integer :: i

    if (xout .lt. xin(1)) then
      print*,"xout ", xout, "xin min ", xin(1)
      stop "Error in interpolation: xout < min(xin)"
    end if
    if (xout .gt. xin(lenxin)) then
    print*,"xout ", xout, "xin max ", xin(lenxin), lenxin
    stop "Error in interpolation: xout > max(xin)"
  end if

    i = 1
    do while (xout .gt. xin(i))
      i = i+1
    end do

    yout = yin(i-1)+(yin(i)-yin(i-1))/(xin(i)-xin(i-1))*(xout-xin(i-1))
    return
  end subroutine interp1



end module akmod


subroutine get_alpha_kappa(nq,nv)
  use akmod
  implicit none
  !v^2nv, q^2nq
  integer, intent(in) :: nq, nv
  integer i
  character*300 afilename, kfilename

  open(99,file="ak_files/mVect.dat")
  do i=1,nlinesinaktable
  read(99,*) muVect(i)
  end do
  close(99)

  if ((nq .eq. 0) .and. (nv .eq. 0)) then
     afilename = "ak_files/aVect_0.dat"
     kfilename = "ak_files/kVect_0.dat"
   else if (nq == 1) then
     afilename = "ak_files/aVect_q2.dat"
     kfilename = "ak_files/kVect_q2.dat"
   else if (nq == 2) then
     afilename = "ak_files/aVect_q4.dat"
     kfilename = "ak_files/kVect_q4.dat"
   else if (nq == -1) then
     afilename = "ak_files/aVect_qm2.dat"
     kfilename = "ak_files/kVect_qm2.dat"
   else if (nv == 1) then
     afilename = "ak_files/aVect_v2.dat"
     kfilename = "ak_files/kVect_v2.dat"
   else if (nv == 2) then
     afilename = "ak_files/aVect_v4.dat"
     kfilename = "ak_files/kVect_v4.dat"
   else if (nv == -1) then
     afilename = "ak_files/aVect_vm2.dat"
     kfilename = "ak_files/kVect_vm2.dat"
   end if


  open(99,file=afilename)
  open(95,file=kfilename)
  do i=1,nlinesinaktable
  read(99,*) alphavect(i)
  read(95,*) kappavect(i)
  end do
  close(99)
  close(95)

  return
end subroutine get_alpha_kappa
