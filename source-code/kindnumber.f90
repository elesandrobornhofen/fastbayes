module kindnumber
!
! DESCRIPTION:
!   specify precision of numbers
!
! REVISIONS:
!   02/06/09 wittenburg
!
  implicit none

  integer,parameter :: ik0=selected_int_kind(1)
  integer,parameter :: ik2=selected_int_kind(4)
  integer,parameter :: ik4=selected_int_kind(9)
  integer,parameter :: ikxl=selected_int_kind(12)
  integer,parameter :: rk4=selected_real_kind(6,37)
  integer,parameter :: rk8=selected_real_kind(12,60)
!  integer,parameter :: rk8=selected_real_kind(15,307)

end module kindnumber
