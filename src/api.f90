!> C API
module arvo_api
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_char, c_null_char
  use arvo_main, only: arvo
  implicit none(type, external)

contains
  subroutine arvo_c(n_atoms, coordinates, radii, probe_radius, V, S, stat, errmsg) &
    bind(c, name="arvo")
    !! Calculate molecular volume and surface
    !> Number of atoms
    integer(c_int), value, intent(in) :: n_atoms
    !> Coordinates (Å)
    real(c_double), intent(in) :: coordinates(3, n_atoms)
    !> vdW radii (Å)
    real(c_double), intent(in) :: radii(n_atoms)
    !> Probe radius
    real(c_double), value, intent(in) :: probe_radius
    !> Molecular volume
    real(c_double), intent(out) :: V
    !> Molecular surface
    real(c_double), intent(out) :: S
    !> Return code
    integer(c_int), intent(out) :: stat
    !> Error message
    character(c_char), intent(out) :: errmsg(*)

    integer :: i
    character(:), allocatable :: errmsg_f

    ! Call subroutine
    call arvo(n_atoms, coordinates, radii, probe_radius, V, S, stat, errmsg_f)

    ! Convert error message to C format.
    if (allocated(errmsg_f)) then
      do i = 1, len(errmsg_f)
        errmsg(i) = errmsg_f(i:i)
      end do
      errmsg(len(errmsg_f) + 1) = c_null_char
    end if

  end subroutine arvo_c
end module arvo_api
