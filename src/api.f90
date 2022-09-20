!> C API
module arvo_api
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_char, c_null_char
  use arvo_main, only: arvo
  implicit none(type, external)

contains
  subroutine arvo_c(n_atoms, coordinates, radii, probe_radius, volume, area, ns_v, ns_a, stat, errmsg) &
    bind(c, name="arvo")
    !! Calculate molecular volume and surface
    !> Number of atoms
    integer(c_int), value, intent(in) :: n_atoms
    !> Coordinates (Å)
    real(c_double), intent(in) :: coordinates(3, n_atoms)
    !> vdW radii (Å)
    real(c_double), intent(in) :: radii(n_atoms)
    !> Probe radius (Å)
    real(c_double), value, intent(in) :: probe_radius
    !> Molecular volume (Å^3)
    real(c_double), intent(out) :: volume
    !> Molecular surface (Å^2)
    real(c_double), intent(out) :: area
    !> Volume per atom (Å^3)
    real(c_double), intent(out) :: ns_v(n_atoms)
    !> Surface per atom (Å^2)
    real(c_double), intent(out) :: ns_a(n_atoms)
    !> Return code
    integer(c_int), intent(out) :: stat
    !> Error message
    character(c_char), intent(out) :: errmsg(*)

    integer :: i
    character(:), allocatable :: errmsg_f

    ! Call subroutine
    call arvo(n_atoms, coordinates, radii, probe_radius, volume, area, ns_v, ns_a, stat, errmsg_f)

    ! Convert error message to C format.
    if (allocated(errmsg_f)) then
      do i = 1, len(errmsg_f)
        errmsg(i) = errmsg_f(i:i)
      end do
      errmsg(len(errmsg_f) + 1) = c_null_char
    end if

  end subroutine arvo_c
end module arvo_api
