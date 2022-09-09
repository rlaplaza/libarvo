module test_integration
  use testdrive, only: new_unittest, unittest_type, error_type, check
  use arvo_env, only: wp
  use arvo_main, only: arvo
  implicit none(type, external)
  private

  public :: collect_suite_integration

contains

!> Collect all exported unit tests
  subroutine collect_suite_integration(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
                new_unittest("One sphere", test_onesphere), &
                new_unittest("Three spheres", test_threespheres), &
                new_unittest("Three overlapping spheres", test_threeospheres) &
                ]

  end subroutine collect_suite_integration

  subroutine test_onesphere(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: coordinates(1, 3), radii(1), probe_radius, V, S
    integer :: stat, n
    character(:), allocatable :: errmsg

    n = 1
    probe_radius = 0._wp
    coordinates = reshape([0._wp, 0._wp, 0._wp], [1, 3])
    radii = [1.7_wp]
    call arvo(n, coordinates, radii, probe_radius, V, S, stat, errmsg)

    call check(error, V, 20.5795_wp, thr=0.001_wp)
    if (allocated(error)) return
    call check(error, S, 36.3168_wp, thr=0.001_wp)
    if (allocated(error)) return

  end subroutine test_onesphere

  subroutine test_threespheres(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: coordinates(3, 3), radii(3), probe_radius, V, S
    integer :: stat, n
    character(:), allocatable :: errmsg

    n = 3
    probe_radius = 0._wp
    coordinates = reshape([0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 3.4_wp, 3.4_wp, 0._wp, 0._wp], [3, 3])
    radii = [1.7_wp, 1.7_wp, 1.7_wp]
    call arvo(n, coordinates, radii, probe_radius, V, S, stat, errmsg)

    call check(error, V, 61.7386_wp, thr=0.001_wp)
    if (allocated(error)) return

    ! Check stat calculation failed.
    call check(error, stat, 0, message=errmsg)
    if (allocated(error)) return

  end subroutine test_threespheres

  subroutine test_threeospheres(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: coordinates(3, 3), radii(3), probe_radius, V, S
    integer :: stat, n
    character(:), allocatable :: errmsg

    n = 3
    probe_radius = 1.20_wp
    coordinates = reshape([0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 1.90_wp, 1.90_wp, 0._wp, 0._wp], [3, 3])
    radii = [1.7_wp, 1.7_wp, 1.7_wp]
    call arvo(n, coordinates, radii, probe_radius, V, S, stat, errmsg)

    ! Check stat calculation failed.
    call check(error, stat, 0, message=errmsg)
    if (allocated(error)) return

  end subroutine test_threeospheres
  
end module test_integration
