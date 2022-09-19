# libarvo

Library for calculating molecular surfaces and volumes using the algorithm implemented by J. Buša and coworkers.[^1] This library is not meant as a standalone application but rather as a resource for other packages. libarvo is (re)written in Fortran90 with a C interface and Python API.

## Installation

### pip

To install the Python API with the embedded libarvo shared library you can use pip.

```shell
pip install libarvo
```

It is also possible to install directly from the GitHub repository.
You will need a Fortran compiler such as GFortran to compile the shared library but the whole process is automated via [scikit-build](https://github.com/scikit-build/scikit-build).

```shell
pip install git+https://github.com/rlaplaza/libarvo
```

### conda

Another option to install the Python API that doesn't require a compiler is with conda.

```shell
conda install -c conda-forge libarvo
```

### cmake

The shared library can be built and installed with cmake. An example worklfow is given below where you need to replace $PREFIX with the desired directory.

```shell
FC=gfortran cmake -B _build -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release
cmake --build _build
cmake --install _build
```

### fpm

To use as a dependency in your Fortran project with [fpm](https://github.com/fortran-lang/fpm), add the following to your `fpm.toml`.

```toml
[dependencies]
libarvo.git = "https://github.com/rlaplaza/libconeangle"
```

## Usage

### Python API

There is only one function: `molecular_vs`. An example is given below for PdCO.

> ⚠️ All atoms are zero-index in the Python API.

```python
>>> from libarvo import molecular_vs
>>> import numpy as np
>>> coordinates =  np.array([[0.0, 0.0, -0.52], [0.0, 0.0, 1.76], [0.0, 0.0, 2.86]])
>>> radii = np.array([2.1, 1.7, 1.52])
>>> index_metal = 0 # Zero-indexed
>>> c_angle, axis, tangent_atoms = cone_angle(coordinates, radii, index_metal)
>>> c_angle
96.4237340645161
>>> axis
array([0., 0., 1.])
>>> tangent_atoms # Also zero-indexed
[1]
```

### Fortran API

The Fortran API exposes the function `arvo` with the follow signature.

> ⚠️ All atoms are one-index in the Fortran API.

```fortran
subroutine arvo(coordinates, radii, index_metal, alpha, axis, tangent_atoms, stat, errmsg)
  !! Calculate cone angle, cone axis and tangent atoms
  !> Coordinates (Å) (shape: 3, number of atoms)
  real(wp), intent(in) :: coordinates(:, :)
  !> vdW radii (Å)
  real(wp), intent(in) :: radii(:)
  !> Index of metal atom
  integer, intent(in) :: index_metal
  !> Cone angle (degrees)
  real(wp), intent(out) :: alpha
  !> Cone axis (Å)
  real(wp), intent(out) :: axis(3)
  !> Indices of atoms tangent to cone
  integer, intent(out) :: tangent_atoms(3)
  !> Return code
  integer, intent(out) :: stat
  !> Error message
  character(:), allocatable, intent(out) :: errmsg
  ...
  
end subroutine arvo
```

Here is one example of how it could be used as given in the demo [program](app/demo.f90).

```fortran
program demo
  use coneangle_main, only: cone_angle

  integer, parameter :: dp = selected_real_kind(15, 307)
  real(dp) :: coordinates(3, 3), radii(3), alpha, axis(3)
  integer :: tangent_atoms(3), stat
  character(:), allocatable :: errmsg

  coordinates = reshape([0._dp, 0._dp, -0.52_dp, 0._dp, 0._dp, 1.76_dp, 0._dp, 0._dp, 2.86_dp], [3, 3])
  radii = [2.1_dp, 1.7_dp, 1.52_dp]
  call cone_angle(coordinates, radii, 1, alpha, axis, tangent_atoms, stat, errmsg)
  write (*, *) "Cone angle:", alpha
  write (*, *) "Cone axis:", axis
  write (*, *) "Tangent atoms:", tangent_atoms
end program demo
```

A minimal [FORD](https://github.com/Fortran-FOSS-Programmers/ford) documentation can be built from `pages.md`

### C API

The C API exposes the subroutine `arvo_c` with the C name `arvo`. It's signature is the same as for the Fortran subroutine. 

> ⚠️ All atoms are zero-index in the C API.

```fortran
subroutine arvo_c(n_atoms, coordinates, radii, index_metal, alpha, axis, tangent_atoms, stat, errmsg) &
  bind(c, name="cone_angle")
  !! Calculate cone angle, cone axis and tangent atoms
  !> Number of atoms
  integer(c_int), value, intent(in) :: n_atoms
  !> Coordinates (Å)
  real(c_double), intent(in) :: coordinates(3, n_atoms)
  !> vdW radii (Å)
  real(c_double), intent(in) :: radii(n_atoms)
  !> Index of metal atom
  integer(c_int), value, intent(in) :: index_metal
  !> Cone angle (degrees)
  real(c_double), intent(out) :: alpha
  !> Cone axis (Å)
  real(c_double), intent(out) :: axis(3)
  !> Indices of atoms tangent to cone
  integer(c_int), intent(out) :: tangent_atoms(3)
  !> Return code
  integer(c_int), intent(out) :: stat
  !> Error message
  character(c_char), intent(out) :: errmsg(*)
  
  ...

end subroutine arvo_c
```

The C header file can be found [here](include/cone_angle.h). An [example](libarvo/lib.py) is given for loading the shared library with ctypes in the Python API.

## Acknowledgements

Any published work derived from the use of libarvo should cite the original publications for the original F77 implementation and algorithm.[^1]

## References

[^1]: J. Buša, J. Džurina, E. Hayryan, S. Hayryan, C.-K. Hu, J. Plavka, I. Pokorný, J. Skřivánek and M.-C. Wu, *Computer Physics Communications*, **2005**, *165*, 59–96. https://doi.org/10.1016/j.cpc.2004.08.002

