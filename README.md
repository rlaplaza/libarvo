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
libarvo.git = "https://github.com/rlaplaza/libarvo"
```

## Usage

### Python API

There is only one function: `molecular_vs`. An example is given below for PdCO.

```python
>>> from libarvo import molecular_vs
>>> import numpy as np
>>> coordinates =  np.array([[0.0, 0.0, -0.52], [0.0, 0.0, 1.76], [0.0, 0.0, 2.86]])
>>> radii = np.array([2.1, 1.7, 1.52])
>>> probe_radius = 0.0
>>> volume, surface = molecular_vs(coordinates, radii, probe_radius)
>>> volume
59.465621235350206
>>> surface
81.70751658028372
```

### Fortran API

The Fortran API exposes the function `arvo` with the follow signature.

```fortran
subroutine arvo(ns, sc, sr, pr, volume, area, stat, errmsg)
  !! Computing surface area and volume of the overlapping spheres
  !> Number of spheres
  integer, intent(in) :: ns
  !> sc : sphere coordinates (Å) (shape: 3, number of atoms)
  !> sr : sphere radii (Å) (shape: number of atoms)
  !> pr : probe radius (Å)
  real(wp), intent(in) :: sc(3, ns), sr(ns), pr
  !> surface area and volume (Å² and Å³ respectively) 
  real(wp), intent(out) ::  area, volume
  !> Return code
  integer, intent(out) :: stat
  !> Error message
  character(len=:), allocatable, intent(out) :: errmsg
  ...
  
end subroutine arvo
```

A minimal [FORD](https://github.com/Fortran-FOSS-Programmers/ford) documentation can be built from `pages.md`

### C API

The C API exposes the subroutine `arvo_c` with the C name `arvo`. It's signature is the same as for the Fortran subroutine. 

```fortran
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
  ...

end subroutine arvo_c
```

The C header file can be found [here](include/arvo.h). An [example](libarvo/lib.py) is given for loading the shared library with ctypes in the Python API.

## Acknowledgements

Any published work derived from the use of libarvo should cite the original publications for the original F77 implementation and algorithm.[^1]

The [libconeangle](https://github.com/kjelljorner/libconeangle) library by K. Jorner is acknowledged for inspiration regarding packaging, distribution and coding style. 

## References

[^1]: J. Buša, J. Džurina, E. Hayryan, S. Hayryan, C.-K. Hu, J. Plavka, I. Pokorný, J. Skřivánek and M.-C. Wu, *Computer Physics Communications*, **2005**, *165*, 59–96. https://doi.org/10.1016/j.cpc.2004.08.002

