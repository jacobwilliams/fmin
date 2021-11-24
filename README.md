Derivative free 1D function minimizer in modern Fortran

![Build Status](https://github.com/jacobwilliams/fmin/actions/workflows/CI.yml/badge.svg)

### Compiling

* A [Fortran Package Manager](https://github.com/fortran-lang/fpm) file is also included, so that the library and tests cases can be compiled with FPM. For example:

```
fpm build --profile release
fpm test --profile release
```

### Example

```fortran
program test

use fmin_module
use iso_fortran_env, only: wp => real64 ! double precision

real(wp) :: ax, bx, xmin, xerr

real(wp),parameter :: ax  = -4.0_wp       ! lower bound
real(wp),parameter :: bx  = 0.0_wp        ! upper bound
real(wp),parameter :: tol = 1.0e-8_wp     ! tolerance
real(wp),parameter :: pi  = acos(-1.0_wp) ! pi
real(wp),parameter :: x  = -pi/2.0_wp     ! true answer

xmin = fmin(func,ax,bx,tol) ! compute the minimum

xerr = xmin - x  ! difference from true value

write(*,*) 'xmin       = ', xmin
write(*,*) 'xmin exact = ', x
write(*,*) 'xmin error = ', xerr

contains

    function func(x) result(f)

    implicit none

    real(wp),intent(in) :: x  !! indep. variable
    real(wp)            :: f  !! function value `f(x)`

    f = sin(x)

    end function func

end program test
```

The output is:

```text
 xmin       =   -1.5707963254967554
 xmin exact =   -1.5707963267948966
 xmin error =    1.2981411501300499E-009
```

### Documentation

 * The API documentation for the current ```master``` branch can be found [here](https://jacobwilliams.github.io/fmin/).  This is generated by processing the source files with [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

### License

 * The Fmin source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/fmin/blob/master/LICENSE) (BSD-3).

### See also

  * Richard Brent, "[Algorithms for Minimization Without Derivatives](https://maths-people.anu.edu.au/~brent/pub/pub011.html)",
    Prentice - Hall, Inc. (1973)
  * [fmin from Netlib](http://www.netlib.org/fmm/fmin.f)
  * [roots-fortran](https://github.com/jacobwilliams/roots-fortran) (1D derivative-free roots solvers)
