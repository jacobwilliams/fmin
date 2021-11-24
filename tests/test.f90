!*****************************************************************************************
!>
!  Tests for [[fmin_module]].

    program test

    use fmin_module
    use iso_fortran_env, only: wp => real64 ! double precision


    real(wp) :: ax, bx, xmin, xerr, x

    real(wp),parameter :: pi  = acos(-1.0_wp)
    real(wp),parameter :: tol = 1.0e-8_wp

    ax = -4.0_wp
    bx = 0.0_wp
    x = -pi/2.0_wp  ! actual answer

    xmin = fmin(func,ax,bx,tol)

    xerr = xmin - x  ! difference from try value

    write(*,*) 'xmin       = ', xmin
    write(*,*) 'xmin exact = ', x
    write(*,*) 'xmin error = ', xerr

    if (abs(xerr) > 10.0_wp * tol) then
        error stop 'test failed'
    end if

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Test function to minimize.

    function func(x) result(f)

    implicit none

    real(wp),intent(in) :: x  !! indep. variable
    real(wp)            :: f  !! function value `f(x)`

    f = sin(x)

    end function func
!*****************************************************************************************

!*****************************************************************************************
    end program test
!*****************************************************************************************