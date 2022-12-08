program main
    implicit none

    integer :: ni, nj
    integer :: i, j
    double precision :: x, y
    double precision :: phimin, phimax
    double precision :: dx, dxinv, dy, dyinv
    double precision :: xl, yl
    double precision :: t, dt
    double precision :: a, b, temperature, kappa
    double precision :: u, v
    double precision, dimension(:, :), allocatable :: phi
    integer :: maxstep, step
    integer :: dataou
    character(32) fname

    dataou = 100
    maxstep = 40000

    ni = 64
    nj = 64

    a = 1.0d0
    b = 1.0d0
    kappa = 0.1d0
    temperature = 0.293d0

    dt = 2.5e-2
    dx = 1.0d0
    dy = 1.0d0

    phimin = 0.265d0
    phimax = 0.405d0

    xl = dx*dble(ni)
    yl = dy*dble(nj)
    dxinv = 1.0d0/dx
    dyinv = 1.0d0/dy

    u = 0.5d0
    v = 0.5d0

    include'allocate.h'
    write (*, '("Courant Number      ",20e20.10)') abs(u*dt/dx)
    call init(ni, nj, dx, dy, phi, phimin, phimax, 6.0d0)
    call bndset(ni, nj, phi)
    include'mkphi.h'

    do step = 1, maxstep
        call calphi(ni, nj, u, v, dxinv, dyinv, phi, dt, a, b, temperature, kappa)

        call bndset(ni, nj, phi)
        if (mod(step, dataou) == 0) then
            include'mkphi.h'
        end if
    end do

end program main
