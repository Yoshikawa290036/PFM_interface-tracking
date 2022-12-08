subroutine calphi(ni, nj, u, v, dxinv, dyinv, phi, dt, a, b, temperature, kappa)
    implicit none
    integer :: ni, nj
    double precision :: u, v, dxinv, dyinv
    double precision :: phi(-6:ni+7, -6:nj+7)
    double precision :: phipx(-6:ni+7, -6:nj+7), dphix, dphipx
    double precision :: phipy(-6:ni+7, -6:nj+7), dphiy, dphipy

    double precision :: ddphix, ddphiy, ddphi(-6:ni+7, -6:nj+7)
    ! double precision :: nphi1(-6:ni+7, -6:nj+7), nphi2(-6:ni+7, -6:nj+7)
    double precision :: dt
    double precision :: a, b, temperature, kappa
    integer :: i, j
    double precision :: advx(-6:ni+7, -6:nj+7), chpx(-6:ni+7, -6:nj+7)
    double precision :: advy(-6:ni+7, -6:nj+7), chpy(-6:ni+7, -6:nj+7)
    double precision :: Jpx(-6:ni+7, -6:nj+7), zetapx
    double precision :: Jpy(-6:ni+7, -6:nj+7), zetapy
    double precision :: Gamma
    double precision :: m1x, m2x, p0x, p1x, p2x
    double precision :: m1y, m2y, p0y, p1y, p2y
    double precision :: dddphipx, dchpx
    double precision :: dddphipy, dchpy

    Gamma = 12.0d0

    ! ddphi = nabla nabla phi
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, m1x, m2x, p0x, p1x, p2x, m1y, m2y, p0y, p1y, p2y, ddphix, ddphiy)
    do j = -4, nj+5
        do i = -4, ni+5
            m2x = -1.0d0/24.0d0*dxinv**2*phi(i-2, j)
            m1x =  1.0d0/24.0d0*dxinv**2*phi(i-1, j)*16.0d0
            p0x = -1.0d0/24.0d0*dxinv**2*phi(i  , j)*30.0d0
            p1x =  1.0d0/24.0d0*dxinv**2*phi(i+1, j)*16.0d0
            p2x = -1.0d0/24.0d0*dxinv**2*phi(i+2, j)
            ddphix = m2x+m1x+p0x+p1x+p2x

            m2y = -1.0d0/24.0d0*dxinv**2*phi(i, j-2)
            m1y =  1.0d0/24.0d0*dxinv**2*phi(i, j-1)*16.0d0
            p0y = -1.0d0/24.0d0*dxinv**2*phi(i, j  )*30.0d0
            p1y =  1.0d0/24.0d0*dxinv**2*phi(i, j+1)*16.0d0
            p2y = -1.0d0/24.0d0*dxinv**2*phi(i, j+2)
            ddphiy = m2y+m1y+p0y+p1y+p2y
            ddphi(i,j) = ddphix + ddphiy
        end do
    end do
!$OMP  END PARALLEL DO


! cal x direction of advection and chemical potential term
! ========================================================
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)
    do j = -6, nj+7
        do i = -5, ni+5
            phipx(i, j) = 1.0d0/16.0d0* &
                    & (-phi(i-1, j)+9.0d0*phi(i, j)+9.0d0*phi(i+1, j)-phi(i+2, j))
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, zetapx, dphipx, dddphipx)
    do j = -4, nj+5
        do i = -2, ni+3
            zetapx = (temperature/((1.0d0-b*phipx(i, j))**2))-2.0d0*a*phipx(i, j)
            call nabla(dxinv, phi(i-1,j)  , phi(i,j)  , phi(i+1,j)  , phi(i+2,j)  , dphipx  )
            call nabla(dxinv, ddphi(i-1,j), ddphi(i,j), ddphi(i+1,j), ddphi(i+2,j), dddphipx)

            Jpx(i, j) = -zetapx*dphipx+kappa*phipx(i, j)*dddphipx
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, dphix, dchpx)
    do j = -1, nj+3
        do i = 0, ni+2
            call nabla(dxinv, phipx(i-2,j), phipx(i-1,j), phipx(i,j), phipx(i+1,j), dphix)
            advx(i, j) = u*dphix

            call nabla(dxinv, phipx(i-2,j)*Jpx(i-2,j), phipx(i-1,j)*Jpx(i-1,j), phipx(i,j)*Jpx(i,j), phipx(i+1,j)*Jpx(i+1,j), dchpx)
            chpx(i,j) = gamma*dchpx

            ! nphi1(i,j) = phi(i,j)-dt*(advx(i,j)+chpx(i,j))
        end do
    end do
!$OMP  END PARALLEL DO

! cal y direction of advection and chemical potential term
! ========================================================
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)
    do j = -5, nj+5
        do i = -6, ni+7
            phipy(i, j) = 1.0d0/16.0d0* &
                    & (-phi(i, j-1)+9.0d0*phi(i, j)+9.0d0*phi(i, j+1)-phi(i, j+2))
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, zetapy, dphipy, dddphipy)
    do j = -2, nj+3
        do i = -4, ni+5
            zetapy = (temperature/((1.0d0-b*phipy(i, j))**2))-2.0d0*a*phipy(i, j)
            call nabla(dyinv, phi(i,j-1)  , phi(i,j)  , phi(i,j+1)  , phi(i,j+2)  , dphipy  )
            call nabla(dyinv, ddphi(i,j-1), ddphi(i,j), ddphi(i,j+1), ddphi(i,j+2), dddphipy)

            Jpy(i, j) = -zetapy*dphipy+kappa*phipy(i, j)*dddphipy
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, dphiy, dchpy)
    do j = 0, nj+2
        do i = -1, ni+3
            call nabla(dyinv, phipy(i,j-2), phipy(i,j-1), phipy(i,j), phipy(i,j+1), dphiy)
            advy(i, j) = v*dphiy

            call nabla(dyinv, phipy(i,j-2)*Jpy(i,j-2), phipy(i,j-1)*Jpy(i,j-1), phipy(i,j)*Jpy(i,j), phipy(i,j+1)*Jpy(i,j+1), dchpy)
            chpy(i,j) = gamma*dchpy

            ! nphi1(i,j) = phi(i,j)-dt*(advy(i,j)+chpy(i,j))
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
    do j = 0, nj+1
        do i = 0, ni+1
            phi(i,j) = phi(i,j)-dt*(advx(i,j)+advy(i,j)+chpx(i,j)+chpy(i,j))
            ! phi(i,j) = nphi1(i,j)
        end do
    end do
!$OMP  END PARALLEL DO

end subroutine calphi
