! ans(i+1/2)

subroutine nabla(dxinv, m2, m1, p1, p2, ans)
    implicit none
    double precision :: dxinv, m2, m1, p1, p2, ans
    double precision :: mm2, mm1, pp1, pp2

    mm2 = 1.0d0/24.0d0*dxinv*m2
    mm1 = 1.0d0/24.0d0*dxinv*m1*27.0d0
    pp1 = 1.0d0/24.0d0*dxinv*p1*27.0d0
    pp2 = 1.0d0/24.0d0*dxinv*p2

    ans = mm2 - mm1 + pp1 - pp2
end subroutine nabla
