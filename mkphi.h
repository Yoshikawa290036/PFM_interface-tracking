write (fname, '("xyphi",i7.7)') step

open (10, file=fname)

do j = -6, nj+7
    do i = -6, ni+7
        x = (dble(i)-0.5d0)*dx
        y = (dble(j)-0.5d0)*dy
        write (10, '(20e20.10)') x, y, phi(i,j)
    end do
    write (10, *)
end do
