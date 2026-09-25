pure subroutine lpalaceian(psi_, gr, psi_copy_)
    !!compute the laplaceian of psi_
    !!write the result to psi_copy_
    !!assumes the gr is in the parent scope
    complex(wp), intent(in) :: psi_
    type(grid_params), intent(in) :: gr
    complex(wp), intent(out) :: psi_copy_

    integer :: ii, jj, kk
    integer :: iend

    !first handle the inside
    do kk = 2, size(psi, 3) - 1
    do jj = 2, size(psi, 2) - 1
    do ii = 2, size(psi, 1) - 1

        psi_copy_(ii, jj, kk) = &
            (psi_(ii + 1, jj, kk) - 2._wp*psi_(ii, jj, kk) + psi_(ii - 1, jj, kk))/gr%dr(1)**2 + &
            (psi_(ii, jj + 1, kk) - 2._wp*psi_(ii, jj, kk) + psi_(ii, jj - 1, kk))/gr%dr(2)**2 + &
            (psi_(ii, jj, kk + 1) - 2._wp*psi_(ii, jj, kk) + psi_(ii, jj, kk - 1))/gr%dr(3)**2
    end do
    end do
    end do

    !TODO: Conider something more inteligent and compact
    !loop sides to be evaluated with forward and backward differences

    !x
    do kk = 2, size(psi, 3) - 1
    do jj = 2, size(psi, 2) - 1
        !forward difference
        psi_copy_(1, jj, kk) = &
            (2._wp*psi_(1, jj, kk) - 5._wp*psi_(2, jj, kk) + &
             4._wp*psi_(3, jj, kk) - psi_(4, jj, kk))/gr%dr(1)**2 + &
            (psi_(1, jj + 1, kk) - 2._wp*psi_(1, jj, kk) + psi_(1, jj - 1, kk))/gr%dr(2)**2 + &
            (psi_(1, jj, kk + 1) - 2._wp*psi_(1, jj, kk) + psi_(1, jj, kk - 1))/gr%dr(3)**2

        iend = size(psi_copy_, 1)

        !backward difference
        psi_copy_(iend, jj, kk) = &
            (2._wp*psi_(iend, jj, kk) - 5._wp*psi_(iend - 1, jj, kk) + &
             4._wp*psi_(iend - 2, jj, kk) - psi_(iend - 3, jj, kk))/gr%dr(1)**2 + &
            (psi_(iend, jj + 1, kk) - 2._wp*psi_(iend, jj, kk) + psi_(iend, jj - 1, kk))/gr%dr(2)**2 + &
            (psi_(iend, jj, kk + 1) - 2._wp*psi_(iend, jj, kk) + psi_(iend, jj, kk - 1))/gr%dr(3)**2
    end do
    end do
    !y
    do kk = 2, size(psi, 3) - 1
    do ii = 2, size(psi, 1) - 1
        psi_copy_(ii, jj, kk) = &
            (psi_(ii + 1, 1, kk) - 2._wp*psi_(ii, 1, kk) + psi_(ii - 1, 1, kk))/gr%dr(1)**2 + &
            (2._wp*psi_(ii, 1, kk) - 5._wp*psi_(ii, 2, kk) + &
             4._wp*psi_(ii, 3, kk) - psi_(ii, 4, kk))/gr%dr(2)**2 + &
            (psi_(ii, 1, kk + 1) - 2._wp*psi_(ii, 1, kk) + psi_(ii, 1, kk - 1))/gr%dr(3)**2

        iend = size(psi_copy_, 2)

        psi_copy_(ii, jj, kk) = &
            (psi_(ii + 1, iend, kk) - 2._wp*psi_(ii, iend, kk) + psi_(ii - 1, iend, kk))/gr%dr(1)**2 + &
            (2._wp*psi_(ii, iend, kk) - 5._wp*psi_(ii, iend - 1, kk) + &
             4._wp*psi_(ii, iend - 2, kk) - psi_(ii, iend - 3, kk))/gr%dr(2)**2 + &
            (psi_(ii, iend, kk + 1) - 2._wp*psi_(ii, iend, kk) + psi_(ii, iend, kk - 1))/gr%dr(3)**2
    end do
    end do
    !z
    do jj = 2, size(psi, 2) - 1
    do ii = 2, size(psi, 1) - 1
        psi_copy_(ii, jj, kk) = &
            (psi_(ii + 1, jj, 1) - 2._wp*psi_(ii, jj, 1) + psi_(ii - 1, jj, 1))/gr%dr(1)**2 + &
            (psi_(ii, jj + 1, 1) - 2._wp*psi_(ii, jj, 1) + psi_(1, jj - 1, 1))/gr%dr(2)**2 + &
            (2._wp*psi_(ii, jj, 1) - 5._wp*psi_(ii, jj, 2) + &
             4._wp*psi_(ii, jj, 3) - psi_(ii, jj, 4))/gr%dr(3)**2

        iend = size(psi_copy_, 3)

        psi_copy_(ii, jj, kk) = &
            (psi_(ii + 1, jj, iend) - 2._wp*psi_(ii, jj, iend) + psi_(ii - 1, jj, iend))/gr%dr(1)**2 + &
            (psi_(ii, jj + 1, iend) - 2._wp*psi_(ii, jj, iend) + psi_(ii, jj - 1, iend))/gr%dr(2)**2 + &
            (2._wp*psi_(ii, jj, iend) - 5._wp*psi_(ii, jj, iend - 1) + &
             4._wp*psi_(ii, jj, iend - 2) - psi_(ii, jj, iend - 3))/gr%dr(3)**2
    end do
    end do
end subroutine laplaceian
