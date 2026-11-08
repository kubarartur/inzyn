submodule(hamiltonian) submod_laplaceian
    implicit none
contains
    module subroutine laplaceian(gr, psi_in, psi_out)
        type(grid_params), intent(in) :: gr
        complex(wp), intent(in) :: psi_in(:, :, :)
        complex(wp), intent(out) :: psi_out(:, :, :)

        integer :: endx
        integer :: endy
        integer :: endz
        integer :: ii, jj, kk

        endx = size(psi_in, 1)
        endy = size(psi_in, 2)
        endz = size(psi_in, 3)

            do kk = 2, endz - 1
            do jj = 2, endy - 1
            do ii = 2, endx - 1
                psi_out(ii, jj, kk) = &
                    (psi_in(ii + 1, jj, kk) - 2._wp*psi_in(ii, jj, kk) + psi_in(ii - 1, jj, kk))/gr%dr(1)**2 + &
                    (psi_in(ii, jj + 1, kk) - 2._wp*psi_in(ii, jj, kk) + psi_in(ii, jj - 1, kk))/gr%dr(2)**2 + &
                    (psi_in(ii, jj, kk + 1) - 2._wp*psi_in(ii, jj, kk) + psi_in(ii, jj, kk - 1))/gr%dr(3)**2
            end do
            end do
            end do
    end subroutine laplaceian
end submodule
