module neutronsrf

implicit none
public :: initial, mc_step_inelastic, deallocate_

integer, public :: n_p, n_r, n_t, n_a, n_th, n_s, target_
real, dimension(:), allocatable, public :: x, y, z, ts, energies, mean_free_path
integer, dimension(33), public :: seed
real, public :: t, energy
logical, public :: store

contains

subroutine initial()
    call random_seed(put = seed)
    allocate(x(0 : n_p - 1), y(0 : n_p - 1), z(0 : n_p - 1), ts(0 : n_p - 1), energies(0 : n_p - 1), mean_free_path(0 : n_p - 1)) 
    return
end subroutine

subroutine mc_step_inelastic()
        real, dimension(3) :: densities, molar_mass, sigma_thermal_s
        real, dimension(3) :: sigma_fast_s, sigma_fast_c, sigma_thermal_c
        real, dimension(3) :: xi, sigma_c_f, sigma_c_t
        character(len = 15) :: file_name
        integer :: i, j, k
        logical :: exit_
        real :: sigma_c, sigma_tot, e, rnd, ran1, ran2, ran3
        real :: lambda_, l, v, p_s, phi, cos_theta, theta, pi
        
        pi = 16.*atan(1./5.) - 4.*atan(1./239.)
        densities = [2.1, 7.9, 2.4]
        molar_mass = [12., 56., 10.]
        sigma_thermal_s = [5., 10., 200.]
        sigma_thermal_c = [0.002, 2., 5.]
        sigma_fast_s = [2., 20., 2.]
        sigma_fast_c = [0.0001, 0.003, 0.4]
        xi = [0.16, 0.04, 0.19]
        
        do i = 0, n_p-1
            sigma_c_t = sigma_thermal_c
            sigma_c_f = sigma_fast_c
            e = energy
            j = target_
            if (store .eqv. .true.) then
                write(unit = file_name, fmt="(i2.2,a)") i, ".txt" 
                file_name = "output"//file_name
                file_name = trim(file_name)
                open (unit = 1, file = file_name, action = "write", status = "replace")
                write(1, *) 0, 0, -0.25
                write(1, *) 0, 0, 0
            end if
            exit_ = .false.
            k = 0
            do while (exit_ .eqv. .false.)
                if (e < 0.025) then
                    n_th = n_th + 1
                    exit_ = .true.
                end if
                if (e > 10.) then
                    sigma_c = sigma_fast_c(j) * (1000./sqrt(e))
                    sigma_tot = sigma_fast_s(j) + sigma_c
                    p_s = sigma_fast_s(j)/sigma_tot
                    lambda_ = molar_mass(j)/(densities(j) * 60.2 * sigma_tot)
                else
                    sigma_c = sigma_thermal_c(j) * (3.2/sqrt(E))
                    sigma_tot = sigma_fast_s(j) + sigma_c
                    p_s = sigma_thermal_s(j)/sigma_tot
                    lambda_ = molar_mass(j)/(densities(j) * 60.2 * sigma_tot)
                end if
                call random_number(rnd)
                if (rnd <= p_s) then
                    n_s = n_s + 1
                    k = k + 1
                    e = (1. - xi(j)) * e
                    call random_number(ran1)
                    call random_number(ran2)
                    call random_number(ran3)
                    l = - lambda_ * log(ran1)
                    v = sqrt(2. * e) * 1000000.
                    ts(i) = ts(i) + l/v
                    mean_free_path(i) = (l + mean_free_path(i))/real(k)
                    phi = 2. * pi * ran2
                    cos_theta = 1. - 2. * ran3
                    theta = acos(cos_theta)
                    z(i) = l * cos_theta + z(i)
                    if (store .eqv. .true.) then
                        x(i) = l * sin(theta) * cos(phi) + x(i)
                        y(i) = l * sin(theta) * sin(phi) + y(i)
                        write(1,*) x(i), y(i), z(i)
                    end if

                    if (z(i) < 0) then
                        n_r = n_r + 1
                        exit_ = .true.
                    end if
                    if (z(i) > t) then
                        n_t = n_t + 1
                        exit_ = .true.
                    end if
                else
                    n_a = n_a + 1
                    exit_ = .true.
                end if
            end do
            energies(i) = e
            if (store .eqv. .true.) then
                close(1)
            end if
        end do
       
        return
end subroutine

subroutine deallocate_()
    deallocate(x, y, z, ts, energies, mean_free_path) 
end subroutine

end module neutronsrf