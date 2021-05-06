program celestial_bodies_simluation_by_Runge_Kutta_4th_order_method

        implicit none
        integer :: io = 0
        integer :: i, j, k
        integer :: nb_bodies = 0
        character(len=1),allocatable,dimension(:) :: body_type
        character(len=300) :: input_file = ""
        double precision,parameter :: G = 1d-1
        double precision :: rij = 0d0, rij3 = 0d0
        double precision :: fij(3) = 0d0
        double precision,allocatable,dimension(:) :: m
        double precision,allocatable,dimension(:,:) :: x, sx, v, sv, f, sf
        double precision,allocatable,dimension(:,:) :: kx1, kx2, kx3, kx4
        double precision,allocatable,dimension(:,:) :: kv1, kv2, kv3, kv4
        double precision :: kinetic_energy=0d0, potential_energy=0d0
        double precision,parameter :: dt = 1d-4
        double precision,parameter :: total_time = 100.0
        integer,parameter :: nb_max_steps = int(total_time / dt)
        integer,parameter :: output_interval = 1.0 / dt

        call get_command_argument(1, input_file)
        open(12,file=trim(input_file),status='old')
        do
         read(12,'()',iostat=io)
         if(io .ne. 0) exit
         nb_bodies = nb_bodies + 1
        end do
        close(12)

        allocate(body_type(nb_bodies))
        allocate(m(nb_bodies))
        allocate(x(nb_bodies,3))
        allocate(sx(nb_bodies,3))
        allocate(v(nb_bodies,3))
        allocate(sv(nb_bodies,3))
        allocate(f(nb_bodies,3))
        allocate(sf(nb_bodies,3))
        allocate(kx1(nb_bodies,3))
        allocate(kx2(nb_bodies,3))
        allocate(kx3(nb_bodies,3))
        allocate(kx4(nb_bodies,3))
        allocate(kv1(nb_bodies,3))
        allocate(kv2(nb_bodies,3))
        allocate(kv3(nb_bodies,3))
        allocate(kv4(nb_bodies,3))

        open(12,file=trim(input_file),status='old')
        do i = 1, nb_bodies
         read(12,*) body_type(i), m(i), x(i,:), v(i,:)
        end do
        close(12)

        do k = 1, nb_max_steps
         call calc_forces
         call calc_k1_k2_k3_k4
         call calc_next_positions
         call calc_velocities
!        if (mod(k,output_interval) .eq. 1) call calc_energy
         if (mod(k,output_interval) .eq. 1) call output
        end do

        deallocate(body_type,m,x,sx,v,sv,f,sf,kx1,kx2,kx3,kx4,kv1,kv2,kv3,kv4)

contains
subroutine calc_forces

        f = 0d0
        do i = 1, nb_bodies
         do j = i+1, nb_bodies

          rij = sqrt((x(i,1)-x(j,1))**2d0+(x(i,2)-x(j,2))**2d0+(x(i,3)-x(j,3))**2d0)

          rij3 = rij**3d0

          fij(1) = G*m(i)*m(j)*(x(j,1)-x(i,1))/rij3
          fij(2) = G*m(i)*m(j)*(x(j,2)-x(i,2))/rij3
          fij(3) = G*m(i)*m(j)*(x(j,3)-x(i,3))/rij3

          f(i,1) = f(i,1) + fij(1)
          f(i,2) = f(i,2) + fij(2)
          f(i,3) = f(i,3) + fij(3)
          f(j,1) = f(j,1) - fij(1)
          f(j,2) = f(j,2) - fij(2)
          f(j,3) = f(j,3) - fij(3)

         end do
        end do

end subroutine

subroutine calc_sforces

        sf = 0d0
        do i = 1, nb_bodies
         do j = i+1, nb_bodies

          rij = sqrt((sx(i,1)-sx(j,1))**2d0+(sx(i,2)-sx(j,2))**2d0+(sx(i,3)-sx(j,3))**2d0)

          rij3 = rij**3d0

          fij(1) = G*m(i)*m(j)*(sx(j,1)-sx(i,1))/rij3
          fij(2) = G*m(i)*m(j)*(sx(j,2)-sx(i,2))/rij3
          fij(3) = G*m(i)*m(j)*(sx(j,3)-sx(i,3))/rij3

          sf(i,1) = sf(i,1) + fij(1)
          sf(i,2) = sf(i,2) + fij(2)
          sf(i,3) = sf(i,3) + fij(3)
          sf(j,1) = sf(j,1) - fij(1)
          sf(j,2) = sf(j,2) - fij(2)
          sf(j,3) = sf(j,3) - fij(3)

         end do
        end do

end subroutine

subroutine calc_next_positions

        x = x + (kx1+2d0*kx2+2d0*kx3+kx4)/6d0*dt

end subroutine

subroutine calc_k1_k2_k3_k4

        kx1 = v
        do i = 1, nb_bodies
         kv1(i,1) = f(i,1) / m(i)
         kv1(i,2) = f(i,2) / m(i)
         kv1(i,3) = f(i,3) / m(i)
        end do

        sx = x + kx1 * 0.5d0 * dt

        call calc_sforces
        call calc_svelocities_2

        kx2 = sv
        do i = 1, nb_bodies
         kv2(i,1) = sf(i,1) / m(i)
         kv2(i,2) = sf(i,2) / m(i)
         kv2(i,3) = sf(i,3) / m(i)
        end do

        sx = x + kx2 * 0.5d0 * dt

        call calc_sforces
        call calc_svelocities_2

        kx3 = sv
        do i = 1, nb_bodies
         kv3(i,1) = sf(i,1) / m(i)
         kv3(i,2) = sf(i,2) / m(i)
         kv3(i,3) = sf(i,3) / m(i)
        end do

        sx = x + kx3 * dt

        call calc_sforces
        call calc_svelocities_1

        kx4 = sv
        do i = 1, nb_bodies
         kv4(i,1) = sf(i,1) / m(i)
         kv4(i,2) = sf(i,2) / m(i)
         kv4(i,3) = sf(i,3) / m(i)
        end do

end subroutine

subroutine calc_velocities

        v = v + (kv1+2d0*kv2+2d0*kv3+kv4)/6d0*dt

end subroutine

subroutine calc_svelocities_1

        do i = 1, nb_bodies
         sv(i,1) = v(i,1) + sf(i,1) * dt / m(i)
         sv(i,2) = v(i,2) + sf(i,2) * dt / m(i)
         sv(i,3) = v(i,3) + sf(i,3) * dt / m(i)
        end do

end subroutine

subroutine calc_svelocities_2

        do i = 1, nb_bodies
         sv(i,1) = v(i,1) + 0.5d0 * sf(i,1) * dt / m(i)
         sv(i,2) = v(i,2) + 0.5d0 * sf(i,2) * dt / m(i)
         sv(i,3) = v(i,3) + 0.5d0 * sf(i,3) * dt / m(i)
        end do

end subroutine

subroutine calc_energy

        kinetic_energy = 0d0
        potential_energy = 0d0

        do i = 1, nb_bodies
         kinetic_energy = kinetic_energy + 0.5d0 * (v(i,1)**2d0+v(i,2)**2d0+v(i,3)**2d0)
        end do

        do i = 1, nb_bodies
         do j = i+1, nb_bodies
          rij = sqrt((x(i,1)-x(j,1))**2d0+(x(i,2)-x(j,2))**2d0+(x(i,3)-x(j,3))**2d0)
          potential_energy = potential_energy - G * m(i) * m(j) / rij
         end do
        end do

        write(*,*) kinetic_energy, potential_energy, kinetic_energy+potential_energy 

end subroutine

subroutine output

        write(*,*) nb_bodies
        write(*,*)

        do i = 1, nb_bodies
        write(*,*) body_type(i), x(i,:)
        end do

end subroutine

end
