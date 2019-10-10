! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE hybrid
#ifdef HYBRID

  USE balance
  USE current_smooth
  USE diagnostics
  USE injectors
  USE particles
  USE particle_migration
  USE random_generator
#ifdef LANDAU_LIFSHITZ
  USE landau_lifshitz
#endif
#ifdef PHOTONS
  USE photons
#endif
#ifdef BREMSSTRAHLUNG
  USE bremsstrahlung
#endif

  IMPLICIT NONE

CONTAINS

  SUBROUTINE run_hybrid_PIC(push, halt, force_dump)

    REAL(num), ALLOCATABLE :: heat_capacity(:,:)
    LOGICAL :: push, halt, force_dump

    ! A copy of the EPOCH PIC loop, modified to run in the hybrid scheme
    DO
      ! Check we have not passed the end condition
      IF ((step >= nsteps .AND. nsteps >= 0) &
          .OR. (time >= t_end) .OR. halt) EXIT

      ! Timing information.
      ! Functions/subroutines found in housekeeping/timer.f90
      IF (timer_collect) THEN
        CALL timer_stop(c_timer_step)
        CALL timer_reset
        timer_first(c_timer_step) = timer_walltime
      END IF

      ! Radiation scripts
#ifdef LANDAU_LIFSHITZ
      ! Landau-Lifshitz classical radiation reaction
      IF (use_landau_lifshitz .AND. time > landau_lifshitz_start_time &
          .AND. push) THEN
        CALL classical_radiation_reaction()
      END IF
#endif
#ifdef PHOTONS
      ! Non-linear Compton scatter/synchrotron radiation calculation (photons
      ! can be generated)
      IF (use_qed .AND. time > qed_start_time .AND. push) THEN
        CALL qed_update_optical_depth()
      END IF
#endif
#ifdef BREMSSTRAHLUNG
      ! Bremsstrahlung radiation calculation (photons can be generated)
      IF (use_bremsstrahlung .AND. time > bremsstrahlung_start_time &
          .AND. push) THEN
        CALL bremsstrahlung_update_optical_depth()
      END IF
#endif

      ! Evaluate fields a half timestep ahead of the particles
      IF (use_hybrid_fields) THEN
        CALL half_B_step
        CALL calculate_E
      END IF

      ! Logical flag set to true when particles can start to move
      push = (time >= particle_push_start_time)

      ! The following scripts will only be executed if particles can move
      IF (push) THEN

        ! Inject particles into the simulation
        CALL run_injectors

        ! .FALSE. this time to use load balancing threshold
        IF (use_balance) CALL balance_workload(.FALSE.)

        ! Move particles, leapfrogging over E and B
        CALL push_particles

        ! Obtain heat capacity to calculate the temperature change in
        ! hybrid_collisions and ohmic_heating
        ALLOCATE(heat_capacity(1-ng:nx+ng, 1-ng:ny+ng))
        CALL get_heat_capacity(heat_capacity)

        ! Calculates collisional drag/scattering, and updates grid temperature
        IF (use_hybrid_collisions) THEN
          CALL hybrid_collisions(heat_capacity)
        END IF

        ! Updates grid temperature due to Ohmic heating
        CALL ohmic_heating(heat_capacity)
        DEALLOCATE(heat_capacity)

        ! Now that temperature has been fully updated, re-evaluate resistivity
        CALL update_resistivity

        ! Migrate particle species if they pass the migration criteria
        IF (use_particle_migration) CALL migrate_particles(step)

        ! For current smoothing, see housekeeping/current_smooth.F90
        CALL current_finish

        ! See housekeeping/partlist.F90
        CALL update_particle_count
      END IF

      CALL check_for_stop_condition(halt, force_dump)
      IF (halt) EXIT
      step = step + 1
      time = time + dt / 2.0_num
      CALL output_routines(step)
      time = time + dt / 2.0_num

      ! Iterate B and E, such that they are evaluated at the same time as the
      ! particles (note: main PIC loop also does this after the output dump)
      IF (use_hybrid_fields) THEN
        CALL half_B_step
        CALL calculate_E
      END IF
    END DO

  END SUBROUTINE run_hybrid_PIC



  SUBROUTINE initialise_hybrid

    ! This subroutine initialises the hybrid arrays, and sets the values of some
    ! constants, to speed up the code

    REAL(num) :: resistivity_init
    INTEGER :: ix, iy

    ! Preset useful constants
    hybrid_const_dx = 1.0_num / (mu0 * dx)
    hybrid_const_dy = 1.0_num / (mu0 * dy)
    hybrid_const_K_to_eV = kb / q0
    hybrid_const_dt_by_dx = 0.5_num * dt / dx
    hybrid_const_dt_by_dy = 0.5_num * dt / dy
    hybrid_D = 0.5_num * hybrid_ne * q0**4 / pi / (epsilon0**2)
    hybrid_ln_s = 4.0_num * epsilon0 * h_planck &
        / (hybrid_Z**(1.0_num/3.0_num) * m0 * q0**2)
    hybrid_const_ZeV = hybrid_Z**(-4.0_num/3.0_num) * hybrid_const_K_to_eV

    ! Allocate additional arrays for running in hybrid mode. These require extra
    ! remapping scripts in balance.F90 (for domain change in load-balance)
    ALLOCATE(resistivity(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(hybrid_Tb(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(ion_charge(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(ion_density(1-ng:nx+ng,1-ng:ny+ng))
    ALLOCATE(ion_temp(1-ng:nx+ng,1-ng:ny+ng))

    resistivity_init = calc_resistivity(hybrid_Tb_init)

    ! Initialise arrays
    DO iy = 1, ny
      DO ix = 1, nx
        ! Set background temperature and resistivity to those matching
        ! hybrid_Tb_init, which is read in from input.deck
        hybrid_Tb(ix,iy) = hybrid_Tb_init
        resistivity(ix,iy) = resistivity_init

        ! These arrays currently serve no purpose, but will be used when
        ! ionisation routines are added
        ion_charge(ix,iy) = 0.0_num
        ion_density(ix,iy) = 0.0_num
        ion_temp(ix,iy) = 0.0_num
      END DO
    END DO

  END SUBROUTINE initialise_hybrid



  FUNCTION calc_resistivity(Tb)

    ! Calculates the resistivity for the background temperature Tb, using the
    ! equation from H. M. Milchberg, et al, 1988. Phys. Rev. Lett., 61(20),
    ! p.2364.

    REAL(num), INTENT(IN) :: Tb
    REAL(num) :: Tb_eV
    REAL(num) :: calc_resistivity

    ! Convert Tb from Kelvin to eV
    Tb_eV = hybrid_const_K_to_eV * Tb

    ! Calculate resistivity for a METAL
    calc_resistivity = Tb_eV &
        / (5.0e6_num + 170.0_num*Tb_eV**2.5_num + 3.0e5_num*Tb_eV)

  END FUNCTION calc_resistivity



  SUBROUTINE half_B_step

    ! This subroutine performs a half-step in the magnetic field, assuming a
    ! constant electric field and using:
    !
    ! dB/dt = -curl(E)

    ! Calculate the initial electric field
    INTEGER :: ix, iy

    ! Update B by half a timestep
    DO iy = 1, ny
      DO ix = 1, nx
        bx(ix, iy) = bx(ix, iy) &
            - hybrid_const_dt_by_dy * (ez(ix  , iy+1) - ez(ix  , iy  ))

        by(ix, iy) = by(ix, iy) &
            + hybrid_const_dt_by_dx * (ez(ix+1, iy  ) - ez(ix  , iy  ))

        bz(ix, iy) = bz(ix, iy) &
            - hybrid_const_dt_by_dx * (ey(ix+1, iy  ) - ey(ix  , iy  )) &
            + hybrid_const_dt_by_dy * (ex(ix  , iy+1) - ex(ix  , iy  ))
      END DO
    END DO

  END SUBROUTINE half_B_step



  SUBROUTINE calculate_E

    ! Calculates the electric field for the current values of global variables
    ! bx, by, bz, j, and resistivity, using the equation:
    !
    ! E = resistivity * (curl(B)/mu_0 - J)
    !
    ! We have precalculated hybrid_const_dx = 1/(mu_0*dx)

    ! Note that B is staggered from E, whereas J is evaulated at the same point
    ! as E. Resistivity is a cell centred variable.

    INTEGER :: ix, iy

    DO iy = 1, ny
      DO ix = 1, nx
        ex(ix,iy) = &
            + 0.5_num * (resistivity(ix+1,iy) + resistivity(ix,iy)) &
            * (+ hybrid_const_dy * (bz(ix,iy) - bz(ix,iy-1)) &
               - jx(ix,iy))

        ey(ix,iy) = &
            + 0.5_num * (resistivity(ix,iy+1) + resistivity(ix,iy)) &
            * (- hybrid_const_dx * (bz(ix,iy) - bz(ix-1,iy)) &
               - jy(ix,iy))

        ez(ix,iy) = &
            + resistivity(ix,iy) &
            * (+ hybrid_const_dx * (by(ix,iy) - by(ix-1,iy)) &
               - hybrid_const_dy * (bx(ix,iy) - bx(ix,iy-1)) &
               - jz(ix,iy))
      END DO
    END DO

  END SUBROUTINE calculate_E



  SUBROUTINE get_heat_capacity(heat_capacity)

    ! Calculates heat capacity as described by Davies, et al, (2002). Phys. Rev.
    ! E, 65(2), 026407
    !
    ! C = 0.3 + 1.2*T'*(2.2 + T')/(1.1 + T')^2
    !
    ! Where T' = Z^(-4/3) * Tb, and Tb is measured in eV

    REAL(num), INTENT(OUT) :: heat_capacity(1-ng:nx+ng, 1-ng:ny+ng)
    REAL(num) :: T_prime, hybrid_const_ZeV
    INTEGER :: ix, iy

    DO iy = 1, ny
      DO ix = 1, nx
        T_prime = hybrid_const_ZeV * hybrid_Tb(ix,iy)
        heat_capacity(ix,iy) = 0.3_num &
            + 1.2_num * T_prime * (2.2_num + T_prime)/(1.1_num + T_prime)**2
      END DO
    END DO

  END SUBROUTINE get_heat_capacity



  SUBROUTINE hybrid_collisions(heat_capacity)

    ! Calculates collisional drag/scattering, and updates grid temperature based
    ! on energy lost by the particles

    REAL(num), INTENT(IN) :: heat_capacity(1-ng:nx+ng,1-ng:ny+ng)
    INTEGER :: ispecies
    INTEGER(i8) :: ipart
    REAL(num) :: ipart_mc, m2c4
    REAL(num) :: px, py, pz, p, gamma, v, ipart_KE, ln_lambda_L, delta_p
    REAL(num) :: rand1, rand2, rand_scatter, ln_lambda_S, delta_theta, delta_phi
    REAL(num) :: frac_p, ux, uy, uz, frac_uz, frac
    REAL(num) :: sin_t, cos_t, cos_p, sin_t_cos_p, sin_t_sin_p
    REAL(num) :: p_new_2, weight, dE, part_x, part_y, part_C, delta_Tb
    TYPE(particle), POINTER :: current, next

    ! Generate a normally distributed random number for the scatter angle.
    ! The argument here refers to the standard deviation. The mean is
    ! assumed zero. In Davies, et al, (2002). Phys. Rev. E, 65(2), 026407, this
    ! random number is a function of time, which might mean constant for all
    ! particles at a given timestep, t.
    rand_scatter = random_box_muller(1.0_num)

    ! Cycle over all non-photon species
    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      IF (species_list(ispecies)%species_type == c_species_id_photon) CYCLE

      ipart_mc  = 1.0_num / (c * species_list(ispecies)%mass)
      m2c4 = species_list(ispecies)%mass**2 * c**4

      ! Cycle over all particles in the species
      DO ipart = 1, species_list(ispecies)%attached_list%count
        next => current%next

        ! Find collisional drag momentum loss
        px = current%part_p(1)
        py = current%part_p(2)
        pz = current%part_p(3)
        p = SQRT(px**2 + py**2 + pz**2)
        gamma = SQRT((p * ipart_mc)**2 + 1.0_num)
        v = p * ipart_mc * c / gamma
        ipart_KE = (gamma - 1.0_num) * species_list(ispecies)%mass * c**2
        ln_lambda_L = LOG(ipart_KE/hybrid_Iex) + 0.5_num*LOG(gamma + 1.0_num) &
            + 0.909_num/gamma**2 - 0.818_num/gamma - 0.246_num
        delta_p = - 0.5_num * hybrid_D / (m0 * v**2) * ln_lambda_L * dt

        ! Calculate scattering angles
        ln_lambda_S = LOG(hybrid_ln_S * p)
        delta_theta = SQRT(hybrid_Z * hybrid_D * ln_lambda_S  * dt / v) &
            * rand_scatter / p
        delta_phi = 2.0_num * pi * random()

        ! Apply scattering angles delta_theta and delta_phi to rotate p
        frac_p = 1.0_num / p
        ux = px * frac_p
        uy = py * frac_p
        uz = pz * frac_p
        sin_t = SIN(delta_theta)
        IF (ABS(1.0_num - uz) < 1.0e-5_num) THEN
          px = p * sin_t * COS(delta_phi)
          py = p * sin_t * SIN(delta_phi)
          pz = p * uz / ABS(uz) * COS(delta_theta)
        ELSE
          frac_uz = 1.0_num/SQRT(1 - uz**2)
          cos_t = COS(delta_theta)
          cos_p = COS(delta_phi)
          sin_t_cos_p = sin_t*cos_p
          sin_t_sin_p = sin_t*SIN(delta_phi)
          px = p*(ux*cos_t + sin_t_cos_p*ux*uz*frac_uz - sin_t_sin_p*uy*frac_uz)
          py = p*(uy*cos_t + sin_t_cos_p*uy*uz*frac_uz + sin_t_sin_p*ux*frac_uz)
          pz = p*(uz*cos_t + cos_p*sin_t*(uz**2 - 1.0_num)*frac_uz)
        END IF

        ! Reduce particle momentum
        frac = MAX(1.0_num + delta_p / ABS(p), 0.0_num)
        current%part_p(1) = px * frac
        current%part_p(2) = py * frac
        current%part_p(3) = pz * frac

        ! Calculate energy change due to collisions
        p_new_2 = current%part_p(1)**2 + current%part_p(2)**2 &
            + current%part_p(3)**2
#ifndef PER_SPECIES_WEIGHT
        weight = current%weight
#else
        IF(rank == 0) THEN
          PRINT*,'*** ERROR ***'
          PRINT*,'Currently, the code uses particle weight to calculate ' &
              'energy loss in the hybrid collisions. This cannot be accessed ' &
              'if the precompiler flag -DPER_SPECIES_WEIGHT is switched on.'
          PRINT*,'Code will terminate'
          CALL abort_code(c_err_bad_value)
        END IF
#endif
        dE = (SQRT(p_new_2*c**2 + m2c4) - SQRT((p*c)**2 + m2c4)) * weight

        ! Get heat capacity at the particle position
        part_x = current%part_pos(1) - x_grid_min_local
        part_y = current%part_pos(2) - y_grid_min_local
        CALL hy_grid_centred_var_at_particle(part_x, part_y, part_C, &
            heat_capacity)

        ! Calculate the temperature increase, and add this to Tb. We use -dE, as
        ! the energy gain for temperature is equal to the energy loss of the
        ! electron. We use T[Kelvin] = T[Joule] / kb. Final term is energy per
        ! unit volume, and dz = 1m in epoch2d
        delta_Tb = 1.0_num/(hybrid_ne * part_C * kb)*(-dE / (dx*dy))

        ! Write temperature change to the grid
        CALL add_particle_var_to_grid(part_x, part_y, delta_Tb, hybrid_Tb)

        current => next

      END DO
    END DO

  END SUBROUTINE hybrid_collisions



  SUBROUTINE ohmic_heating(heat_capacity)

    ! Calculates the Ohmic heating of the simulation grid, as described by
    ! Davies, et al, (2002). Phys. Rev. E, 65(2), 026407. At this point in the
    ! simulation, we are at the end of a timestep, with E evaluated in the
    ! middle of this timestep. Assuming this is the average electric field over
    ! the timestep, we have a power per unit volume of E**2/resistivity

    REAL(num), INTENT(IN) :: heat_capacity(1-ng:nx+ng, 1-ng:ny+ng)
    REAL(num) :: E2, fac
    INTEGER :: ix, iy

    fac = dt/(hybrid_ne*kb)

    DO iy = 1, ny
      DO ix = 1, nx
        ! Tb is a cell-centred variable, but E has stagger - need to average
        E2 = (0.5_num*(ex(ix,iy) + ex(ix-1,iy)))**2 &
            +(0.5_num*(ey(ix,iy) + ey(ix,iy-1)))**2 &
            + ez(ix,iy)**2

        ! Calculate Ohmic heating, avoiding the 0/0 NaN
        IF (resistivity(ix,iy) > 0.0_num) THEN
          hybrid_Tb(ix,iy) = hybrid_Tb(ix,iy) &
              + E2*fac/(resistivity(ix,iy) * heat_capacity(ix,iy))
        END IF
      END DO
    END DO

  END SUBROUTINE ohmic_heating



  SUBROUTINE update_resistivity

    ! Loops over the local resistivity grid and updates each point

    INTEGER :: ix, iy

    DO iy = 1, ny
      DO ix = 1, nx
         resistivity(ix,iy) = calc_resistivity(hybrid_Tb(ix,iy))
      END DO
    END DO

  END SUBROUTINE update_resistivity



  SUBROUTINE hy_grid_centred_var_at_particle(part_x, part_y, part_var, grid_var)

    ! Calculates the value of a grid-centred variable "part_var" stored in the
    ! grid "grid_var", averaged over the particle shape for a particle at
    ! position (part_x, part_y)

    REAL(num), INTENT(in) :: part_x, part_y
    REAL(num), INTENT(in) :: grid_var(1-ng:nx+ng, 1-ng:ny+ng)
    REAL(num), INTENT(out) :: part_var
    INTEGER :: cell_x1, cell_y1
    REAL(num) :: cell_x_r, cell_y_r
    REAL(num) :: cell_frac_x, cell_frac_y
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = (1.0_num)**c_ndims
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    ! The following method is lifted from photons.F90 (field_at_particle), for
    ! the cell-centered fields, taking into account the various particle shapes
#ifdef PARTICLE_SHAPE_TOPHAT
    cell_x_r = part_x / dx - 0.5_num
    cell_y_r = part_y / dy - 0.5_num
#else
    cell_x_r = part_x / dx
    cell_y_r = part_y / dy
#endif
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1
    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"

    part_var = &
          gy(-2) * (gx(-2) * grid_var(cell_x1-2,cell_y1-2) &
        +           gx(-1) * grid_var(cell_x1-1,cell_y1-2) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1-2) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1-2) &
        +           gx( 2) * grid_var(cell_x1+2,cell_y1-2)) &
        + gy(-1) * (gx(-2) * grid_var(cell_x1-2,cell_y1-1) &
        +           gx(-1) * grid_var(cell_x1-1,cell_y1-1) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1-1) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1-1) &
        +           gx( 2) * grid_var(cell_x1+2,cell_y1-1)) &
        + gy( 0) * (gx(-2) * grid_var(cell_x1-2,cell_y1  ) &
        +           gx(-1) * grid_var(cell_x1-1,cell_y1  ) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1  ) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1  ) &
        +           gx( 2) * grid_var(cell_x1+2,cell_y1  )) &
        + gy( 1) * (gx(-2) * grid_var(cell_x1-2,cell_y1+1) &
        +           gx(-1) * grid_var(cell_x1-1,cell_y1+1) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1+1) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1+1) &
        +           gx( 2) * grid_var(cell_x1+2,cell_y1+1)) &
        + gy( 2) * (gx(-2) * grid_var(cell_x1-2,cell_y1+2) &
        +           gx(-1) * grid_var(cell_x1-1,cell_y1+2) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1+2) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1+2) &
        +           gx( 2) * grid_var(cell_x1+2,cell_y1+2))
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
    part_var = &
          gy( 0) * (gx( 0) * grid_var(cell_x1  ,cell_y1  ) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1  )) &
        + gy( 1) * (gx( 0) * grid_var(cell_x1  ,cell_y1+1) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1+1))
#else
#include "triangle/gx.inc"
    part_var = &
          gy(-1) * (gx(-1) * grid_var(cell_x1-1,cell_y1-1) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1-1) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1-1)) &
        + gy( 0) * (gx(-1) * grid_var(cell_x1-1,cell_y1  ) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1  ) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1  )) &
        + gy( 1) * (gx(-1) * grid_var(cell_x1-1,cell_y1+1) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1+1) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1+1))
#endif
    part_var = fac*part_var

  END SUBROUTINE hy_grid_centred_var_at_particle



  SUBROUTINE add_particle_var_to_grid(part_x, part_y, part_var, grid)

    ! Adds the value of a variable (part_var), which is evaluated on a particle
    ! at position (part_x, part_y), to the grid (grid). Particle weighting
    ! assumes we are summing to a cell-centred variable (no stagger). This uses
    ! the same cell indices and weighting as "hy_grid_centred_var_at_particle".

    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(INOUT) :: grid(1-ng:nx+ng, 1-ng:ny+ng)
    REAL(num), INTENT(IN) :: part_var
    INTEGER :: cell_x1, cell_y1
    REAL(num) :: cell_x_r, cell_y_r
    REAL(num) :: cell_frac_x, cell_frac_y
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
    REAL(num) :: scaled_part_var
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = (1.0_num)**c_ndims
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

#ifdef PARTICLE_SHAPE_TOPHAT
    cell_x_r = part_x / dx - 0.5_num
    cell_y_r = part_y / dy - 0.5_num
#else
    cell_x_r = part_x / dx
    cell_y_r = part_y / dy
#endif
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1
    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1

    scaled_part_var = fac*part_var

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#include "bspline3/part_to_grid.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#include "tophat/part_to_grid.inc"
#else
#include "triangle/gx.inc"
#include "triangle/part_to_grid.inc"
#endif

  END SUBROUTINE add_particle_var_to_grid

#endif
END MODULE hybrid
