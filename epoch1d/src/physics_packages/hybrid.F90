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
  USE boundary
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
        CALL hybrid_bremsstrahlung_update_optical_depth()
      END IF
#endif

      ! Evaluate fields a half timestep ahead of the particles
      IF (use_hybrid_fields) THEN
        CALL half_B_step
        CALL hybrid_B_bc
        CALL calculate_E
        CALL hybrid_E_bc
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
        CALL get_heat_capacity

        ! Calculates collisional drag/scattering, and updates grid temperature
        IF (use_hybrid_collisions) THEN
          CALL hybrid_collisions
        END IF

        ! Updates grid temperature due to Ohmic heating
        CALL ohmic_heating
        CALL clear_heat_capacity

        ! Now that temperature has been fully updated, re-evaluate resistivity
        CALL update_resistivity

        ! Migrate particle species if they pass the migration criteria
        IF (use_particle_migration) CALL migrate_particles(step)

        ! Pass current to neighbouring ranks (housekeeping/current_smooth.F90)
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
        CALL hybrid_B_bc
        CALL calculate_E
        CALL hybrid_E_bc
      END IF
    END DO

  END SUBROUTINE run_hybrid_PIC



  SUBROUTINE initialise_hybrid

    ! This subroutine initialises the hybrid arrays, and sets the values of some
    ! constants, to speed up the code

    REAL(num) :: resistivity_init, max_ne
    REAL(num) :: sum_ne(1-ng:nx+ng)
    INTEGER :: ix, i_sol
    INTEGER :: io, iu

    IF (use_hybrid) THEN
      ! Preset useful constants
      hybrid_const_dx = 1.0_num / (mu0 * dx)
      hybrid_const_K_to_eV = kb / q0
      hybrid_const_dt_by_dx = 0.5_num * dt / dx

      ! Preset useful constants for solids
      ALLOCATE(hybrid_const_heat(1-ng:nx+ng))
      sum_ne = 0
      DO i_sol = 1, solid_count
        solid_array(i_sol)%el_density = solid_array(i_sol)%hybrid_Z &
            * solid_array(i_sol)%ion_density
        solid_array(i_sol)%theta_fac = solid_array(i_sol)%hybrid_Z * q0**4 &
            / (2.0_num * pi * epsilon0**2 )
        solid_array(i_sol)%hybrid_ln_s = 4.0_num * epsilon0 * h_planck &
            / (solid_array(i_sol)%hybrid_Z**(1.0_num/3.0_num) * m0 * q0**2)
        solid_array(i_sol)%hybrid_const_ZeV = &
            solid_array(i_sol)%hybrid_Z**(-4.0_num/3.0_num) &
            * hybrid_const_K_to_eV

        sum_ne = sum_ne + solid_array(i_sol)%el_density
      END DO
      hybrid_const_heat = 1.0_num/(kb * sum_ne**2)

      ! Allocate additional arrays for running in hybrid mode. These require
      ! extra remapping scripts in balance.F90 (for domain change in
      ! load-balance)
      ALLOCATE(resistivity(1-ng:nx+ng))
      ALLOCATE(resistivity_model(1-ng:nx+ng))
      ALLOCATE(hybrid_Tb(1-ng:nx+ng))
      ALLOCATE(ion_charge(1-ng:nx+ng))
      ALLOCATE(ion_density(1-ng:nx+ng))
      ALLOCATE(ion_temp(1-ng:nx+ng))

      ! Assign a resistivity model to each cell based on present solids
      ! (default to vacuum)
      resistivity_model = c_resist_vacuum
      DO ix = 1-ng, nx+ng
        max_ne = 0.0_num
        DO i_sol = 1, solid_count
          ! Use the same resistivity model as the solid with the highest local
          ! electron number density
          IF (solid_array(i_sol)%el_density(ix) > max_ne) THEN
            max_ne = solid_array(i_sol)%el_density(ix)
            resistivity_model(ix) = solid_array(i_sol)%material
          END IF
        END DO
      END DO

      ! Initialise arrays
      DO ix = 1-ng, nx+ng
        ! Set background temperature and resistivity to those matching
        ! hybrid_Tb_init, which is read in from input.deck
        hybrid_Tb(ix) = hybrid_Tb_init

        ! These arrays currently serve no purpose, but will be used when
        ! ionisation routines are added
        ion_charge(ix) = 0.0_num
        ion_density(ix) = 0.0_num
        ion_temp(ix) = 0.0_num
      END DO
      CALL update_resistivity

      ! Do we have a solid background?
      IF (solid_count == 0) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No solids specified! Please set all solid ion ', &
                'number densities to positive values.'
            WRITE(io,*) 'Code will terminate.'
          END DO
        END IF
        errcode = c_err_bad_value + c_err_terminate
      END IF

      ! Check the solid parameters have been entered correctly
      DO i_sol = 1, solid_count

        ! Check ion number density is positive
        IF (MINVAL(solid_array(i_sol)%ion_density) < 0.0_num) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Please ensure the solid ion number densities are ', &
                  'positive values in all cells.'
              WRITE(io,*) 'Code will terminate.'
            END DO
          END IF
          errcode = c_err_bad_value + c_err_terminate
        END IF

        ! Check the atomic number is positive
        IF (solid_array(i_sol)%hybrid_Z < 0.0_num) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Please set all solid atomic numbers to positive ', &
                  'values.'
              WRITE(io,*) 'Code will terminate.'
            END DO
          END IF
          errcode = c_err_bad_value + c_err_terminate
        END IF

        ! Check the mean excitation energy is positive (only used in collisions)
        IF (solid_array(i_sol)%hybrid_Iex < 0.0_num &
            .AND. use_hybrid_collisions) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Please set all solid mean excitation energies to ', &
                  'positive values for hybrid collisions'
              WRITE(io,*) 'Code will terminate.'
            END DO
          END IF
          errcode = c_err_bad_value + c_err_terminate
        END IF
      END DO

      IF (rank == 0) PRINT*, 'Code is running in hybrid mode'

    ELSE
      ! Do not try to output hybrid variables if we aren't running in hybrid
      ! mode
      IF (rank == 0) THEN
        PRINT*, ''
        PRINT*, 'Code is not running in hybrid mode'
        PRINT*, 'Any attempts to output hybrid-only variables (like ', &
           'resistivity) will be ignored'
        PRINT*, 'Switch off the -DHYBRID pre-processor flag to stop this ', &
            'message from printing'
        PRINT*, ''
      END IF
    END IF

  END SUBROUTINE initialise_hybrid



  FUNCTION calc_resistivity_conductor(Tb)

    ! Calculates the resistivity for the background temperature Tb, using the
    ! equation from H. M. Milchberg, et al, 1988. Phys. Rev. Lett., 61(20),
    ! p.2364.

    REAL(num), INTENT(IN) :: Tb
    REAL(num) :: Tb_eV
    REAL(num) :: calc_resistivity_conductor

    ! Convert Tb from Kelvin to eV
    Tb_eV = hybrid_const_K_to_eV * Tb

    ! Calculate resistivity
    calc_resistivity_conductor = Tb_eV &
        / (5.0e6_num + 170.0_num*Tb_eV**2.5_num + 3.0e5_num*Tb_eV)

  END FUNCTION calc_resistivity_conductor



  FUNCTION calc_resistivity_insulator(Tb)

    ! Calculates the resistivity for the background temperature Tb, using the
    ! heuristic model for plastic from Davies, et al, (2002). Phys. Rev. E,
    ! 65(2), 026407

    REAL(num), INTENT(IN) :: Tb
    REAL(num) :: Tb_eV
    REAL(num) :: calc_resistivity_insulator

    ! Convert Tb from Kelvin to eV
    Tb_eV = hybrid_const_K_to_eV * Tb

    ! Calculate resistivity
    calc_resistivity_insulator =  1.0_num / (4.3e5_num + 1.3e3_num*Tb**1.5_num)

  END FUNCTION calc_resistivity_insulator



  SUBROUTINE half_B_step

    ! This subroutine performs a half-step in the magnetic field, assuming a
    ! constant electric field and using:
    !
    ! dB/dt = -curl(E)
    !
    ! We calculate into 1 ghost cell to allow non-zero curls across simulation
    ! boundaries

    INTEGER :: ix

    ! Update B by half a timestep
    DO ix = 0, nx+1
      by(ix) = by(ix) &
          + hybrid_const_dt_by_dx * (ez(ix+1) - ez(ix))

      bz(ix) = bz(ix) &
          - hybrid_const_dt_by_dx * (ey(ix+1) - ey(ix))
    END DO

  END SUBROUTINE half_B_step



  SUBROUTINE hybrid_B_bc

    ! This script updates the ghost cells on the local processor, by either
    ! passing over those from the neighbouring processor, or by ensuring the
    ! curls beyond the boundary are zero

    INTEGER :: i

    ! Pass neighbouring ghost cells
    CALL field_bc(bx, ng)
    CALL field_bc(by, ng)
    CALL field_bc(bz, ng)

    ! Special cases for boundary processors
    DO i = 1, 2*c_ndims
      CALL field_zero_curl(bx, i)
      CALL field_zero_curl(by, i)
      CALL field_zero_curl(bz, i)
    END DO

  END SUBROUTINE hybrid_B_bc



  SUBROUTINE calculate_E

    ! Calculates the electric field for the current values of global variables
    ! bx, by, bz, j, and resistivity, using the equation:
    !
    ! E = resistivity * (curl(B)/mu_0 - J)
    !
    ! We have precalculated hybrid_const_dx = 1/(mu_0*dx)

    ! Note that B is staggered from E, whereas J is evaulated at the same point
    ! as E. Resistivity is a cell centred variable. We calculate into 1 ghost
    ! cell to allow non-zero curls across simulation boundaries

    INTEGER :: ix

    DO ix = 0, nx+1
      ex(ix) = &
          - 0.5_num * (resistivity(ix+1) + resistivity(ix)) &
          * jx(ix)

      ey(ix) = &
          + resistivity(ix) &
          * (- hybrid_const_dx * (bz(ix) - bz(ix-1)) &
          - jy(ix))

      ez(ix) = &
          + resistivity(ix) &
          * (+ hybrid_const_dx * (by(ix) - by(ix-1)) &
          - jz(ix))
    END DO

  END SUBROUTINE calculate_E



  SUBROUTINE hybrid_E_bc

    ! This script updates the ghost cells on the local processor, by either
    ! passing over those from the neighbouring processor, or by ensuring the
    ! curls beyond the boundary are zero

    INTEGER :: i

    ! Pass neighbouring ghost cells
    CALL field_bc(ex, ng)
    CALL field_bc(ey, ng)
    CALL field_bc(ez, ng)

    ! Special cases for boundary processors
    DO i = 1, 2*c_ndims
      CALL field_zero_curl(ex, i)
      CALL field_zero_curl(ey, i)
      CALL field_zero_curl(ez, i)
    END DO

  END SUBROUTINE hybrid_E_bc



  SUBROUTINE get_heat_capacity

    ! Calculates heat capacity as described by Davies, et al, (2002). Phys. Rev.
    ! E, 65(2), 026407, for each element of solid_array
    !
    ! C = 0.3 + 1.2*T'*(2.2 + T')/(1.1 + T')^2
    !
    ! Where T' = Z^(-4/3) * Tb, and Tb is measured in eV

    REAL(num) :: T_prime
    INTEGER :: ix, i_sol

    DO i_sol = 1, solid_count
      ALLOCATE(solid_array(i_sol)%heat_capacity(1-ng:nx+ng))

      DO ix = 1-ng, nx+ng
        T_prime = solid_array(i_sol)%hybrid_const_ZeV * hybrid_Tb(ix)
        solid_array(i_sol)%heat_capacity(ix) = 0.3_num &
            + 1.2_num * T_prime * (2.2_num + T_prime)/(1.1_num + T_prime)**2
      END DO
    END DO

  END SUBROUTINE get_heat_capacity



  SUBROUTINE hybrid_collisions

    ! Calculates collisional drag/scattering, and updates grid temperature based
    ! on energy lost by the particles

    INTEGER :: ispecies
    INTEGER(i8) :: ipart
    REAL(num) :: heat_const(1-ng:nx+ng)
    REAL(num) :: part_heat_const
    REAL(num) :: ipart_mc, m2c4, sum_dp, sum_dtheta, part_ne
    REAL(num) :: px, py, pz, p, gamma, v, ipart_KE, ln_lambda_L_terms, delta_p
    REAL(num) :: rand1, rand2, rand_scatter, ln_lambda_S, delta_theta, delta_phi
    REAL(num) :: frac_p, ux, uy, uz, frac_uz, frac
    REAL(num) :: sin_t, cos_t, cos_p, sin_t_cos_p, sin_t_sin_p
    REAL(num) :: p_new_2, weight, dE, part_x, part_C, delta_Tb
    REAL(num) :: idx
    INTEGER :: ix, i_sol
    TYPE(particle), POINTER :: current, next

    idx = 1.0_num/dx

    ! Generate a normally distributed random number for the scatter angle.
    ! The argument here refers to the standard deviation. The mean is
    ! assumed zero. In Davies, et al, (2002). Phys. Rev. E, 65(2), 026407, this
    ! random number is a function of time, which might mean constant for all
    ! particles at a given timestep, t.
    rand_scatter = random_box_muller(1.0_num)

    ! 1 / (kb * sum(ne) * dx) - last bit is per unit volume, but dy and dz are
    ! 1m in epoch1d
    heat_const = hybrid_const_heat / dx

    ! Loop over all non-photon species
    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      IF (species_list(ispecies)%species_type == c_species_id_photon) CYCLE

      ipart_mc  = 1.0_num / (c * species_list(ispecies)%mass)
      m2c4 = species_list(ispecies)%mass**2 * c**4

      ! Loop over all particles in the species
      DO ipart = 1, species_list(ispecies)%attached_list%count
        next => current%next

        ! Extract particle variables
        part_x = current%part_pos - x_grid_min_local
        px = current%part_p(1)
        py = current%part_p(2)
        pz = current%part_p(3)
        p = SQRT(px**2 + py**2 + pz**2)
        gamma = SQRT((p * ipart_mc)**2 + 1.0_num)
        v = p * ipart_mc * c / gamma
        ipart_KE = (gamma - 1.0_num) * species_list(ispecies)%mass * c**2
        ln_lambda_L_terms = 0.5_num*LOG(gamma + 1.0_num) + 0.909_num/gamma**2 &
            - 0.818_num/gamma - 0.246_num

        ! Extract solid variables averaged over all solids present in this cell
        sum_dp = 0.0_num
        sum_dtheta = 0.0_num
        DO i_sol = 1, solid_count
          CALL hy_grid_centred_var_at_particle(part_x, part_ne, &
              solid_array(i_sol)%el_density)

          ! If ne is zero, then the solid is not in this cell, so don't
          ! calculate the constants
          IF (part_ne < TINY(1.0_num)) CYCLE

          ! Very low energy particles (~50 eV) will make ln_lambda_S go negative
          ! causing a NaN in delta_theta. The hybrid mode isn't designed for low
          ! energy electrons, so delta_theta will be ignored in this case
          sum_dp = sum_dp + part_ne * ( ln_lambda_L_terms &
              + LOG(MAX(1.0_num, ipart_KE/solid_array(i_sol)%hybrid_Iex)))
          sum_dtheta = sum_dtheta + solid_array(i_sol)%theta_fac * part_ne &
              * LOG(MAX(1.0_num, solid_array(i_sol)%hybrid_ln_s * p)) * dt / v
        END DO

        ! Collisional changes to momentum and direction
        delta_p =  hybrid_const_dp * sum_dp * dt &
            / (species_list(ispecies)%mass * v**2)
        delta_theta = SQRT(sum_dtheta) / p * rand_scatter
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
          px = p*(ux*cos_t + sin_t_cos_p*ux*uz*frac_uz &
              - sin_t_sin_p*uy*frac_uz)
          py = p*(uy*cos_t + sin_t_cos_p*uy*uz*frac_uz &
              + sin_t_sin_p*ux*frac_uz)
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
        weight = species_list(ispecies)%weight
#endif
        dE = (SQRT(p_new_2*c**2 + m2c4) - SQRT((p*c)**2 + m2c4)) * weight

        ! Get effective heat capacity at the particle position
        CALL hy_grid_centred_var_at_particle(part_x, part_heat_const, &
            heat_const)
        CALL get_effective_heat_capacity(part_x, part_C)

        ! Calculate the temperature increase, and add this to Tb. We use -dE,
        ! as the energy gain for temperature is equal to the energy loss of
        ! the electron.
        delta_Tb = - dE * part_C * part_heat_const

        ! Write temperature change to the grid (ignores particle shape)
        ! Impose maximum temperature of 1e30 to prevent temperature rising
        ! unphysically
        ix = MAX(1,CEILING(part_x*idx))
        hybrid_Tb(ix) = MIN(hybrid_Tb(ix) + delta_Tb, 1.0e30_num)

        current => next

      END DO
    END DO

    ! Pass new temperature values to ghost cells of neighbouring processors
    CALL field_bc(hybrid_Tb, ng)

  END SUBROUTINE hybrid_collisions



  SUBROUTINE get_effective_heat_capacity(part_x, part_C)

    ! Temperature increase for one background species would be:
    !
    ! (energy change per unit vol.) * (vol. of e-, 1/ne) / (heat cap. of e-)
    !
    ! This implementation of the Davies model can use multiple backgrounds, so
    ! we define an effective heat capacity to find the total temperature rise.
    ! This assumes the energy is split evenly between all e- in the cell, which
    ! convert the energy into a temperature increase using different C values
    !
    ! (effective heat capacity) = sum(ne/C)

    REAL(num), INTENT(IN) :: part_x
    REAL(num), INTENT(OUT) :: part_C
    INTEGER :: i_sol
    REAL(num) :: C_sol, part_ne

    part_C = 0.0_num

    ! Loop over all solids
    DO i_sol = 1, solid_count
      ! Get heat capacity at particle position for current solid
      CALL hy_grid_centred_var_at_particle(part_x, C_sol, &
          solid_array(i_sol)%heat_capacity)
      CALL hy_grid_centred_var_at_particle(part_x, part_ne, &
          solid_array(i_sol)%el_density)

      ! Contribution to the effective heat capacity
      part_C = part_C + part_ne / C_sol
    END DO

  END SUBROUTINE get_effective_heat_capacity



  SUBROUTINE ohmic_heating

    ! Calculates the Ohmic heating of the simulation grid, as described by
    ! Davies, et al, (2002). Phys. Rev. E, 65(2), 026407. At this point in the
    ! simulation, we are at the end of a timestep, with E evaluated in the
    ! middle of this timestep. Assuming this is the average electric field over
    ! the timestep, we have a power per unit volume of E**2/resistivity

    REAL(num) :: E2
    REAL(num) :: eff_heat_capacity(1-ng:nx+ng)
    REAL(num) :: fac(1-ng:nx+ng)
    INTEGER :: ix, i_sol

    fac = dt * hybrid_const_heat

    ! Obtain the effective heat capacity: sum(ne/C)
    eff_heat_capacity = 0.0_num
    DO i_sol = 1, solid_count
      DO ix = 1, nx
        eff_heat_capacity(ix) = eff_heat_capacity(ix) &
            + solid_array(i_sol)%el_density(ix) &
            / solid_array(i_sol)%heat_capacity(ix)
      END DO
    END DO

    ! Loop over all grid points to find temperature change
    DO ix = 1, nx
      ! Tb is a cell-centred variable, but E has stagger - need to average
      E2 = (0.5_num*(ex(ix) + ex(ix-1)))**2 + ey(ix)**2 + ez(ix)**2

      ! Calculate Ohmic heating, avoiding the 0/0 NaN. Cap temperature at 1.0e30
      IF (resistivity(ix) > 0.0_num .AND. &
          eff_heat_capacity(ix) > 0.0_num) THEN
        hybrid_Tb(ix) = MIN(hybrid_Tb(ix) &
            + E2*fac(ix)*eff_heat_capacity(ix)/resistivity(ix), 1.0e30_num)
      END IF
    END DO

    ! Pass new temperature values to ghost cells of neighbouring processors
    CALL field_bc(hybrid_Tb, ng)

  END SUBROUTINE ohmic_heating



  SUBROUTINE clear_heat_capacity

    ! De-allocate the heat capacity in each solid. This array is recalculated
    ! each time-step anyway, so there's no need to save and resize it in the
    ! load balancer

    INTEGER :: i_sol

    DO i_sol = 1, solid_count
      DEALLOCATE(solid_array(i_sol)%heat_capacity)
    END DO

  END SUBROUTINE clear_heat_capacity



  SUBROUTINE update_resistivity

    ! Loops over the resistivity grid and updates each point based on the
    ! resistivity model in that cell

    INTEGER :: ix

    DO ix = 1-ng, nx+ng
      SELECT CASE (resistivity_model(ix))
      CASE(c_resist_vacuum)
        resistivity(ix) = 0.0_num
      CASE(c_resist_conductor)
        resistivity(ix) = calc_resistivity_conductor(hybrid_Tb(ix))
      CASE(c_resist_insulator)
        resistivity(ix) = calc_resistivity_insulator(hybrid_Tb(ix))
      END SELECT
    END DO

  END SUBROUTINE update_resistivity



  SUBROUTINE field_zero_curl(field, boundary)

    ! If we are on the simulation boundary, let the ghost cells match those of
    ! the first ghost cell (e.g. by(nx+2:nx+ng) = by(nx+1)). This forces
    ! the gradient across the first two ghost cells to be zero

    INTEGER, INTENT(IN) :: boundary
    REAL(num), DIMENSION(1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, nn

    IF (boundary == c_bd_x_min .AND. x_min_boundary) THEN
      DO i = 1, ng-1
        field(-i) = field(0)
      END DO
    ELSE IF (boundary == c_bd_x_max .AND. x_max_boundary) THEN
      DO i = 1, ng-1
        field(nx+1+i) = field(nx+1)
      END DO
    END IF

  END SUBROUTINE field_zero_curl



  SUBROUTINE hy_grid_centred_var_at_particle(part_x, part_var, grid_var)

    ! Calculates the value of a grid-centred variable "part_var" stored in the
    ! grid "grid_var", averaged over the particle shape for a particle at part_x

    REAL(num), INTENT(IN) :: part_x
    REAL(num), INTENT(IN) :: grid_var(1-ng:nx+ng)
    REAL(num), INTENT(OUT) :: part_var
    INTEGER :: cell_x1
    REAL(num) :: cell_x_r
    REAL(num) :: cell_frac_x
    REAL(num), DIMENSION(sf_min:sf_max) :: gx
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
#else
    cell_x_r = part_x / dx
#endif
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
    part_var = &
          gx(-2) * grid_var(cell_x1-2) &
        + gx(-1) * grid_var(cell_x1-1) &
        + gx( 0) * grid_var(cell_x1  ) &
        + gx( 1) * grid_var(cell_x1+1) &
        + gx( 2) * grid_var(cell_x1+2)
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
    part_var = &
          gx(0) * grid_var(cell_x1  ) &
        + gx(1) * grid_var(cell_x1+1)
#else
#include "triangle/gx.inc"
    part_var = &
          gx(-1) * grid_var(cell_x1-1) &
        + gx( 0) * grid_var(cell_x1  ) &
        + gx( 1) * grid_var(cell_x1+1)
#endif
    part_var = fac*part_var

  END SUBROUTINE hy_grid_centred_var_at_particle

#endif
END MODULE hybrid
