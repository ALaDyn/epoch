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
  USE hy_elastic_scatter
  USE hy_laser
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
        CALL run_hybrid_lasers

        ! .FALSE. this time to use load balancing threshold
        IF (use_balance) CALL balance_workload(.FALSE.)

        ! Move particles, leapfrogging over E and B
        CALL push_particles

        ! Calculate scattering from elastic collisions
        IF (use_hybrid_scatter) THEN
          CALL run_elastic_scatter
        END IF

        ! Obtain heat capacity to calculate the temperature change in
        ! hybrid_collisions and ohmic_heating
        CALL get_heat_capacity

        ! Calculates ionisational energy loss, and updates grid temperature
        IF (use_hybrid_collisions) THEN
          CALL run_ionisation_loss
        END IF

        ! Removes slow particles from the simulation
        CALL remove_slow_particles

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
    INTEGER :: ix, iy, iz, i_sol
    INTEGER :: max_id
    INTEGER :: io, iu

    IF (use_hybrid) THEN

      IF (use_hybrid_collisions) THEN
        CALL setup_hy_ionisation_loss
      END IF

      IF (use_hybrid_scatter) THEN
        CALL setup_hy_elastic_scatter
      END IF

      ! Preset useful constants
      hybrid_const_dx = 1.0_num / (mu0 * dx)
      hybrid_const_dy = 1.0_num / (mu0 * dy)
      hybrid_const_dz = 1.0_num / (mu0 * dz)
      hybrid_const_K_to_eV = kb / q0
      hybrid_const_dt_by_dx = 0.5_num * dt / dx
      hybrid_const_dt_by_dy = 0.5_num * dt / dy
      hybrid_const_dt_by_dz = 0.5_num * dt / dz
      hybrid_idx = 1.0_num/dx
      hybrid_idy = 1.0_num/dy
      hybrid_idz = 1.0_num/dz

      ! Preset useful constants for solids
      ALLOCATE(hybrid_const_heat(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(hybrid_const_sum_ne(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      hybrid_const_sum_ne = 0
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
        solid_array(i_sol)%Iex_term = 2.0_num &
            / (solid_array(i_sol)%hybrid_Iex/mc2)**2
        solid_array(i_sol)%dEdx_C = 1.0_num + 2.0_num &
            * LOG(solid_array(i_sol)%hybrid_Iex/(h_bar*q0)*SQRT(epsilon0 * m0))

        hybrid_const_sum_ne = hybrid_const_sum_ne &
            + solid_array(i_sol)%el_density
      END DO
      hybrid_const_heat = 0.0_num
      DO iz = 1-ng,nz+ng
        DO iy = 1-ng,ny+ng
          DO ix = 1-ng,nx+ng
            IF (hybrid_const_sum_ne(ix,iy,iz) > 0.0_num) THEN
              hybrid_const_heat(ix,iy,iz) = 1.0_num &
                / (kb * hybrid_const_sum_ne(ix,iy,iz)**2)
            END IF
          END DO
        END DO
      END DO

      ! Allocate additional arrays for running in hybrid mode. These require
      ! extra remapping scripts in balance.F90 (for domain change in
      ! load-balance)
      ALLOCATE(resistivity(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(resistivity_model(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(comp_id(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(hybrid_Tb(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(ion_charge(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(ion_density(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
      ALLOCATE(ion_temp(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))

      ! Cycle through all solids and ensure they have a compound ID
      ! First check how many compounds already have IDs
      max_id = 0
      DO i_sol = 1,solid_count
        IF (solid_array(i_sol)%compound_id > max_id) THEN
          max_id = solid_array(i_sol)%compound_id
        END IF
      END DO
      ! Assign compound IDs to the remaining solids, assuming they are all
      ! separate
      DO i_sol = 1,solid_count
        IF (solid_array(i_sol)%compound_id < 0) THEN
          max_id = max_id + 1
          solid_array(i_sol)%compound_id = max_id
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*)
              WRITE(io,*) '*** WARNING ***'
              WRITE(io,*) 'A solid has not been given a compound_id value.'
              WRITE(io,*) 'Code will assume it has no other elements present.'
              WRITE(io,*)
            END DO
          END IF
        END IF
      END DO

      ! Assign a resistivity model to each cell based on present solids
      ! (default to vacuum). Also write the compound_id to the comp_id grid
      resistivity_model = c_resist_vacuum
      comp_id = 0
      DO iz = 1-ng, nz+ng
        DO iy = 1-ng, ny+ng
          DO ix = 1-ng, nx+ng
            max_ne = 0.0_num
            DO i_sol = 1, solid_count
              ! Use the same resistivity model and comp_id as the solid with the
              ! highest local electron number density
              IF (solid_array(i_sol)%el_density(ix,iy,iz) > max_ne) THEN
                max_ne = solid_array(i_sol)%el_density(ix,iy,iz)
                resistivity_model(ix,iy,iz) = solid_array(i_sol)%material
                comp_id(ix,iy,iz) = solid_array(i_sol)%compound_id
              END IF
            END DO
          END DO
        END DO
      END DO

      ! Initialise arrays
      DO iz = 1-ng, nz+ng
        DO iy = 1-ng, ny+ng
          DO ix = 1-ng, nx+ng
            ! Set background temperature and resistivity to those matching
            ! hybrid_Tb_init, which is read in from input.deck
            hybrid_Tb(ix,iy,iz) = hybrid_Tb_init

            ! These arrays currently serve no purpose, but will be used when
            ! ionisation routines are added
            ion_charge(ix,iy,iz) = 0.0_num
            ion_density(ix,iy,iz) = 0.0_num
            ion_temp(ix,iy,iz) = 0.0_num
          END DO
        END DO
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

        ! Check the radiation length is positive (only used in elastic scatter)
        IF (solid_array(i_sol)%rad_len < 0.0_num &
            .AND. use_hybrid_scatter) THEN
          IF (rank == 0) THEN
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Please set all solid radiation lengths to ', &
                  'positive values for hybrid scatter'
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
    calc_resistivity_insulator =  1.0_num &
        / (4.3e5_num + 1.3e3_num*Tb_eV**1.5_num)

  END FUNCTION calc_resistivity_insulator



  SUBROUTINE half_B_step

    ! This subroutine performs a half-step in the magnetic field, assuming a
    ! constant electric field and using:
    !
    ! dB/dt = -curl(E)
    !
    ! We calculate into 1 ghost cell to allow non-zero curls across simulation
    ! boundaries
    !
    ! If two different compounds occur over a boundary, then their resistivities
    ! will be discontinuous, as will their corresponding fields. To prevent
    ! these discontinuities from causing the fields to change non-physically, we
    ! will ignore the gradients over these boundaries (like we do over the
    ! simulation edge)

    INTEGER :: ix, iy, iz
    REAL(num) :: dt_dey_dx, dt_dez_dx, dt_dex_dy, dt_dez_dy
    REAL(num) :: dt_dex_dz, dt_dey_dz

    ! Update B by half a timestep
    DO iz = 0, nz+1
      DO iy = 0, ny+1
        DO ix = 0, nx+1

          ! Calculate derivatives across cells of the same compound
          dt_dey_dx = 0.0_num
          dt_dez_dx = 0.0_num
          dt_dex_dy = 0.0_num
          dt_dez_dy = 0.0_num
          dt_dex_dz = 0.0_num
          dt_dey_dz = 0.0_num
          IF (comp_id(ix,iy,iz) == comp_id(ix+1,iy,iz)) THEN
            dt_dey_dx = (ey(ix+1,iy,iz) - ey(ix,iy,iz)) * hybrid_const_dt_by_dx
            dt_dez_dx = (ez(ix+1,iy,iz) - ez(ix,iy,iz)) * hybrid_const_dt_by_dx
          END IF
          IF (comp_id(ix,iy,iz) == comp_id(ix,iy+1,iz)) THEN
            dt_dex_dy = (ex(ix,iy+1,iz) - ex(ix,iy,iz)) * hybrid_const_dt_by_dy
            dt_dez_dy = (ez(ix,iy+1,iz) - ez(ix,iy,iz)) * hybrid_const_dt_by_dy
          END IF
          IF (comp_id(ix,iy,iz) == comp_id(ix,iy,iz+1)) THEN
            dt_dex_dz = (ex(ix,iy,iz+1) - ex(ix,iy,iz)) * hybrid_const_dt_by_dz
            dt_dey_dz = (ez(ix,iy,iz+1) - ez(ix,iy,iz)) * hybrid_const_dt_by_dz
          END IF

          ! Update fields
          bx(ix, iy, iz) = bx(ix, iy, iz) - dt_dez_dy + dt_dey_dz
          by(ix, iy, iz) = by(ix, iy, iz) - dt_dex_dz + dt_dez_dx
          bz(ix, iy, iz) = bz(ix, iy, iz) - dt_dey_dx + dt_dex_dy
        END DO
      END DO
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

    ! If two different compounds occur over a boundary, then their resistivities
    ! will be discontinuous, as will their corresponding fields. To prevent
    ! these discontinuities from causing the fields to change non-physically, we
    ! will ignore the gradients over these boundaries (like we do over the
    ! simulation edge)

    ! Note that B is staggered from E, whereas J is evaulated at the same point
    ! as E. Resistivity is a cell centred variable. We calculate into 1 ghost
    ! cell to allow non-zero curls across simulation boundaries

    INTEGER :: ix, iy, iz
    REAL(num) :: dby_dx, dbz_dx, dbx_dy, dbz_dy, dbx_dz, dby_dz
    REAL(num) :: res_x, res_y, res_z

    DO iz = 0, nz+1
      DO iy = 0, ny+1
        DO ix = 0, nx+1

          ! Calculate derivatives across cells of the same compound
          dby_dx = 0.0_num
          dbz_dx = 0.0_num
          dbx_dy = 0.0_num
          dbz_dy = 0.0_num
          dbx_dz = 0.0_num
          dby_dz = 0.0_num
          IF (comp_id(ix,iy,iz) == comp_id(ix-1,iy,iz)) THEN
            dby_dx = (by(ix,iy,iz) - by(ix-1,iy,iz)) * hybrid_const_dx
            dbz_dx = (bz(ix,iy,iz) - bz(ix-1,iy,iz)) * hybrid_const_dx
          END IF
          IF (comp_id(ix,iy,iz) == comp_id(ix,iy-1,iz)) THEN
            dbx_dy = (bx(ix,iy,iz) - bx(ix,iy-1,iz)) * hybrid_const_dy
            dbz_dy = (bz(ix,iy,iz) - bz(ix,iy-1,iz)) * hybrid_const_dy
          END IF
          IF (comp_id(ix,iy,iz) == comp_id(ix,iy,iz-1)) THEN
            dbx_dz = (bx(ix,iy,iz) - bx(ix,iy,iz-1)) * hybrid_const_dz
            dby_dz = (by(ix,iy,iz) - by(ix,iy,iz-1)) * hybrid_const_dz
          END IF

          ! Calculate average resistivity at the electric field position (do not
          ! average over two different compounds)
          res_x = resistivity(ix,iy,iz)
          res_y = res_x
          res_z = res_x
          IF (comp_id(ix,iy,iz) == comp_id(ix+1,iy,iz)) THEN
            res_x = 0.5_num * (resistivity(ix+1,iy,iz) + res_x)
          END IF
          IF (comp_id(ix,iy,iz) == comp_id(ix,iy+1,iz)) THEN
            res_y = 0.5_num * (resistivity(ix,iy+1,iz) + res_y)
          END IF
          IF (comp_id(ix,iy,iz) == comp_id(ix,iy,iz+1)) THEN
            res_z = 0.5_num * (resistivity(ix,iy,iz+1) + res_z)
          END IF

          ! Update fields
          ex(ix,iy,iz) = res_x * (dbz_dy - dby_dz - jx(ix,iy,iz))
          ey(ix,iy,iz) = res_y * (dbx_dz - dbz_dx - jy(ix,iy,iz))
          ez(ix,iy,iz) = res_z * (dby_dx - dbx_dy - jz(ix,iy,iz))
        END DO
      END DO
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
    INTEGER :: ix, iy, iz, i_sol

    DO i_sol = 1, solid_count
      ALLOCATE(solid_array(i_sol)%heat_capacity(1-ng:nx+ng,1-ng:ny+ng, &
          1-ng:nz+ng))

      DO iz = 1-ng, nz+ng
        DO iy = 1-ng, ny+ng
          DO ix = 1-ng, nx+ng
            T_prime = solid_array(i_sol)%hybrid_const_ZeV * hybrid_Tb(ix,iy,iz)
            solid_array(i_sol)%heat_capacity(ix,iy,iz) = 0.3_num &
                + 1.2_num * T_prime * (2.2_num + T_prime)/(1.1_num + T_prime)**2
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE get_heat_capacity



  SUBROUTINE run_elastic_scatter

    ! Runs the user-requested elastic scatter model

    IF (elastic_scatter_model == c_hy_davies) THEN
      CALL Davies_elastic_scatter
    ELSE IF (elastic_scatter_model == c_hy_urban) THEN
      CALL Urban_elastic_scatter
    END IF

  END SUBROUTINE run_elastic_scatter



  SUBROUTINE Davies_elastic_scatter

    ! Calculates collisional scattering using the method discussed in  Davies,
    ! et al, (2002). Phys. Rev. E, 65(2), 026407

    INTEGER :: ispecies
    INTEGER(i8) :: ipart
    REAL(num) :: ipart_mc, sum_dtheta, part_ne
    REAL(num) :: p, gamma, v
    REAL(num) :: rand_scatter, delta_theta, delta_phi
    REAL(num) :: part_x, part_y, part_z
    INTEGER :: i_sol
    TYPE(particle), POINTER :: current, next

    ! Precalculate repeated variables
    ipart_mc  = 1.0_num / mc0

    ! Generate a normally distributed random number for the scatter angle.
    ! The argument here refers to the standard deviation. The mean is
    ! assumed zero.
    rand_scatter = random_box_muller(1.0_num)

    ! Loop over all electron species
    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      IF (.NOT. species_list(ispecies)%species_type == c_species_id_electron) &
          CYCLE

      ! Loop over all particles in the species
      DO ipart = 1, species_list(ispecies)%attached_list%count
        next => current%next

        ! Extract particle variables
        part_x = current%part_pos(1) - x_grid_min_local
        part_y = current%part_pos(2) - y_grid_min_local
        part_z = current%part_pos(3) - z_grid_min_local
        p = SQRT(current%part_p(1)**2 + current%part_p(2)**2 &
            + current%part_p(3)**2)
        ! No scatter for immobile electrons
        IF (p < c_tiny) THEN
          current => next
          CYCLE
        END IF
        gamma = SQRT((p * ipart_mc)**2 + 1.0_num)
        v = p * ipart_mc * c / gamma

        ! Extract solid variables averaged over all solids present in this cell
        sum_dtheta = 0.0_num
        DO i_sol = 1, solid_count
          CALL hy_grid_centred_var_at_particle(part_x, part_y, part_z, &
              part_ne, solid_array(i_sol)%el_density)

          ! If ne is zero, then the solid is not in this cell, so don't
          ! calculate the constants
          IF (part_ne < TINY(1.0_num)) CYCLE

          ! Very low energy particles (~50 eV) will make ln_lambda_S go negative
          ! causing a NaN in delta_theta. The hybrid mode isn't designed for low
          ! energy electrons, so delta_theta will be ignored in this case
          sum_dtheta = sum_dtheta + solid_array(i_sol)%theta_fac * part_ne &
              * LOG(MAX(1.0_num, solid_array(i_sol)%hybrid_ln_s * p)) * dt / v
        END DO

        ! Collisional changes to direction
        delta_theta = SQRT(sum_dtheta) / p * rand_scatter
        delta_phi = 2.0_num * pi * random()

        ! Apply rotation
        CALL rotate_p(current, COS(delta_theta), delta_phi, p)

        current => next
      END DO
    END DO

  END SUBROUTINE Davies_elastic_scatter



  SUBROUTINE remove_slow_particles

    ! Removes particles from the simulation if their energy drops below
    ! min_hybrid_energy. If running in PIC-hybrid mode, this removal only occurs
    ! in pure hybrid regions. Ionisation energy loss will reduce the KE of
    ! particles to zero if their total energy is below min_hybrid_energy,
    ! dumping the KE locally as a temperature increase. These particles will
    ! have zero KE when this subroutine is called, so remove all particles with
    ! zero KE

    INTEGER :: ispecies
    TYPE(particle), POINTER :: current, next

    ! Check each electron species
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type /= c_species_id_electron) CYCLE

      ! Cycle through all electrons in this species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        next=>current%next

        ! If the particle has lost all its KE, remove it from the simulation
        IF (current%particle_energy - mc2 < 1.0e-10_num*mc2) THEN
          CALL remove_particle_from_partlist(&
              species_list(ispecies)%attached_list, current)
          IF (track_ejected_particles) THEN
            CALL add_particle_to_partlist(&
                ejected_list(ispecies)%attached_list, current)
          ELSE
            DEALLOCATE(current)
          END IF
        END IF

        current=>next
      END DO
    END DO

  END SUBROUTINE remove_slow_particles



  SUBROUTINE ohmic_heating

    ! Calculates the Ohmic heating of the simulation grid, as described by
    ! Davies, et al, (2002). Phys. Rev. E, 65(2), 026407. At this point in the
    ! simulation, we are at the end of a timestep, with E evaluated in the
    ! middle of this timestep. Assuming this is the average electric field over
    ! the timestep, we have a power per unit volume of E**2/resistivity

    REAL(num) :: E2
    REAL(num) :: eff_heat_capacity(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng)
    REAL(num) :: fac(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng)
    INTEGER :: ix, iy, iz, i_sol

    fac = dt * hybrid_const_heat

    ! Obtain the effective heat capacity: sum(ne/C)
    eff_heat_capacity = 0.0_num
    DO i_sol = 1, solid_count
      DO iz = 1-ng, nz+ng
        DO iy = 1-ng, ny+ng
          DO ix = 1-ng, nx+ng
            eff_heat_capacity(ix,iy,iz) = eff_heat_capacity(ix,iy,iz) &
                + solid_array(i_sol)%el_density(ix,iy,iz) &
                / solid_array(i_sol)%heat_capacity(ix,iy,iz)
          END DO
        END DO
      END DO
    END DO

    ! Loop over all grid points to find temperature change
    DO iz = 1-ng, nz+ng
      DO iy = 1-ng, ny+ng
        DO ix = 1-ng, nx+ng

          ! Sum contributions from squares of electric fields
          E2 = 0.0_num
          IF (comp_id(ix,iy,iz) == comp_id(ix-1,iy,iz)) THEN
            ! Tb is a cell-centred variable, but E has stagger - need to average
            E2 = E2 + (0.5_num*(ex(ix,iy,iz) + ex(ix-1,iy,iz)))**2
          ELSE
            ! Don't average electric fields in different regimes
            E2 = E2 + ex(ix,iy,iz)**2
          END IF
          IF (comp_id(ix,iy,iz) == comp_id(ix,iy-1,iz)) THEN
            E2 = E2 + (0.5_num*(ey(ix,iy,iz) + ey(ix,iy-1,iz)))**2
          ELSE
            E2 = E2 + ey(ix,iy,iz)**2
          END IF
          IF (comp_id(ix,iy,iz) == comp_id(ix,iy,iz-1)) THEN
            E2 = E2 + (0.5_num*(ez(ix,iy,iz) + ez(ix,iy,iz-1)))**2
          ELSE
            E2 = E2 + ez(ix,iy,iz)**2
          END IF

          ! Calculate Ohmic heating, avoiding the 0/0 NaN. Cap temp at 1.0e30
          IF (resistivity(ix,iy,iz) > 0.0_num .AND. &
              eff_heat_capacity(ix,iy,iz) > 0.0_num) THEN
            hybrid_Tb(ix,iy,iz) =   MIN(hybrid_Tb(ix,iy,iz) &
                + E2*fac(ix,iy,iz)*eff_heat_capacity(ix,iy,iz) &
                /resistivity(ix,iy,iz), 1.0e30_num)
          END IF
        END DO
      END DO
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

    INTEGER :: ix, iy, iz

    DO iz = 1-ng, nz+ng
      DO iy = 1-ng, ny+ng
        DO ix = 1-ng, nx+ng
          SELECT CASE (resistivity_model(ix,iy,iz))
          CASE(c_resist_vacuum)
            resistivity(ix,iy,iz) = 0.0_num
          CASE(c_resist_conductor)
            resistivity(ix,iy,iz) = &
                calc_resistivity_conductor(hybrid_Tb(ix,iy,iz))
          CASE(c_resist_insulator)
            resistivity(ix,iy,iz) = &
                calc_resistivity_insulator(hybrid_Tb(ix,iy,iz))
          END SELECT
        END DO
      END DO
    END DO

  END SUBROUTINE update_resistivity



  SUBROUTINE field_zero_curl(field, boundary)

    ! If we are on the simulation boundary, let the ghost cells match those of
    ! the first ghost cell (e.g. by(nx+2:nx+ng,1,1) = by(nx+1,1,1)). This forces
    ! the gradient across the first two ghost cells to be zero

    INTEGER, INTENT(IN) :: boundary
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, nn

    IF (boundary == c_bd_x_min .AND. x_min_boundary) THEN
      DO i = 1, ng-1
        field(-i,:,:) = field(0,:,:)
      END DO
    ELSE IF (boundary == c_bd_x_max .AND. x_max_boundary) THEN
      DO i = 1, ng
        field(nx+i,:,:) = field(nx,:,:)
      END DO

    ELSE IF (boundary == c_bd_y_min .AND. y_min_boundary) THEN
      DO i = 1, ng-1
        field(:,-i,:) = field(:,0,:)
      END DO
    ELSE IF (boundary == c_bd_y_max .AND. y_max_boundary) THEN
      DO i = 1, ng
        field(:,ny+i,:) = field(:,ny,:)
      END DO

    ELSE IF (boundary == c_bd_z_min .AND. z_min_boundary) THEN
      DO i = 1, ng-1
        field(:,:,-i) = field(:,:,0)
      END DO
    ELSE IF (boundary == c_bd_z_max .AND. z_max_boundary) THEN
      DO i = 1, ng
        field(:,:,nz+i) = field(:,:,nz)
      END DO
    END IF

  END SUBROUTINE field_zero_curl

#endif
END MODULE hybrid
