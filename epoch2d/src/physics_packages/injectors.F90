! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2017 Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE injectors

  USE shared_data
  USE partlist
  USE particle_temperature
  USE evaluator
  USE random_generator

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_injector(boundary, injector)

    INTEGER, INTENT(IN) :: boundary
    TYPE(injector_block), INTENT(INOUT) :: injector

    injector%npart_per_cell = 0
    injector%species = -1
    injector%boundary = boundary
    injector%t_start = 0.0_num
    injector%t_end = t_end
    injector%has_t_end = .FALSE.
    injector%density_min = 0.0_num
    injector%use_flux_injector = .FALSE.

    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      ALLOCATE(injector%dt_inject(1-ng:ny+ng))
      ALLOCATE(injector%depth(1-ng:ny+ng))
    END IF

    IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
      ALLOCATE(injector%dt_inject(1-ng:nx+ng))
      ALLOCATE(injector%depth(1-ng:nx+ng))
    END IF

    injector%depth = 1.0_num
    injector%dt_inject = -1.0_num
    NULLIFY(injector%next)

    injector%inject_from_file = .FALSE.
    injector%file_finished = .FALSE.
    injector%x_data_given = .FALSE.
    injector%y_data_given = .FALSE.
    injector%px_data_given = .FALSE.
    injector%py_data_given = .FALSE.
    injector%pz_data_given = .FALSE.
    injector%t_data_given = .FALSE.
#ifndef PER_SPECIES_WEIGHT
    injector%w_data_given = .FALSE.
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
    injector%id_data_given = .FALSE.
#endif

  END SUBROUTINE init_injector



  SUBROUTINE attach_injector(injector)

    TYPE(injector_block), POINTER :: injector
    INTEGER :: boundary

    boundary = injector%boundary

    IF (boundary == c_bd_x_min) THEN
      CALL attach_injector_to_list(injector_x_min, injector)
    ELSE IF (boundary == c_bd_x_max) THEN
      CALL attach_injector_to_list(injector_x_max, injector)
    ELSE IF (boundary == c_bd_y_min) THEN
      CALL attach_injector_to_list(injector_y_min, injector)
    ELSE IF (boundary == c_bd_y_max) THEN
      CALL attach_injector_to_list(injector_y_max, injector)
    END IF

  END SUBROUTINE attach_injector



  ! Actually does the attaching of the injector to the correct list
  SUBROUTINE attach_injector_to_list(list, injector)

    TYPE(injector_block), POINTER :: list
    TYPE(injector_block), POINTER :: injector
    TYPE(injector_block), POINTER :: current

    IF (ASSOCIATED(list)) THEN
      current => list
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => injector
    ELSE
      list => injector
    END IF

  END SUBROUTINE attach_injector_to_list



  SUBROUTINE deallocate_injectors

    CALL deallocate_injector_list(injector_x_min)
    CALL deallocate_injector_list(injector_x_max)
    CALL deallocate_injector_list(injector_y_min)
    CALL deallocate_injector_list(injector_y_max)

  END SUBROUTINE deallocate_injectors



  SUBROUTINE deallocate_injector_list(list)

    TYPE(injector_block), POINTER :: list
    TYPE(injector_block), POINTER :: current, next
    INTEGER :: i

    current => list
    DO WHILE(ASSOCIATED(current))
      next => current%next
      IF (current%density_function%init) &
          CALL deallocate_stack(current%density_function)
      IF (ASSOCIATED(current%dt_inject)) DEALLOCATE(current%dt_inject)
      IF (ASSOCIATED(current%depth)) DEALLOCATE(current%depth)
      DO i = 1, 3
        IF (current%temperature_function(i)%init) &
            CALL deallocate_stack(current%temperature_function(i))
        IF (current%drift_function(i)%init) &
            CALL deallocate_stack(current%drift_function(i))
      END DO
      DEALLOCATE(current)
      current => next
    END DO

  END SUBROUTINE deallocate_injector_list



  SUBROUTINE open_injector_files(injector)

    ! Called in deck_injector_block after we have read the injector variables,
    ! and only for injectors with the inject_from_file flag. The file units are
    ! chosen such that each variable for each injector has a unique unit.

    TYPE(injector_block), POINTER :: injector

    ! Injection times (mandatory)
    OPEN(UNIT=custom_base_unit + (injector%custom_id-1)*9 + unit_t, &
        FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%t_data))

#ifndef PER_SPECIES_WEIGHT
    ! Particle weights (mandatory)
    OPEN(UNIT=custom_base_unit + (injector%custom_id-1)*9 + unit_w, &
        FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%w_data))
#endif

    ! Position (positions that aren't in boundary dimension are mandatory)
    IF (injector%x_data_given) THEN
    OPEN(UNIT=custom_base_unit + (injector%custom_id-1)*9 + unit_x, &
          FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%x_data))
    END IF
    IF (injector%y_data_given) THEN
    OPEN(UNIT=custom_base_unit + (injector%custom_id-1)*9 + unit_y, &
          FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%y_data))
    END IF

    ! Momentum, optional - if missing, will be set to 0
    IF (injector%px_data_given) THEN
    OPEN(UNIT=custom_base_unit + (injector%custom_id-1)*9 + unit_px, &
          FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%px_data))
    END IF
    IF (injector%py_data_given) THEN
      OPEN(UNIT=custom_base_unit + (injector%custom_id-1)*9 + unit_py, &
          FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%py_data))
    END IF
    IF (injector%pz_data_given) THEN
      OPEN(UNIT=custom_base_unit + (injector%custom_id-1)*9 + unit_pz, &
          FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%pz_data))
    END IF

#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
    ! Particle ID, optional - if missing, particle ID will be left empty
    IF (injector%id_data_given) THEN
      OPEN(UNIT=custom_base_unit + (injector%custom_id-1)*9 + unit_id, &
          FILE=TRIM(data_dir) // '/' // TRIM(injector_filenames%id_data))
    END IF
#endif

  END SUBROUTINE open_injector_files



  SUBROUTINE run_injectors

    TYPE(injector_block), POINTER :: current

    IF (x_min_boundary) THEN
      current => injector_x_min
      DO WHILE(ASSOCIATED(current))
        IF (.NOT. current%inject_from_file) THEN
          CALL run_single_injector(current, c_bd_x_min)
        END IF
        current => current%next
      END DO
    END IF

    IF (x_max_boundary) THEN
      current => injector_x_max
      DO WHILE(ASSOCIATED(current))
        IF (.NOT. current%inject_from_file) THEN
          CALL run_single_injector(current, c_bd_x_max)
        END IF
        current => current%next
      END DO
    END IF

    IF (y_min_boundary) THEN
      current => injector_y_min
      DO WHILE(ASSOCIATED(current))
        IF (.NOT. current%inject_from_file) THEN
          CALL run_single_injector(current, c_bd_y_min)
        END IF
        current => current%next
      END DO
    END IF

    IF (y_max_boundary) THEN
      current => injector_y_max
      DO WHILE(ASSOCIATED(current))
        IF (.NOT. current%inject_from_file) THEN
          CALL run_single_injector(current, c_bd_y_max)
        END IF
        current => current%next
      END DO
    END IF

    ! Injectors exist on all ranks. Above, the IF (..._boundary) lines prevent
    ! non-boundary processors from generating particles. When injecting from
    ! file, all processors must know which particles have been added, so that if
    ! the load balancer swaps injectors around they can continue where the other
    ! left off. Hence, all processors must read the files (only boundary
    ! processors will add the particles though)
    current => injector_x_max
    DO WHILE(ASSOCIATED(current))
      IF (current%inject_from_file) THEN
        CALL run_file_injection(current)
      END IF
      current => current%next
    END DO

    current => injector_x_min
    DO WHILE(ASSOCIATED(current))
      IF (current%inject_from_file) THEN
        CALL run_file_injection(current)
      END IF
      current => current%next
    END DO

    current => injector_y_max
    DO WHILE(ASSOCIATED(current))
      IF (current%inject_from_file) THEN
        CALL run_file_injection(current)
      END IF
      current => current%next
    END DO

    current => injector_y_min
    DO WHILE(ASSOCIATED(current))
      IF (current%inject_from_file) THEN
        CALL run_file_injection(current)
      END IF
      current => current%next
    END DO

  END SUBROUTINE run_injectors



  SUBROUTINE run_single_injector(injector, direction)

    TYPE(injector_block), POINTER :: injector
    INTEGER, INTENT(IN) :: direction
    REAL(num) :: bdy_pos, bdy_space
    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    REAL(num) :: mass, typical_mc2, p_therm, p_inject_drift, density_grid
    REAL(num) :: gamma_mass, v_inject, density, vol
    REAL(num) :: npart_ideal, itemp, v_inject_s
    REAL(num), DIMENSION(3) :: temperature, drift
    INTEGER :: parts_this_time, ipart, idir, dir_index, ii
    INTEGER, DIMENSION(c_ndims-1) :: perp_dir_index, nel
    REAL(num), DIMENSION(c_ndims-1) :: perp_cell_size, cur_cell
    TYPE(parameter_pack) :: parameters
    REAL(num), DIMENSION(3) :: dir_mult
    LOGICAL :: first_inject, flux_fn

    IF (time < injector%t_start .OR. time > injector%t_end) RETURN

    ! If you have a moving window that has started moving then unless you
    ! EXPLICITLY give a t_end value to the injector stop the injector
    IF (move_window .AND. window_started .AND. .NOT. injector%has_t_end) &
        RETURN

    flux_fn = .FALSE.
    dir_mult = 1.0_num

    IF (direction == c_bd_x_min) THEN
      parameters%pack_ix = 0
      nel = (/ny/)
      perp_cell_size = (/dy/)
      perp_dir_index = (/2/)
      dir_index = 1
      bdy_pos = x_min
      bdy_space = -dx
      IF (injector%use_flux_injector) THEN
        flux_fn = .TRUE.
        dir_mult(dir_index) = 1.0_num
      END IF
    ELSE IF (direction == c_bd_x_max) THEN
      parameters%pack_ix = nx
      nel = (/ny/)
      perp_cell_size = (/dy/)
      perp_dir_index = (/2/)
      dir_index = 1
      bdy_pos = x_max
      bdy_space = dx
      IF (injector%use_flux_injector) THEN
        flux_fn = .TRUE.
        dir_mult(dir_index) = -1.0_num
      END IF
    ELSE IF (direction == c_bd_y_min) THEN
      parameters%pack_iy = 0
      nel = (/nx/)
      perp_cell_size = (/dx/)
      perp_dir_index = (/1/)
      dir_index = 2
      bdy_pos = y_min
      bdy_space = -dy
      IF (injector%use_flux_injector) THEN
        flux_fn = .TRUE.
        dir_mult(dir_index) = 1.0_num
      END IF
    ELSE IF (direction == c_bd_y_max) THEN
      parameters%pack_iy = ny
      nel = (/nx/)
      perp_cell_size = (/dx/)
      perp_dir_index = (/1/)
      dir_index = 2
      bdy_pos = y_max
      bdy_space = dy
      IF (injector%use_flux_injector) THEN
        flux_fn = .TRUE.
        dir_mult(dir_index) = -1.0_num
      END IF
    ELSE
      RETURN
    END IF

    vol = ABS(bdy_space)
    DO idir = 1, c_ndims-1
      vol = vol * perp_cell_size(idir)
    END DO

    mass = species_list(injector%species)%mass
    typical_mc2 = (mass * c)**2
    cur_cell = 0.0_num

    CALL create_empty_partlist(plist)
    DO ii = 1, nel(1)
      DO idir = 1, c_ndims-1
        IF (perp_dir_index(idir) == 1) cur_cell(idir) = x(ii)
        IF (perp_dir_index(idir) == 2) cur_cell(idir) = y(ii)
      END DO

      parameters%use_grid_position = .TRUE.
      CALL assign_pack_value(parameters, perp_dir_index(1), ii)

      IF (injector%dt_inject(ii) > 0.0_num) THEN
        npart_ideal = dt / injector%dt_inject(ii)
        itemp = random_box_muller(0.5_num * SQRT(npart_ideal &
            * (1.0_num - npart_ideal / REAL(injector%npart_per_cell, num)))) &
            + npart_ideal
        injector%depth(ii) = injector%depth(ii) - itemp
        first_inject = .FALSE.

        IF (injector%depth(ii) >= 0.0_num) CYCLE
      ELSE
        first_inject = .TRUE.
      END IF

      CALL populate_injector_properties(injector, parameters, density_grid, &
          temperature, drift)

      IF (density_grid < injector%density_min) CYCLE

      ! Assume agressive maximum thermal momentum, all components
      ! like hottest component
      p_therm = SQRT(mass * kb * MAXVAL(temperature))
      p_inject_drift = drift(dir_index)
      gamma_mass = SQRT((p_therm + p_inject_drift)**2 + typical_mc2) / c
      v_inject_s = p_inject_drift / gamma_mass
      v_inject = ABS(v_inject_s)

      injector%dt_inject(ii) = ABS(bdy_space) &
          / MAX(injector%npart_per_cell * v_inject, c_tiny)
      IF (first_inject) THEN
        ! On the first run of the injectors it isn't possible to decrement
        ! the optical depth until this point
        npart_ideal = dt / injector%dt_inject(ii)
        itemp = random_box_muller(0.5_num * SQRT(npart_ideal &
            * (1.0_num - npart_ideal / REAL(injector%npart_per_cell, num)))) &
            + npart_ideal
        injector%depth(ii) = injector%depth(ii) - itemp
      END IF

      parts_this_time = FLOOR(ABS(injector%depth(ii) - 1.0_num))
      injector%depth(ii) = injector%depth(ii) + REAL(parts_this_time, num)

      DO ipart = 1, parts_this_time
        CALL create_particle(new)

        new%part_pos = 0.0_num
        DO idir = 1, c_ndims-1
          new%part_pos(perp_dir_index(idir)) = &
              (random() - 0.5_num) * perp_cell_size(idir) + cur_cell(idir)
        END DO

        new%part_pos(dir_index) = bdy_pos + 0.5_num * bdy_space * png &
            - random() * v_inject_s * dt
        parameters%pack_pos = new%part_pos
        parameters%use_grid_position = .FALSE.

        CALL populate_injector_properties(injector, parameters, density, &
            temperature, drift)

        DO idir = 1, 3
          IF (flux_fn) THEN
            new%part_p(idir) = flux_momentum_from_temperature(mass, &
                temperature(idir), drift(idir)) * dir_mult(idir)
          ELSE
            new%part_p(idir) = momentum_from_temperature(mass, &
                temperature(idir), drift(idir))
          END IF
        END DO
#ifdef PER_PARTICLE_CHARGE_MASS
        new%charge = species_list(injector%species)%charge
        new%mass = mass
#endif
#ifndef PER_SPECIES_WEIGHT
        new%weight = vol * density / REAL(injector%npart_per_cell, num)
#endif
        CALL add_particle_to_partlist(plist, new)
      END DO
    END DO

    CALL append_partlist(species_list(injector%species)%attached_list, plist)

  END SUBROUTINE run_single_injector



  SUBROUTINE populate_injector_properties(injector, parameters, density, &
      temperature, drift)

    TYPE(injector_block), POINTER :: injector
    TYPE(parameter_pack), INTENT(IN) :: parameters
    REAL(num), INTENT(OUT) :: density
    REAL(num), DIMENSION(3), INTENT(OUT) :: temperature, drift
    INTEGER :: errcode, i

    errcode = 0
    density = MAX(evaluate_with_parameters(injector%density_function, &
        parameters, errcode), 0.0_num)

    ! Stack can only be time varying if valid. Change if this isn't true
    DO i = 1, 3
      IF (injector%temperature_function(i)%init) THEN
        temperature(i) = &
            MAX(evaluate_with_parameters(injector%temperature_function(i), &
                parameters, errcode), 0.0_num)
      ELSE
        temperature(i) = 0.0_num
      END IF
      IF (injector%drift_function(i)%init) THEN
        drift(i) = &
            evaluate_with_parameters(injector%drift_function(i), &
                                     parameters, errcode)
      ELSE
        drift(i) = 0.0_num
      END IF
    END DO

    IF (errcode /= c_err_none) CALL abort_code(errcode)

  END SUBROUTINE populate_injector_properties



  SUBROUTINE assign_pack_value(parameters, dir_index, p_value)

    TYPE(parameter_pack), INTENT(INOUT) :: parameters
    INTEGER, INTENT(IN) :: dir_index
    INTEGER, INTENT(IN) :: p_value

    IF (dir_index == 1) THEN
      parameters%pack_ix = p_value
    ELSE IF (dir_index == 2) THEN
      parameters%pack_iy = p_value
    END IF

  END SUBROUTINE assign_pack_value



  SUBROUTINE run_file_injection(injector)

    ! This subroutine reads particles from files (opened after the full injector
    ! block has been read, in a call from deck_injector_block.f90). We read from
    ! the injection time file (t_data) until we hit a particle injected after
    ! the NEXT timestep (time + dt). We position injected particles such that
    ! they pass the injection boundary at the time specified in t_data
    ! (assuming constant velocity), only injecting with processors on the
    ! appropriate boundary

    TYPE(injector_block), POINTER :: injector
    REAL(num) :: mass, inv_m2c2
    REAL(num) :: x_in, y_in, px_in, py_in, pz_in
#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: w_in
#endif
#ifdef PARTICLE_ID4
    INTEGER :: id_in
#elif PARTICLE_ID
    INTEGER(i8) :: id_in
#endif
    INTEGER :: boundary
    REAL(num) :: next_time
    REAL(num) :: vx, vy, gamma, inv_gamma_mass, x_start, y_start, time_to_bdy
    TYPE(particle), POINTER :: new
    TYPE(particle_list) :: plist
    LOGICAL :: no_particles_added, skip_processor

    ! Set to true if any of the associated files have reached the end of file
    IF (injector%file_finished) RETURN

    mass = species_list(injector%species)%mass
    inv_m2c2 = 1.0_num/(mass*c)**2
    no_particles_added = .TRUE.

    ! Add particles until we reach an injection time greater than the next
    ! timestep (particles must pass the injection boundary in the following
    ! particle push), or until there are no more particles to add
    DO
      ! We always start with injector%next_time known. Global time is a half
      ! timestep ahead of particle time when this is called
      IF (.NOT. injector%next_time < time + 0.5_num*dt ) EXIT

      ! If on x boundary read the y position, and vice versa
      boundary = injector%boundary
      IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
        CALL read_injector_real(unit_y, y_start, injector)
      ELSE
        CALL read_injector_real(unit_x, x_start, injector)
      END IF

#ifndef PER_SPECIES_WEIGHT
      ! Read weight data
      CALL read_injector_real(unit_w, w_in, injector)
#endif

      ! Read momentum data
      IF (injector%px_data_given) THEN
        CALL read_injector_real(unit_px, px_in, injector)
      ELSE
        px_in = 0.0_num
      END IF

      IF (injector%py_data_given) THEN
        CALL read_injector_real(unit_py, py_in, injector)
      ELSE
        py_in = 0.0_num
      END IF

      IF (injector%pz_data_given) THEN
        CALL read_injector_real(unit_pz, pz_in, injector)
      ELSE
        pz_in = 0.0_num
      END IF

      ! Read id data
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
      IF (injector%id_data_given) THEN
#ifdef PARTICLE_ID4
        CALL read_injector_int4(unit_id, id_in, injector)
#elif PARTICLE_ID
        CALL read_injector_int(unit_id, id_in, injector)
#endif
      ELSE
        id_in = 0
      END IF
#endif

      ! Ensure we still have values for each variable here
      IF (injector%file_finished) EXIT

      ! Identify processors which aren't adding this particle to the simulation
      boundary = injector%boundary
      skip_processor = .FALSE.
      IF (boundary == c_bd_x_min) THEN
        IF (.NOT. x_min_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_x_max) THEN
        IF (.NOT. x_max_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_y_min) THEN
        IF (.NOT. y_min_boundary) skip_processor = .TRUE.
      ELSE IF (boundary == c_bd_y_max) THEN
        IF (.NOT. y_max_boundary) skip_processor = .TRUE.
      END IF

      ! Skip the following phases if this rank isn't on the right boundary
      IF (skip_processor) THEN
        ! Calculate time for the next particle to be injected
        CALL read_injector_real(unit_t, next_time, injector)
        ! If there are no more particles to add, exit the loop, otherwise save
        ! the next injection time
        IF (injector%file_finished) EXIT
        injector%next_time = next_time
        CYCLE
      END IF

      ! If code has been restarted from an output dump, the position in our file
      ! will be lost. All particles with injection time < simulation time were
      ! injected in the previous run, and shouldn't be added again
      IF (injector%next_time < time - 0.5_num*dt) THEN
        CALL read_injector_real(unit_t, next_time, injector)
        IF (.NOT. injector%file_finished) injector%next_time = next_time
        CYCLE
      END IF

      ! Only ranks on the same boundary as the particle can reach here
      ! Calculate particle velocity
      gamma = SQRT(1.0_num + (px_in**2 + py_in**2 + pz_in**2)*inv_m2c2)
      inv_gamma_mass = 1.0_num/(gamma*mass)
      vx = px_in*inv_gamma_mass
      vy = py_in*inv_gamma_mass

      ! Calculate position of injection such that paritlces reach the boundary
      ! at next_time. Note that global time is a half timestep ahead of the time
      ! our particles are at
      time_to_bdy = (injector%next_time - (time-0.5_num*dt))
      IF (boundary == c_bd_x_min) THEN
        x_in = x_min - time_to_bdy * vx
        y_in = y_start - time_to_bdy * vy
      ELSE IF (boundary == c_bd_x_max) THEN
        x_in = x_max - time_to_bdy * vx
        y_in = y_start - time_to_bdy * vy
      ELSE IF (boundary == c_bd_y_min) THEN
        x_in = x_start - time_to_bdy * vx
        y_in = y_min - time_to_bdy * vy
      ELSE IF (boundary == c_bd_y_max) THEN
        x_in = x_start - time_to_bdy * vx
        y_in = y_max - time_to_bdy * vy
      END IF

      ! Now we have start position, find the processor which will be adding the
      ! particle
      IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
        ! Skip all processors which are at the wrong y position
        IF (y_in <= y_min) THEN
          IF (.NOT. y_min_boundary) THEN
            skip_processor = .TRUE.
          END IF
        ELSE IF (y_in > y_max) THEN
          IF (.NOT. y_max_boundary) THEN
            skip_processor = .TRUE.
          END IF
        ELSE
          IF (y_in <= y_min_local .OR. y_in > y_max_local) THEN
            skip_processor = .TRUE.
          END IF
        END IF

      ELSE IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
        ! Skip all processors which are at the wrong x position
        IF (x_in <= x_min) THEN
          IF (.NOT. x_min_boundary) THEN
            skip_processor = .TRUE.
          END IF
        ELSE IF (x_in > x_max) THEN
          IF (.NOT. x_max_boundary) THEN
            skip_processor = .TRUE.
          END IF
        ELSE
          IF (x_in <= x_min_local .OR. x_in > x_max_local) THEN
            skip_processor = .TRUE.
          END IF
        END IF
      END IF

      ! Skip the particle adding phase if this rank isn't adding the particle
      IF (skip_processor) THEN
        ! Calculate time for the next particle to be injected
        CALL read_injector_real(unit_t, next_time, injector)
        ! If there are no more particles to add, exit the loop, otherwise save
        ! the next injection time
        IF (injector%file_finished) EXIT
        injector%next_time = next_time
        CYCLE
      END IF

      ! Create the particle and assign properties
      CALL create_particle(new)
      new%part_pos(1) = x_in
      new%part_pos(2) = y_in
      new%part_p(1) = px_in
      new%part_p(2) = py_in
      new%part_p(3) = pz_in
#ifdef PER_PARTICLE_CHARGE_MASS
      new%charge = species_list(injector%species)%charge
      new%mass = mass
#endif
#ifndef PER_SPECIES_WEIGHT
      new%weight = w_in
#endif
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
      new%id = id_in
#endif

      ! Add the particle to our dummy particle list plist
      IF (no_particles_added) THEN
        ! Create an empty particle list the first time a particle is added
        CALL create_empty_partlist(plist)
        no_particles_added = .FALSE.
      END IF
      CALL add_particle_to_partlist(plist, new)

      ! Calculate the next time
      CALL read_injector_real(unit_t, next_time, injector)
      IF (injector%file_finished) THEN
        EXIT
      ELSE
        injector%next_time = next_time
      END IF
    END DO

    ! Append particles to the main particle list for this species on this
    ! processor if we have added particles
    IF (.NOT. no_particles_added) THEN
      CALL append_partlist(species_list(injector%species)%attached_list, plist)
    END IF

  END SUBROUTINE run_file_injection



  SUBROUTINE read_injector_real(unit_code, value, injector)

    ! Opens the file specified by the unit_code, reads a real variable and
    ! changes value to match. Use the unit codes defined in shared_data.F90 -
    ! e.g. unit_t for the t_data file. Also checks for the end of file condition

    TYPE(injector_block), POINTER :: injector
    REAL(num) :: value
    INTEGER :: unit_code, eof

    ! Read file, checking a value has been assigned to time_in, and that we have
    ! not reached the end of the file
    READ(custom_base_unit + (injector%custom_id-1)*9 + unit_code, *, &
        IOSTAT = eof) value

    ! File read successfully
    IF (eof == 0) RETURN

    ! Triggered if there are no more values to read
    IF (eof < 0) THEN
      injector%file_finished = .TRUE.

    ! Triggered if illegal value entered
    ELSE
      IF(rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Illegal value found in read_injector_real'
        PRINT*,'Injected particle will not behave as expected'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

  END SUBROUTINE



#ifdef PARTICLE_ID
  SUBROUTINE read_injector_int(unit_code, value, injector)

    ! Opens the file specified by the unit_code, reads an integer variable and
    ! changes value to match. Use the unit codes defined in shared_data.F90 -
    ! e.g. unit_t for the t_data file. Also checks for the end of file condition

    TYPE(injector_block), POINTER :: injector
    INTEGER(i8) :: value
    INTEGER :: unit_code, eof

    ! Read file, checking a value has been assigned to time_in, and that we have
    ! not reached the end of the file
    READ(custom_base_unit + (injector%custom_id-1)*9 + unit_code, *, &
        IOSTAT = eof) value

    ! File read successfully
    IF (eof == 0) RETURN

    ! Triggered if there are no more values to read
    IF (eof < 0) THEN
      injector%file_finished = .TRUE.

    ! Triggered if illegal value entered
    ELSE
      IF(rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Illegal value found in read_injector_real'
        PRINT*,'Injected particle will not behave as expected'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

  END SUBROUTINE
#endif



#ifdef PARTICLE_ID4
  SUBROUTINE read_injector_int4(unit_code, value, injector)

    ! Opens the file specified by the unit_code, reads an integer variable and
    ! changes value to match. Use the unit codes defined in shared_data.F90 -
    ! e.g. unit_t for the t_data file. Also checks for the end of file condition

    TYPE(injector_block), POINTER :: injector
    INTEGER :: value
    INTEGER :: unit_code, eof

    ! Read file, checking a value has been assigned to time_in, and that we have
    ! not reached the end of the file
    READ(custom_base_unit + (injector%custom_id-1)*9 + unit_code, *, &
        IOSTAT = eof) value

    ! File read successfully
    IF (eof == 0) RETURN

    ! Triggered if there are no more values to read
    IF (eof < 0) THEN
      injector%file_finished = .TRUE.

    ! Triggered if illegal value entered
    ELSE
      IF(rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Illegal value found in read_injector_real'
        PRINT*,'Injected particle will not behave as expected'
        CALL abort_code(c_err_bad_value)
      END IF
    END IF

  END SUBROUTINE
#endif

END MODULE injectors
