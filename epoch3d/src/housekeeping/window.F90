! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009      Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE window

  USE boundary
  USE partlist
  USE evaluator

  IMPLICIT NONE

  LOGICAL, SAVE :: window_started
  REAL(num), ALLOCATABLE :: density(:,:), temperature(:,:,:), drift(:,:,:)
  REAL(num), SAVE :: window_shift_fraction

CONTAINS

  SUBROUTINE initialise_window

#ifdef PER_SPECIES_WEIGHT
    INTEGER :: ierr
#endif

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    ALLOCATE(density(-2:ny+3,-2:nz+3))
    ALLOCATE(temperature(-2:ny+3,-2:nz+3, 1:3))
    ALLOCATE(drift(-2:ny+3,-2:nz+3, 1:3))
    window_started = .FALSE.
#else
    IF (rank == 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    ENDIF
    CALL abort_code(c_err_pp_options_missing)
#endif

  END SUBROUTINE initialise_window



  SUBROUTINE deallocate_window

    IF (ALLOCATED(density)) DEALLOCATE(density)
    IF (ALLOCATED(temperature)) DEALLOCATE(temperature)
    IF (ALLOCATED(drift)) DEALLOCATE(drift)

  END SUBROUTINE deallocate_window



#ifndef PER_SPECIES_WEIGHT
  SUBROUTINE shift_window(window_shift_cells)

    INTEGER, INTENT(IN) :: window_shift_cells
    INTEGER :: iwindow

    ! Shift the window round one cell at a time.
    ! Inefficient, but it works
    DO iwindow = 1, window_shift_cells
      CALL insert_particles

      ! Shift the box around
      x_grid_min = x_grid_min + dx
      x_grid_max = x_grid_max + dx
      x_grid_mins = x_grid_mins + dx
      x_grid_maxs = x_grid_maxs + dx
      x_grid_min_local = x_grid_min_local + dx
      x_grid_max_local = x_grid_max_local + dx
      x_min = x_min + dx
      x_max = x_max + dx
      x_min_local = x_min_local + dx
      x_max_local = x_max_local + dx

      x = x + dx
      x_global = x_global + dx
      xb_global = xb_global + dx

      CALL remove_particles

      ! Shift fields around
      CALL shift_fields
    ENDDO

  END SUBROUTINE shift_window



  SUBROUTINE shift_fields

    INTEGER :: j, k

    CALL shift_field(ex, ng)
    CALL shift_field(ey, ng)
    CALL shift_field(ez, ng)

    CALL shift_field(bx, ng)
    CALL shift_field(by, ng)
    CALL shift_field(bz, ng)

    CALL shift_field(jx, jng)
    CALL shift_field(jy, jng)
    CALL shift_field(jz, jng)

    IF (x_max_boundary) THEN
      DO k = -2, nz+3
        DO j = -2, ny+3
          ! Fix incoming field cell.
          ex(nx,j,k)   = ex_x_max(j,k)
          ex(nx+1,j,k) = ex_x_max(j,k)
          ey(nx+1,j,k) = ey_x_max(j,k)
          ez(nx+1,j,k) = ez_x_max(j,k)
          ex(nx-1,j,k) = 0.5_num * (ex(nx-2,j,k) + ex(nx,j,k))
          ey(nx,j,k)   = 0.5_num * (ey(nx-1,j,k) + ey(nx+1,j,k))
          ez(nx,j,k)   = 0.5_num * (ez(nx-1,j,k) + ez(nx+1,j,k))
          bx(nx+1,j,k) = bx_x_max(j,k)
          by(nx,j,k)   = by_x_max(j,k)
          bz(nx,j,k)   = bz_x_max(j,k)
          bx(nx,j,k)   = 0.5_num * (bx(nx-1,j,k) + bx(nx+1,j,k))
          by(nx-1,j,k) = 0.5_num * (by(nx-2,j,k) + by(nx,j,k))
          bz(nx-1,j,k) = 0.5_num * (bz(nx-2,j,k) + bz(nx,j,k))
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE shift_fields



  SUBROUTINE shift_field(field, ng)

    INTEGER, INTENT(IN) :: ng
    REAL(num), DIMENSION(1-ng:,1-ng:,1-ng:), INTENT(INOUT) :: field
    INTEGER :: i, j, k

    ! Shift field to the left by one cell
    DO k = 1-ng, nz+ng
    DO j = 1-ng, ny+ng
    DO i = 1-ng, nx+ng-1
      field(i,j,k) = field(i+1,j,k)
    ENDDO
    ENDDO
    ENDDO

    CALL field_bc(field, ng)

  END SUBROUTINE shift_field



  SUBROUTINE insert_particles

    TYPE(particle), POINTER :: current
    TYPE(particle_list), POINTER :: append_list
    INTEGER :: ispecies, i, iy, iz
    INTEGER(i8) :: ipart, npart_per_cell, n0
    REAL(num), DIMENSION(3) :: temp_local, drift_local
    REAL(num) :: weight_local, density_local, mass, npart_frac
    TYPE(parameter_pack) :: parameters
    TYPE(particle_species), POINTER :: species

    ! This subroutine injects particles at the right hand edge of the box

    ! Only processors on the right need do anything
    IF (.NOT.x_max_boundary) RETURN

    errcode = c_err_none

    DO ispecies = 1, n_species
      ALLOCATE(append_list)
      CALL create_empty_partlist(append_list)
      npart_per_cell = FLOOR(species_list(ispecies)%npart_per_cell, KIND=i8)
      npart_frac = species_list(ispecies)%npart_per_cell - npart_per_cell
      IF (npart_frac > 0) THEN
        n0 = 0
      ELSE
        n0 = 1
      ENDIF

      parameters%pack_ix = nx
      DO iz = -2, nz+3
        parameters%pack_iz = iz
        DO iy = -2, ny+3
          parameters%pack_iy = iy
          density_local = evaluate_with_parameters( &
              species_list(ispecies)%density_function, parameters, errcode)
          IF (density_local < species_list(ispecies)%&
              initial_conditions%density_min) CYCLE
          DO ipart = n0, npart_per_cell
            ! Place extra particle based on probability
            IF (ipart == 0) THEN
              IF (npart_frac < random()) CYCLE
            ENDIF
            CALL create_particle(current)
            current%part_pos(1) = x_grid_max + dx + (random() - 0.5_num) * dx
            current%part_pos(2) = y(iy) + (random() - 0.5_num) * dy
            current%part_pos(3) = z(iz) + (random() - 0.5_num) * dz

          parameters%use_grid_position=.FALSE.
          parameters%pack_pos = current%part_pos
          !Setup particle density at the exact position of the new particle
          density_local = evaluate_with_parameters( &
            species_list(ispecies)%density_function, parameters, errcode)
          IF (density_local > species_list(ispecies)%&
              initial_conditions%density_max) THEN
            density_local = species_list(ispecies)%initial_conditions%&
                density_max
          ENDIF
          weight_local = dx * dy * dz * density_local/species_list(ispecies)&
              %npart_per_cell

          IF (.NOT. species_list(ispecies)%dist_fn_set) THEN
            DO i = 1, 3
              temp_local(i) = evaluate_with_parameters( &
                  species_list(ispecies)%temperature_function(i), &
                  parameters, errcode)
              drift_local(i) = evaluate_with_parameters( &
                  species_list(ispecies)%drift_function(i), parameters, errcode)
            ENDDO
            IF (species_list(ispecies)%use_maxwell_juettner) THEN
              current%part_p = momentum_from_temperature_relativistic(&
                  mass, temp_local, &
                  species_list(ispecies)%fractional_tail_cutoff)
            ELSE
              DO i = 1, 3
                current%part_p(i) = momentum_from_temperature(&
                    mass, temp_local(i))
              ENDDO
            ENDIF
            CALL particle_drift_lorentz_transform(current, mass, drift_local)
          ENDIF
          current%weight = weight_local

#ifdef PARTICLE_DEBUG
            current%processor = rank
            current%processor_at_t0 = rank
#endif
            IF (species_list(ispecies)%dist_fn_set) THEN
              species=>species_list(ispecies)
              CALL sample_partlist_from_distfn(species, append_list)
            ENDIF
            CALL add_particle_to_partlist(append_list, current)
          ENDDO
        ENDDO
      ENDDO

      CALL append_partlist(species_list(ispecies)%attached_list, append_list)
      DEALLOCATE(append_list)
    ENDDO

  END SUBROUTINE insert_particles



  SUBROUTINE remove_particles

    TYPE(particle), POINTER :: current, next
    INTEGER :: ispecies

    ! Only processors on the left need do anything
    IF (x_min_boundary) THEN
      DO ispecies = 1, n_species
        current => species_list(ispecies)%attached_list%head
        DO WHILE(ASSOCIATED(current))
          next => current%next
          IF (current%part_pos(1) < x_min) THEN
            CALL remove_particle_from_partlist(&
                species_list(ispecies)%attached_list, current)
            DEALLOCATE(current)
          ENDIF
          current => next
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE remove_particles
#endif



  SUBROUTINE moving_window

#ifndef PER_SPECIES_WEIGHT
    REAL(num) :: window_shift_real
    INTEGER :: window_shift_cells
#else
    INTEGER :: ierr
#endif

    IF (.NOT. move_window) RETURN

#ifndef PER_SPECIES_WEIGHT
    IF (.NOT. window_started) THEN
      IF (time >= window_start_time) THEN
        bc_field(c_bd_x_min) = bc_x_min_after_move
        bc_field(c_bd_x_max) = bc_x_max_after_move
        CALL setup_particle_boundaries
        IF (.NOT.ic_from_restart) window_shift_fraction = 0.0_num
        window_started = .TRUE.
      ENDIF
    ENDIF

    ! If we have a moving window then update the window position
    IF (window_started) THEN
      window_shift_fraction = window_shift_fraction + dt * window_v_x / dx
      window_shift_cells = FLOOR(window_shift_fraction)
      ! Allow for posibility of having jumped two cells at once
      IF (window_shift_cells > 0) THEN
        window_shift_real = REAL(window_shift_cells, num)
        IF (use_offset_grid) THEN
          window_shift(1) = window_shift(1) + window_shift_real * dx
        ENDIF
        CALL shift_window(window_shift_cells)
        CALL particle_bcs
        window_shift_fraction = window_shift_fraction - window_shift_real
      ENDIF
    ENDIF
#else
    IF (rank == 0) THEN
      WRITE(*,*) 'moving windows only available when using', &
          ' per particle weighting'
    ENDIF
    CALL abort_code(c_err_pp_options_missing)
#endif

  END SUBROUTINE moving_window

END MODULE window
