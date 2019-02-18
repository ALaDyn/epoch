MODULE return_boundary

  USE shared_data
  USE injectors

  IMPLICIT NONE

  PUBLIC :: setup_return_boundaries, create_return_injector, &
      update_return_injector

CONTAINS


  SUBROUTINE setup_return_boundaries

    INTEGER :: i
    TYPE(particle_species), POINTER :: current

    DO i = 1, n_species
      current => species_list(i)
      IF (ANY(current%bc_particle == c_bc_return)) THEN
        !Create and setup injector for "return current"
        CALL create_return_injector(current%bc_particle, i)
      END IF

    END DO

  END SUBROUTINE setup_return_boundaries


  !Create an injector for a return-boundary species
  SUBROUTINE create_return_injector(bcs, ispecies)

    INTEGER, INTENT(IN) :: ispecies
    INTEGER, DIMENSION(c_ndims*2), INTENT(IN) :: bcs
    TYPE(injector_block), POINTER :: working_injector
    INTEGER :: i, j
    TYPE(parameter_pack) :: parameters

    DO i = 1, c_ndims*2
      IF (bcs(i) /= c_bc_return) CYCLE
      use_injectors = .TRUE.
      need_random_state = .TRUE.
      ALLOCATE(working_injector)
      CALL init_injector(i,  working_injector)
      working_injector%species = ispecies
      !TODO makes sense to pass species here, rather than n_species
      working_injector%npart_per_cell = 1 !Tenporary, fixed after load
      IF (i == 1) THEN
        species_list(ispecies)%injector_x_min => working_injector
        parameters%pack_ix = 1
      ELSE IF (i == 2) THEN
        species_list(ispecies)%injector_x_max => working_injector
        parameters%pack_ix = nx
      ENDIF

      DO j = 1, 3
        working_injector%temperature(j) = evaluate_with_parameters( &
            species_list(ispecies)%temperature_function(j), parameters, errcode)
        working_injector%drift(j) = evaluate_with_parameters( &
            species_list(ispecies)%drift_function(j), parameters, errcode)
      END DO
      working_injector%density = evaluate_with_parameters( &
          species_list(ispecies)%density_function, parameters, errcode)

      CALL attach_injector(working_injector)

    END DO

  END SUBROUTINE

  !Fix up the ppc now we've loaded the particles
  SUBROUTINE finish_return_boundaries

    TYPE(injector_block), POINTER :: working_injector
    INTEGER :: i

    DO i = 1, n_species
      IF (ANY(species_list(i)%bc_particle == c_bc_return)) THEN
        working_injector => species_list(i)%injector_x_min
        working_injector%npart_per_cell = &
            FLOOR(species_list(i)%npart_per_cell)
        CALL update_return_injector(working_injector)

        working_injector => species_list(i)%injector_x_max
        working_injector%npart_per_cell = &
            FLOOR(species_list(i)%npart_per_cell)
        CALL update_return_injector(working_injector)
      END IF
    END DO

  END SUBROUTINE

  !Update temperature and drift for a return-injector species
  SUBROUTINE update_return_injector(injector)

    TYPE(injector_block), INTENT(IN), POINTER :: injector
    TYPE(particle_species), POINTER :: species
    !TODO Need to verfy this really is a return injector?

    species => species_list(injector%species)

    IF (injector%boundary == c_bd_x_min) THEN
      injector%temperature = species%ext_temp_x_min
      injector%drift(1) = species%ext_drift_x_min
    ELSE IF (injector%boundary == c_bd_x_max) THEN
      injector%temperature = species%ext_temp_x_max
      injector%drift(1) = species%ext_drift_x_max
    END IF
    CALL update_dt_inject(injector)

  END SUBROUTINE


END MODULE
