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

MODULE deck_injector_block

  USE strings_advanced
  USE shunt
  USE evaluator
  USE injectors
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: injector_deck_initialise, injector_deck_finalise
  PUBLIC :: injector_block_start, injector_block_end
  PUBLIC :: injector_block_handle_element, injector_block_check

  TYPE(injector_block), POINTER :: working_injector
  LOGICAL :: boundary_set = .FALSE.
  INTEGER :: boundary

CONTAINS

  SUBROUTINE injector_deck_initialise

    NULLIFY(injector_x_min)
    NULLIFY(injector_x_max)
    NULLIFY(injector_y_min)
    NULLIFY(injector_y_max)
    NULLIFY(injector_z_min)
    NULLIFY(injector_z_max)

  END SUBROUTINE injector_deck_initialise



  SUBROUTINE injector_deck_finalise

  END SUBROUTINE injector_deck_finalise



  SUBROUTINE injector_block_start

    IF (deck_state == c_ds_first) RETURN

    ! Every new laser uses the internal time function
    ALLOCATE(working_injector)

  END SUBROUTINE injector_block_start



  SUBROUTINE injector_block_end

    REAL(num) :: first_time

    IF (deck_state == c_ds_first) RETURN

    ! To grant each custom injector a unique set of units for I/O, we assign
    ! each one an injector ID
    IF (working_injector%inject_from_file) THEN
      custom_injector_count = custom_injector_count + 1
      working_injector%custom_id = custom_injector_count
      CALL open_injector_files(working_injector)

      ! Our file injecting routines require the first injection time
      CALL read_injector_real(unit_t, first_time, working_injector)
      IF (.NOT. working_injector%file_finished) &
          working_injector%next_time = first_time
    END IF

    CALL attach_injector(working_injector)
    boundary_set = .FALSE.

  END SUBROUTINE injector_block_end



  FUNCTION injector_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, i, filename_error_ignore
    INTEGER :: io, iu
    CHARACTER(LEN=string_length) :: filename
    LOGICAL :: got_filename

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'boundary')) THEN
      ! If the boundary has already been set, simply ignore further calls to it
      IF (boundary_set) RETURN
      boundary = as_boundary_print(value, element, errcode)
      boundary_set = .TRUE.
      CALL init_injector(boundary, working_injector)
      RETURN
    END IF

    IF (.NOT. boundary_set) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Cannot set injector properties before boundary is set'
        END DO
        CALL abort_code(c_err_required_element_not_set)
      END IF
      extended_error_string = 'boundary'
      errcode = c_err_required_element_not_set
      RETURN
    END IF

    IF (str_cmp(element, 'species')) THEN
      working_injector%species = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_start')) THEN
      working_injector%t_start = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_end')) THEN
      working_injector%t_end = as_time_print(value, element, errcode)
      working_injector%has_t_end = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'use_flux_maxwellian')) THEN
      working_injector%use_flux_injector = as_logical_print(value, element, &
          errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'npart_per_cell')) THEN
      working_injector%npart_per_cell = as_integer_print(value, element, &
          errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'density_min')) THEN
      working_injector%density_min = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'density')) THEN
      CALL initialise_stack(working_injector%density_function)
      CALL tokenize(value, working_injector%density_function, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'temp_x')) THEN
      i = 1
      CALL initialise_stack(working_injector%temperature_function(i))
      CALL tokenize(value, working_injector%temperature_function(i), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'temp_y')) THEN
      i = 2
      CALL initialise_stack(working_injector%temperature_function(i))
      CALL tokenize(value, working_injector%temperature_function(i), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'temp_z')) THEN
      i = 3
      CALL initialise_stack(working_injector%temperature_function(i))
      CALL tokenize(value, working_injector%temperature_function(i), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'temp')) THEN
      DO i = 1, 3
        CALL initialise_stack(working_injector%temperature_function(i))
        CALL tokenize(value, working_injector%temperature_function(i), errcode)
      END DO
      RETURN
    END IF

    IF (str_cmp(element, 'drift_x')) THEN
      i = 1
      CALL initialise_stack(working_injector%drift_function(i))
      CALL tokenize(value, working_injector%drift_function(i), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_y')) THEN
      i = 2
      CALL initialise_stack(working_injector%drift_function(i))
      CALL tokenize(value, working_injector%drift_function(i), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'drift_z')) THEN
      i = 3
      CALL initialise_stack(working_injector%drift_function(i))
      CALL tokenize(value, working_injector%drift_function(i), errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'inject_from_file')) THEN
      working_injector%inject_from_file = as_logical_print(value, element, &
          errcode)
      RETURN
    END IF

    ! Below here only expect to be dealing with filenames
    CALL get_filename(value, filename, got_filename, filename_error_ignore)

    IF (str_cmp(element, 'x_data')) THEN
      IF (got_filename) THEN
        injector_filenames%x_data = TRIM(filename)
        working_injector%x_data_given = .TRUE.
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 'y_data')) THEN
      IF (got_filename) THEN
        injector_filenames%y_data = TRIM(filename)
        working_injector%y_data_given = .TRUE.
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 'z_data')) THEN
      IF (got_filename) THEN
        injector_filenames%z_data = TRIM(filename)
        working_injector%z_data_given = .TRUE.
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 'px_data')) THEN
      IF (got_filename) THEN
        injector_filenames%px_data = TRIM(filename)
        working_injector%px_data_given = .TRUE.
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 'py_data')) THEN
      IF (got_filename) THEN
        injector_filenames%py_data = TRIM(filename)
        working_injector%py_data_given = .TRUE.
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 'pz_data')) THEN
      IF (got_filename) THEN
        injector_filenames%pz_data = TRIM(filename)
        working_injector%pz_data_given = .TRUE.
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 't_data')) THEN
      IF (got_filename) THEN
        injector_filenames%t_data = TRIM(filename)
        working_injector%t_data_given = .TRUE.
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
    END IF

    IF (str_cmp(element, 'w_data')) THEN
#ifndef PER_SPECIES_WEIGHT
      IF (got_filename) THEN
        injector_filenames%w_data = TRIM(filename)
        working_injector%w_data_given = .TRUE.
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
#else
      IF (rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Cannot assign weights from ', filename, ' as compiler flag ', &
            '-DPER_SPECIES_WEIGHT is used.'
        PRINT*,'Code will terminate.'
      END IF
      errcode = c_err_bad_value
#endif
    END IF

    IF (str_cmp(element, 'id_data')) THEN
#if defined(PARTICLE_ID4) || defined(PARTICLE_ID)
      IF (got_filename) THEN
        injector_filenames%id_data = TRIM(filename)
        working_injector%id_data_given = .TRUE.
        RETURN
      END IF
      errcode = c_err_bad_value
      RETURN
#else
      IF (rank == 0) THEN
        PRINT*,'*** ERROR ***'
        PRINT*,'Cannot assign weights from ', filename, ' as compiler flag ', &
            '-DPER_SPECIES_WEIGHT is used.'
        PRINT*,'Code will terminate.'
      END IF
      errcode = c_err_bad_value
#endif
    END IF

    errcode = c_err_unknown_element

  END FUNCTION injector_block_handle_element



  FUNCTION injector_block_check() RESULT(errcode)

    INTEGER :: errcode
    TYPE(injector_block), POINTER :: current
    INTEGER :: error, io, iu
    LOGICAL :: custom_error

    custom_error = .FALSE.
    use_injectors = .FALSE.

    errcode = c_err_none
    error = 0
    current => injector_x_min
    DO WHILE(ASSOCIATED(current))
      IF (current%species == -1) error = IOR(error, 1)
      IF (.NOT. current%density_function%init &
          .AND. .NOT. current%inject_from_file) error = IOR(error, 2)
      IF (error == 0) use_injectors = .TRUE.

      ! Specifying momentum and ID is optional for custom injectors, but weight,
      ! injection time (and position in 2D and 3D codes) must be given.
      IF (current%inject_from_file) THEN
        IF (.NOT. current%y_data_given) custom_error = .TRUE.
        IF (.NOT. current%z_data_given) custom_error = .TRUE.
        IF (.NOT. current%t_data_given) custom_error = .TRUE.
#ifndef PER_SPECIES_WEIGHT
        IF (.NOT. current%w_data_given) custom_error = .TRUE.
#endif
      END IF
      current => current%next
    END DO

    current => injector_x_max
    DO WHILE(ASSOCIATED(current))
      IF (current%species == -1) error = IOR(error, 1)
      IF (.NOT. current%density_function%init &
          .AND. .NOT. current%inject_from_file) error = IOR(error, 2)
      IF (error == 0) use_injectors = .TRUE.

      ! Specifying momentum and ID is optional for custom injectors, but weight,
      ! injection time (and position in 2D and 3D codes) must be given.
      IF (current%inject_from_file) THEN
        IF (.NOT. current%y_data_given) custom_error = .TRUE.
        IF (.NOT. current%z_data_given) custom_error = .TRUE.
        IF (.NOT. current%t_data_given) custom_error = .TRUE.
#ifndef PER_SPECIES_WEIGHT
        IF (.NOT. current%w_data_given) custom_error = .TRUE.
#endif
      END IF
      current => current%next
    END DO

    current => injector_y_min
    DO WHILE(ASSOCIATED(current))
      IF (current%species == -1) error = IOR(error, 1)
      IF (.NOT. current%density_function%init &
          .AND. .NOT. current%inject_from_file) error = IOR(error, 2)
      IF (error == 0) use_injectors = .TRUE.

      ! Specifying momentum and ID is optional for custom injectors, but weight,
      ! injection time (and position in 2D and 3D codes) must be given.
      IF (current%inject_from_file) THEN
        IF (.NOT. current%x_data_given) custom_error = .TRUE.
        IF (.NOT. current%z_data_given) custom_error = .TRUE.
        IF (.NOT. current%t_data_given) custom_error = .TRUE.
#ifndef PER_SPECIES_WEIGHT
        IF (.NOT. current%w_data_given) custom_error = .TRUE.
#endif
      END IF
      current => current%next
    END DO

    current => injector_y_max
    DO WHILE(ASSOCIATED(current))
      IF (current%species == -1) error = IOR(error, 1)
      IF (.NOT. current%density_function%init &
          .AND. .NOT. current%inject_from_file) error = IOR(error, 2)
      IF (error == 0) use_injectors = .TRUE.

      ! Specifying momentum and ID is optional for custom injectors, but weight,
      ! injection time (and position in 2D and 3D codes) must be given.
      IF (current%inject_from_file) THEN
        IF (.NOT. current%x_data_given) custom_error = .TRUE.
        IF (.NOT. current%z_data_given) custom_error = .TRUE.
        IF (.NOT. current%t_data_given) custom_error = .TRUE.
#ifndef PER_SPECIES_WEIGHT
        IF (.NOT. current%w_data_given) custom_error = .TRUE.
#endif
      END IF
      current => current%next
    END DO

    current => injector_z_min
    DO WHILE(ASSOCIATED(current))
      IF (current%species == -1) error = IOR(error, 1)
      IF (.NOT. current%density_function%init &
          .AND. .NOT. current%inject_from_file) error = IOR(error, 2)
      IF (error == 0) use_injectors = .TRUE.

      ! Specifying momentum and ID is optional for custom injectors, but weight,
      ! injection time (and position in 2D and 3D codes) must be given.
      IF (current%inject_from_file) THEN
        IF (.NOT. current%x_data_given) custom_error = .TRUE.
        IF (.NOT. current%y_data_given) custom_error = .TRUE.
        IF (.NOT. current%t_data_given) custom_error = .TRUE.
#ifndef PER_SPECIES_WEIGHT
        IF (.NOT. current%w_data_given) custom_error = .TRUE.
#endif
      END IF
      current => current%next
    END DO

    current => injector_z_max
    DO WHILE(ASSOCIATED(current))
      IF (current%species == -1) error = IOR(error, 1)
      IF (.NOT. current%density_function%init &
          .AND. .NOT. current%inject_from_file) error = IOR(error, 2)
      IF (error == 0) use_injectors = .TRUE.

      ! Specifying momentum and ID is optional for custom injectors, but weight,
      ! injection time (and position in 2D and 3D codes) must be given.
      IF (current%inject_from_file) THEN
        IF (.NOT. current%x_data_given) custom_error = .TRUE.
        IF (.NOT. current%y_data_given) custom_error = .TRUE.
        IF (.NOT. current%t_data_given) custom_error = .TRUE.
#ifndef PER_SPECIES_WEIGHT
        IF (.NOT. current%w_data_given) custom_error = .TRUE.
#endif
      END IF
      current => current%next
    END DO

    IF (error /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          IF (IAND(error, 1) /= 1) &
              WRITE(io,*) 'Must define a species for every injector'
          IF (IAND(error, 2) /= 0) &
              WRITE(io,*) 'Must define a density for every injector'
#ifndef PER_SPECIES_WEIGHT
          IF (custom_error) WRITE(io,*) 'Must include injection time, ', &
              'position, and particle weights'
#else
          IF (custom_error) WRITE(io,*) 'Must include injection time and ', &
              'particle position'
#endif
        END DO
      END IF
      errcode = c_err_missing_elements
      RETURN
    END IF

    ! To avoid using up all the unit integers for file reading and writing,
    ! enforce an upper limit on the number of custom injectors (currently 1e6)
    IF (custom_injector_count > custom_injector_limit) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Maximum custom injector limit exceeded!'
          WRITE(io,*) 'Cannot specify more than ', custom_injector_limit, &
              ' injectors.'
        END DO
      END IF
      errcode = c_err_bad_setup
      RETURN
    END IF

  END FUNCTION injector_block_check

END MODULE deck_injector_block
