! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2010 Chris Brady <C.S.Brady@warwick.ac.uk>
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

MODULE deck_hy_laser_block

  USE strings_advanced
  USE hy_laser
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: hy_laser_deck_initialise, hy_laser_deck_finalise
  PUBLIC :: hy_laser_block_start, hy_laser_block_end
  PUBLIC :: hy_laser_block_handle_element, hy_laser_block_check

  TYPE(hy_laser_block), POINTER :: working_laser
  LOGICAL :: hy_boundary_set = .FALSE.
  INTEGER :: hy_boundary

CONTAINS

  SUBROUTINE hy_laser_deck_initialise

  END SUBROUTINE hy_laser_deck_initialise



  SUBROUTINE hy_laser_deck_finalise

#ifndef HYBRID
    IF (use_hybrid) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to add hy_laser!'
          WRITE(io,*) 'Please recompile with the -DHYBRID preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif

  END SUBROUTINE hy_laser_deck_finalise



  SUBROUTINE hy_laser_block_start

#ifdef HYBRID
    IF (deck_state == c_ds_first) RETURN

    ! Every new laser uses the internal time function
    ALLOCATE(working_laser)
    working_laser%use_time_function = .FALSE.
    working_laser%use_profile_function = .FALSE.
    working_laser%use_omega_function = .FALSE.
#endif

  END SUBROUTINE hy_laser_block_start



  SUBROUTINE hy_laser_block_end

#ifdef HYBRID
    IF (deck_state == c_ds_first) RETURN

    CALL attach_hy_laser(working_laser)
    hy_boundary_set = .FALSE.
#endif

  END SUBROUTINE hy_laser_block_end



  FUNCTION hy_laser_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    REAL(num) :: dummy
    INTEGER :: io, iu

#ifdef HYBRID
    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'boundary')) THEN
      ! If the hy_boundary has already been set, simply ignore further calls to it
      IF (hy_boundary_set) RETURN
      hy_boundary = as_boundary_print(value, element, errcode)
      hy_boundary_set = .TRUE.
      CALL init_hy_laser(hy_boundary, working_laser)
      RETURN
    END IF

    IF (.NOT. hy_boundary_set) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Cannot set hy_laser properties before boundary is set'
        END DO
        CALL abort_code(c_err_required_element_not_set)
      END IF
      extended_error_string = 'hy_boundary'
      errcode = c_err_required_element_not_set
      RETURN
    END IF

    IF (str_cmp(element, 'ppc')) THEN
      working_laser%ppc = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'profile')) THEN
      CALL initialise_stack(working_laser%profile_function)
      CALL tokenize(value, working_laser%profile_function, errcode)
      working_laser%profile = 0.0_num
      working_laser%use_profile_function = .TRUE.
      CALL hy_laser_update_profile(working_laser)
      IF (.NOT. working_laser%profile_function%is_time_varying) THEN
        CALL deallocate_stack(working_laser%profile_function)
        working_laser%use_profile_function = .FALSE.
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'omega')) THEN
      CALL initialise_stack(working_laser%omega_function)
      CALL tokenize(value, working_laser%omega_function, errcode)
      working_laser%omega = 0.0_num
      working_laser%omega_func_type = c_of_omega
      working_laser%use_omega_function = .TRUE.
      CALL hy_laser_update_omega(working_laser)
      IF (.NOT. working_laser%omega_function%is_time_varying) THEN
        CALL deallocate_stack(working_laser%omega_function)
        working_laser%use_omega_function = .FALSE.
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'frequency')) THEN
      CALL initialise_stack(working_laser%omega_function)
      CALL tokenize(value, working_laser%omega_function, errcode)
      working_laser%omega = 0.0_num
      working_laser%omega_func_type = c_of_freq
      working_laser%use_omega_function = .TRUE.
      CALL hy_laser_update_omega(working_laser)
      IF (.NOT. working_laser%omega_function%is_time_varying) THEN
        CALL deallocate_stack(working_laser%omega_function)
        working_laser%use_omega_function = .FALSE.
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'lambda') .OR. str_cmp(element, 'wavelength')) THEN
      CALL initialise_stack(working_laser%omega_function)
      CALL tokenize(value, working_laser%omega_function, errcode)
      working_laser%omega = 0.0_num
      working_laser%omega_func_type = c_of_lambda
      working_laser%use_omega_function = .TRUE.
      CALL hy_laser_update_omega(working_laser)
      IF (.NOT. working_laser%omega_function%is_time_varying) THEN
        CALL deallocate_stack(working_laser%omega_function)
        working_laser%use_omega_function = .FALSE.
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 't_profile')) THEN
      working_laser%use_time_function = .TRUE.
      CALL initialise_stack(working_laser%time_function)
      CALL tokenize(value, working_laser%time_function, errcode)
      ! evaluate it once to check that it's a valid block
      dummy = evaluate(working_laser%time_function, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'id')) THEN
      working_laser%id = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'intensity')) THEN
      ! [W/m²]
      working_laser%intensity = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_start')) THEN
      working_laser%t_start = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't_end')) THEN
      working_laser%t_end = as_time_print(value, element, errcode)
      working_laser%has_t_end = .TRUE.
      RETURN
    END IF

    IF (str_cmp(element, 'profile_min')) THEN
      working_laser%profile_min = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'efficiency')) THEN
      working_laser%efficiency = as_time_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'species')) THEN
      working_laser%species = as_integer_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'mean_energy')) THEN
      IF (str_cmp(value, 'a0')) THEN
        working_laser%mean = c_mean_a0
      ELSE IF (str_cmp(value, 'Wilks') .OR. str_cmp(value, 'wilks')) THEN
        working_laser%mean = c_mean_wilks
      ELSE
        extended_error_string = 'Unrecognised mean energy model'
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'energy_dist')) THEN
      IF (str_cmp(value, 'exp')) THEN
        working_laser%e_dist = e_dist_exp
      ELSE IF (str_cmp(value, 'mono')) THEN
        working_laser%e_dist = e_dist_mono
      ELSE IF (str_cmp(value, 'tophat') .OR. str_cmp(value, 'top_hat')) THEN
        working_laser%e_dist = e_dist_tophat
      ELSE IF (str_cmp(value, 'exp_weight')) THEN
        working_laser%e_dist = e_dist_exp_weight
      ELSE
        extended_error_string = 'Unrecognised energy distribution model'
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'angular_dist')) THEN
      IF (str_cmp(value, 'uniform')) THEN
        working_laser%ang_dist = c_ang_uniform
      ELSE IF (str_cmp(value, 'cos')) THEN
        working_laser%ang_dist = c_ang_cos
      ELSE IF (str_cmp(value, 'beam')) THEN
        working_laser%ang_dist = c_ang_beam
      ELSE
        extended_error_string = 'Unrecognised angular distribution model'
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'theta_max')) THEN
      working_laser%user_theta_max = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'cos_n_power')) THEN
      working_laser%cos_n_power = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'top_hat_L')) THEN
      working_laser%top_hat_L = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'mean_mult')) THEN
      working_laser%mean_mult = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_sheng_dir')) THEN
      working_laser%use_sheng_dir = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sheng_angle')) THEN
      working_laser%sheng_angle = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_moore_max')) THEN
      working_laser%use_moore_max = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'theta')) THEN
      working_laser%theta_mean = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'phi')) THEN
      working_laser%phi_mean = as_real_print(value, element, errcode)
      RETURN
    END IF
#endif

    errcode = c_err_unknown_element

  END FUNCTION hy_laser_block_handle_element



  FUNCTION hy_laser_block_check() RESULT(errcode)

    INTEGER :: errcode
    TYPE(hy_laser_block), POINTER :: current
    INTEGER :: error, io, iu

    errcode = c_err_none

#ifdef HYBRID
    error = 0
    current => hy_laser_x_min
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%intensity < 0.0_num) error = IOR(error, 2)
      IF (current%ppc < 0) error = IOR(error, 4)
      IF (current%species < 0) error = IOR(error, 8)
      IF (current%efficiency < 0.0_num) error = IOR(error, 16)
      IF (current%mean < 0) error = IOR(error, 32)
      IF (current%e_dist < 0) error = IOR(error, 64)
      IF (current%ang_dist < 0) error = IOR(error, 128)
      IF (current%efficiency > 1.0_num) error = IOR(error, 256)
      IF (ABS(current%top_hat_L-0.5_num) >= 0.5_num) error = IOR(error, 512)
      IF (ABS(current%sheng_angle) > 0.5_num*pi) error = IOR(error, 1024)
      IF (current%mean_mult < 0.0_num) error = IOR(error, 2048)
      current => current%next
    END DO

    current => hy_laser_x_max
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%intensity < 0.0_num) error = IOR(error, 2)
      IF (current%ppc < 0) error = IOR(error, 4)
      IF (current%species < 0) error = IOR(error, 8)
      IF (current%efficiency < 0.0_num) error = IOR(error, 16)
      IF (current%mean < 0) error = IOR(error, 32)
      IF (current%e_dist < 0) error = IOR(error, 64)
      IF (current%ang_dist < 0) error = IOR(error, 128)
      IF (current%efficiency > 1.0_num) error = IOR(error, 256)
      IF (ABS(current%top_hat_L-0.5_num) >= 0.5_num) error = IOR(error, 512)
      IF (ABS(current%sheng_angle) > 0.5_num*pi) error = IOR(error, 1024)
      IF (current%mean_mult < 0.0_num) error = IOR(error, 2048)
      current => current%next
    END DO

    current => hy_laser_y_min
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%intensity < 0.0_num) error = IOR(error, 2)
      IF (current%ppc < 0) error = IOR(error, 4)
      IF (current%species < 0) error = IOR(error, 8)
      IF (current%efficiency < 0.0_num) error = IOR(error, 16)
      IF (current%mean < 0) error = IOR(error, 32)
      IF (current%e_dist < 0) error = IOR(error, 64)
      IF (current%ang_dist < 0) error = IOR(error, 128)
      IF (current%efficiency > 1.0_num) error = IOR(error, 256)
      IF (ABS(current%top_hat_L-0.5_num) >= 0.5_num) error = IOR(error, 512)
      IF (ABS(current%sheng_angle) > 0.5_num*pi) error = IOR(error, 1024)
      IF (current%mean_mult < 0.0_num) error = IOR(error, 2048)
      current => current%next
    END DO

    current => hy_laser_y_max
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%intensity < 0.0_num) error = IOR(error, 2)
      IF (current%ppc < 0) error = IOR(error, 4)
      IF (current%species < 0) error = IOR(error, 8)
      IF (current%efficiency < 0.0_num) error = IOR(error, 16)
      IF (current%mean < 0) error = IOR(error, 32)
      IF (current%e_dist < 0) error = IOR(error, 64)
      IF (current%ang_dist < 0) error = IOR(error, 128)
      IF (current%efficiency > 1.0_num) error = IOR(error, 256)
      IF (ABS(current%top_hat_L-0.5_num) >= 0.5_num) error = IOR(error, 512)
      IF (ABS(current%sheng_angle) > 0.5_num*pi) error = IOR(error, 1024)
      IF (current%mean_mult < 0.0_num) error = IOR(error, 2048)
      current => current%next
    END DO

    current => hy_laser_z_min
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%intensity < 0.0_num) error = IOR(error, 2)
      IF (current%ppc < 0) error = IOR(error, 4)
      IF (current%species < 0) error = IOR(error, 8)
      IF (current%efficiency < 0.0_num) error = IOR(error, 16)
      IF (current%mean < 0) error = IOR(error, 32)
      IF (current%e_dist < 0) error = IOR(error, 64)
      IF (current%ang_dist < 0) error = IOR(error, 128)
      IF (current%efficiency > 1.0_num) error = IOR(error, 256)
      IF (ABS(current%top_hat_L-0.5_num) >= 0.5_num) error = IOR(error, 512)
      IF (ABS(current%sheng_angle) > 0.5_num*pi) error = IOR(error, 1024)
      IF (current%mean_mult < 0.0_num) error = IOR(error, 2048)
      current => current%next
    END DO

    current => hy_laser_z_max
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%intensity < 0.0_num) error = IOR(error, 2)
      IF (current%ppc < 0) error = IOR(error, 4)
      IF (current%species < 0) error = IOR(error, 8)
      IF (current%efficiency < 0.0_num) error = IOR(error, 16)
      IF (current%mean < 0) error = IOR(error, 32)
      IF (current%e_dist < 0) error = IOR(error, 64)
      IF (current%ang_dist < 0) error = IOR(error, 128)
      IF (current%efficiency > 1.0_num) error = IOR(error, 256)
      IF (ABS(current%top_hat_L-0.5_num) >= 0.5_num) error = IOR(error, 512)
      IF (ABS(current%sheng_angle) > 0.5_num*pi) error = IOR(error, 1024)
      IF (current%mean_mult < 0.0_num) error = IOR(error, 2048)
      current => current%next
    END DO

    IF (IAND(error, 1) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "lambda" or "omega" for every hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define an "intensity" for every hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 4) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "ppc" for every hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 8) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "species" for every hy_laser to inject to.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 16) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define an "efficiency" (laser energy to electron', &
              ' energy) for every'
          WRITE(io,*)'hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 32) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must set a model for mean electron energy with', &
              ' "mean_energy" for every'
          WRITE(io,*) 'hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 64) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must set a model for the energy distribution with', &
              ' "energy_dist" for every'
          WRITE(io,*) 'hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 128) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must set a model for the angular distribution with', &
              ' "angular_dist" for every'
          WRITE(io,*) 'hy_laser.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 256) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Efficiency cannot exceed 1.0 (100% conversion from', &
              ' laser energy to e- energy).'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF

    IF (IAND(error, 512) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Top-hat fractional-width L must satisfy 0 < L < 1.'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF

    IF (IAND(error, 1024) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Incident (Sheng) angle cannot exceed pi/2 radians.'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF

    IF (IAND(error, 2048) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'mean_mult must be greater than zero.'
        END DO
      END IF
      errcode = c_err_bad_value
    END IF
#endif

  END FUNCTION hy_laser_block_check

END MODULE deck_hy_laser_block
