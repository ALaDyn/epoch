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

MODULE deck_hybrid_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: hybrid_deck_initialise, hybrid_deck_finalise
  PUBLIC :: hybrid_block_start, hybrid_block_end
  PUBLIC :: hybrid_block_handle_element, hybrid_block_check

CONTAINS

  SUBROUTINE hybrid_deck_initialise

  END SUBROUTINE hybrid_deck_initialise



  SUBROUTINE hybrid_deck_finalise

    INTEGER :: io, iu
#ifdef HYBRID
    IF (deck_state == c_ds_first) RETURN
    IF (use_hybrid .AND. use_hybrid_collisions) need_random_state = .TRUE.

#else
    IF (use_hybrid) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_hybrid=T" in the', &
              ' "hybrid" block.'
          WRITE(io,*) 'Please recompile with the -DHYBRID preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif

  END SUBROUTINE hybrid_deck_finalise



  SUBROUTINE hybrid_block_start

  END SUBROUTINE hybrid_block_start



  SUBROUTINE hybrid_block_end

  END SUBROUTINE hybrid_block_end



  FUNCTION hybrid_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode, io, iu

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'use_hybrid') &
        .OR. str_cmp(element, 'hybrid')) THEN
      use_hybrid = as_logical_print(value, element, errcode)
      RETURN
    END IF

#ifdef HYBRID
    IF (str_cmp(element, 'use_hybrid_fields') &
        .OR. str_cmp(element, 'use_fields')) THEN
      use_hybrid_fields = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_hybrid_collisions') &
        .OR. str_cmp(element, 'use_collisions')) THEN
      use_hybrid_collisions = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_hybrid_scatter') &
        .OR. str_cmp(element, 'use_scatter')) THEN
      use_hybrid_scatter = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'produce_delta_rays') &
        .OR. str_cmp(element, 'produce_delta')) THEN
      produce_delta_rays = as_logical_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'temperature') &
        .OR. str_cmp(element, 'background_temp') &
        .OR. str_cmp(element, 'hybrid_temp') &
        .OR. str_cmp(element, 'start_temp')) THEN
      hybrid_Tb_init = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'min_delta_energy')) THEN
      min_delta_energy = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'min_delta_KE')) THEN
      min_delta_energy = as_real_print(value, element, errcode) + mc2
      RETURN
    END IF

    IF (str_cmp(element, 'min_hybrid_energy')) THEN
      min_hybrid_energy = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'min_hybrid_KE')) THEN
      min_hybrid_energy = as_real_print(value, element, errcode) + mc2
      RETURN
    END IF

    IF (str_cmp(element, 'elastic_model')) THEN
      IF (str_cmp(value, 'urban') .OR. str_cmp(value, 'Urban')) THEN
        elastic_scatter_model = c_hy_urban
      ELSE IF (str_cmp(value, 'davies') .OR. str_cmp(value, 'Davies')) THEN
        elastic_scatter_model = c_hy_davies
      ELSE
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Unrecognised elastic scatter model: ', value
            WRITE(io,*) 'This will be ignored.'
          END DO
        END IF
      END IF
      RETURN
    END IF

    errcode = c_err_unknown_element
#endif

  END FUNCTION hybrid_block_handle_element



  FUNCTION hybrid_block_check() RESULT(errcode)

    INTEGER :: errcode
#ifdef HYBRID
    INTEGER :: io, iu
#endif

    errcode = c_err_none

#ifdef HYBRID
    IF (hybrid_Tb_init < 0.0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Please set the background temperature to a positive &
              value.'
          WRITE(io,*) 'Code will terminate.'
        END DO
      END IF

      errcode = c_err_bad_value + c_err_terminate
    END IF

    ! Delta-ray emission with kinetic energy below 1 keV is treated as
    ! continuous energy loss
    IF (produce_delta_rays .AND. min_delta_energy < 1.0e3_num * q0 + mc2 ) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Delta-rays with under 1 keV kinetic energy are ', &
              'treated as a continuous energy loss.'
          WRITE(io,*) 'This code will not add delta rays below 1 keV ', &
              'kinetic energy to the simulation.'
        END DO
      END IF
      min_delta_energy = 1.0e3*q0 + mc2
    END IF

    ! We shouldn't be able to generate delta-rays with lower energy than the
    ! minimum allowed hybrid energy
    IF (produce_delta_rays .AND. min_delta_energy < min_hybrid_energy) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'min_delta_energy < min_hybrid_energy, meaning some', &
              'delta-rays would be added to the simulation, and immediately ', &
              'removed. This code will run with: '
          WRITE(io,*) 'min_delta_energy = ', min_hybrid_energy/q0, ' eV'
          WRITE(io,*) 'min_hybrid_energy = ', min_hybrid_energy/q0, 'eV'
        END DO
      END IF
      min_delta_energy = min_hybrid_energy
    END IF
#endif

  END FUNCTION hybrid_block_check

END MODULE deck_hybrid_block
