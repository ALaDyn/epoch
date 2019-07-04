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

MODULE deck_landau_lifshitz_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: landau_lifshitz_deck_initialise, landau_lifshitz_deck_finalise
  PUBLIC :: landau_lifshitz_block_start, landau_lifshitz_block_end
  PUBLIC :: landau_lifshitz_block_handle_element

CONTAINS

  SUBROUTINE landau_lifshitz_deck_initialise

  END SUBROUTINE landau_lifshitz_deck_initialise



  SUBROUTINE landau_lifshitz_deck_finalise

    INTEGER :: io, iu

#ifdef LANDAU_LIFSHITZ

    IF (deck_state == c_ds_first) RETURN

    IF (rank == 0) THEN
#ifdef PHOTONS
      IF (use_qed .AND. use_landau_lifshitz) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to run with both qed and and Landau Lifshitz', &
              ' simultaneously.'
        END DO
        CALL abort_code(c_err_io_error)
      END IF
#endif
    END IF

#else
    IF (use_landau_lifshitz) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Unable to set "use_landau_lifshitz=T" in the', &
              ' "landau_lifshitz" block.'
          WRITE(io,*) 'Please recompile with the -DLANDAU_LIFSHITZ', &
              'preprocessor flag.'
        END DO
      END IF
      CALL abort_code(c_err_pp_options_missing)
    END IF
#endif

END SUBROUTINE landau_lifshitz_deck_finalise



  SUBROUTINE landau_lifshitz_block_start

  END SUBROUTINE landau_lifshitz_block_start



  SUBROUTINE landau_lifshitz_block_end

  END SUBROUTINE landau_lifshitz_block_end



  FUNCTION landau_lifshitz_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode

    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'use_landau_lifshitz') &
        .OR. str_cmp(element, 'landau_lifshitz')) THEN
      use_landau_lifshitz = as_logical_print(value, element, errcode)
      RETURN
    ENDIF

#ifdef LANDAU_LIFSHITZ
    IF (str_cmp(element, 'landau_lifshitz_start_time') &
        .OR. str_cmp(element, 'start_time')) THEN
      landau_lifshitz_start_time = as_real_print(value, element, errcode)
      RETURN
    ENDIF

    errcode = c_err_unknown_element
#endif

  END FUNCTION landau_lifshitz_block_handle_element

END MODULE deck_landau_lifshitz_block
