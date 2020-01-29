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

MODULE deck_boundaries_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: boundary_deck_initialise, boundary_deck_finalise
  PUBLIC :: boundary_block_start, boundary_block_end
  PUBLIC :: boundary_block_handle_element, boundary_block_check

  INTEGER, PARAMETER :: boundary_block_nbase = 2 * c_ndims
  INTEGER, PARAMETER :: boundary_block_elements = 3 * boundary_block_nbase + 4
  LOGICAL, DIMENSION(boundary_block_elements) :: boundary_block_done
  CHARACTER(LEN=string_length), DIMENSION(boundary_block_elements) :: &
      boundary_block_name = (/ &
          'bc_x_min         ', &
          'bc_x_max         ', &
          'bc_y_min         ', &
          'bc_y_max         ', &
          'bc_z_min         ', &
          'bc_z_max         ', &
          'bc_x_min_field   ', &
          'bc_x_max_field   ', &
          'bc_y_min_field   ', &
          'bc_y_max_field   ', &
          'bc_z_min_field   ', &
          'bc_z_max_field   ', &
          'bc_x_min_particle', &
          'bc_x_max_particle', &
          'bc_y_min_particle', &
          'bc_y_max_particle', &
          'bc_z_min_particle', &
          'bc_z_max_particle', &
          'cpml_thickness   ', &
          'cpml_kappa_max   ', &
          'cpml_a_max       ', &
          'cpml_sigma_max   ' /)
  CHARACTER(LEN=string_length), DIMENSION(boundary_block_elements) :: &
      alternate_name = (/ &
          'xbc_left          ', &
          'xbc_right         ', &
          'ybc_down          ', &
          'ybc_up            ', &
          'zbc_back          ', &
          'zbc_front         ', &
          'xbc_left_field    ', &
          'xbc_right_field   ', &
          'ybc_down_field    ', &
          'ybc_up_field      ', &
          'zbc_back_field    ', &
          'zbc_front_field   ', &
          'xbc_left_particle ', &
          'xbc_right_particle', &
          'ybc_down_particle ', &
          'ybc_up_particle   ', &
          'zbc_back_particle ', &
          'zbc_front_particle', &
          'cpml_thickness    ', &
          'cpml_kappa_max    ', &
          'cpml_a_max        ', &
          'cpml_sigma_max    ' /)

CONTAINS

  SUBROUTINE boundary_deck_initialise

    IF (deck_state /= c_ds_first) RETURN
    boundary_block_done = .FALSE.

  END SUBROUTINE boundary_deck_initialise



  SUBROUTINE boundary_deck_finalise

  END SUBROUTINE boundary_deck_finalise



  SUBROUTINE boundary_block_start

  END SUBROUTINE boundary_block_start



  SUBROUTINE boundary_block_end

  END SUBROUTINE boundary_block_end



  FUNCTION boundary_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    INTEGER :: loop, elementselected, itmp
    INTEGER, PARAMETER :: nbase = boundary_block_nbase
    CHARACTER(LEN=string_length) :: filename
    LOGICAL :: got_filename
    INTEGER :: filename_error_ignore

    errcode = c_err_none
    IF (deck_state /= c_ds_first) RETURN

    ! Generate tables for TNSA boundaries if appropriate

    IF (str_cmp(element, 'tnsa_escape')) THEN
      CALL get_filename(value, filename, got_filename, filename_error_ignore)
      IF (got_filename) THEN
        CALL make_table(filename, c_tnsa_escape)
      ELSE
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    IF (str_cmp(element, 'tnsa_reflect')) THEN
      CALL get_filename(value, filename, got_filename, filename_error_ignore)
      IF (got_filename) THEN
        CALL make_table(filename, c_tnsa_reflect)
      ELSE
        errcode = c_err_bad_value
      END IF
      RETURN
    END IF

    errcode = c_err_unknown_element

    ! Assign boundary types

    elementselected = 0

    DO loop = 1, boundary_block_elements
      IF (str_cmp(element, TRIM(ADJUSTL(boundary_block_name(loop)))) &
          .OR. str_cmp(element, TRIM(ADJUSTL(alternate_name(loop))))) THEN
        elementselected = loop
        EXIT
      END IF
    END DO

    IF (elementselected == 0) RETURN

    IF (boundary_block_done(elementselected)) THEN
      errcode = c_err_preset_element
      RETURN
    END IF
    boundary_block_done(elementselected) = .TRUE.
    errcode = c_err_none

    IF (elementselected <= nbase) THEN
      boundary_block_done(elementselected+  nbase) = .TRUE.
      boundary_block_done(elementselected+2*nbase) = .TRUE.
    END IF

    SELECT CASE (elementselected)
    CASE(1)
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_x_min) = itmp
      bc_particle(c_bd_x_min) = itmp
    CASE(2)
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_x_max) = itmp
      bc_particle(c_bd_x_max) = itmp
    CASE(3)
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_y_min) = itmp
      bc_particle(c_bd_y_min) = itmp
    CASE(4)
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_y_max) = itmp
      bc_particle(c_bd_y_max) = itmp
    CASE(5)
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_z_min) = itmp
      bc_particle(c_bd_z_min) = itmp
    CASE(6)
      itmp = as_bc_print(value, element, errcode)
      bc_field(c_bd_z_max) = itmp
      bc_particle(c_bd_z_max) = itmp
    CASE(nbase+1)
      bc_field(c_bd_x_min) = as_bc_print(value, element, errcode)
      boundary_block_done(1)  = .TRUE.
    CASE(nbase+2)
      bc_field(c_bd_x_max) = as_bc_print(value, element, errcode)
      boundary_block_done(2)  = .TRUE.
    CASE(nbase+3)
      bc_field(c_bd_y_min) = as_bc_print(value, element, errcode)
      boundary_block_done(3)  = .TRUE.
    CASE(nbase+4)
      bc_field(c_bd_y_max) = as_bc_print(value, element, errcode)
      boundary_block_done(4)  = .TRUE.
    CASE(nbase+5)
      bc_field(c_bd_z_min) = as_bc_print(value, element, errcode)
      boundary_block_done(5)  = .TRUE.
    CASE(nbase+6)
      bc_field(c_bd_z_max) = as_bc_print(value, element, errcode)
      boundary_block_done(6)  = .TRUE.
    CASE(2*nbase+1)
      bc_particle(c_bd_x_min) = as_bc_print(value, element, errcode)
      boundary_block_done(1)  = .TRUE.
    CASE(2*nbase+2)
      bc_particle(c_bd_x_max) = as_bc_print(value, element, errcode)
      boundary_block_done(2)  = .TRUE.
    CASE(2*nbase+3)
      bc_particle(c_bd_y_min) = as_bc_print(value, element, errcode)
      boundary_block_done(3)  = .TRUE.
    CASE(2*nbase+4)
      bc_particle(c_bd_y_max) = as_bc_print(value, element, errcode)
      boundary_block_done(4)  = .TRUE.
    CASE(2*nbase+5)
      bc_particle(c_bd_z_min) = as_bc_print(value, element, errcode)
      boundary_block_done(5)  = .TRUE.
    CASE(2*nbase+6)
      bc_particle(c_bd_z_max) = as_bc_print(value, element, errcode)
      boundary_block_done(6)  = .TRUE.
    CASE(3*nbase+1)
      cpml_thickness = as_integer_print(value, element, errcode)
    CASE(3*nbase+2)
      cpml_kappa_max = as_real_print(value, element, errcode)
    CASE(3*nbase+3)
      cpml_a_max = as_real_print(value, element, errcode)
    CASE(3*nbase+4)
      cpml_sigma_max = as_real_print(value, element, errcode)
    END SELECT

  END FUNCTION boundary_block_handle_element



  FUNCTION boundary_block_check() RESULT(errcode)

    INTEGER :: errcode
    INTEGER :: idx, io, iu
    INTEGER, PARAMETER :: nbase = boundary_block_nbase
    LOGICAL :: error

    errcode = c_err_none

    DO idx = 1, nbase
      IF (.NOT.boundary_block_done(idx+nbase) &
          .AND. .NOT.boundary_block_done(idx+2*nbase)) THEN
        boundary_block_done(idx+  nbase) = .TRUE.
        boundary_block_done(idx+2*nbase) = .TRUE.
      END IF
    END DO

    DO idx = 1, boundary_block_elements - 4
      IF (.NOT. boundary_block_done(idx)) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Required boundary block element "' &
                // TRIM(ADJUSTL(boundary_block_name(idx))) // '" absent.'
            WRITE(io,*) 'Please create this entry in the input deck'
          END DO
        END IF
        errcode = c_err_missing_elements
      END IF
    END DO

    ! TNSA boundaries only apply to particles, replace field TNSA with reflect
    DO idx = 1, 2*c_ndims
      IF (bc_field(idx) == c_bc_tnsa) bc_field(idx) = c_bc_reflect
    END DO

    ! Sanity check on periodic boundaries
    error = .FALSE.
    DO idx = 1, c_ndims
      IF (bc_field(2*idx-1) == c_bc_periodic &
          .AND. bc_field(2*idx) /= c_bc_periodic) error = .TRUE.
      IF (bc_field(2*idx-1) /= c_bc_periodic &
          .AND. bc_field(2*idx) == c_bc_periodic) error = .TRUE.
      IF (bc_particle(2*idx-1) == c_bc_periodic &
          .AND. bc_particle(2*idx) /= c_bc_periodic) error = .TRUE.
      IF (bc_particle(2*idx-1) /= c_bc_periodic &
          .AND. bc_particle(2*idx) == c_bc_periodic) error = .TRUE.
    END DO

    IF (error) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Periodic boundaries must be specified on both sides', &
              ' of the domain.'
        END DO
      END IF
      CALL abort_code(c_err_bad_value)
    END IF

  END FUNCTION boundary_block_check



  ! Allocates the global arrays tnsa_escape and tnsa_reflect
  SUBROUTINE make_table(file_name, table_id)

    CHARACTER(LEN=string_length) :: file_name
    INTEGER :: table_id
    INTEGER :: i_err, eof, line_count, iline, io, iu
    REAL(num) :: temp(2)
    LOGICAL :: pre_allocated

    ! Open the file
    OPEN(unit=lu, file=TRIM(file_name), STATUS='OLD')

    ! First pass through the file to find the number of data values
    eof = 0
    line_count = 0
    DO WHILE (eof == 0)
      READ(lu, *, IOSTAT = eof) temp(1:2)
      ! Successfully read line
      IF (eof == 0) THEN
        line_count = line_count + 1
      ! Line contains illegal values
      ELSE IF (eof > 1) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Illegal value present in: ', file_name
            WRITE(io,*) 'Code will continue to run but behaviour is unknown'
          END DO
        END IF
        eof = 0
      END IF
    END DO

    ! Allocate tables
    pre_allocated = .TRUE.

    ! Escape table
    IF (table_id == c_tnsa_escape) THEN
      IF (.NOT. ALLOCATED(tnsa_escape)) THEN
        pre_allocated = .FALSE.
        ALLOCATE(tnsa_escape(2,line_count))
      END IF
    ! Reflect table
    ELSE IF (table_id == c_tnsa_reflect) THEN
      IF (.NOT. ALLOCATED(tnsa_reflect)) THEN
        pre_allocated = .FALSE.
        ALLOCATE(tnsa_reflect(2,line_count))
      END IF
    END IF

    ! Has the user already written this table?
    IF (pre_allocated) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Multiple attempts to write TNSA tables will be ignored'
        END DO
      END IF
    END IF

    ! Prepare file for second pass
    CLOSE(UNIT=lu)
    OPEN(unit=lu, file=TRIM(file_name), STATUS='OLD')

    ! Write table
    DO iline = 1, line_count
      IF (table_id==c_tnsa_escape) READ(lu,*,IOSTAT=eof) tnsa_escape(:,iline)
      IF (table_id==c_tnsa_reflect) READ(lu,*,IOSTAT=eof) tnsa_reflect(:,iline)
    END DO

    ! Close file
    CLOSE(UNIT=lu)

  END SUBROUTINE make_table

END MODULE deck_boundaries_block
