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

MODULE deck_pic_hybrid_block

  USE strings_advanced
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pic_hybrid_deck_initialise, pic_hybrid_deck_finalise
  PUBLIC :: pic_hybrid_block_start, pic_hybrid_block_end
  PUBLIC :: pic_hybrid_block_handle_element, pic_hybrid_block_check

CONTAINS

  SUBROUTINE pic_hybrid_deck_initialise

  END SUBROUTINE pic_hybrid_deck_initialise



  SUBROUTINE pic_hybrid_deck_finalise

    IF (deck_state == c_ds_first) RETURN

  END SUBROUTINE pic_hybrid_deck_finalise



  SUBROUTINE pic_hybrid_block_start

  END SUBROUTINE pic_hybrid_block_start



  SUBROUTINE pic_hybrid_block_end

  END SUBROUTINE pic_hybrid_block_end



  FUNCTION pic_hybrid_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
    errcode = c_err_none
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'use_pic_hybrid') &
        .OR. str_cmp(element, 'use_hybrid_pic')) THEN
      use_pic_hybrid = as_logical_print(value, element, errcode)
      RETURN
    END IF

#ifdef PIC_HYBRID
    IF (str_cmp(element, 'overlap_depth') &
        .OR. str_cmp(element, 'pic_hybrid_depth')) THEN
      overlap_depth = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'use_fields') &
        .OR. str_cmp(element, 'use_pic_hybrid_fields')) THEN
      use_pic_hybrid_fields = as_logical_print(value, element, errcode)
      RETURN
    END IF
#endif

    errcode = c_err_unknown_element

  END FUNCTION pic_hybrid_block_handle_element



  FUNCTION pic_hybrid_block_check() RESULT(errcode)

    INTEGER :: errcode
    errcode = c_err_none

  END FUNCTION pic_hybrid_block_check

END MODULE deck_pic_hybrid_block
