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

MODULE bremsstrahlung
#ifdef BREMSSTRAHLUNG

  USE partlist
  USE calc_df

  IMPLICIT NONE

CONTAINS

  !Function to ensure all required species are present
  FUNCTION check_bremsstrahlung_variables()
    INTEGER :: check_bremsstrahlung_variables
    INTEGER :: io, iu
    INTEGER :: ispecies
    INTEGER :: first_electron = -1

    check_bremsstrahlung_variables = c_err_none

    ! No special species required if bremsstrahlung is turned off
    IF (.NOT.use_bremsstrahlung) RETURN

    ! No special species required if we only do radiation reaction
    IF (.NOT.produce_bremsstrahlung_photons) RETURN

    ! Identify if there exists any electron species
    DO ispecies = 1, n_species
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .AND. first_electron == -1) THEN
        first_electron = ispecies
      ENDIF
    ENDDO

    ! Print warning if there is no electron species
    IF (first_electron < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No electron species specified.'
          WRITE(io,*) 'Specify using "identify:electron".'
          WRITE(io,*) 'Bremsstrahlung routines require at least one species' , &
              ' of electrons.'
        ENDDO
      ENDIF
      check_bremsstrahlung_variables = c_err_missing_elements
      RETURN
    ENDIF

#ifdef PHOTONS
    ! photon_species can act as bremsstrahlung_photon_species if no
    ! bremsstrahlung species is defined
    IF (bremsstrahlung_photon_species == -1 &
        .AND. .NOT. photon_species == -1 ) THEN
      bremsstrahlung_photon_species = photon_species
    END IF
#endif

    !Check if there exists a species to populate with bremsstrahlung photons
    IF (bremsstrahlung_photon_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No photon species specified. Specify using ', &
              '"identify:brem_photon"'
        ENDDO
      ENDIF
      check_bremsstrahlung_variables = c_err_missing_elements
      RETURN
    ENDIF

  END FUNCTION check_bremsstrahlung_variables

#endif
END MODULE bremsstrahlung
