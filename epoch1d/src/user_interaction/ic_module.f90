MODULE ic_module

  USE shared_data
  USE helper
  USE particles

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: manual_load

CONTAINS

#ifdef DELTAF_METHOD
  ! Find local particle temperature: at the moment, just copied and
  ! pasted from particle_temperature.F90

  SUBROUTINE params_local(current, temperature, drift, temp_local, drift_local)

    TYPE(particle), POINTER, INTENT(IN) :: current
    REAL(num), DIMENSION(-2:), INTENT(IN) :: temperature, drift
    REAL(num), INTENT(INOUT) :: temp_local, drift_local
    INTEGER :: ix, iy, iz

#include "particle_head.inc"

    ! Assume that temperature is cell centred
#include "particle_to_grid.inc"

    temp_local = 0.0_num
    drift_local = 0.0_num
    DO ix = sf_min, sf_max
       temp_local = temp_local + gx(ix) &
            * temperature(cell_x+ix)
       drift_local = drift_local + gx(ix) &
            * drift(cell_x+ix)
    ENDDO

  END SUBROUTINE params_local
#endif



  SUBROUTINE manual_load

#ifdef DELTAF_METHOD
    REAL(num) :: Tx, Ty, Tz, driftx, drifty, driftz
    REAL(num) :: f0_exponent, distribution, mass, npart_per_cell, idx
    REAL(num) :: two_kb_mass, two_pi_kb_mass3, part_weight
    REAL(num), PARAMETER :: two_kb = 2.0_num * kb
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current
    TYPE(particle_species), POINTER :: species
    INTEGER :: ipart, ispecies
#if DELTAF_DEBUG
    REAL(num) :: weight_back, f0_back
#endif

    ! f0 calculation: mainly, we need to calculate the phase space volumes.
    ! Calculate this based on the loading parameters. Easy to check
    ! that this is OK for a Maxwellian load by setting f0 = f0_back,
    ! and making sure the weights cancel.

    idx = 1.0_num / dx

    DO ispecies = 1, n_species
      species => species_list(ispecies)
      partlist => species%attached_list
      current => partlist%head

      mass = species%mass
      two_kb_mass = two_kb * mass
      two_pi_kb_mass3 = (pi * two_kb_mass)**3
      part_weight = species_list(ispecies)%weight

      ipart = 0
      DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
        mass = current%mass
        two_kb_mass = two_kb * mass
        two_pi_kb_mass3 = (pi * two_kb_mass)**3
#endif
        CALL params_local(current, initial_conditions(ispecies)%temp(:,1), &
            initial_conditions(ispecies)%drift(:,1), Tx, driftx)
        CALL params_local(current, initial_conditions(ispecies)%temp(:,2), &
            initial_conditions(ispecies)%drift(:,2), Ty, drifty)
        CALL params_local(current, initial_conditions(ispecies)%temp(:,3), &
            initial_conditions(ispecies)%drift(:,3), Tz, driftz)

        f0_exponent = ((current%part_p(1) - driftx)**2 / Tx &
                     + (current%part_p(2) - drifty)**2 / Ty &
                     + (current%part_p(3) - driftz)**2 / Tz) / two_kb_mass

        npart_per_cell = current%pvol

        ! We want to calculate the distribution of markers.
        distribution = EXP(-f0_exponent) * npart_per_cell * idx &
            / SQRT(two_pi_kb_mass3 * Tx * Ty * Tz)
        current%pvol = 1.0_num / distribution

#if DELTAF_DEBUG
        f0_back = f0(ispecies, mass, current%part_p)

        ! Checks for correct particle weight calculation.
        weight_back = f0_back * current%pvol
#ifndef PER_SPECIES_WEIGHT
        part_weight = current%weight
#endif
        WRITE(*,*) ipart, distribution, f0_exponent, npart_per_cell, &
            SQRT(two_pi_kb_mass3 * Tx * Ty * Tz), kb, mass
        WRITE(*,*) ipart, 'R', EXP(-f0_exponent), EXP(-f0_exponent) &
            * npart_per_cell / SQRT(two_pi_kb_mass3 * Tx * Ty * Tz)
        WRITE(*,*) ipart, 'Q', distribution, f0_back, weight_back, part_weight
#endif

        current => current%next
        ipart = ipart + 1
      ENDDO
    ENDDO
#endif

  END SUBROUTINE manual_load

END MODULE ic_module
