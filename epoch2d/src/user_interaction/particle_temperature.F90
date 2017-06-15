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

MODULE particle_temperature

  USE shared_data
  USE random_generator
  USE evaluator

  IMPLICIT NONE

CONTAINS

  ! Subroutine to initialise a thermal particle distribution
  SUBROUTINE setup_particle_temperature(part_species)

    TYPE(particle_species), POINTER :: part_species
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, temp_local
    REAL(num), DIMENSION(3) :: drift_local    
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    INTEGER :: idir
    TYPE(parameter_pack) :: parameters
#include "particle_head.inc"

    partlist => part_species%attached_list
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = part_species%mass
#endif

      ! Assume that temperature is cell centred
#include "particle_to_grid.inc"

      parameters%use_grid_position = .FALSE.
      parameters%pack_ix = cell_x
      parameters%pack_iy = cell_y
      parameters%pack_pos = current%part_pos
      DO idir = 1, 3
        temp_local = evaluate_with_parameters( &
            part_species%temperature_function(idir), parameters, errcode)
        current%part_p(idir) = momentum_from_temperature(mass,temp_local)
        drift_local(idir) = evaluate_with_parameters( &
            part_species%drift_function(idir), parameters, errcode)
      ENDDO
      CALL particle_drift_lorentz_transform(current, mass, drift_local)
      current => current%next
      ipart = ipart + 1
    ENDDO

  END SUBROUTINE setup_particle_temperature



  ! Subroutine to initialise a thermal particle distribution
  SUBROUTINE setup_particle_temperature_relativistic(part_species)

    TYPE(particle_species), POINTER :: part_species
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass
    REAL(num), DIMENSION(3) :: drift_local, temp_local
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    INTEGER :: ix, idir
    TYPE(parameter_pack) :: parameters
#include "particle_head.inc"

    partlist => part_species%attached_list
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = part_species%mass
#endif

#include "particle_to_grid.inc"
      parameters%use_grid_position = .FALSE.
      parameters%pack_ix = cell_x
      parameters%pack_iy = cell_y
      parameters%pack_pos = current%part_pos
      DO idir = 1, 3
        temp_local(idir) = evaluate_with_parameters( &
            part_species%temperature_function(idir), parameters, errcode)
        drift_local(idir) = evaluate_with_parameters( &
            part_species%drift_function(idir), parameters, errcode)
      ENDDO

      current%part_p = momentum_from_temperature_relativistic(mass, &
          temp_local, part_species%fractional_tail_cutoff)

      CALL particle_drift_lorentz_transform(current, mass, drift_local)
      !Advanced the linked list
      current=>current%next
      ipart = ipart + 1
    END DO

  END SUBROUTINE setup_particle_temperature_relativistic


  !Subroutine takes a particle and a drift momentum and Lorentz transforms
  !The particle momentum subject to the specified drift
  SUBROUTINE particle_drift_lorentz_transform(part, mass, drift)

    TYPE(particle), POINTER :: part
    REAL(num), INTENT(IN) :: mass
    REAL(num), DIMENSION(3), INTENT(IN) :: drift
    REAL(num) :: gamma_drift, gamma_part, e_prime
    INTEGER :: idir

    gamma_drift = SQRT(1.0_num + DOT_PRODUCT(drift/(mass*c),drift/(mass*c)))
    gamma_part = SQRT(1.0_num + DOT_PRODUCT(part%part_p/(mass*c),&
        part%part_p/(mass*c)))
    e_prime = gamma_part * mass * c**2
    DO idir = 1, 3
      part%part_p(idir) = part%part_p(idir) * gamma_drift + &
          drift(idir)/(mass*c**2) * e_prime
    ENDDO

  END SUBROUTINE particle_drift_lorentz_transform



  FUNCTION momentum_from_temperature(mass, temperature)

    REAL(num), INTENT(IN) :: mass, temperature
    REAL(num) :: momentum_from_temperature

    REAL(num) :: stdev
    REAL(num) :: rand1, rand2, w
    REAL(num), SAVE :: val
    LOGICAL, SAVE :: cached = .FALSE.

    ! This is a basic polar Box-Muller transform
    ! It generates gaussian distributed random numbers
    ! The standard deviation (stdev) is related to temperature

    stdev = SQRT(temperature * kb * mass)

    IF (cached) THEN
      cached = .FALSE.
      momentum_from_temperature = val * stdev
    ELSE
      cached = .TRUE.

      DO
        rand1 = random()
        rand2 = random()

        rand1 = 2.0_num * rand1 - 1.0_num
        rand2 = 2.0_num * rand2 - 1.0_num

        w = rand1**2 + rand2**2

        IF (w > c_tiny .AND. w < 1.0_num) EXIT
      ENDDO

      w = SQRT((-2.0_num * LOG(w)) / w)

      momentum_from_temperature = rand1 * w * stdev
      val = rand2 * w
    ENDIF

  END FUNCTION momentum_from_temperature



  FUNCTION momentum_from_temperature_relativistic(mass, temperature, cutoff)
    REAL(num), INTENT(IN) :: mass
    REAL(num), DIMENSION(3), INTENT(IN) :: temperature
    REAL(num), INTENT(IN) :: cutoff
    REAL(num), DIMENSION(3) :: momentum_from_temperature_relativistic

    !Three parameters for calculating the range of momenta
    !Includes different combinations of physical constants
    REAL(num), PARAMETER :: param1 = -3.07236e-40_num
    REAL(num), PARAMETER :: param2 = 2.35985e-80_num
    !c^2/kb
    REAL(num), PARAMETER :: c2ok = 6.509658203714208e39_num
    REAL(num) :: random_number, probability
    REAL(num) :: momentum_x, momentum_y, momentum_z, momentum
    REAL(num) :: temp, temp_max, p_max_x, p_max_y, p_max_z, p_max
    REAL(num) :: dof
    REAL(num) :: inter1, inter2, inter3

    temp = SUM(temperature)
    temp_max = MAXVAL(temperature)
    dof=COUNT(temperature > c_tiny)
    !If there are no degrees of freedom them sampling is unnecessary
    !plasma is cold
    IF (dof == 0) THEN
      momentum_from_temperature_relativistic = 0.0_num
      RETURN
    ENDIF

    p_max = SQRT(param1 * mass * temp * LOG(cutoff) + &
        param2 * temp**2 * LOG(cutoff)**2)/mass

    p_max_x = p_max * SQRT(temperature(1)/temp)
    p_max_y = p_max * SQRT(temperature(2)/temp)
    p_max_z = p_max * SQRT(temperature(3)/temp)

    !Loop around until a momentum is accepted for this particle
    DO
      !Generate random x and y momenta between p_min and p_max
      momentum_x = (random() * 2.0_num * p_max_x - p_max_x)
      momentum_y = (random() * 2.0_num * p_max_y - p_max_y)
      momentum_z = (random() * 2.0_num * p_max_z - p_max_z)
      momentum = SQRT(momentum_x**2+momentum_y**2+momentum_z**2)
      !From that value, have to generate the probability that a particle
      !with that momentum should be accepted.
      !This is just the particle distribution function scaled to have
      !a maximum of 1 (or lower).
      !In general you will have to work this out yourself
      inter1 = temperature(1)/temp
      inter2 = temperature(2)/temp
      inter3 = temperature(3)/temp
      probability = EXP(-c2ok * mass / temp * (SQRT(1.0_num + &
          1.0_num/MAX(inter1,c_tiny) * momentum_x**2 + &
          1.0_num/MAX(inter2,c_tiny) * momentum_y**2 + &
          1.0_num/MAX(inter3,c_tiny) * momentum_z**2)-1.0_num))
      !Once you know your probability you just generate a random number
      !between 0 and 1 and if the generated number is less than the
      !probability then accept the particle and exit this loop.
      random_number=random()
      IF (random_number <= probability) EXIT
    END DO

    momentum_from_temperature_relativistic = &
        (/ momentum_x , momentum_y , momentum_z /) * m0 * c

  END FUNCTION momentum_from_temperature_relativistic


  ! Function for generating momenta of thermal particles in a particular
  ! direction, e.g. the +x direction.
  ! These satisfy a Rayleigh distribution, formed by combining two
  ! normally-distributed (~N(0,sigma)) random variables as follows:
  ! Z = SQRT(X**2 + Y**2)
  FUNCTION flux_momentum_from_temperature(mass, temperature)

    REAL(num), INTENT(IN) :: mass, temperature
    REAL(num) :: flux_momentum_from_temperature
    REAL(num) :: mom1, mom2

    mom1 = momentum_from_temperature(mass, temperature)
    mom2 = momentum_from_temperature(mass, temperature)

    flux_momentum_from_temperature = SQRT(mom1**2 + mom2**2)

  END FUNCTION flux_momentum_from_temperature
 

 
  ! Subroutine to initialise a non-thermal particle distribution
  SUBROUTINE setup_particle_distfn(part_species)

    TYPE(particle_species), POINTER :: part_species
    TYPE(particle_list), POINTER :: partlist
    INTEGER(i8) :: ipart, iit, ipart_global, iit_global
    REAL(num) :: average_its, mass
    CHARACTER(LEN=25) :: string
#include "particle_head.inc"

    partlist => part_species%attached_list
    CALL sample_partlist_from_distfn(part_species, partlist, iit, ipart)
    CALL MPI_REDUCE(ipart, ipart_global, 1, MPI_INTEGER8, MPI_SUM, 0, comm, &
                    errcode)
    CALL MPI_REDUCE(iit, iit_global, 1, MPI_INTEGER8, MPI_SUM, 0, comm, errcode)
    average_its = REAL(iit_global, num) / REAL(ipart_global, num)

    IF (rank == 0) THEN
      WRITE(string,'(F8.1)') average_its
      WRITE(*,*) 'Setup distribution for particles of species ', '"' &
          // TRIM(part_species%name) // '"', ' taking ' // TRIM(string) &
          // ' iteratations per particle on average'
      IF (average_its >= 20.0_num) THEN
        WRITE(*,*) '***WARNING***'
        WRITE(*,*) 'Average iterations is high. ', &
            'Possibly try smaller momentum range'
      ENDIF
#ifndef NO_IO
      WRITE(stat_unit,*) 'Setup distribution for particles of species ', '"' &
          // TRIM(part_species%name) // '"', ' taking ' // TRIM(string) &
          // ' iteratations per particle on average'
      IF (average_its >= 20.0_num) THEN
        WRITE(stat_unit,*) '***WARNING***'
        WRITE(stat_unit,*) 'Average iterations is high. ', &
            'Possibly try smaller momentum range'
      ENDIF
#endif
    ENDIF

  END SUBROUTINE setup_particle_distfn



  !> Routine to sample a distribution function for an entire
  !> particle list
  SUBROUTINE sample_partlist_from_distfn(part_species, partlist, iit_r, ipart_r)

    REAL(num), DIMENSION(3,2) :: ranges
    TYPE(particle_species), POINTER :: part_species
    TYPE(particle_list), POINTER :: partlist
    INTEGER(i8), INTENT(OUT), OPTIONAL :: iit_r, ipart_r
    INTEGER(i8) :: iit, ipart
    TYPE(particle), POINTER :: current
    REAL(num), DIMENSION(3) :: drift_local
    INTEGER :: err, idir
    REAL(num) :: mass
    TYPE(parameter_pack) :: parameters
#include "particle_head.inc"

    current => partlist%head
    ipart = 0
    iit = 0

    DO WHILE(ASSOCIATED(current))
#include "particle_to_grid.inc"
#ifdef PER_PARTICLE_CHARGE_MASS
      mass = current%mass
#else
      mass = part_species%mass
#endif

      parameters%use_grid_position = .FALSE.
      parameters%pack_ix = cell_x
      parameters%pack_pos = current%part_pos
      err = 0
      CALL evaluate_with_parameters_to_array(part_species%dist_fn_range(1), &
          parameters, 2, ranges(1,:), err)
      CALL evaluate_with_parameters_to_array(part_species%dist_fn_range(2), &
          parameters, 2, ranges(2,:), err)
      CALL evaluate_with_parameters_to_array(part_species%dist_fn_range(3), &
          parameters, 2, ranges(3,:), err)
      IF (err /= c_err_none .AND. rank == 0) THEN
        PRINT*, 'Unable to evaluate distribution function ranges'
        CALL abort_code(c_err_bad_setup)
      ENDIF

      CALL sample_from_deck_expression(part_species%dist_fn, parameters, &
          ranges, iit)

      current%part_p = parameters%pack_p * c * mass

      DO idir = 1, 3
        drift_local(idir) = evaluate_with_parameters( &
            part_species%drift_function(idir), parameters, errcode)
      ENDDO
      CALL particle_drift_lorentz_transform(current, mass, drift_local)
      current => current%next
      ipart = ipart + 1
    ENDDO

    IF (PRESENT(iit_r)) iit_r=iit
    IF (PRESENT(ipart_r)) ipart_r = ipart

  END SUBROUTINE sample_partlist_from_distfn



  SUBROUTINE sample_from_deck_expression(stack, parameters, ranges, iit_r)

    TYPE(primitive_stack), INTENT(INOUT) :: stack
    TYPE(parameter_pack), INTENT(INOUT) :: parameters
    REAL(num), DIMENSION(3,2), INTENT(IN) :: ranges
    INTEGER(i8), INTENT(INOUT), OPTIONAL :: iit_r
    REAL(num) :: setlevel
    INTEGER :: err
    INTEGER(i8) :: iit

    err = c_err_none
    IF (PRESENT(iit_r)) iit = iit_r
    DO
      ! These lines are setting global variables that are later used by
      ! the deck parser
      parameters%pack_p(1) = random() * (ranges(1,2) - ranges(1,1)) &
          + ranges(1,1)
      parameters%pack_p(2) = random() * (ranges(2,2) - ranges(2,1)) &
          + ranges(2,1)
      parameters%pack_p(3) = random() * (ranges(3,2) - ranges(3,1)) &
          + ranges(3,1)

      !pack spatial information has already been set before calling
      CALL basic_evaluate(stack, parameters, err)
      setlevel = pop_off_eval()
      IF (err /= c_err_none .AND. rank == 0) THEN
        PRINT*, 'Unable to evaluate distribution function'
        CALL abort_code(c_err_bad_setup)
      ENDIF

      iit = iit + 1

      IF (random() <= setlevel) EXIT
    ENDDO
    IF (PRESENT(iit_r)) iit_r = iit

  END SUBROUTINE sample_from_deck_expression

END MODULE particle_temperature
