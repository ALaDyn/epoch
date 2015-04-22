MODULE ic_module

  USE shared_data
  USE helper
  USE particles

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: manual_load

CONTAINS

  ! Find local particle temperature: at the moment, just copied and
  ! pasted from particle_temperature.F90
  SUBROUTINE params_local(current,temperature,drift, temp_local,drift_local)
    TYPE(particle), POINTER, INTENT(IN) :: current
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: temperature
    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: drift
    INTEGER :: ix, iy, iz
    REAL(num), INTENT(INOUT) :: temp_local, drift_local
#include "particle_head.inc"

    temp_local = 0.0_num
    drift_local = 0.0_num

#include "particle_to_grid.inc"

    DO iz = sf_min, sf_max
       DO iy = sf_min, sf_max
          DO ix = sf_min, sf_max
             temp_local = temp_local + gx(ix) * gy(iy) * gz(iz) &
                  * temperature(cell_x+ix, cell_y+iy, cell_z+iz)
             drift_local = drift_local + gx(ix) * gy(iy) * gz(iz) &
                  * drift(cell_x+ix, cell_y+iy, cell_z+iz)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE params_local

  SUBROUTINE manual_load
    REAL(num) :: Tx,Ty,Tz,driftx,drifty,driftz
    REAL(num) :: f0_exponent, distribution, mass, npart_per_cell
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current
    TYPE(particle_species), POINTER :: species
    INTEGER :: ipart, ispecies
    REAL(num) :: part_weight,weight_back, f0_back
#ifdef DELTAF_METHOD
    ! f0 calculation: mainly, we need to calculate the phase space volumes.
    ! Calculate this based on the loading parameters. Easy to check
    ! that this is OK for a Maxwellian load by setting f0 = f0_back,
    ! and making sure the weights cancel.

    DO ispecies = 1, n_species
       species => species_list(ispecies)
       partlist => species%attached_list
       current => partlist%head
       ipart = 0
       DO WHILE(ipart < partlist%count)
#ifdef PER_PARTICLE_CHARGE_MASS
          mass = current%mass
#else
          mass = species%mass
#endif
          CALL params_local(current, initial_conditions(ispecies)%temp(:,:,:,1),   &
               &    initial_conditions(ispecies)%drift(:,:,:,1),          &
               &    Tx,driftx)
          CALL params_local(current, initial_conditions(ispecies)%temp(:,:,:,2),   &
               &    initial_conditions(ispecies)%drift(:,:,:,2),          &
               &    Ty,drifty)
          CALL params_local(current, initial_conditions(ispecies)%temp(:,:,:,3),   &
               &    initial_conditions(ispecies)%drift(:,:,:,3),          &
               &    Tz,driftz)

          f0_exponent = (  (current%part_p(1)-driftx)**2/(2*kb*Tx*mass)    & 
               &         + (current%part_p(2)-drifty)**2/(2*kb*Ty*mass)    &
               &         + (current%part_p(3)-driftz)**2/(2*kb*Tz*mass) )
          !
          npart_per_cell = current%pvol
           ! We want to calculate the distribution of markers. 
          distribution = exp(-f0_exponent)*(npart_per_cell/(dx*dy*dz))/sqrt( (2*pi*kb*mass)**3*Tx*Ty*Tz )
          current%pvol = 1.0/distribution

          f0_back = f0(ispecies,mass,current%part_p(1),current%part_p(2),current%part_p(3))

          ! Checks for correct particle weight calculation.
          weight_back = f0_back * current%pvol
#ifdef PER_SPECIES_WEIGHT
          part_weight = species_list(ispecies)%weight
#endif
#ifndef PER_SPECIES_WEIGHT
          part_weight = current%weight
#endif
        !  WRITE (*,*) ipart,distribution, f0_exponent, npart_per_cell, sqrt( (2*pi*kb*mass)**3*Tx*Ty*Tz ),kb,mass
        !  WRITE (*,*) ipart,'R',exp(-f0_exponent),exp(-f0_exponent)*npart_per_cell/sqrt( (2*pi*kb*mass)**3*Tx*Ty*Tz )
        !  WRITE (*,*) ipart,'Q',distribution,f0_back,weight_back,part_weight

          current => current%next
          ipart = ipart + 1
       ENDDO
    ENDDO

#endif

  END SUBROUTINE manual_load

END MODULE ic_module
