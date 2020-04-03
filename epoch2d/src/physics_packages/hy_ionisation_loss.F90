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
!
! This module contains a PIC implementation of the ionisation energy loss
! functions in Geant4. This conists of a continuous energy loss for electrons
! as they pass through a solid, and a discrete electron-emission from the
! background, known as Moller scattering or "delta rays". The equivalent process
! for e+ ionisation of background e- (Bhahba scatter) has not yet been
! implemented.

MODULE hy_ionisation_loss
#ifdef HYBRID

  USE partlist
  USE calc_df

  IMPLICIT NONE

  ! Variables declared here may be seen and changed by all functions and
  ! subroutines below the contains line, but by nothing outside this module
  INTEGER, PRIVATE :: ispecies
  REAL(num), PRIVATE :: min_p2

CONTAINS

  SUBROUTINE setup_hy_ionisation_loss

    ! Initialises the optical depths for all electron species

    INTEGER :: iu, io
    LOGICAL :: found
    TYPE(particle_species) :: current_species
    TYPE(particle), POINTER :: current
    REAL(num) :: p_tau

    ! Keep optical depths unchanged if loading from restart dump
    IF (ic_from_restart) RETURN

    ! Randomly sample the optical depth from an exponential distribution for
    ! every particle in every species
    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        p_tau = random()
        current%optical_depth_delta = -LOG(1.0_num - p_tau)
        current => current%next
      END DO
    END DO

    min_p2 = (min_hybrid_energy**2 - mc2**2)/c**2

  END SUBROUTINE setup_hy_ionisation_loss



  FUNCTION check_ionisation_loss_variables()

    ! Function to ensure en electron species is present for ionisation delta
    ! rays.

    INTEGER :: check_ionisation_loss_variables
    INTEGER :: io, iu

    check_ionisation_loss_variables = c_err_none

    ! No special species required if hybrid collisions are turned off
    IF (.NOT. use_hybrid_collisions) RETURN

    ! No special species required if we only do energy loss and scatter
    IF (.NOT. produce_delta_rays) RETURN

    ! If the user hasn't chosen a delta_electron species, then assign the first
    ! electron species
    IF (delta_electron_species < 0) THEN
      DO ispecies = 1,n_species
        IF (species_list(ispecies)%species_type == c_species_id_electron) THEN
          delta_electron_species = ispecies
          ! Warn user about possibility of unexpected electrons in this species
          IF (rank == 0) THEN
            DO iu = 1, nio_units
              io = io_units(iu)
              WRITE(io,*) ''
              WRITE(io,*) '*** Warning ***'
              WRITE(io,*) 'No delta-ray electron species specified'
              WRITE(io,*) 'Delta-ray electrons will be written to: '
              WRITE(io,*) TRIM(species_list(ispecies)%name)
              WRITE(io,*) 'Delta-ray species can be specified with &
                  "identify:delta_electron" in the'
              WRITE(io,*) 'species block'
              WRITE(io,*) ''
            END DO
          END IF
          EXIT
        END IF
      END DO
    END IF

    ! Check if there exists a species to populate with delta e-
    IF (delta_electron_species < 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'No delta-ray electron species specified. Specify &
              using "identify:delta_electron"'
          END DO
        END IF
        check_ionisation_loss_variables = c_err_missing_elements
        RETURN
    END IF

  END FUNCTION check_ionisation_loss_variables



  SUBROUTINE run_ionisation_loss

    ! Updates the optical depth associated with the production of delta-rays
    ! (high energy electrons pulled from the solid background via ionisation),
    ! also applies continuous energy loss. This is the subroutine called by the
    ! hybrid loop.

    INTEGER :: i_sol, part_count, i_part
    TYPE(particle), POINTER :: current
    REAL(num) :: heat_const(1-ng:nx+ng,1-ng:ny+ng)
    REAL(num) :: part_p2, gamma_rel, part_KE
    REAL(num) :: part_x, part_y, part_ne, part_ne_test, part_v, dE
    REAL(num) :: max_ne, part_Iex_2, part_Iex_C

    ! 1 / (kb * sum(ne) * dx * dy) - last bit is per unit volume, but dz is 1m
    ! in epoch2d. Sum(ne) is the sum over all solids in a given cell (an array),
    ! not the sum over all cells
    heat_const = hybrid_const_heat / (dx * dy)

    ! Update the optical depth for each electron species
    DO ispecies = 1, n_species

      ! Only update the optical depths for the electron species
      IF (species_list(ispecies)%species_type /= c_species_id_electron) CYCLE

      ! Only update particles present at the start of the step (don't remove
      ! KE from delta-rays just created)
      part_count = species_list(ispecies)%attached_list%count
      current => species_list(ispecies)%attached_list%head
      DO i_part = 1, part_count

        ! Only update particles in the PIC-Hybrid region
#ifdef PIC_HYBRID
        IF (use_pic_hybrid) THEN
          IF (.NOT. current%in_hybrid) THEN
            current=>current%next
            CYCLE
          END IF
        END IF
#endif

        ! Only consider particles with non-zero energy
        part_p2 = current%part_p(1)**2 + current%part_p(2)**2 &
            + current%part_p(3)**2
        IF (part_p2 < c_tiny) THEN
          current=>current%next
          CYCLE
        END IF

        ! Get total number density at electron position
        part_x = current%part_pos(1) - x_grid_min_local
        part_y = current%part_pos(2) - y_grid_min_local
        CALL hy_grid_centred_var_at_particle(part_x, part_y, part_ne, &
            hybrid_const_sum_ne)

        ! No background solid, so no ionisation contribution
        IF (part_ne <= c_tiny) THEN
          current => current%next
          CYCLE
        END IF

        ! The hybrid model requires a single mean excitation energy for all
        ! solids at a given position. Assume the mean excitation energy felt by
        ! the electron belongs to the solid of the highest ne at the particle
        ! position
        max_ne = 0.0_num
        DO i_sol = 1, solid_count
          CALL hy_grid_centred_var_at_particle(part_x, part_y, part_ne_test, &
              solid_array(i_sol)%el_density)
          IF (part_ne_test > max_ne) THEN
            max_ne = part_ne_test
            part_Iex_2 = solid_array(i_sol)%Iex_term
            part_Iex_C = solid_array(i_sol)%dEdx_C
          END IF
        END DO

        ! Energy loss due to creation of delta-rays with energy below discrete
        ! treatment threshold
        CALL continuous_energy_loss(current, part_ne, part_Iex_2, part_Iex_C, &
            dE)

        ! Convert energy loss to increase in grid temperature
        CALL calculate_heating(dE, part_x, part_y, heat_const)

        ! Recalculate particle variables after continuous energy loss
        part_p2 = current%part_p(1)**2 + current%part_p(2)**2 &
            + current%part_p(3)**2
        gamma_rel = SQRT(part_p2/mc0**2 + 1.0_num)
        part_KE = (gamma_rel-1.0_num)*mc2
        current%particle_energy = part_KE + mc2
        part_v = SQRT(part_p2) * c**2 / current%particle_energy

        ! Don't update the optical depth if the electron kinetic energy is
        ! below 2*KE_cut_delta (maximum delta-ray energy from an electron is
        ! half its current KE)
        IF (part_KE < two_KE_cut_delta) THEN
          current => current%next
          CYCLE
        END IF

        ! Update the optical depth
        current%optical_depth_delta = current%optical_depth_delta &
            - cross_sec_delta(part_v, gamma_rel, part_KE) * part_ne * dt

        ! If optical depth dropped below zero generate electron and reset
        ! optical depth
        IF (current%optical_depth_delta <= 0.0_num) THEN
          CALL generate_delta(current, delta_electron_species, gamma_rel, &
              part_v, part_p2, part_KE, dE)
          current%optical_depth_delta = reset_optical_depth()

          ! Convert energy loss to increase in grid temperature if delta ray
          ! isn't added to the simulation
          IF (dE > c_tiny) THEN
            CALL calculate_heating(dE, part_x, part_y, heat_const)
          END IF
        END IF

        current => current%next
      END DO
    END DO

    ! Pass new temperature values to ghost cells of neighbouring processors
    CALL field_bc(hybrid_Tb, ng)

  END SUBROUTINE run_ionisation_loss



  SUBROUTINE continuous_energy_loss(current, part_ne, Iex_2, Iex_C, dE)

    ! Calculate the ionisation energy loss of the current particle due to the
    ! creation of delta-rays below the kinetic energy KE_cut_delta (shared_data)
    !
    ! The energy loss dE is saved, to be later used in heating the grid
    !
    ! This uses the Berger-Seltzer formula, using the method described in the
    ! Geant4 Physics Reference Manual, Section 10.1.2
    !
    ! part_ne: Total electron number density from all solids at e- position
    ! Iex_2 = 2/(Iex/(mc2))^2: To speed up dEdx calculation
    ! Iex_C = 1 + 2*LOG(Iex*SQRT(eps0*m0)/(hbar*q0)): speed up C calc for
    !                                                 dens_corr switch

    TYPE(particle), POINTER :: current
    REAL(num), INTENT(IN) :: part_ne, Iex_2, Iex_C
    REAL(num), INTENT(OUT) :: dE
    REAL(num) :: part_p2, gamma, part_KE, part_v
    REAL(num) :: inv_gam2, tau, tau_up, dens_corr, p_new, p_frac
    REAL(num) :: x_val, dEdx_C, dEdx_x0, dEdx_x1, dEdx_a

    ! Calculate electron variables
    part_p2 = current%part_p(1)**2 + current%part_p(2)**2 &
        + current%part_p(3)**2
    gamma = SQRT(part_p2/mc0**2 + 1.0_num)
    part_KE = (gamma-1.0_num)*mc2
    current%particle_energy = part_KE + mc2
    part_v = SQRT(part_p2) * c**2 / current%particle_energy

    ! Calculate energy change (low energy particles lose all their energy)
    IF (current%particle_energy < min_hybrid_energy) THEN
      dE = part_KE
    ELSE
      ! Precalculate repeated terms
      inv_gam2 = 1.0_num/gamma**2
      tau = gamma - 1.0_num
      tau_up = MIN(c_hybrid_tau_min, 0.5_num*tau)

      ! Obtain density correction factor delta
      IF (part_KE < c_hybrid_dens_lim) THEN
        ! Density correction is always zero if x_val < 0.2, which is equivalent
        ! to a kinetic energy of ~445 keV
        dens_corr = 0.0_num
      ELSE
        ! Obtain sampling variables
        x_val = LOG(gamma**2 - 1.0_num)*c_hybrid_dEdx_x
        dEdx_C = Iex_C - LOG(part_ne)

        ! Switch for x0 and x1 based on C. First statement is equivalent to
        ! IF (I < 100eV)
        IF (Iex_2 > c_hybrid_C_switch) THEN
          IF (dEdx_C <= 3.681_num) THEN
            dEdx_x0 = 0.2_num
          ELSE
            dEdX_x0 = 0.326_num*dEdx_C - 1.0_num
          END IF
          dEdx_x1 = 2.0_num
        ELSE
          IF (dEdx_C <= 5.215_num) THEN
            dEdx_x0 = 0.2_num
          ELSE
            dEdX_x0 = 0.326_num*dEdx_C - 1.5_num
          END IF
          dEdx_x1 = 3.0_num
        END IF

        ! Switch for the density correction based on x, x0 and x1
        IF (x_val <= dEdx_x0) THEN
          ! x <= x0
          dens_corr = 0.0_num
        ELSE IF (x_val < dEdx_x1) THEN
          ! x0 < x < x1
          dEdx_a = (dEdx_C  - 4.606_num*dEdx_x0)/(dEdx_x1 - dEdx_x0)**3
          dens_corr = 4.606_num*x_val - dEdx_C + dEdx_a*(dEdx_x1 - x_val)**3
        ELSE
          ! x => x1
          dens_corr = 4.606_num*x_val - dEdx_C
        END IF
      END IF

      ! Calculate energy loss
      dE = c_hybrid_ion_dedx*part_ne/part_v &
          * (LOG(Iex_2 * (gamma + 1.0_num) &
          * (tau + tau_up) * tau_up &
          * (1.0_num - 0.5_num*tau_up)**((2.0_num*tau + 1.0_num)*inv_gam2)) &
          - 1.0_num - (part_v/c)**2 + tau/(tau - tau_up) &
          + 0.5_num*tau_up**2*inv_gam2 - dens_corr) * dt
    END IF

    ! Apply energy loss
    IF (dE >= part_KE) THEN
      ! Particle cannot lose more energy than its current kinetic energy
      current%part_p(:) = 0.0_num
      dE = part_KE
    ELSE
      p_new = SQRT((current%particle_energy-dE)**2 - mc2**2)/c
      p_frac = p_new/SQRT(part_p2)
      current%part_p(:) = p_frac * current%part_p(:)
    END IF
    current%particle_energy = current%particle_energy - dE

    ! True energy loss includes weight
#ifdef PER_SPECIES_WEIGHT
    dE = dE * species_list(ispecies)%weight
#else
    dE = dE * current%weight
#endif

  END SUBROUTINE continuous_energy_loss



  SUBROUTINE calculate_heating(dE, part_x, part_y, heat_const)

    ! The current particle has deposited energy dE into the solid, at position
    ! (part_x, part_y), due to ionisation energy loss. By finding the heat
    ! capacity at that point, this is converted to a temperature increase.
    !
    ! The temperature rise of a compound solid will be approximated to:
    !
    ! dT = (Energy change per unit vol.) * Sum(ne/C) / [Sum(ne)]² / kb
    !
    ! Where Sum refers to a sum over solids. We do this so that if we have two
    ! equivalent solids of half their usual electron density, this returns:
    !
    ! dT = (Energy change per unit vol.) / (C * ne_Total) / kb
    !
    ! which is the result for a single solid of normal density.
    !
    ! heat_const = 1 / kb / sum(ne²) / (cell vol.)
    ! effective heat capacity = sum(ne/C)

    REAL(num), INTENT(IN) :: dE, part_x, part_y
    REAL(num), INTENT(IN) :: heat_const(1-ng:nx+ng, 1-ng:ny+ng)
    INTEGER :: ix, iy
    REAL(num) :: part_C, part_heat_const, delta_Tb

    CALL hy_grid_centred_var_at_particle(part_x, part_y, part_heat_const, &
        heat_const)
    CALL get_effective_heat_capacity(part_x, part_y, part_C)

    ! Calculate the temperature increase, and add this to Tb
    delta_Tb = dE * part_C * part_heat_const

    ! Write temperature change to the grid (ignores particle shape)
    ix = MAX(1,CEILING(part_x*hybrid_idx))
    iy = MAX(1,CEILING(part_y*hybrid_idy))
    hybrid_Tb(ix,iy) = hybrid_Tb(ix,iy) + delta_Tb

  END SUBROUTINE calculate_heating



  SUBROUTINE get_effective_heat_capacity(part_x, part_y, part_C)

    ! Temperature increase for one background species would be:
    !
    ! (energy change per unit vol.) * (vol. of e-, 1/ne) / (heat cap. of e-)
    !
    ! This implementation of the Davies model can use multiple backgrounds, so
    ! we define an effective heat capacity to find the total temperature rise.
    ! This assumes the energy is split evenly between all e- in the cell, which
    ! convert the energy into a temperature increase using different C values
    !
    ! (effective heat capacity) = sum(ne/C)

    REAL(num), INTENT(IN) :: part_x, part_y
    REAL(num), INTENT(OUT) :: part_C
    INTEGER :: i_sol
    REAL(num) :: C_sol, part_ne

    part_C = 0.0_num

    ! Loop over all solids
    DO i_sol = 1, solid_count
      ! Get heat capacity at particle position for current solid
      CALL hy_grid_centred_var_at_particle(part_x, part_y, C_sol, &
          solid_array(i_sol)%heat_capacity)
      CALL hy_grid_centred_var_at_particle(part_x, part_y, part_ne, &
          solid_array(i_sol)%el_density)

      ! Contribution to the effective heat capacity
      part_C = part_C + part_ne / C_sol
    END DO

  END SUBROUTINE get_effective_heat_capacity



  FUNCTION cross_sec_delta(part_v, gamma_rel, part_KE)

    ! Calculate the cross section of creating a delta-ray with energy over
    ! KE_cut_delta. This uses the same formula as can be found in the Geant4
    ! physics reference manual for Moller scatter.
    !
    ! The proper form here should be 1/v², however we must multiply the
    ! cross section by ne*v*dt to get a change in optical depth, so we have
    ! cancelled out one of the 1/v inside the cross section with the v outside

    REAL(num), INTENT(IN) :: part_v, gamma_rel, part_KE
    REAL(num) :: gam_minus_1, inv_gam2, KE_cut_frac, one_minus_KEcf, inv_KEcf
    REAL(num) :: cross_sec_delta

    ! Precalculate repeated terms
    gam_minus_1 = gamma_rel - 1.0_num
    inv_gam2 = 1.0_num/(gamma_rel**2)
    KE_cut_frac = KE_cut_delta/part_KE
    one_minus_KEcf = 1.0_num - KE_cut_frac
    inv_KEcf = 1.0_num/KE_cut_frac

    ! Get cross section
    cross_sec_delta = c_hybrid_ion_sig /(part_v * gam_minus_1) &
        * ((gam_minus_1**2 * inv_gam2)*(0.5_num - KE_cut_frac) &
        + inv_KEcf - 1.0_num/one_minus_KEcf &
        - (2.0_num*gamma_rel - 1.0_num)*inv_gam2*LOG(one_minus_KEcf*inv_KEcf))

  END FUNCTION cross_sec_delta



  FUNCTION reset_optical_depth()

    ! Draws a new random number for the exponentially distributed optical depths

    REAL(num) :: reset_optical_depth
    REAL(num) :: p_tau

    p_tau = random()
    reset_optical_depth = -LOG(1.0_num - p_tau)

  END FUNCTION reset_optical_depth



  SUBROUTINE generate_delta(electron, idelta, gamma, part_v, part_p2, part_KE, &
      dE)

    ! Generates a delta-ray (electron ionised from background) with properties
    ! derived from an incident electron, and writes the delta-electron to the
    ! species idelta. We use the sampling method outlined in the Geant4 Physics
    ! Reference Manual

    TYPE(particle), POINTER :: electron
    INTEGER, INTENT(IN) :: idelta
    REAL(num), INTENT(IN) :: gamma, part_v, part_p2, part_KE
    REAL(num), INTENT(OUT) :: dE
    TYPE(particle), POINTER :: delta_ray
    REAL(num) :: E0_delta, E0_delta_term, gam2, rej1, rej2, rej3
    REAL(num) :: test_E, inv_minus_tE, accept, inv_c
    REAL(num) :: E_el, E_delta, p_el, p_init, p_delta, p_delta_frac, p_el_frac
    REAL(num) :: cos_th_el, cos_th_de, phi_el, phi_de

    ! Pre-calculate terms for the accept/reject sampling of the delta-ray energy
    E0_delta = KE_cut_delta/part_KE
    E0_delta_term = 1.0_num - 2.0_num*E0_delta
    gam2 = gamma**2
    rej1 = (gamma - 1.0_num)**2
    rej2 = -2.0_num*(gam2 + gamma) + 1.0_num
    rej3 = 2.25_num*gam2 - 2.5_num*gamma + 1.25_num

    ! Sample delta-ray energy using method in Geant4 Physics Reference Manual,
    ! Section 10.1.4
    DO
      ! Generate test energy fraction
      test_E = E0_delta / (1.0_num - E0_delta_term*rand())

      ! Acceptance check
      inv_minus_tE = 1.0_num/(1.0_num - test_E)
      accept = rej1*test_E**2 + rej2*test_E*inv_minus_tE + gam2*inv_minus_tE**2
      IF (rej3 * accept >= rand()) EXIT
    END DO
    E_delta = (test_E * part_KE + mc2)

    ! Calculate the momentum magnitudes of the delta-ray electron p_delta, and
    ! the incident electron after recoil p_el. The latter comes from
    ! conservation of energy
    inv_c = 1.0_num/c
    p_delta = SQRT(E_delta**2 - mc2**2) * inv_c
    E_el = electron%particle_energy + mc2 - E_delta
    p_el = SQRT(E_el**2 - mc2**2) * inv_c

    ! Scatter angle theta and azimuthal rotation phi from conservation of
    ! momentum
    p_init = SQRT(part_p2)
    cos_th_el = (p_el**2 - p_delta**2 + part_p2)/(2.0_num*p_init*p_el)
    phi_el = 2*pi*rand()

    ! Only add delta-rays to the delta-ray species if the user has requested,
    ! "produce_delta_rays", AND the delta-ray energy is over the cut-off
    ! min_delta_energy (this is always greater than or equal to
    ! min_hybrid_energy)
    IF (E_delta > min_delta_energy .AND. produce_delta_rays) THEN

      ! Create a delta_ray (new electron) at the incident electron position.
      ! This function also creates new optical depths for all relevant processes
      CALL create_particle(delta_ray)
      delta_ray%part_pos = electron%part_pos

      ! Set momentum to the direction of the incident particle, apply scatter
      ! later
      p_delta_frac = p_delta/p_init
      delta_ray%part_p(1) = electron%part_p(1) * p_delta_frac
      delta_ray%part_p(2) = electron%part_p(2) * p_delta_frac
      delta_ray%part_p(3) = electron%part_p(3) * p_delta_frac

      ! Remaining variables
      delta_ray%particle_energy = E_delta
      delta_ray%weight = electron%weight

      ! Calculate the scatter angle from conservation of momentum
      cos_th_de = (p_delta**2 - p_el**2 + part_p2)/(2.0_num*p_init*p_delta)
      phi_de = phi_el - pi

      ! Apply scatter angle
      CALL rotate_p(delta_ray, cos_th_de, phi_de, p_delta)

      ! Add particle to list
      CALL add_particle_to_partlist(species_list(idelta)%attached_list, &
          delta_ray)

      ! No energy needs adding to the grid, energy is conserved
      dE = 0.0_num

    ELSE
      ! Particle has created a discrete delta-ray, but this hasn't been added to
      ! the simulation. To conserve energy, we will assume the delta ray
      ! deposits its energy locally, and calculate the temperature increase
      dE = E_delta - mc2

      ! True energy loss includes weight
#ifdef PER_SPECIES_WEIGHT
      dE = dE * species_list(ispecies)%weight
#else
      dE = dE * electron%weight
#endif
    END IF

    ! Apply momentum reduction to the generating electron
    p_el_frac = p_el/p_init
    electron%part_p(1) = electron%part_p(1) * p_el_frac
    electron%part_p(2) = electron%part_p(2) * p_el_frac
    electron%part_p(3) = electron%part_p(3) * p_el_frac
    electron%particle_energy = electron%particle_energy - (E_delta - mc2)

    ! Apply rotation to the generating electron
    CALL rotate_p(electron, cos_th_el, phi_el, p_el)

  END SUBROUTINE generate_delta



  SUBROUTINE rotate_p(part, cos_theta, phi, part_p)

    ! Let the polar direction be defined as the initial momentum direction of
    ! the particle, part. This subroutine rotates that momentum direction by
    ! theta in the polar direction (we read in cos(theta)), and phi in the
    ! azimuthal direction, without changing the magnitude.
    !
    ! If we have already calculated the magnitude of the particle's momentum,
    ! this can also be fed into the subroutine to speed up the calculation

    TYPE(particle), POINTER :: part
    REAL(num), INTENT(IN) :: cos_theta, phi
    REAL(num), OPTIONAL :: part_p
    REAL(num) :: p, frac_p, pcos_theta, sin_theta, psin_theta
    REAL(num) :: ux, uy, uz, pfrac_uz, cos_phi, term_1, term_2

    ! Extract particle momentum
    IF (PRESENT(part_p)) THEN
      p = part_p
    ELSE
      p = SQRT(part%part_p(1)**2 + part%part_p(2)**2 + part%part_p(3)**2)
    END IF
    frac_p = 1.0_num / p

    ! Precalculate repeated terms
    pcos_theta = p * cos_theta
    sin_theta = SQRT(1.0_num - cos_theta**2)
    psin_theta = p*sin_theta
    uz = part%part_p(3) * frac_p

    IF (ABS(1.0_num - uz) < 1.0e-5_num) THEN
      ! Special case if the polar direction points along z
      part%part_p(1) = psin_theta * COS(phi)
      part%part_p(2) = psin_theta * SIN(phi)
      part%part_p(3) = pcos_theta * SIGN(1.0_num, uz)
    ELSE
      ! Precalculate repeated terms
      ux = part%part_p(1) * frac_p
      uy = part%part_p(2) * frac_p
      pfrac_uz = p/SQRT(1.0_num - uz**2)
      cos_phi = COS(phi)
      term_1 = sin_theta*cos_phi*uz*pfrac_uz + pcos_theta
      term_2 = sin_theta*SIN(phi)*pfrac_uz

      part%part_p(1) = ux * term_1 - uy * term_2
      part%part_p(2) = uy * term_1 + ux * term_2
      part%part_p(3) = uz * pcos_theta &
          + cos_phi * sin_theta * (uz**2 - 1.0_num) * pfrac_uz
    END IF

  END SUBROUTINE rotate_p



  SUBROUTINE hy_grid_centred_var_at_particle(part_x, part_y, part_var, grid_var)

    ! Calculates the value of a grid-centred variable "part_var" stored in the
    ! grid "grid_var", averaged over the particle shape for a particle at
    ! position (part_x, part_y)

    REAL(num), INTENT(in) :: part_x, part_y
    REAL(num), INTENT(in) :: grid_var(1-ng:nx+ng, 1-ng:ny+ng)
    REAL(num), INTENT(out) :: part_var
    INTEGER :: cell_x1, cell_y1
    REAL(num) :: cell_x_r, cell_y_r
    REAL(num) :: cell_frac_x, cell_frac_y
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = (1.0_num)**c_ndims
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    ! The following method is lifted from photons.F90 (field_at_particle), for
    ! the cell-centered fields, taking into account the various particle shapes
#ifdef PARTICLE_SHAPE_TOPHAT
    cell_x_r = part_x / dx - 0.5_num
    cell_y_r = part_y / dy - 0.5_num
#else
    cell_x_r = part_x / dx
    cell_y_r = part_y / dy
#endif
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1
    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1

#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"

    part_var = &
          gy(-2) * (gx(-2) * grid_var(cell_x1-2,cell_y1-2) &
        +           gx(-1) * grid_var(cell_x1-1,cell_y1-2) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1-2) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1-2) &
        +           gx( 2) * grid_var(cell_x1+2,cell_y1-2)) &
        + gy(-1) * (gx(-2) * grid_var(cell_x1-2,cell_y1-1) &
        +           gx(-1) * grid_var(cell_x1-1,cell_y1-1) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1-1) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1-1) &
        +           gx( 2) * grid_var(cell_x1+2,cell_y1-1)) &
        + gy( 0) * (gx(-2) * grid_var(cell_x1-2,cell_y1  ) &
        +           gx(-1) * grid_var(cell_x1-1,cell_y1  ) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1  ) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1  ) &
        +           gx( 2) * grid_var(cell_x1+2,cell_y1  )) &
        + gy( 1) * (gx(-2) * grid_var(cell_x1-2,cell_y1+1) &
        +           gx(-1) * grid_var(cell_x1-1,cell_y1+1) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1+1) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1+1) &
        +           gx( 2) * grid_var(cell_x1+2,cell_y1+1)) &
        + gy( 2) * (gx(-2) * grid_var(cell_x1-2,cell_y1+2) &
        +           gx(-1) * grid_var(cell_x1-1,cell_y1+2) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1+2) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1+2) &
        +           gx( 2) * grid_var(cell_x1+2,cell_y1+2))
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
    part_var = &
          gy( 0) * (gx( 0) * grid_var(cell_x1  ,cell_y1  ) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1  )) &
        + gy( 1) * (gx( 0) * grid_var(cell_x1  ,cell_y1+1) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1+1))
#else
#include "triangle/gx.inc"
    part_var = &
          gy(-1) * (gx(-1) * grid_var(cell_x1-1,cell_y1-1) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1-1) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1-1)) &
        + gy( 0) * (gx(-1) * grid_var(cell_x1-1,cell_y1  ) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1  ) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1  )) &
        + gy( 1) * (gx(-1) * grid_var(cell_x1-1,cell_y1+1) &
        +           gx( 0) * grid_var(cell_x1  ,cell_y1+1) &
        +           gx( 1) * grid_var(cell_x1+1,cell_y1+1))
#endif
    part_var = fac*part_var

  END SUBROUTINE hy_grid_centred_var_at_particle

#endif
END MODULE hy_ionisation_loss
