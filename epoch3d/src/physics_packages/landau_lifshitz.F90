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

MODULE landau_lifshitz
#ifdef LANDAU_LIFSHITZ

  USE partlist

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_landau_lifshitz()

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies

    ! Write out the start position of every particle
    DO ispecies = 1,n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        current%start_x = current%part_pos(1)
        current => current%next
      END DO
    END DO

  END SUBROUTINE setup_landau_lifshitz



  SUBROUTINE classical_radiation_reaction()

    TYPE(particle), POINTER :: current
    INTEGER :: ispecies
    REAL(num) :: part_x, part_y, part_z, part_ux, part_uy, part_uz, gamma_rel
    REAL(num) :: eta, e_at_part(3), b_at_part(3)
    REAL(num) :: power_loss, Gaunt, total_p, energy_loss, new_energy, new_p

    DO ispecies=1,n_species
      current => species_list(ispecies)%attached_list%head

      ! Landau Lifshitz does not apply to uncharged particles
      IF (species_list(ispecies)%species_type == c_species_id_photon) CYCLE

      DO WHILE(ASSOCIATED(current))

        ! Find eta at the particle position
        part_x  = current%part_pos(1) - x_grid_min_local
        part_y  = current%part_pos(2) - y_grid_min_local
        part_z  = current%part_pos(3) - z_grid_min_local
        part_ux = current%part_p(1) / mc0
        part_uy = current%part_p(2) / mc0
        part_uz = current%part_p(3) / mc0
        gamma_rel = SQRT(part_ux**2 + part_uy**2 + part_uz**2 + 1.0_num)
        CALL field_at_particle(part_x, part_y, part_z, e_at_part, b_at_part)
        current%ex_at_part = e_at_part(1)
        eta = calculate_eta(part_ux, part_uy, part_uz, gamma_rel, e_at_part, &
            b_at_part)
        current%power_eta = eta

        ! Only calculate RR for moving partilces
        IF (gamma_rel - 1.0_num > c_tiny) THEN
          ! Calculate momentum change
          Gaunt = (1.0_num + 4.8_num*(1.0_num + eta) &
              *log(1.0_num + 1.7_num*eta) + 2.4_num*eta**2)**(-2.0_num/3.0_num)
          ! tau_c is (reduced Compton wavelength / c)
          power_loss = 2.0_num/3.0_num/tau_c*alpha_f*m0*c**2*eta**2*Gaunt
          energy_loss = power_loss*dt
          new_energy = gamma_rel*m0*c**2 - energy_loss
          new_p = SQRT(new_energy**2 - m0**2*c**4)/c

          ! Change particle momentum
          total_p = SQRT(current%part_p(1)**2 + current%part_p(2)**2 &
              + current%part_p(3)**2)
          IF (new_p < total_p) THEN
            current%part_p(1) = new_p/total_p * current%part_p(1)
            current%part_p(2) = new_p/total_p * current%part_p(2)
            current%part_p(3) = new_p/total_p * current%part_p(3)
          ELSE
            current%part_p(1) = 0
            current%part_p(2) = 0
            current%part_p(3) = 0
          END IF
        END IF

        current => current%next
      END DO
    END DO

  END SUBROUTINE classical_radiation_reaction



  FUNCTION calculate_eta(part_ux, part_uy, part_uz, gamma_rel, e_at_part, &
      b_at_part)

    REAL(num) :: calculate_eta
    REAL(num), INTENT(IN) :: part_ux, part_uy, part_uz, gamma_rel
    REAL(num), INTENT(IN) :: e_at_part(3), b_at_part(3)
    REAL(num) :: beta_x, beta_y, beta_z
    REAL(num) :: flperp(3), i_e, tau0, roland_eta
    REAL(num) :: lambdac, coeff_eta, moduclip, moduclip2, u_dot_e

    moduclip2 = MAX(part_ux**2 + part_uy**2 + part_uz**2, c_tiny)
    moduclip = SQRT(moduclip2)

    beta_x = part_ux / gamma_rel
    beta_y = part_uy / gamma_rel
    beta_z = part_uz / gamma_rel

    lambdac = h_bar / mc0

    coeff_eta = SQRT(3.0_num * lambdac / (2.0_num * alpha_f * m0 * c**3))

    u_dot_e = (part_ux * e_at_part(1) + part_uy * e_at_part(2) &
        + part_uz * e_at_part(3)) / moduclip2

    flperp(1) = q0 * (e_at_part(1) - u_dot_e * part_ux &
        + c * (beta_y * b_at_part(3) - beta_z * b_at_part(2)))

    flperp(2) = q0 * (e_at_part(2) - u_dot_e * part_uy &
        + c * (beta_z * b_at_part(1) - beta_x * b_at_part(3)))

    flperp(3) = q0 * (e_at_part(3) - u_dot_e * part_uz &
        + c * (beta_x * b_at_part(2) - beta_y * b_at_part(1)))

    ! Dipole emission intensity

    tau0 = q0**2 / (6.0_num * pi * epsilon0 * m0 * c**3)

    i_e = tau0 * gamma_rel**2 * (flperp(1)**2 + flperp(2)**2 + flperp(3)**2 &
        + (q0 * (beta_x * e_at_part(1) + beta_y * e_at_part(2) &
        + beta_z * e_at_part(3)) / moduclip)**2) / m0

    roland_eta = coeff_eta * SQRT(i_e)

    ! Determine eta from fields
    calculate_eta = roland_eta

  END FUNCTION calculate_eta



  SUBROUTINE field_at_particle(part_x, part_y, part_z, e_at_part, b_at_part)

    REAL(num), INTENT(IN) :: part_x, part_y, part_z
    REAL(num), INTENT(OUT) :: e_at_part(3), b_at_part(3)

    ! Contains the integer cell position of the particle in x, y, z
    INTEGER :: cell_x1, cell_x2, cell_y1, cell_y2, cell_z1, cell_z2

    ! Contains the floating point version of the cell number (never actually
    ! used)
    REAL(num) :: cell_x_r, cell_y_r, cell_z_r

    ! The fraction of a cell between the particle position and cell boundary
    REAL(num) :: cell_frac_x, cell_frac_y, cell_frac_z

    ! Weighting factors as Eqn 4.77 page 25 of manual
    ! Eqn 4.77 would be written as
    ! F(j-1) * gmx + F(j) * g0x + F(j+1) * gpx
    ! Defined at the particle position
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz

    ! Defined at the particle position - 0.5 grid cell in each direction
    ! This is to deal with the grid stagger
    REAL(num), DIMENSION(sf_min:sf_max) :: hx, hy, hz
    ! Temporary variables
    INTEGER :: dcellx, dcelly, dcellz
    REAL(num) :: ex_part, ey_part, ez_part
    REAL(num) :: bx_part, by_part, bz_part

    ! Particle weighting multiplication factor
#ifdef PARTICLE_SHAPE_BSPLINE3
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (1.0_num / 24.0_num)**c_ndims
#elif  PARTICLE_SHAPE_TOPHAT
    REAL(num), PARAMETER :: fac = (1.0_num)**c_ndims
#else
    REAL(num) :: cf2
    REAL(num), PARAMETER :: fac = (0.5_num)**c_ndims
#endif

    ! Grid cell position as a fraction.
#ifdef PARTICLE_SHAPE_TOPHAT
    cell_x_r = part_x / dx - 0.5_num
    cell_y_r = part_y / dy - 0.5_num
    cell_z_r = part_z / dz - 0.5_num
#else
    cell_x_r = part_x / dx
    cell_y_r = part_y / dy
    cell_z_r = part_z / dz
#endif
    ! Round cell position to nearest cell
    cell_x1 = FLOOR(cell_x_r + 0.5_num)
    ! Calculate fraction of cell between nearest cell boundary and particle
    cell_frac_x = REAL(cell_x1, num) - cell_x_r
    cell_x1 = cell_x1 + 1

    cell_y1 = FLOOR(cell_y_r + 0.5_num)
    cell_frac_y = REAL(cell_y1, num) - cell_y_r
    cell_y1 = cell_y1 + 1

    cell_z1 = FLOOR(cell_z_r + 0.5_num)
    cell_frac_z = REAL(cell_z1, num) - cell_z_r
    cell_z1 = cell_z1 + 1

    ! Particle weight factors as described in the manual, page25
    ! These weight grid properties onto particles
    ! Also used to weight particle properties onto grid, used later
    ! to calculate J
    ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/gx.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/gx.inc"
#else
#include "triangle/gx.inc"
#endif

    ! Now redo shifted by half a cell due to grid stagger.
    ! Use shifted version for ex in X, ey in Y, ez in Z
    ! And in Y&Z for bx, X&Z for by, X&Y for bz
    cell_x2 = FLOOR(cell_x_r)
    cell_frac_x = REAL(cell_x2, num) - cell_x_r + 0.5_num
    cell_x2 = cell_x2 + 1

    cell_y2 = FLOOR(cell_y_r)
    cell_frac_y = REAL(cell_y2, num) - cell_y_r + 0.5_num
    cell_y2 = cell_y2 + 1

    cell_z2 = FLOOR(cell_z_r)
    cell_frac_z = REAL(cell_z2, num) - cell_z_r + 0.5_num
    cell_z2 = cell_z2 + 1

    dcellx = 0
    dcelly = 0
    dcellz = 0
    ! NOTE: These weights require an additional multiplication factor!
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/hx_dcell.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/hx_dcell.inc"
#else
#include "triangle/hx_dcell.inc"
#endif

    ! These are the electric and magnetic fields interpolated to the
    ! particle position. They have been checked and are correct.
    ! Actually checking this is messy.
#ifdef PARTICLE_SHAPE_BSPLINE3
#include "bspline3/e_part.inc"
#include "bspline3/b_part.inc"
#elif  PARTICLE_SHAPE_TOPHAT
#include "tophat/e_part.inc"
#include "tophat/b_part.inc"
#else
#include "triangle/e_part.inc"
#include "triangle/b_part.inc"
#endif

    ! update particle momenta using weighted fields
    ! ex_part etc are NOT fields at particle, but fac times
    ! field

    e_at_part(1) = fac * ex_part
    e_at_part(2) = fac * ey_part
    e_at_part(3) = fac * ez_part

    b_at_part(1) = fac * bx_part
    b_at_part(2) = fac * by_part
    b_at_part(3) = fac * bz_part

  END SUBROUTINE field_at_particle


#endif
END MODULE landau_lifshitz
