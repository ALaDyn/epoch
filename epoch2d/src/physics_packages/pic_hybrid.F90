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

MODULE pic_hybrid
#ifdef PIC_HYBRID

! We require hybrid to be pre-compiled if using PIC_HYBRID
#ifndef HYBRID
#define(HYBRID)
#endif

  USE hybrid
  USE fields

  IMPLICIT NONE

CONTAINS

  SUBROUTINE initialise_pic_hybrid

    IF (use_pic_hybrid) THEN
      ALLOCATE(is_hybrid(1-ng:nx+ng,1-ng:ny+ng))
      ALLOCATE(field_frac(1-ng:nx+ng,1-ng:ny+ng))
      CALL set_is_hybrid
      CALL set_field_frac

      IF (use_pic_hybrid_fields) THEN
        ALLOCATE(ex_pic(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(ey_pic(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(ez_pic(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(bx_pic(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(by_pic(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(bz_pic(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(ex_hy(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(ey_hy(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(ez_hy(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(bx_hy(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(by_hy(1-ng:nx+ng,1-ng:ny+ng))
        ALLOCATE(bz_hy(1-ng:nx+ng,1-ng:ny+ng))

        ex_pic = ex
        ey_pic = ey
        ez_pic = ez
        bx_pic = bx
        by_pic = by
        bz_pic = bz
        ex_hy = ex
        ey_hy = ey
        ez_hy = ez
        bx_hy = bx
        by_hy = by
        bz_hy = bz
      END IF
    END IF

  END SUBROUTINE initialise_pic_hybrid



  SUBROUTINE set_is_hybrid

    ! This subroutine will loop through all cells, and assign each one a logical
    ! is_hybrid value. This is set to true if any solid species has a non-zero
    ! number-density in this cell

    INTEGER :: ix, iy
    INTEGER :: i_sol

    is_hybrid = .FALSE.
    DO i_sol = 1,solid_count
      DO iy = 1-ng,ny+ng
        DO ix = 1-ng,nx+ng
          IF (solid_array(i_sol)%ion_density(ix,iy) > c_tiny) THEN
            is_hybrid(ix,iy) = .TRUE.
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE set_is_hybrid



  SUBROUTINE set_field_frac

    ! This subroutine will loop through all cells, setting all PIC cells to 1,
    ! and hybrid cells to 0. We will then average each cell with the neighbours
    ! and itself to smooth out the boundaries. This weighting will be used when
    ! running the field updater to allow for a smoother transition between
    ! algorithms. Only cells with field_frac < 1 will be updated in the average
    ! to avoid a discontinuous rise in resistivity on the boundary.

    INTEGER :: ix, iy, i_smooth
    INTEGER :: smooth_count_x, smooth_count_y
    REAL(num) :: x_up(1-ng:nx+ng,1-ng:ny+ng), x_down(1-ng:nx+ng,1-ng:ny+ng)
    REAL(num) :: y_up(1-ng:nx+ng,1-ng:ny+ng), y_down(1-ng:nx+ng,1-ng:ny+ng)
    REAL(num) :: temp(1-ng:nx+ng,1-ng:ny+ng)
    REAL(num) :: interp_x, interp_y

    ! Initialise field_frac to is_hybrid
    field_frac = 1.0_num
    DO iy = 1-ng,ny+ng
      DO ix = 1-ng,nx+ng
        IF (is_hybrid(ix,iy)) field_frac(ix,iy) = 0.0_num
      END DO
    END DO

    ! How many cells will be averaged over?
    smooth_count_x = CEILING(overlap_depth/dx)
    smooth_count_y = CEILING(overlap_depth/dy)

    ! What is the difference between consecutive field_frac values
    interp_x = 1.0_num/REAL(smooth_count_x, num)
    interp_y = 1.0_num/REAL(smooth_count_y, num)

    ! No smoothing will be performed if overlap_depth = 0
    IF (MAX(smooth_count_x, smooth_count_y) == 0) RETURN

    ! 1D interpolations: loop over all cells, each time a 0 touches a non-zero
    ! in a given direction, replace that zero with the next value along in the
    ! interpolation
    !
    ! e.g. first pass, say 1 sits below 0. If smooth count is 10, interp_x = 0.1
    ! so replace the 0 with 0.9 (1-0.1). The next pass, another 0 is above 0.9,
    ! replace with 0.8 (0.9-0.1), etc. Update ghost cells each pass to allow
    ! interpolation over boundaries

    ! x_up
    x_up = field_frac
    DO i_smooth = 1,smooth_count_x
      temp = x_up
      DO iy = 1,ny
        DO ix = 1,nx
          ! Test if zero is touching a non-zero in the "x up" direction
          IF (x_up(ix,iy) < c_tiny) THEN
            IF (x_up(ix-1,iy) > 0.0_num) THEN
              temp(ix,iy) = temp(ix-1,iy) - interp_x
            END IF
          END IF
        END DO
      END DO
      x_up = temp
      ! Pass neighbouring ghost cells
      CALL field_bc(x_up, ng)
    END DO

    ! x_down
    x_down = field_frac
    DO i_smooth = 1,smooth_count_x
      temp = x_down
      DO iy = 1,ny
        DO ix = 1,nx
          ! Test if zero is touching a non-zero in the "x up" direction
          IF (x_down(ix,iy) < c_tiny) THEN
            IF (x_down(ix+1,iy) > 0.0_num) THEN
              temp(ix,iy) = temp(ix+1,iy) - interp_x
            END IF
          END IF
        END DO
      END DO
      x_down = temp
      ! Pass neighbouring ghost cells
      CALL field_bc(x_down, ng)
    END DO

    ! y_up
    y_up = field_frac
    DO i_smooth = 1,smooth_count_y
      temp = y_up
      DO iy = 1,ny
        DO ix = 1,nx
          ! Test if zero is touching a non-zero in the "x up" direction
          IF (y_up(ix,iy) < c_tiny) THEN
            IF (y_up(ix,iy-1) > 0.0_num) THEN
              temp(ix,iy) = temp(ix,iy-1) - interp_y
            END IF
          END IF
        END DO
      END DO
      y_up = temp
      ! Pass neighbouring ghost cells
      CALL field_bc(y_up, ng)
    END DO

    ! y_down
    y_down = field_frac
    DO i_smooth = 1,smooth_count_y
      temp = y_down
      DO iy = 1,ny
        DO ix = 1,nx
          ! Test if zero is touching a non-zero in the "x up" direction
          IF (y_down(ix,iy) < c_tiny) THEN
            IF (y_down(ix,iy+1) > 0.0_num) THEN
              temp(ix,iy) = temp(ix,iy+1) - interp_y
            END IF
          END IF
        END DO
      END DO
      y_down = temp
      ! Pass neighbouring ghost cells
      CALL field_bc(y_down, ng)
    END DO

    ! These interpolations can overlap with each other. To combine them all, we
    ! will add the fractions together, weighted by what remains when the
    ! previous fractions are subtracted from one
    field_frac = x_up + (1.0_num - x_up)*x_down
    field_frac = field_frac + (1.0_num - field_frac)*y_up
    field_frac = field_frac + (1.0_num - field_frac)*y_down
    CALL field_bc(y_down, ng)

  END SUBROUTINE set_field_frac



  SUBROUTINE set_in_hybrid

    ! Loops over all the particles and checks which algorithm region they are in

    INTEGER :: ispecies, ipart, ix, iy
    INTEGER :: idx, idy
    REAL(num) :: part_x, part_y
    TYPE(particle), POINTER :: current, next

    idx = 1.0_num/dx
    idy = 1.0_num/dy

    ! Loop over all particle species
    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head

      ! Loop over all particles in the species
      DO ipart = 1, species_list(ispecies)%attached_list%count
        next => current%next

        ! Check type of current cell
        part_x = current%part_pos(1) - x_grid_min_local
        part_y = current%part_pos(2) - y_grid_min_local
        ix = MAX(1,CEILING(part_x*idx))
        iy = MAX(1,CEILING(part_y*idy))
        IF (is_hybrid(ix,iy)) THEN
          current%in_hybrid = .TRUE.
        ELSE
          current%in_hybrid = .FALSE.
        END IF

        current => next
      END DO
    END DO

  END SUBROUTINE set_in_hybrid



  SUBROUTINE electron_transport

    ! Run the hybrid routines which don't update the fields (collisions and
    ! heating)

    CALL set_in_hybrid
    IF (use_hybrid_scatter) THEN
      CALL run_elastic_scatter
    END IF
    CALL get_heat_capacity
    IF (use_hybrid_collisions) THEN
      CALL run_ionisation_loss
  !    CALL ion_bethe_bloch   ! Continuous approximation only at present
    END IF
    CALL remove_slow_particles
    CALL pic_hybrid_ohmic_heating
    CALL clear_heat_capacity
    CALL update_resistivity

  END SUBROUTINE electron_transport



  SUBROUTINE pic_hybrid_fields_half

    ! Calculate the field update with both traditional PIC and hybrid PIC, and
    ! only apply the update relevant to the code region. In the overlap region,
    ! the fields will be interpolated using the field_frac values.
    !
    ! Iterates E and B from t to t+dt/2

    INTEGER :: ix, iy
    REAL(num) :: hy_frac

    ! Option to deactivate field solver
    IF (use_pic_hybrid_fields) THEN
      ! PIC update
      CALL update_pic_E
      CALL update_pic_B(.TRUE.)

      ! Hybrid update
      CALL update_hybrid_B
      CALL update_hybrid_E

      ! Overlap solutions
      DO iy = 1-ng,ny+ng
        DO ix = 1-ng,nx+ng
          hy_frac = 1.0_num - field_frac(ix,iy)
          ex(ix,iy) = field_frac(ix,iy)*ex_pic(ix,iy) + hy_frac*ex_hy(ix,iy)
          ey(ix,iy) = field_frac(ix,iy)*ey_pic(ix,iy) + hy_frac*ey_hy(ix,iy)
          ez(ix,iy) = field_frac(ix,iy)*ez_pic(ix,iy) + hy_frac*ez_hy(ix,iy)
          bx(ix,iy) = field_frac(ix,iy)*bx_pic(ix,iy) + hy_frac*bx_hy(ix,iy)
          by(ix,iy) = field_frac(ix,iy)*by_pic(ix,iy) + hy_frac*by_hy(ix,iy)
          bz(ix,iy) = field_frac(ix,iy)*bz_pic(ix,iy) + hy_frac*bz_hy(ix,iy)
        END DO
      END DO
    END IF

  END SUBROUTINE pic_hybrid_fields_half



  SUBROUTINE pic_hybrid_fields_final

    ! Calculate the field update with both traditional PIC and hybrid PIC, and
    ! only apply the update relevant to the code region. In the overlap region,
    ! the fields will be interpolated using the field_frac values.
    !
    ! Iterates E and B from t+dt/2 to t+dt

    INTEGER :: ix, iy
    REAL(num) :: hy_frac

    ! Option to deactivate field solver
    IF (use_pic_hybrid_fields) THEN
      ! PIC update
      CALL update_pic_B(.FALSE.)
      CALL update_pic_E

      ! Hybrid update
      CALL update_hybrid_B
      CALL update_hybrid_E

      ! Overlap solutions
      DO iy = 1-ng,ny+ng
        DO ix = 1-ng,nx+ng
          hy_frac = 1.0_num - field_frac(ix,iy)
          ex(ix,iy) = field_frac(ix,iy)*ex_pic(ix,iy) + hy_frac*ex_hy(ix,iy)
          ey(ix,iy) = field_frac(ix,iy)*ey_pic(ix,iy) + hy_frac*ey_hy(ix,iy)
          ez(ix,iy) = field_frac(ix,iy)*ez_pic(ix,iy) + hy_frac*ez_hy(ix,iy)
          bx(ix,iy) = field_frac(ix,iy)*bx_pic(ix,iy) + hy_frac*bx_hy(ix,iy)
          by(ix,iy) = field_frac(ix,iy)*by_pic(ix,iy) + hy_frac*by_hy(ix,iy)
          bz(ix,iy) = field_frac(ix,iy)*bz_pic(ix,iy) + hy_frac*bz_hy(ix,iy)
        END DO
      END DO
    END IF

  END SUBROUTINE pic_hybrid_fields_final



  SUBROUTINE update_pic_E

    ! Updates the PIC electric fields. We do not update this in cells which have
    ! a field_frac < 1e-10 (pure hybrid cells have field_frac = 0.0)
    !
    ! Note, initial electric fields are zero everywhere and will remain zero if
    ! we don't update them. This clamps the fields to zero on the PIC-Hybrid
    ! boundary, which corresponds to reflecting boundary conditions

    INTEGER :: ix, iy
    REAL(num) :: cx1, cy1

    ! Pre-derive constants
    hdt  = 0.5_num * dt
    hdtx = hdt / dx
    hdty = hdt / dy
    cx1 = hdtx * c**2
    cy1 = hdty * c**2
    fac = hdt / epsilon0

    ! Loop over all cells
    DO iy = 1, ny
      DO ix = 1, nx

        ! Ignore pure hybrid cells
        IF (field_frac(ix,iy) < 1.0e-10) CYCLE

        ! Update PIC fields
        ex_pic(ix, iy) = ex_pic(ix, iy) &
            + cy1 * (bz_pic(ix  , iy  ) - bz_pic(ix  , iy-1)) &
            - fac * jx(ix, iy)
        ey_pic(ix, iy) = ey_pic(ix, iy) &
            - cx1 * (bz_pic(ix  , iy  ) - bz_pic(ix-1, iy  )) &
            - fac * jy(ix, iy)
        ez_pic(ix, iy) = ez_pic(ix, iy) &
            + cx1 * (by_pic(ix  , iy  ) - by_pic(ix-1, iy  )) &
            - cy1 * (bx_pic(ix  , iy  ) - bx_pic(ix  , iy-1)) &
            - fac * jz(ix, iy)
      END DO
    END DO

    ! Address boundary conditions (note these use the global field values)
    ex = ex_pic
    ey = ey_pic
    ez = ez_pic
    CALL efield_bcs
    ex_pic = ex
    ey_pic = ey
    ez_pic = ez

  END SUBROUTINE



  SUBROUTINE update_pic_B(bc_switch)

    ! Updates the PIC magnetic fields fields. We do not update this in cells
    ! which have a field_frac < 1e-10 (pure hybrid cells have field_frac = 0.0)
    !
    ! Note, initial magnetic fields are zero everywhere and will remain zero if
    ! we don't update them. This clamps the fields to zero on the PIC-Hybrid
    ! boundary, which corresponds to reflecting boundary conditions

    LOGICAL :: bc_switch
    INTEGER :: ix, iy
    REAL(num) :: cx1, cy1

    DO iy = 1, ny
      DO ix = 1, nx

        ! Ignore pure hybrid cells
        IF (field_frac(ix,iy) < 1.0e-10) CYCLE

        bx_pic(ix, iy) = bx_pic(ix, iy) &
            - hdty * (ez_pic(ix  , iy+1) - ez_pic(ix  , iy  ))
        by_pic(ix, iy) = by_pic(ix, iy) &
            + hdtx * (ez_pic(ix+1, iy  ) - ez_pic(ix  , iy  ))
        bz_pic(ix, iy) = bz_pic(ix, iy) &
            - hdtx * (ey_pic(ix+1, iy  ) - ey_pic(ix  , iy  )) &
            + hdty * (ex_pic(ix  , iy+1) - ex_pic(ix  , iy  ))
      END DO
    END DO

    ! Address boundary conditions (note these use the global field values)
    bx = bx_pic
    by = by_pic
    bz = bz_pic
    IF (bc_switch) THEN
      CALL bfield_bcs(.TRUE.)
    ELSE
      CALL bfield_final_bcs
    END IF
    bx_pic = bx
    by_pic = by
    bz_pic = bz

  END SUBROUTINE



  SUBROUTINE update_hybrid_B

    ! This subroutine performs a half-step in the magnetic field, assuming a
    ! constant electric field and using:
    !
    ! dB/dt = -curl(E)
    !
    ! We calculate into 1 ghost cell to allow non-zero curls across simulation
    ! boundaries
    !
    ! This routine only updates hybrid fields in hybrid cells - any gradients
    ! across the PIC-hybrid divide will be ignored

    INTEGER :: ix, iy, i

    ! Update B by half a timestep, only in the hybrid regime
    DO iy = 0, ny+1
      DO ix = 0, nx+1

        ! Don't update the hybrid fields for non-hybrid cells
        IF (.NOT. is_hybrid(ix,iy)) CYCLE

        ! Ensure all required cells are hybrid cells
        IF (is_hybrid(ix,iy+1)) THEN
          bx_hy(ix, iy) = bx_hy(ix, iy) &
              - hybrid_const_dt_by_dy * (ez_hy(ix  , iy+1) - ez_hy(ix  , iy  ))
          bz_hy(ix, iy) = bz_hy(ix, iy) &
              + hybrid_const_dt_by_dy * (ex_hy(ix  , iy+1) - ex_hy(ix  , iy  ))
        END IF

        IF (is_hybrid(ix+1,iy)) THEN
          by_hy(ix, iy) = by_hy(ix, iy) &
              + hybrid_const_dt_by_dx * (ez_hy(ix+1, iy  ) - ez_hy(ix  , iy  ))
          bz_hy(ix, iy) = bz_hy(ix, iy) &
              - hybrid_const_dt_by_dx * (ey_hy(ix+1, iy  ) - ey_hy(ix  , iy  ))
        END IF

      END DO
    END DO

    ! Pass neighbouring ghost cells
    CALL field_bc(bx_hy, ng)
    CALL field_bc(by_hy, ng)
    CALL field_bc(bz_hy, ng)

    ! Special cases for boundary processors
    DO i = 1, 2*c_ndims
      CALL field_zero_curl(bx_hy, i)
      CALL field_zero_curl(by_hy, i)
      CALL field_zero_curl(bz_hy, i)
    END DO

  END SUBROUTINE update_hybrid_B



  SUBROUTINE update_hybrid_E

    ! Calculates the electric field for the current values of global variables
    ! bx, by, bz, j, and resistivity, using the equation:
    !
    ! E = resistivity * (curl(B)/mu_0 - J)
    !
    ! We have precalculated hybrid_const_dx = 1/(mu_0*dx)

    ! Note that B is staggered from E, whereas J is evaulated at the same point
    ! as E. Resistivity is a cell centred variable. We calculate into 1 ghost
    ! cell to allow non-zero curls across simulation boundaries
    !
    ! We do not calculate E values in non-hybrid cells, and assume zero curl
    ! over PIC-Hybrid boundaries

    INTEGER :: ix, iy, i
    REAL(num) :: mean_res, bdif1, bdif2

    DO iy = 0, ny+1
      DO ix = 0, nx+1

        ! Don't update the hybrid fields for non-hybrid cells
        IF (.NOT. is_hybrid(ix,iy)) CYCLE
                                                                ! NOTE: IF CURRENT DENSITY IS CAUSING THE PROBLEM DUE TO
                                                                ! PARTICLE SHAPE EXTENDING INTO BOTH REGIONS, THEN TURN
                                                                ! OFF CURRENT DENSITY CALCULATION IN OVERLAP REGION:
                                                                ! IF(FIELD_FRAC(IX,IY) > 0.5), IGNORE J
        ! ex
        ! Resistivity at ex(ix,iy)
        IF (is_hybrid(ix+1,iy)) THEN
          mean_res = 0.5_num * (resistivity(ix+1,iy) + resistivity(ix,iy))
        ELSE
          mean_res = resistivity(ix,iy)
        END IF
        ! B curl at ex(ix,iy)
        IF (is_hybrid(ix, iy-1)) THEN
          bdif1 = hybrid_const_dy * (bz_hy(ix,iy) - bz_hy(ix,iy-1))
        ELSE
          bdif1 = 0.0_num
        END IF
        ! Update ex(ix,iy)
        ex_hy(ix,iy) = mean_res * (bdif1 - jx(ix,iy))

        ! ey
        ! Resistivity at ey(ix,iy)
        IF (is_hybrid(ix,iy+1)) THEN
          mean_res = 0.5_num * (resistivity(ix,iy+1) + resistivity(ix,iy))
        ELSE
          mean_res = resistivity(ix,iy)
        END IF
        ! B curl at ey(ix,iy)
        IF (is_hybrid(ix-1, iy)) THEN
          bdif1 = - hybrid_const_dx * (bz_hy(ix,iy) - bz_hy(ix-1,iy))
        ELSE
          bdif1 = 0.0_num
        END IF
        ! Update ey(ix,iy)
        ey_hy(ix,iy) = mean_res * (bdif1 - jy(ix,iy))

        ! ez
        ! Part of B curl at ez(ix,iy)
        IF (is_hybrid(ix-1, iy)) THEN
          bdif1 = hybrid_const_dx * (by_hy(ix,iy) - by_hy(ix-1,iy))
        ELSE
          bdif1 = 0.0_num
        END IF
        ! Remaining B curl at ez(ix,iy)
        IF (is_hybrid(ix, iy-1)) THEN
          bdif2 = - hybrid_const_dy * (bx_hy(ix,iy) - bx_hy(ix,iy-1))
        ELSE
          bdif2 = 0.0_num
        END IF
        ! Update ez(ix,iy)
        ez_hy(ix,iy) = resistivity(ix,iy) * (bdif1 + bdif2 - jz(ix,iy))
      END DO
    END DO

    ! Pass neighbouring ghost cells
    CALL field_bc(ex_hy, ng)
    CALL field_bc(ey_hy, ng)
    CALL field_bc(ez_hy, ng)

    ! Special cases for boundary processors
    DO i = 1, 2*c_ndims
      CALL field_zero_curl(ex_hy, i)
      CALL field_zero_curl(ey_hy, i)
      CALL field_zero_curl(ez_hy, i)
    END DO

  END SUBROUTINE update_hybrid_E



  SUBROUTINE pic_hybrid_ohmic_heating

    ! Calculates the Ohmic heating using the hybrid fields ONLY. This way the
    ! much higher laser fields cannot contribute to heating the solid. Same form
    ! as ohmic_heating in hybrid.F90

    REAL(num) :: E2
    REAL(num) :: eff_heat_capacity(1-ng:nx+ng,1-ng:ny+ng)
    REAL(num) :: fac(1-ng:nx+ng,1-ng:ny+ng)
    INTEGER :: ix, iy, i_sol

    fac = dt * hybrid_const_heat

    ! Obtain the effective heat capacity: sum(ne/C)
    eff_heat_capacity = 0.0_num
    DO i_sol = 1, solid_count
      DO iy = 1, ny
        DO ix = 1, nx
          ! No Ohmic heating in PIC region
          IF (.NOT. is_hybrid(ix,iy)) CYCLE

          eff_heat_capacity(ix,iy) = eff_heat_capacity(ix,iy) &
              + solid_array(i_sol)%el_density(ix,iy) &
              / solid_array(i_sol)%heat_capacity(ix,iy)
        END DO
      END DO
    END DO

    ! Loop over all grid points to find temperature change
    DO iy = 1, ny
      DO ix = 1, nx
        ! No Ohmic heating in PIC region
        IF (.NOT. is_hybrid(ix,iy)) CYCLE

        ! Tb is a cell-centred variable, but E has stagger - need to average
        ! Use the HYBRID fields (ignore files outside hybrid region)
        E2 = ez_hy(ix,iy)**2
        IF (is_hybrid(ix-1,iy)) THEN
          E2 = E2 + (0.5_num*(ex_hy(ix,iy) + ex_hy(ix-1,iy)))**2
        ELSE
          E2 = E2 + (0.5_num*ex_hy(ix,iy))**2
        END IF
        IF (is_hybrid(ix,iy-1)) THEN
          E2 = E2 + (0.5_num*(ey_hy(ix,iy) + ey_hy(ix,iy-1)))**2
        ELSE
          E2 = E2 + (0.5_num*ey_hy(ix,iy))**2
        END IF

        ! Calculate Ohmic heating, avoiding the 0/0 NaN. Cap temp at 1.0e30
        IF (resistivity(ix,iy) > 0.0_num .AND. &
            eff_heat_capacity(ix,iy) > 0.0_num) THEN
          hybrid_Tb(ix,iy) = MIN(hybrid_Tb(ix,iy) &
              + E2*fac(ix,iy)*eff_heat_capacity(ix,iy)/resistivity(ix,iy), &
              1.0e30_num)
        END IF

        ! IF (HYBRID_TB(IX,IY) > 1E10_NUM) THEN
        !   PRINT*, 'STUPID TEMPERATURE'
        !   PRINT*, 'IX,IY:', IX, IY
        !   PRINT*, 'E: ', ex_hy(ix,iy), ey_hy(ix,iy), ez_hy(ix,iy)
        !   PRINT*, 'E_STAG: ', ex_hy(ix-1,iy), ey_hy(ix,iy-1)
        !   PRINT*, 'IS_HYBRID: ', is_hybrid(ix,iy), is_hybrid(ix-1,iy), is_hybrid(ix,iy-1)
        !   PRINT*, 'HEATS: ', fac(ix,iy), eff_heat_capacity(ix,iy), resistivity(ix,iy)
        !   PRINT*, 'TEMP: ', HYBRID_TB(IX, IY)
        !   STOP
        ! ENDIF
      END DO
    END DO

    ! Pass new temperature values to ghost cells of neighbouring processors
    CALL field_bc(hybrid_Tb, ng)

  END SUBROUTINE pic_hybrid_ohmic_heating



  SUBROUTINE ion_bethe_bloch

    ! The hybrid collisions are only designed for electrons, but in PIC-hybrid,
    ! ions may enter the hybrid region. This subroutine calculates the
    ! ionisation energy loss of ions (Bethe-Bloch), using the approximate form
    ! given in eq. (2.8) in "Berger et al (1984), Report 37, Journal of the
    ! International Commission on Radiation Units and Measurements"

    INTEGER :: ispecies
    INTEGER(i8) :: ipart
    REAL(num) :: heat_const(1-ng:nx+ng,1-ng:ny+ng)
    REAL(num) :: part_heat_const, c_bb_1
    REAL(num) :: ipart_mc, m2c4, sum_dp, part_ne, frac
    REAL(num) :: px, py, pz, p, gamma, v, ipart_KE, ln_terms, delta_p, beta2
    REAL(num) :: sin_t, cos_t, cos_p, sin_t_cos_p, sin_t_sin_p
    REAL(num) :: p_new_2, weight, dE, part_x, part_y, part_C, delta_Tb
    REAL(num) :: idx, idy
    INTEGER :: ix, iy, i_sol
    TYPE(particle), POINTER :: current, next

    idx = 1.0_num/dx
    idy = 1.0_num/dy

    c_bb_1 = 2.0_num*m0*c**2

    ! 1 / (kb * sum(ne) * dx * dy) - last bit is per unit volume, but dz is 1m
    ! in epoch2d
    heat_const = hybrid_const_heat / (dx * dy)

    ! Loop over all ion species
    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      IF (species_list(ispecies)%species_type == c_species_id_electron &
          .OR. species_list(ispecies)%species_type == c_species_id_photon) &
          CYCLE

      ipart_mc  = 1.0_num / (c * species_list(ispecies)%mass)
      m2c4 = species_list(ispecies)%mass**2 * c**4

      ! Loop over all particles in the species
      DO ipart = 1, species_list(ispecies)%attached_list%count
        next => current%next

        ! Ignore collisional energy loss in PIC region
        IF (.NOT. current%in_hybrid) THEN
          next => current%next
          current => next
          CYCLE
        END IF

        ! Extract particle variables
        part_x = current%part_pos(1) - x_grid_min_local
        part_y = current%part_pos(2) - y_grid_min_local
        px = current%part_p(1)
        py = current%part_p(2)
        pz = current%part_p(3)
        p = SQRT(px**2 + py**2 + pz**2)
        gamma = SQRT((p * ipart_mc)**2 + 1.0_num)
        v = p * ipart_mc * c / gamma
        ipart_KE = (gamma - 1.0_num) * species_list(ispecies)%mass * c**2
        ln_terms = c_bb_1*(gamma**2 - 1.0_num)
        beta2 = (v/c)**2

        ! Extract solid variables averaged over all solids present in this cell
        sum_dp = 0.0_num
        DO i_sol = 1, solid_count
          CALL hy_grid_centred_var_at_particle(part_x, part_y, part_ne, &
              solid_array(i_sol)%el_density)

          ! If ne is zero, then the solid is not in this cell, so no energy loss
          IF (part_ne < TINY(1.0_num)) CYCLE

          ! Very low energy particles (~50 eV) will make lthe log term go
          ! negative causing a NaN in dp. The hybrid mode isn't designed for low
          ! energy ions, so dp will be ignored in this case
          sum_dp = sum_dp + part_ne * &
              ( LOG(MAX(ln_terms/solid_array(i_sol)%hybrid_Iex, 1.0_num)) &
              - beta2)
        END DO

        ! Collisional changes to momentum
        delta_p =  hybrid_const_dp_ion * sum_dp * dt * &
            species_list(ispecies)%charge**2 / v**2

        ! Reduce particle momentum
        frac = MAX(1.0_num + delta_p / ABS(p), 0.0_num)
        current%part_p(1) = px * frac
        current%part_p(2) = py * frac
        current%part_p(3) = pz * frac

        ! Calculate energy change due to collisions
        p_new_2 = current%part_p(1)**2 + current%part_p(2)**2 &
            + current%part_p(3)**2
#ifndef PER_SPECIES_WEIGHT
        weight = current%weight
#else
        weight = species_list(ispecies)%weight
#endif
        dE = (SQRT(p_new_2*c**2 + m2c4) - SQRT((p*c)**2 + m2c4)) * weight

        ! Get effective heat capacity at the particle position
        CALL hy_grid_centred_var_at_particle(part_x, part_y, part_heat_const, &
            heat_const)
        CALL get_effective_heat_capacity(part_x, part_y, part_C)

        ! Calculate the temperature increase, and add this to Tb. We use -dE,
        ! as the energy gain for temperature is equal to the energy loss of
        ! the electron.
        delta_Tb = - dE * part_C * part_heat_const

        ! Write temperature change to the grid (ignores particle shape)
        ! Impose maximum temperature of 1e30 to prevent temperature rising
        ! unphysically
        ix = MAX(1,CEILING(part_x*idx))
        iy = MAX(1,CEILING(part_y*idy))
        hybrid_Tb(ix,iy) = MIN(hybrid_Tb(ix,iy) + delta_Tb, 1.0e30_num)

        current => next

      END DO
    END DO

    ! Pass new temperature values to ghost cells of neighbouring processors
    CALL field_bc(hybrid_Tb, ng)

  END SUBROUTINE ion_bethe_bloch

#endif
END MODULE pic_hybrid
