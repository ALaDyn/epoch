! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2012      Martin Ramsay <M.G.Ramsay@warwick.ac.uk>
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

MODULE laser

  USE custom_laser
  USE lorentz
  USE evaluator

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_laser(boundary, laser)

    INTEGER, INTENT(IN) :: boundary
    TYPE(laser_block), INTENT(INOUT) :: laser
    INTEGER, DIMENSION(2) :: array_range

    laser%boundary = boundary
    laser%id = -1
    laser%use_time_function = .FALSE.
    laser%use_phase_function = .FALSE.
    laser%use_profile_function = .FALSE.
    laser%use_omega_function = .FALSE.
    laser%amp = -1.0_num
    laser%pol_angle = 0.0_num
    laser%t_start = 0.0_num
    laser%t_end = t_end

    laser%kx_mult = 0.0_num

    NULLIFY(laser%profile)
    NULLIFY(laser%phase)
    NULLIFY(laser%next)
    NULLIFY(laser%omega)
    NULLIFY(laser%k)
    NULLIFY(laser%current_integral_phase)

    CALL allocate_with_boundary(laser%profile, boundary)
    CALL allocate_with_boundary(laser%phase, boundary)
    CALL allocate_with_boundary(laser%omega, boundary)
    CALL allocate_with_boundary(laser%current_integral_phase, boundary)
    array_range = get_boundary_range(boundary)
    ALLOCATE(laser%k(array_range(1):array_range(2), c_ndims))
    laser%profile = 1.0_num
    laser%phase = 0.0_num
    laser%k = 0.0_num
    laser%omega = 0.0_num
    laser%current_integral_phase = 0.0_num

  END SUBROUTINE init_laser



  SUBROUTINE setup_laser_phases(laser_init, phases)

    TYPE(laser_block), POINTER :: laser_init
    REAL(num), DIMENSION(1-ng:,:), INTENT(IN) :: phases
    TYPE(laser_block), POINTER :: laser
    INTEGER, DIMENSION(2) :: sizes, array_range, section
    INTEGER :: ilas

    ilas = 1
    laser => laser_init
    sizes = SHAPE(phases)
    DO WHILE(ASSOCIATED(laser))
      IF (sizes(1) == 1) THEN
        laser%current_integral_phase = phases(:,ilas)
      ELSE
        array_range = get_boundary_range(laser%boundary, &
            global_range = section, ghosts = .FALSE.)
        laser%current_integral_phase(array_range(1):array_range(2)) = &
            phases(section(1):section(2),ilas)
      END IF
      ilas = ilas + 1
      laser => laser%next
    END DO

  END SUBROUTINE setup_laser_phases



  SUBROUTINE deallocate_laser(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: idir

    IF (ASSOCIATED(laser%profile)) DEALLOCATE(laser%profile)
    IF (ASSOCIATED(laser%phase)) DEALLOCATE(laser%phase)
    IF (ASSOCIATED(laser%omega)) DEALLOCATE(laser%omega)
    IF (ASSOCIATED(laser%k)) DEALLOCATE(laser%k)
    IF (laser%use_profile_function) &
        CALL deallocate_stack(laser%profile_function)
    IF (laser%use_phase_function) &
        CALL deallocate_stack(laser%phase_function)
    IF (laser%use_time_function) &
        CALL deallocate_stack(laser%time_function)
    IF (laser%use_omega_function) &
        CALL deallocate_stack(laser%omega_function)
    DO idir = 1, c_ndims
      IF (laser%use_k_function(idir)) &
          CALL deallocate_stack(laser%k_function(idir))
    END DO
    DEALLOCATE(laser)

  END SUBROUTINE deallocate_laser



  SUBROUTINE deallocate_lasers

    TYPE(laser_block), POINTER :: current, next

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

    current => laser_y_min
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

    current => laser_y_max
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_laser(current)
      current => next
    END DO

  END SUBROUTINE deallocate_lasers



  ! Subroutine to attach a created laser object to the correct boundary
  SUBROUTINE attach_laser(laser)

    INTEGER :: boundary
    TYPE(laser_block), POINTER :: laser

    boundary = laser%boundary

    IF (boundary == c_bd_x_min) THEN
      n_laser_x_min = n_laser_x_min + 1
      laser%kx_mult = 1.0_num
      laser%initial_pos = [x_global(1), 0.0_num]
      CALL attach_laser_to_list(laser_x_min, laser)
    ELSE IF (boundary == c_bd_x_max) THEN
      n_laser_x_max = n_laser_x_max + 1
      laser%kx_mult = -1.0_num
      laser%initial_pos = [x_global(nx_global), 0.0_num]
      CALL attach_laser_to_list(laser_x_max, laser)
    ELSE IF (boundary == c_bd_y_min) THEN
      n_laser_y_min = n_laser_y_min + 1
      laser%kx_mult = 0.0_num
      laser%initial_pos = [0.0_num, y_global(1)]
      CALL attach_laser_to_list(laser_y_min, laser)
    ELSE IF (boundary == c_bd_y_max) THEN
      n_laser_y_max = n_laser_y_max + 1
      laser%kx_mult = 0.0_num
      laser%initial_pos = [0.0_num, y_global(ny_global)]
      CALL attach_laser_to_list(laser_y_max, laser)
    END IF

  END SUBROUTINE attach_laser



  ! This routine populates the constant elements of a parameter pack
  ! from a laser

  SUBROUTINE populate_pack_from_laser(laser, parameters)

    TYPE(laser_block) :: laser
    TYPE(parameter_pack), INTENT(INOUT) :: parameters

    parameters%pack_ix = 0
    parameters%pack_iy = 0

#ifdef BOOSTED_FRAME
    IF (in_boosted_frame) THEN
      parameters%v_prop = c * laser%kx_mult
    END IF
#endif

    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min)
        parameters%pack_ix = 0
      CASE(c_bd_x_max)
        parameters%pack_ix = nx
      CASE(c_bd_y_min)
        parameters%pack_iy = 0
      CASE(c_bd_y_max)
        parameters%pack_iy = ny
    END SELECT

  END SUBROUTINE populate_pack_from_laser



  SUBROUTINE populate_array_from_stack_1d(laser, array, stack)

    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num), DIMENSION(:), POINTER :: array
    TYPE(primitive_stack), INTENT(INOUT) :: stack
    TYPE(parameter_pack) :: parameters
    INTEGER :: err, i

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        DO i = 1,ny
          parameters%pack_iy = i
          array(i) = &
              evaluate_with_parameters(stack, parameters, err)
        END DO
      CASE(c_bd_y_min, c_bd_y_max)
        DO i = 1,nx
          parameters%pack_ix = i
          array(i) = &
              evaluate_with_parameters(stack, parameters, err)
        END DO
    END SELECT

  END SUBROUTINE populate_array_from_stack_1d



  SUBROUTINE populate_array_from_stack_ndims(laser, array, stacks)

    TYPE(laser_block), INTENT(IN) :: laser
    REAL(num), DIMENSION(:,:), POINTER :: array
    TYPE(primitive_stack), DIMENSION(:), INTENT(INOUT) :: stacks
    TYPE(parameter_pack) :: parameters
    INTEGER :: err, i, idims

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    SELECT CASE(laser%boundary)
      CASE(c_bd_x_min, c_bd_x_max)
        DO idims = 1, c_ndims
          IF (.NOT. stacks(idims)%init) CYCLE
          DO i = 1, ny
            parameters%pack_iy = i
            array(i, idims) = &
                evaluate_with_parameters(stacks(idims), parameters, err)
          END DO
        END DO
      CASE(c_bd_y_min, c_bd_y_max)
        DO idims = 1, c_ndims
          IF (.NOT. stacks(idims)%init) CYCLE
          DO i = 1, nx
            parameters%pack_ix = i
            array(i, idims) = &
                evaluate_with_parameters(stacks(idims), parameters, err)
          END DO
      END DO
    END SELECT

  END SUBROUTINE populate_array_from_stack_ndims



  FUNCTION laser_time_profile(laser)

    TYPE(laser_block), POINTER :: laser
    REAL(num) :: laser_time_profile
    INTEGER :: err
    TYPE(parameter_pack) :: parameters

    err = 0
    IF (laser%use_time_function) THEN
      CALL populate_pack_from_laser(laser, parameters)
      laser_time_profile = evaluate_with_parameters(laser%time_function, &
          parameters, err)
      RETURN
    END IF

    laser_time_profile = custom_laser_time_profile(laser)

  END FUNCTION laser_time_profile



  SUBROUTINE laser_update_phase(laser)

    TYPE(laser_block), POINTER :: laser

    CALL populate_array_from_stack_1d(laser, laser%phase, laser%phase_function)

  END SUBROUTINE laser_update_phase



  SUBROUTINE laser_update_profile(laser)

    TYPE(laser_block), POINTER :: laser

    CALL populate_array_from_stack_1d(laser, laser%profile, &
        laser%profile_function)

  END SUBROUTINE laser_update_profile



  SUBROUTINE laser_update_omega(laser)

    TYPE(laser_block), POINTER :: laser
#ifdef BOOSTED_FRAME
    INTEGER :: iy
#endif

    CALL populate_array_from_stack_1d(laser, laser%omega, laser%omega_function)
    IF (laser%omega_func_type == c_of_freq) &
        laser%omega = 2.0_num * pi * laser%omega
    IF (laser%omega_func_type == c_of_lambda) &
        laser%omega = 2.0_num * pi * c / laser%omega

    laser%k = 0.0_num

#ifdef BOOSTED_FRAME
    IF (in_boosted_frame .AND. (laser%boundary == c_bd_x_min &
        .OR. laser%boundary == c_bd_x_max)) THEN
      DO iy = 1, ny
        laser%omega(iy) = transform_frequency(global_boost_info, &
            laser%omega(iy), laser%omega(iy) / c * laser%kx_mult)
      END DO
    END IF
#endif

  END SUBROUTINE laser_update_omega



  SUBROUTINE laser_update_k(laser)

    TYPE(laser_block), POINTER :: laser
    INTEGER :: err, i
    INTEGER, DIMENSION(2) :: array_range
    TYPE(parameter_pack) :: parameters
#ifdef BOOSTED_FRAME
    REAL(num) :: kboost
#endif

    array_range = get_boundary_range(laser%boundary)

    err = 0
    CALL populate_pack_from_laser(laser, parameters)
    CALL populate_array_from_stack_ndims(laser, laser%k, laser%k_function)

    DO i = array_range(1), array_range(2)
      laser%omega(i) = SQRT(DOT_PRODUCT(laser%k(i,:), laser%k(i,:)) &
          * c**2)
#ifdef BOOSTED_FRAME
      IF (in_boosted_frame) THEN
        laser%omega(i) = transform_frequency(global_boost_info, &
            laser%omega(i), laser%k(i, 1), k_boost = kboost)
        IF (laser%k_function(1)%init) THEN
          laser%k(i,1) = kboost
        END IF
      END IF
#endif
    END DO

  END SUBROUTINE laser_update_k



  SUBROUTINE update_laser_omegas

    TYPE(laser_block), POINTER :: current

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE IF (ANY(current%use_k_function)) THEN
        CALL laser_update_k(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE
        current%current_integral_phase = current%omega * time
      END IF
      current => current%next
    END DO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE IF (ANY(current%use_k_function)) THEN
        CALL laser_update_k(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE
        current%current_integral_phase = current%omega * time
      END IF
      current => current%next
    END DO

    current => laser_y_min
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE IF (ANY(current%use_k_function)) THEN
        CALL laser_update_k(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE
        current%current_integral_phase = current%omega * time
      END IF
      current => current%next
    END DO

    current => laser_y_max
    DO WHILE(ASSOCIATED(current))
      IF (current%use_omega_function) THEN
        CALL laser_update_omega(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE IF (ANY(current%use_k_function)) THEN
        CALL laser_update_k(current)
        current%current_integral_phase = current%current_integral_phase &
            + current%omega * dt
      ELSE
        current%current_integral_phase = current%omega * time
      END IF
      current => current%next
    END DO

  END SUBROUTINE update_laser_omegas



  ! Actually does the attaching of the laser to the correct list
  SUBROUTINE attach_laser_to_list(list, laser)

    TYPE(laser_block), POINTER :: list
    TYPE(laser_block), POINTER :: laser
    TYPE(laser_block), POINTER :: current

    IF (ASSOCIATED(list)) THEN
      current => list
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => laser
    ELSE
      list => laser
    END IF

  END SUBROUTINE attach_laser_to_list



  FUNCTION get_boundary_range(boundary, global_range, ghosts)
    INTEGER, INTENT(IN) :: boundary
    INTEGER, DIMENSION(2), INTENT(OUT), OPTIONAL :: global_range
    LOGICAL, INTENT(IN), OPTIONAL :: ghosts
    INTEGER, DIMENSION(2) :: get_boundary_range
    INTEGER :: ngl

    ngl = ng
    IF (PRESENT(ghosts)) THEN
      IF (.NOT. ghosts) ngl = 0
    END IF

    IF (boundary == c_bd_x_min .OR. boundary == c_bd_x_max) THEN
      get_boundary_range = [1-ngl, ny+ngl]
      IF (PRESENT(global_range)) &
          global_range = [ny_global_min-ngl, ny_global_max+ngl]
    ELSE IF (boundary == c_bd_y_min .OR. boundary == c_bd_y_max) THEN
      get_boundary_range = [1-ngl, nx+ngl]
      IF (PRESENT(global_range)) &
          global_range = [nx_global_min-ngl, nx_global_max+ngl]
    END IF

  END FUNCTION get_boundary_range



  SUBROUTINE allocate_with_boundary(array, boundary)

    REAL(num), DIMENSION(:), POINTER :: array
    INTEGER, INTENT(IN) :: boundary
    INTEGER, DIMENSION(2) :: array_range

    array_range = get_boundary_range(boundary)
    ALLOCATE(array(array_range(1):array_range(2)))

  END SUBROUTINE allocate_with_boundary



  SUBROUTINE set_laser_dt

    REAL(num) :: dt_local
    TYPE(laser_block), POINTER :: current

    dt_laser = HUGE(1.0_num)

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / MAXVAL(current%omega)
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / MAXVAL(current%omega)
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    current => laser_y_min
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / MAXVAL(current%omega)
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    current => laser_y_max
    DO WHILE(ASSOCIATED(current))
      dt_local = 2.0_num * pi / MAXVAL(current%omega)
      dt_laser = MIN(dt_laser, dt_local)
      current => current%next
    END DO

    ! Need at least two iterations per laser period
    ! (Nyquist)
    dt_laser = dt_laser / 2.0_num

  END SUBROUTINE set_laser_dt



  SUBROUTINE outflow_bcs_x_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_min

    laserpos = 1
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_x_min_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(ny))
    ALLOCATE(source2(ny))
    source1 = 0.0_num
    source2 = 0.0_num

    bx(laserpos-1, 1:ny) = bx_x_min(1:ny)

    IF (add_laser(n)) THEN
      current => laser_x_min
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time >= current%t_start .AND. time <= current%t_end) THEN
          IF (current%use_phase_function) CALL laser_update_phase(current)
          IF (current%use_profile_function) CALL laser_update_profile(current)
          t_env = laser_time_profile(current) * current%amp
          DO i = 1,ny
            base = t_env * current%profile(i) &
              * SIN(current%current_integral_phase(i) + current%phase(i) &
              - current%k(i, 1) * (x(laserpos-1)-current%initial_pos(1)) &
              - current%k(i, 2) * (y(i)-current%initial_pos(2)))
            source1(i) = source1(i) + base * COS(current%pol_angle)
            source2(i) = source2(i) + base * SIN(current%pol_angle)
          END DO
        END IF
        current => current%next
      END DO
    END IF

#ifdef BOOSTED_FRAME
    IF (use_boosted_frame) THEN
      source1 = source1 * global_boost_info%lorentz_gamma &
          * (1.0_num - global_boost_info%beta)
      source2 = source2 * global_boost_info%lorentz_gamma &
          * (1.0_num + global_boost_info%beta)
    END IF
#endif

    bz(laserpos-1, 1:ny) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ey_x_min(1:ny) + c * bz_x_min(1:ny)) &
        - 2.0_num * ey(laserpos, 1:ny) &
        + dt_eps * jy(laserpos, 1:ny) &
        + diff * bz(laserpos, 1:ny))

    by(laserpos-1, 1:ny) = sum * (-4.0_num * source2 &
        - 2.0_num * (ez_x_min(1:ny) - c * by_x_min(1:ny)) &
        + 2.0_num * ez(laserpos, 1:ny) &
        - ly * (bx(laserpos, 1:ny) - bx(laserpos, 0:ny-1)) &
        - dt_eps * jz(laserpos, 1:ny) &
        + diff * by(laserpos, 1:ny))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_min, lasers=laser_x_min)
      ELSE
        CALL calc_absorption(c_bd_x_min)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_x_min



  SUBROUTINE outflow_bcs_x_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_x_max

    laserpos = nx
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_x_max_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (lx + c)
    diff = lx - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(ny))
    ALLOCATE(source2(ny))
    source1 = 0.0_num
    source2 = 0.0_num

    bx(laserpos+1, 1:ny) = bx_x_max(1:ny)

    IF (add_laser(n)) THEN
      current => laser_x_max
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time >= current%t_start .AND. time <= current%t_end) THEN
          IF (current%use_phase_function) CALL laser_update_phase(current)
          IF (current%use_profile_function) CALL laser_update_profile(current)
          t_env = laser_time_profile(current) * current%amp
          DO i = 1,ny
            base = t_env * current%profile(i) &
              * SIN(current%current_integral_phase(i) + current%phase(i) &
              - current%k(i, 1) * (x(laserpos-1)-current%initial_pos(1)) &
              - current%k(i, 2) * (y(i)-current%initial_pos(2)))
            source1(i) = source1(i) + base * COS(current%pol_angle)
            source2(i) = source2(i) + base * SIN(current%pol_angle)
          END DO
        END IF
        current => current%next
      END DO
    END IF

#ifdef BOOSTED_FRAME
    IF (use_boosted_frame) THEN
      source1 = source1 * global_boost_info%lorentz_gamma &
          * (1.0_num - global_boost_info%beta)
      source2 = source2 * global_boost_info%lorentz_gamma &
          * (1.0_num + global_boost_info%beta)
    END IF
#endif

    bz(laserpos, 1:ny) = sum * (-4.0_num * source1 &
        - 2.0_num * (ey_x_max(1:ny) - c * bz_x_max(1:ny)) &
        + 2.0_num * ey(laserpos, 1:ny) &
        - dt_eps * jy(laserpos, 1:ny) &
        + diff * bz(laserpos-1, 1:ny))

    by(laserpos, 1:ny) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ez_x_max(1:ny) + c * by_x_max(1:ny)) &
        - 2.0_num * ez(laserpos, 1:ny) &
        + ly * (bx(laserpos, 1:ny) - bx(laserpos, 0:ny-1)) &
        + dt_eps * jz(laserpos, 1:ny) &
        + diff * by(laserpos-1, 1:ny))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_x_max, lasers=laser_x_max)
      ELSE
        CALL calc_absorption(c_bd_x_max)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_x_max



  SUBROUTINE outflow_bcs_y_min

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_y_min

    laserpos = 1
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_y_min_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(nx))
    ALLOCATE(source2(nx))
    source1 = 0.0_num
    source2 = 0.0_num

    by(1:nx, laserpos-1) = by_y_min(1:nx)

    IF (add_laser(n)) THEN
      current => laser_y_min
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time >= current%t_start .AND. time <= current%t_end) THEN
          IF (current%use_phase_function) CALL laser_update_phase(current)
          IF (current%use_profile_function) CALL laser_update_profile(current)
          t_env = laser_time_profile(current) * current%amp
          DO i = 1,nx
            base = t_env * current%profile(i) &
              * SIN(current%current_integral_phase(i) + current%phase(i) &
              - current%k(i, 1) * (x(i)-current%initial_pos(1)) &
              - current%k(i, 2) * (y(laserpos-1)-current%initial_pos(2)))
            source1(i) = source1(i) + base * COS(current%pol_angle)
            source2(i) = source2(i) + base * SIN(current%pol_angle)
          END DO
        END IF
        current => current%next
      END DO
    END IF

#ifdef BOOSTED_FRAME
    IF (use_boosted_frame) THEN
      source2 = source2 * global_boost_info%lorentz_gamma &
          * (1.0_num - global_boost_info%beta)
    END IF
#endif

    bx(1:nx, laserpos-1) = sum * ( 4.0_num * source1 &
        + 2.0_num * (ez_y_min(1:nx) + c * bx_y_min(1:nx)) &
        - 2.0_num * ez(1:nx, laserpos) &
        - lx * (by(1:nx, laserpos) - by(0:nx-1, laserpos)) &
        + dt_eps * jz(1:nx, laserpos) &
        + diff * bx(1:nx, laserpos))

    bz(1:nx, laserpos-1) = sum * (-4.0_num * source2 &
        - 2.0_num * (ex_y_min(1:nx) - c * bz_y_min(1:nx)) &
        + 2.0_num * ex(1:nx, laserpos) &
        - dt_eps * jx(1:nx, laserpos) &
        + diff * bz(1:nx, laserpos))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_y_min, lasers=laser_y_min)
      ELSE
        CALL calc_absorption(c_bd_y_min)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_y_min



  SUBROUTINE outflow_bcs_y_max

    REAL(num) :: t_env
    REAL(num) :: dtc2, lx, ly, sum, diff, dt_eps, base
    REAL(num), DIMENSION(:), ALLOCATABLE :: source1, source2
    INTEGER :: laserpos, n, i
    TYPE(laser_block), POINTER :: current

    n = c_bd_y_max

    laserpos = ny
    IF (bc_field(n) == c_bc_cpml_laser) THEN
      laserpos = cpml_y_max_laser_idx
    END IF
    dtc2 = dt * c**2
    lx = dtc2 / dx
    ly = dtc2 / dy
    sum = 1.0_num / (ly + c)
    diff = ly - c
    dt_eps = dt / epsilon0

    ALLOCATE(source1(nx))
    ALLOCATE(source2(nx))
    source1 = 0.0_num
    source2 = 0.0_num

    by(1:nx, laserpos+1) = by_y_max(1:nx)

    IF (add_laser(n)) THEN
      current => laser_y_max
      DO WHILE(ASSOCIATED(current))
        ! evaluate the temporal evolution of the laser
        IF (time >= current%t_start .AND. time <= current%t_end) THEN
          IF (current%use_phase_function) CALL laser_update_phase(current)
          IF (current%use_profile_function) CALL laser_update_profile(current)
          t_env = laser_time_profile(current) * current%amp
          DO i = 1,nx
            base = t_env * current%profile(i) &
              * SIN(current%current_integral_phase(i) + current%phase(i) &
              - current%k(i, 1) * (x(i)-current%initial_pos(1)) &
              - current%k(i, 2) * (y(laserpos-1)-current%initial_pos(2)))
            source1(i) = source1(i) + base * COS(current%pol_angle)
            source2(i) = source2(i) + base * SIN(current%pol_angle)
          END DO
        END IF
        current => current%next
      END DO
    END IF

#ifdef BOOSTED_FRAME
    IF (use_boosted_frame) THEN
      source2 = source2 * global_boost_info%lorentz_gamma &
          * (1.0_num - global_boost_info%beta)
    END IF
#endif

    bx(1:nx, laserpos) = sum * (-4.0_num * source1 &
        - 2.0_num * (ez_y_max(1:nx) - c * bx_y_max(1:nx)) &
        + 2.0_num * ez(1:nx, laserpos) &
        + lx * (by(1:nx, laserpos) - by(0:nx-1, laserpos)) &
        - dt_eps * jz(1:nx, laserpos) &
        + diff * bx(1:nx, laserpos-1))

    bz(1:nx, laserpos) = sum * ( 4.0_num * source2 &
        + 2.0_num * (ex_y_max(1:nx) + c * bz_y_max(1:nx)) &
        - 2.0_num * ex(1:nx, laserpos) &
        + dt_eps * jx(1:nx, laserpos) &
        + diff * bz(1:nx, laserpos-1))

    DEALLOCATE(source1, source2)

    IF (dump_absorption) THEN
      IF (add_laser(n)) THEN
        CALL calc_absorption(c_bd_y_max, lasers=laser_y_max)
      ELSE
        CALL calc_absorption(c_bd_y_max)
      END IF
    END IF

  END SUBROUTINE outflow_bcs_y_max



  SUBROUTINE calc_absorption(bd, lasers)

    TYPE(laser_block), POINTER, OPTIONAL :: lasers
    INTEGER, INTENT(IN) :: bd
    TYPE(laser_block), POINTER :: current
    REAL(num) :: t_env, dir, dd, factor, lfactor, laser_inject_sum
    REAL(num), DIMENSION(:), ALLOCATABLE :: e1, e2, b1, b2
    INTEGER :: mm, ibc, icell

    ! Note: ideally e1, e2, b1, b2 should be face-centred. However, this is not
    ! possible with 'open' boundaries since E-fields are not defined in the
    ! ghost cell, so we use the cell-centred quantities in the first cell.

    dir = 1.0_num
    mm = 1

    SELECT CASE(bd)
      CASE(c_bd_x_min, c_bd_x_max)
        dd = dy
        mm = ny
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        ibc = 1
        IF (bd == c_bd_x_max) THEN
          dir = -1.0_num
          ibc = nx
        END IF

        e1 = 0.5_num  * (ey(ibc  , 0:ny-1) + ey(ibc, 1:ny  ))
        e2 = ez(ibc, 1:ny)
        b1 = 0.25_num * (bz(ibc-1, 0:ny-1) + bz(ibc, 0:ny-1) &
                       + bz(ibc-1, 1:ny  ) + bz(ibc, 1:ny  ))
        b2 = 0.5_num  * (by(ibc-1, 1:ny  ) + by(ibc, 1:ny  ))

      CASE(c_bd_y_min, c_bd_y_max)
        dd = dx
        mm = nx
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        ibc = 1
        IF (bd == c_bd_y_max) THEN
          dir = -1.0_num
          ibc = ny
        END IF

        e1 = ez(1:nx, ibc)
        e2 = 0.5_num  * (ex(0:nx-1, ibc  ) + ex(1:nx  , ibc))
        b1 = 0.5_num  * (bx(1:nx  , ibc-1) + bx(1:nx  , ibc))
        b2 = 0.25_num * (bz(0:nx-1, ibc-1) + bz(0:nx-1, ibc) &
                       + bz(1:nx  , ibc-1) + bz(1:nx  , ibc))

      CASE DEFAULT
        dd = 0.0_num
        ALLOCATE(e1(mm), e2(mm), b1(mm), b2(mm))

        e1 = 0.0_num
        e2 = 0.0_num
        b1 = 0.0_num
        b2 = 0.0_num
    END SELECT

    factor = dt * dd * dir
    laser_absorb_local = laser_absorb_local &
        + (factor / mu0) * SUM(e1 * b1 - e2 * b2)

    IF (PRESENT(lasers)) THEN
      current => lasers
      DO WHILE(ASSOCIATED(current))
        laser_inject_sum = 0.0_num
        DO icell = 1, mm
          laser_inject_sum = laser_inject_sum + current%profile(icell)**2
        END DO
        t_env = laser_time_profile(current)
        lfactor = 0.5_num * epsilon0 * c * factor * (t_env * current%amp)**2
        laser_inject_local = laser_inject_local + lfactor * laser_inject_sum
        current => current%next
      END DO
    END IF

    DEALLOCATE(e1, e2, b1, b2)

  END SUBROUTINE calc_absorption

END MODULE laser
