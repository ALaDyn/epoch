MODULE boosted_frame

  USE shared_data
  USE utilities
  USE lorentz
  USE partlist
  IMPLICIT NONE

#ifdef BOOSTED_FRAME
  PRIVATE
  PUBLIC :: new_simulation_boost, restart_boost
  PUBLIC :: record_particle, record_fields

  PROCEDURE(transform_position_fn), POINTER :: trans_pos

  REAL(num) :: dt_lab
  CONTAINS

!  FUNCTION test_valid_packages() RESULT(all_ok)
!    LOGICAL :: all_ok

!    all_ok = .TRUE.

#ifdef WORK_DONE_INTEGRATED
!    all_ok = .FALSE.
!    IF (rank == 0) THEN
!      PRINT *,'Work done diagnostic does not function in boosted frame'
!    END IF
#endif

!  END FUNCTION test_valid_packages



  SUBROUTINE new_simulation_boost

    IF (use_boosted_frame) THEN
      trans_pos => transform_position

      CALL transform_domain
      CALL transform_species
      CALL transform_lasers
      CALL transform_io

      in_boosted_frame = .TRUE.
    ELSE
      CALL setup_null_io
    END IF

  END SUBROUTINE new_simulation_boost



  SUBROUTINE restart_boost

  END SUBROUTINE restart_boost



  SUBROUTINE transform_domain

    INTEGER :: ix
    REAL(num) :: time_min

    !Transform time
    dt_lab = dt
    dt = dt * SQRT(1.0_num - global_boost_info%beta**2)

    t_end = t_end * SQRT(1.0_num - global_boost_info%beta**2)

    x_grid_min = trans_pos(global_boost_info, x_grid_min)
    x_grid_max = trans_pos(global_boost_info, x_grid_max)
    x_min = trans_pos(global_boost_info, x_min, time_boost = time_min)
    x_max = trans_pos(global_boost_info, x_max)

    length_x = x_max - x_min

    !Transform local space domain
    x_grid_min_local = trans_pos(global_boost_info, x_grid_min_local)
    x_grid_max_local = trans_pos(global_boost_info, x_grid_max_local)
    x_min_local = trans_pos(global_boost_info, x_min_local)
    x_max_local = trans_pos(global_boost_info, x_max_local)

    !Transform the global domain.
    DO ix = LBOUND(x_global, 1), UBOUND(x_global, 1)
      x_global(ix) = trans_pos(global_boost_info, x_global(ix))
    END DO
    DO ix = LBOUND(xb_global, 1), UBOUND(xb_global, 1)
      xb_global(ix) = trans_pos(global_boost_info, xb_global(ix))
      xb_offset_global(ix) = trans_pos(global_boost_info, xb_offset_global(ix))
    END DO

    !Transform the local domain
    DO ix = LBOUND(x, 1), UBOUND(x, 1)
      x(ix) = trans_pos(global_boost_info, x(ix))
    END DO

    DO ix = LBOUND(xb, 1), UBOUND(xb, 1)
      xb(ix) = trans_pos(global_boost_info, xb(ix))
    END DO

    !Transform the start time in the boosted frame so that there are no
    !negative times in the lab frame
    time = - global_boost_info%beta * x_min / c
    t_end = t_end + time

    !Set dt domain to be the time difference in the boosted frame between the
    !bottom and top of the domain at constant lab time
    dt_domain = global_boost_info%lorentz_gamma * global_boost_info%beta &
        * length_x / c

    dx = x(1) - x(0)

  END SUBROUTINE transform_domain



  SUBROUTINE transform_species(inverse)
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species_in

    DO ispecies = 1, n_species
      species_in => species_list(ispecies)
      current=>species_in%attached_list%head
      DO WHILE (ASSOCIATED(current))
        CALL transform_particle(current, species_in%mass, &
            inverse)
        current => current%next
      END DO
    END DO

  END SUBROUTINE transform_species



  SUBROUTINE transform_particle(current, mass, inverse)
    TYPE(particle), INTENT(INOUT) :: current
    REAL(num), INTENT(IN) :: mass
    LOGICAL, INTENT(IN), OPTIONAL :: inverse

    current%part_p = transform_momentum(global_boost_info, current%part_p, &
        mass = mass)
    current%part_pos = transform_position(global_boost_info, current%part_pos)

  END SUBROUTINE transform_particle



  SUBROUTINE transform_lasers

    TYPE(laser_block), POINTER :: current

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      CALL transform_laser(current)
      current => current%next
    END DO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      CALL transform_laser(current)
      current => current%next
    END DO

  END SUBROUTINE transform_lasers



  SUBROUTINE transform_laser(laser)
    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num) :: kboost

    IF (.NOT. laser%use_omega_function) THEN
      laser%omega = transform_frequency(global_boost_info, laser%omega, &
          laser%omega/c * laser%kx_mult, k_boost = kboost)
      laser%t_start = transform_time(global_boost_info, laser%t_start, &
          laser%x_at_t0)
      laser%t_end = transform_time(global_boost_info, laser%t_end, &
          laser%x_at_t0)

      IF (pi / laser%omega < dt) THEN
        dt = pi / laser%omega
        IF (rank == 0) THEN
          PRINT *,'***WARNING***'
          PRINT *,'Timestep reducing to resolve laser after boost'
          PRINT *,'New timestep is ', dt
        END IF
      END IF

      IF (pi / kboost < dx .AND. rank == 0) THEN
          PRINT *,'***ERROR***'
          PRINT *,'After frame transform a laser has a wavelength below the'
          PRINT *,'grid spacing. Results will be incorrect'
          PRINT *,'Code will terminate'
          CALL abort_code(c_err_bad_setup)
      END IF

      IF (pi / kboost < 10.0_num * dx .AND. rank == 0) THEN
          PRINT *,'***WARNING***'
          PRINT *,'After frame transform a laser is now in a part of the EM'
          PRINT *,'dispersion curve where result quality is likely to be poor.'
          PRINT *,'Code will continute but results should be treated as suspect'
      END IF
    END IF

  END SUBROUTINE transform_laser



  SUBROUTINE setup_null_io
    ALLOCATE(prefix_boosts(1:SIZE(file_prefixes)))
  END SUBROUTINE setup_null_io
 


  SUBROUTINE transform_io

    INTEGER :: iio, ispecies, iprefix, nprefix
    REAL(num) :: dt_dump
    LOGICAL :: error

    nprefix = SIZE(file_prefixes)
    ALLOCATE(prefix_boosts(nprefix))

    dt_dump = HUGE(1.0_num)
    error = .FALSE.
    DO iio = 1, n_io_blocks

      !Transform start and end times into the boosted frame
      io_block_list(iio)%time_start = io_block_list(iio)%time_start &
          * SQRT(1.0_num - global_boost_info%beta**2)
      io_block_list(iio)%time_stop = io_block_list(iio)%time_stop &
          * SQRT(1.0_num - global_boost_info%beta**2)

      !If using time step between snapshots then transform that
      IF (io_block_list(iio)%dt_snapshot > 0.0_num) THEN
        io_block_list(iio)%dt_snapshot_lab = io_block_list(iio)%dt_snapshot
        io_block_list(iio)%time_prev_lab = io_block_list(iio)%dt_snapshot_lab
        io_block_list(iio)%time_prev = time
        io_block_list(iio)%dt_snapshot = io_block_list(iio)%dt_snapshot &
            * SQRT(1.0_num - global_boost_info%beta**2)
        IF (io_block_list(iio)%frame == c_frame_lab) THEN
          dt_dump = MIN(dt_dump, io_block_list(iio)%dt_snapshot)
        END IF
      ELSE IF (io_block_list(iio)%frame == c_frame_lab) THEN
        IF (rank == 0) THEN
          PRINT *,'*** ERROR ***'
          PRINT *,'Output at lab frame time must be specified using'
          PRINT *,'dt_snapshot'
        END IF
        error = .TRUE.
      END IF

      IF (io_block_list(iio)%frame == c_frame_lab) THEN
        !If this output is in the lab frame then current should be output
        !from particles directly rather than using the core current

        io_block_list(iio)%dumpmask(c_dump_jx) = &
            IEOR(io_block_list(iio)%dumpmask(c_dump_jx), c_io_field)
        io_block_list(iio)%dumpmask(c_dump_jy) = &
            IEOR(io_block_list(iio)%dumpmask(c_dump_jy), c_io_field)
        io_block_list(iio)%dumpmask(c_dump_jz) = &
            IEOR(io_block_list(iio)%dumpmask(c_dump_jz), c_io_field)
      END IF

        !Check which prefix this block is associated with and set it up
        !Remember that only one block per prefix is allowed in the lab frame
        DO iprefix = 1, nprefix
          IF (io_block_list(iio)%prefix_index == iprefix) THEN
            IF (prefix_boosts(iprefix)%io_assigned) THEN
              IF (rank == 0) THEN
                PRINT *, '*** ERROR ***'
                PRINT *, 'Can only have one IO block per file prefix in'
                PRINT *, 'lab frame in a boosted frame simulation'
              END IF
              error = .TRUE.
            ELSE
              prefix_boosts(iprefix)%io_assigned = &
                  (io_block_list(iio)%frame == c_frame_lab)
            END IF
            prefix_boosts(iprefix)%next_dump = &
                io_block_list(iio)%dt_snapshot_lab + dt_domain
            prefix_boosts(iprefix)%frame = io_block_list(iio)%frame
          END IF
        END DO

    END DO

    !Go through all prefixes and set up the recording infrastructure
    !if this prefix is in the lab frame
    DO iprefix = 1, nprefix
      IF (prefix_boosts(iprefix)%frame == c_frame_boost) CYCLE
      ALLOCATE(prefix_boosts(iprefix)%particle_lists(1:n_species))
      ALLOCATE(prefix_boosts(iprefix)%field_lists(1-ng:nx+ng, 6))
      DO ispecies = 1, n_species
        prefix_boosts(iprefix)%particle_lists(ispecies) = species_list(ispecies)
        prefix_boosts(iprefix)%particle_lists(ispecies)%count = 0
        CALL create_empty_partlist(&
            prefix_boosts(iprefix)%particle_lists(ispecies)% &
            attached_list, .TRUE.)
      END DO
    END DO

    IF (dt_dump < dt_domain) THEN
      PRINT *, '*** ERROR ***'
      PRINT *, 'Output timescales smaller than a single domain time detected in'
      PRINT *, 'a lab frame output. This will not work correctly' 
      PRINT *, 'Code will abort'
      CALL abort_code(c_err_bad_setup)
    END IF

  END SUBROUTINE transform_io



  SUBROUTINE record_particle(part, init_x, final_x, ispecies)

    TYPE(particle), INTENT(IN) :: part
    REAL(num), INTENT(IN) :: init_x, final_x
    INTEGER, INTENT(IN) :: ispecies
    REAL(num) :: init_t, final_t, mass
    TYPE(particle), POINTER :: new_part
    INTEGER :: iprefix

    init_t = transform_time(global_boost_info, time, init_x, &
        inverse = .TRUE.)
    final_t = transform_time(global_boost_info, time + dt, final_x, &
        inverse = .TRUE.)
#ifndef PER_PARTICLE_CHARGE_MASS
    mass = species_list(ispecies)%mass
#else
    mass = part%mass
#endif

    DO iprefix = 1, SIZE(file_prefixes)
      IF (prefix_boosts(iprefix)%frame == c_frame_boost) CYCLE
      IF (init_t <= prefix_boosts(iprefix)%next_dump .AND. final_t &
          >= prefix_boosts(iprefix)%next_dump) THEN
        ALLOCATE(new_part, SOURCE = part)
        !Transform the particle's momentum and energy to lab frame now but
        !keep position in boosted frame. 
        new_part%part_p = transform_momentum(global_boost_info, &
            new_part%part_p, mass = mass, inverse = .TRUE.)
        CALL add_particle_to_partlist(prefix_boosts(iprefix)%&
            particle_lists(ispecies)%attached_list, new_part)
        prefix_boosts(iprefix)%particle_lists(ispecies)%count = &
            prefix_boosts(iprefix)%particle_lists(ispecies)%count + 1
      END IF
    END DO

  END SUBROUTINE record_particle



  SUBROUTINE record_fields

    REAL(num) :: t_local0, t_local1
    REAL(num), DIMENSION(6) :: fields_in, fields_out
    INTEGER :: iprefix, ix

    DO iprefix = 1, SIZE(file_prefixes)
      IF (prefix_boosts(iprefix)%frame == c_frame_boost) CYCLE
      DO ix = 1, nx
        t_local0 = transform_time(global_boost_info, time - dt, x(ix), &
            inverse = .TRUE.)
        t_local1 = transform_time(global_boost_info, time, x(ix), &
            inverse = .TRUE.)
        IF (t_local0 <= prefix_boosts(iprefix)%next_dump .AND. t_local1 &
            >= prefix_boosts(iprefix)%next_dump) THEN
          fields_in = [ex(ix), ey(ix), ez(ix), bx(ix), by(ix), bz(ix)]
          CALL transform_em_fields(global_boost_info, fields_in, &
              fields_out, inverse = .TRUE.)
          prefix_boosts(iprefix)%field_lists(ix, :) = fields_out
        END IF
      END DO
    END DO

  END SUBROUTINE record_fields

#endif

END MODULE boosted_frame
