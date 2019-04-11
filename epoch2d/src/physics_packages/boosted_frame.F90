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

  FUNCTION test_valid_packages() RESULT(any_errors)
    LOGICAL :: any_errors

    any_errors = .FALSE.

#ifdef WORK_DONE_INTEGRATED
    any_errors = .TRUE.
    IF (rank == 0) THEN
      PRINT *,'*** ERROR ***'
      PRINT *,'Work done diagnostic does not function in boosted frame'
    END IF
#endif

    IF (move_window) THEN
      any_errors = .TRUE.
      IF (rank == 0) THEN
        PRINT *,'*** ERROR ***'
        PRINT *,'Cannot have a moving window inside a boosted frame'
      END IF
    END IF

  END FUNCTION test_valid_packages



  SUBROUTINE new_simulation_boost

    LOGICAL :: any_errors, ret_err

    IF (use_boosted_frame) THEN

      IF (rank == 0) THEN
        PRINT *,'*******************************************'
        PRINT *,'Transforming to boosted frame'
        PRINT *,'*******************************************'
        PRINT '(A, F16.0)', ' Frame vx(m/s) = ', global_boost_info%vx
        PRINT '(A, F16.8)', ' Frame beta    = ', global_boost_info%beta
        PRINT '(A, F16.8)', ' Frame gamma   = ', global_boost_info%lorentz_gamma
        PRINT *,'*******************************************'
      END IF
      trans_pos => transform_position

      !Do not replace these with any_errors = any_errors .OR. transform_*
      !There is no guaranteed short circuit behaviour in Fortran and this will
      !definitely fail on gfortran
      any_errors = .FALSE.

      ret_err = test_valid_packages()
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'Packages OK'

      IF (rank == 0) &
          PRINT *,'*******************************************'
      ret_err = transform_domain(.FALSE.)
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'Domain OK'
      IF (rank == 0) PRINT *,'t_end in boosted frame is ', t_end

      IF (rank == 0) &
          PRINT *,'*******************************************'
      ret_err = transform_species()
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'Particles OK'

      IF (rank == 0) &
          PRINT *,'*******************************************'
      ret_err = test_injectors()
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'Injectors OK'

      IF (rank == 0) &
          PRINT *,'*******************************************'
      ret_err = transform_lasers()
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'Lasers OK'

      IF (rank == 0) &
          PRINT *,'*******************************************'
      ret_err = transform_io()
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'IO OK'

      IF (any_errors) THEN
        IF (rank == 0) THEN
          PRINT *,'*******************************************'
          PRINT *,'THERE ARE ERRORS IN TRANSFORMING TO THE BOOSTED FRAME'
          PRINT *,'FIX THESE ERRORS AND RERUN EPOCH'
          PRINT *,'CODE WILL TERMINATE'
          PRINT *,'*******************************************'
        END IF
        CALL abort_code(c_err_bad_value)
      ELSE
        IF (rank == 0) THEN
          PRINT *,'*******************************************'
          PRINT *,'Transform to boosted frame succeeded'
          PRINT *,'*******************************************'
        END IF
      END IF
      in_boosted_frame = .TRUE.
    ELSE
      CALL setup_null_io
    END IF

  END SUBROUTINE new_simulation_boost



  SUBROUTINE restart_boost

    LOGICAL :: any_errors, ret_err

    IF (use_boosted_frame) THEN

      !If you are restarting from a non boosted frame output boost like an
      !initial condition
      IF (.NOT. restart_in_boosted_frame) THEN
        CALL new_simulation_boost
        RETURN
      END IF

      IF (rank == 0) THEN
        PRINT *,'*******************************************'
        PRINT *,'Transforming restart to boosted frame'
        PRINT *,'*******************************************'
        PRINT '(A, F16.0)', ' Frame vx(m/s) = ', global_boost_info%vx
        PRINT '(A, F16.8)', ' Frame beta    = ', global_boost_info%beta
        PRINT '(A, F16.8)', ' Frame gamma   = ', global_boost_info%lorentz_gamma
        PRINT *,'*******************************************'
      END IF
      trans_pos => transform_position

      !Do not replace these with any_errors = any_errors .OR. transform_*
      !There is no guaranteed short circuit behaviour in Fortran and this will
      !definitely fail on gfortran
      any_errors = .FALSE.

      ret_err = test_valid_packages()
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'Packages OK'

      IF (rank == 0) &
          PRINT *,'*******************************************'
      ret_err = transform_domain(time_only = .TRUE.)
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'Domain OK'
      IF (rank == 0) PRINT *,'t_end in boosted frame is ', t_end
      IF (rank == 0) &
          PRINT *,'*******************************************'
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'Particles not transformed'
      IF (rank == 0) &
          PRINT *,'*******************************************'
      ret_err = test_injectors()
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'Injectors OK'
      IF (rank == 0) &
          PRINT *,'*******************************************'
      ret_err = transform_lasers()
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'Lasers OK'

      IF (rank == 0) &
          PRINT *,'*******************************************'
      ret_err = transform_io()
      any_errors = any_errors .OR. ret_err
      IF (.NOT. ret_err .AND. rank == 0) PRINT *,'IO OK'

      IF (any_errors) THEN
        IF (rank == 0) THEN
          PRINT *,'*******************************************'
          PRINT *,'THERE ARE ERRORS IN TRANSFORMING TO THE BOOSTED FRAME'
          PRINT *,'FIX THESE ERRORS AND RERUN EPOCH'
          PRINT *,'CODE WILL TERMINATE'
          PRINT *,'*******************************************'
        END IF
        CALL abort_code(c_err_bad_value)
      ELSE
        IF (rank == 0) THEN
          PRINT *,'*******************************************'
          PRINT *,'Transform to boosted frame succeeded'
          PRINT *,'*******************************************'
        END IF
      END IF
      in_boosted_frame = .TRUE.
    ELSE
      CALL setup_null_io
    END IF

  END SUBROUTINE restart_boost



  FUNCTION test_injectors() RESULT(any_errors)
    LOGICAL :: any_errors
    TYPE(injector_block), POINTER :: current

    any_errors = .FALSE.

    IF (x_min_boundary) THEN
      current => injector_x_min
      DO WHILE(ASSOCIATED(current))
        IF (current%drift_function(1)%is_time_varying) THEN
          IF (rank == 0) THEN
            PRINT *,'*** ERROR ***'
            PRINT *,'Cannot have time varying momenta in x direction'
            PRINT *,'for injectors in boosted frame'
          END IF
          any_errors = .TRUE.
        END IF
        current => current%next
      END DO
    END IF

    IF (x_max_boundary) THEN
      current => injector_x_max
      DO WHILE(ASSOCIATED(current))
        IF (current%drift_function(1)%is_time_varying) THEN
          IF (rank == 0) THEN
            PRINT *,'*** ERROR ***'       
            PRINT *,'Cannot have time varying momenta in x direction'
            PRINT *,'for injectors in boosted frame'
          END IF
          any_errors = .TRUE.
        END IF
        current => current%next
      END DO
    END IF

  END FUNCTION test_injectors



  FUNCTION transform_domain(time_only) RESULT(any_errors)

    LOGICAL, INTENT(IN) :: time_only
    INTEGER :: ix
    REAL(num) :: time_min
    LOGICAL :: any_errors

    any_errors = .FALSE.

    !Transform time
    dt_lab = dt
    dt = transform_interval_in_prime(global_boost_info, dt)

    t_end = transform_interval_in_prime(global_boost_info, t_end)
    x_min_deck = trans_pos(global_boost_info, x_min_deck, time_val = 0.0_num)

    IF (.NOT. time_only) THEN

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
        xb_offset_global(ix) = trans_pos(global_boost_info, &
            xb_offset_global(ix))
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

    END IF

    t_end = t_end - global_boost_info%beta * x_min_deck / c

    !Set dt domain to be the time difference in the boosted frame between the
    !bottom and top of the domain at constant lab time
    dt_domain = global_boost_info%lorentz_gamma * global_boost_info%beta &
        * length_x / c

    IF (rank == 0) THEN
      PRINT '(A, E16.7)', ' Lab time difference between ends of domain(s)   = '&
          , dt_domain
    END IF

    dx = x(1) - x(0)

  END FUNCTION transform_domain



  FUNCTION get_staggered_field(field, ix, iy, stagger_origin, &
      stagger_destination)
    REAL(num), DIMENSION(1-ng:,1-ng:), INTENT(IN) :: field
    INTEGER, INTENT(IN) :: ix, iy
    INTEGER, DIMENSION(c_ndims), INTENT(IN) :: stagger_origin, &
        stagger_destination
    INTEGER, DIMENSION(c_ndims) :: dstagger
    REAL(num) :: get_staggered_field, field_accum
    REAL(num), DIMENSION(-1:1,-1:1) :: use_cell

    IF (ALL(stagger_origin == stagger_destination)) THEN
      get_staggered_field = field(ix,iy)
      RETURN
    END IF

    use_cell = 1.0_num
    dstagger = stagger_destination - stagger_origin
    IF (dstagger(1) == 1) THEN
      use_cell(-1,:) = 0.0_num
    ELSE IF (dstagger(1) == 0) THEN
      use_cell(-1,:) = 0.0_num
      use_cell( 1,:) = 0.0_num
    ELSE IF (dstagger(1) == -1) THEN
      use_cell( 1,:) = 0.0_num
    END IF

    IF (dstagger(2) == 1) THEN
      use_cell(:,-1) = 0.0_num
    ELSE IF (dstagger(2) == 0) THEN
      use_cell(:,-1) = 0.0_num
      use_cell(:, 1) = 0.0_num
    ELSE IF (dstagger(2) == -1) THEN
      use_cell(:, 1) = 0.0_num
    END IF

    get_staggered_field = SUM(field(ix-1:ix+1,iy-1:iy+1)*use_cell) &
        / SUM(use_cell)

  END FUNCTION get_staggered_field



  FUNCTION get_staggered_fields(ix, iy, stagger_destination)
    INTEGER, INTENT(IN) :: ix, iy
    INTEGER, DIMENSION(c_ndims), INTENT(IN) :: stagger_destination
    REAL(num), DIMENSION(6) :: get_staggered_fields

    get_staggered_fields(1) = get_staggered_field(ex, ix, iy, [1,0], &
        stagger_destination)
    get_staggered_fields(2) = get_staggered_field(ey, ix, iy, [0,1], &
        stagger_destination)
    get_staggered_fields(3) = get_staggered_field(ez, ix, iy, [0,0], &
        stagger_destination)

    get_staggered_fields(4) = get_staggered_field(bx, ix, iy, [0,1], &
        stagger_destination)
    get_staggered_fields(5) = get_staggered_field(by, ix, iy, [1,0], &
        stagger_destination)
    get_staggered_fields(6) = get_staggered_field(bz, ix, iy, [1,1], &
        stagger_destination)
  END FUNCTION get_staggered_fields



  FUNCTION transform_fields() RESULT(any_errors)

    INTEGER :: ix, iy
    REAL(num) :: stagger_ex, stagger_ey, stagger_ez, stagger_bx, stagger_by, &
        stagger_bz
    REAL(num), DIMENSION(6) :: fields_in, fields_out
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: exn, eyn, ezn, bxn, byn, bzn
    LOGICAL :: any_errors

    any_errors = .FALSE.

    ALLOCATE(exn(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(eyn(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(ezn(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(bxn(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(byn(1-ng:nx+ng, 1-ng:ny+ng))
    ALLOCATE(bzn(1-ng:nx+ng, 1-ng:ny+ng))

    DO iy = 1, ny
      DO ix = 1, nx
        fields_in = get_staggered_fields(ix, iy, [1, 0])
        CALL transform_em_fields(global_boost_info, fields_in, fields_out)
        exn(ix,iy) = fields_out(1)

        fields_in = get_staggered_fields(ix, iy, [0, 1])
        CALL transform_em_fields(global_boost_info, fields_in, fields_out)
        eyn(ix,iy) = fields_out(2)

        fields_in = get_staggered_fields(ix, iy, [0, 0])
        CALL transform_em_fields(global_boost_info, fields_in, fields_out)
        ezn(ix,iy) = fields_out(3)

        fields_in = get_staggered_fields(ix, iy, [0, 1])
        CALL transform_em_fields(global_boost_info, fields_in, fields_out)
        bxn(ix,iy) = fields_out(4)

        fields_in = get_staggered_fields(ix, iy, [1, 0])
        CALL transform_em_fields(global_boost_info, fields_in, fields_out)
        byn(ix,iy) = fields_out(5)

        fields_in = get_staggered_fields(ix, iy, [1, 1])
        CALL transform_em_fields(global_boost_info, fields_in, fields_out)
        bzn(ix,iy) = fields_out(6)
      END DO
    END DO

    DEALLOCATE(ex, ey, ez, bx, by, bz)
    CALL MOVE_ALLOC(exn, ex)
    CALL MOVE_ALLOC(eyn, ey)
    CALL MOVE_ALLOC(ezn, ez)
    CALL MOVE_ALLOC(bxn, bx)
    CALL MOVE_ALLOC(byn, by)
    CALL MOVE_ALLOC(bzn, bz)

  END FUNCTION transform_fields



  FUNCTION transform_species(inverse) RESULT(any_errors)
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species_in
    LOGICAL :: any_errors, retval

    any_errors = .FALSE.

    DO ispecies = 1, n_species
      species_in => species_list(ispecies)
      current=>species_in%attached_list%head
      DO WHILE (ASSOCIATED(current))
        retval = transform_particle(current, species_in%mass, inverse)
        any_errors = any_errors .OR. retval
        current => current%next
      END DO
    END DO

  END FUNCTION transform_species



  FUNCTION transform_particle(current, mass, inverse) RESULT(any_errors)
    TYPE(particle), INTENT(INOUT) :: current
    REAL(num), INTENT(IN) :: mass
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    LOGICAL :: any_errors

    any_errors = .FALSE.

    current%part_p = transform_momentum(global_boost_info, current%part_p, &
        mass = mass)
    current%part_pos(1) = transform_position(global_boost_info, &
        current%part_pos(1))

  END FUNCTION transform_particle



  FUNCTION transform_lasers() RESULT(any_errors)

    TYPE(laser_block), POINTER :: current
    LOGICAL :: any_errors, retval

    any_errors = .FALSE.

    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      retval = transform_laser(current)
      any_errors = any_errors .OR. retval
      current => current%next
    END DO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      retval = transform_laser(current)
      any_errors = any_errors .OR. retval
      current => current%next
    END DO

  END FUNCTION transform_lasers



  FUNCTION transform_laser(laser) RESULT(any_errors)
    TYPE(laser_block), INTENT(INOUT) :: laser
    REAL(num) :: kboost, kx, minkboost
    INTEGER :: i
    LOGICAL :: any_errors

    any_errors = .FALSE.

    laser%initial_pos(1) = trans_pos(global_boost_info, laser%initial_pos(1))

    IF (.NOT. laser%use_omega_function .AND. &
        .NOT. ANY(laser%use_k_function)) THEN

      minkboost = c_largest_number
      DO i = LBOUND(laser%omega,1), UBOUND(laser%omega,1)
        IF (laser%omega_func_type == c_of_k) THEN
          kx = laser%k(i, 1)
        ELSE
          kx = laser%omega(i)/c * laser%kx_mult
        END IF

        laser%omega(i) = transform_frequency(global_boost_info, &
            laser%omega(i), kx, k_boost = kboost)
        IF (kboost < minkboost) minkboost=kboost
        IF (laser%omega_func_type == c_of_k) THEN
          laser%k(i, 1) = kboost
        END IF
      END DO

      laser%t_start = transform_interval_in_prime(global_boost_info, &
          laser%t_start) + time
      laser%t_end = transform_interval_in_prime(global_boost_info, &
          laser%t_end) + time

      IF (pi / MAXVAL(laser%omega) < dt) THEN
        dt = pi / MAXVAL(laser%omega)
        IF (rank == 0) THEN
          PRINT *,'***WARNING***'
          PRINT *,'Timestep reducing to resolve laser after boost'
          PRINT *,'New timestep is ', dt
        END IF
      END IF

      IF (pi / MAX(minkboost, c_tiny) < dx .AND. rank == 0) THEN
          PRINT *,'***ERROR***'
          PRINT *,'After frame transform a laser has a wavelength below the'
          PRINT *,'grid spacing. Results will be incorrect'
          PRINT *,'Code will terminate'
          any_errors = .TRUE.
      END IF

      IF (pi / MAX(minkboost, c_tiny) < 10.0_num * dx .AND. rank == 0) THEN
          PRINT *,'***WARNING***'
          PRINT *,'After frame transform a laser is now in a part of the EM'
          PRINT *,'dispersion curve where result quality is likely to be poor.'
          PRINT *,'Code will continute but results should be treated as suspect'
      END IF
    END IF

  END FUNCTION transform_laser



  SUBROUTINE setup_null_io
    ALLOCATE(prefix_boosts(1:SIZE(file_prefixes)))
  END SUBROUTINE setup_null_io
 


  FUNCTION transform_io() RESULT(any_errors)

    INTEGER :: iio, ispecies, iprefix, nprefix, irec
    REAL(num), DIMENSION(:), ALLOCATABLE :: dt_dump
    INTEGER :: n_recorders, max_n_recorders, total_recorders
    LOGICAL :: any_errors

    nprefix = SIZE(file_prefixes)
    ALLOCATE(prefix_boosts(nprefix))
    ALLOCATE(dt_dump(nprefix))
    max_n_recorders = 0
    total_recorders = 0

    dt_dump = c_largest_number
    any_errors = .FALSE.
    DO iio = 1, n_io_blocks

      !Transform start and end times into the boosted frame
      io_block_list(iio)%time_start = &
          transform_interval_in_prime(global_boost_info, &
          io_block_list(iio)%time_start)
      io_block_list(iio)%time_stop = &
          transform_interval_in_prime(global_boost_info, &
          io_block_list(iio)%time_stop)

      !If using time step between snapshots then transform that
      IF (io_block_list(iio)%dt_snapshot > 0.0_num) THEN
        io_block_list(iio)%dt_snapshot_lab = io_block_list(iio)%dt_snapshot
        io_block_list(iio)%time_prev_lab = io_block_list(iio)%dt_snapshot_lab
        io_block_list(iio)%time_prev = time
        io_block_list(iio)%dt_snapshot = transform_interval_in_prime(&
            global_boost_info, io_block_list(iio)%dt_snapshot)
        IF (io_block_list(iio)%frame == c_frame_lab) THEN
          dt_dump(io_block_list(iio)%prefix_index) = &
              MIN(dt_dump(io_block_list(iio)%prefix_index), &
              io_block_list(iio)%dt_snapshot)
        END IF
      ELSE IF (io_block_list(iio)%frame == c_frame_lab) THEN
        IF (rank == 0) THEN
          PRINT *,'*** ERROR ***'
          PRINT *,'Output in lab frame must be specified using'
          PRINT *,'dt_snapshot'
        END IF
        any_errors = .TRUE.
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
              any_errors = .TRUE.
            ELSE
              prefix_boosts(iprefix)%io_assigned = &
                  (io_block_list(iio)%frame == c_frame_lab)
            END IF
            prefix_boosts(iprefix)%dt_snapshot_lab = &
                io_block_list(iio)%dt_snapshot_lab
            prefix_boosts(iprefix)%frame = io_block_list(iio)%frame
          END IF
        END DO

    END DO

    !Go through all prefixes and set up the recording infrastructure
    !if this prefix is in the lab frame
    DO iprefix = 1, nprefix
      IF (prefix_boosts(iprefix)%frame == c_frame_boost) CYCLE
      n_recorders = MAX(CEILING(dt_domain/dt_dump(iprefix)), 1)
      total_recorders = total_recorders + n_recorders
      IF (max_n_recorders < n_recorders) max_n_recorders = n_recorders
      prefix_boosts(iprefix)%n_recorders = n_recorders
      ALLOCATE(prefix_boosts(iprefix)%field_lists(1-ng:nx+ng, 1-ng:ny+ng, 6, &
          n_recorders))
      ALLOCATE(prefix_boosts(iprefix)%particle_recorders(n_recorders))
      ALLOCATE(prefix_boosts(iprefix)%next_dump(n_recorders))
      prefix_boosts%current_recorder = 1
      DO irec = 1, n_recorders
        ALLOCATE(prefix_boosts(iprefix)%particle_recorders(irec)% &
            particle_lists(1:n_species))
        prefix_boosts(iprefix)%next_dump(irec) &
            = REAL(irec-1, num) * prefix_boosts(iprefix)%dt_snapshot_lab
        DO ispecies = 1, n_species
          prefix_boosts(iprefix)%particle_recorders(irec)%&
              particle_lists(ispecies) = species_list(ispecies)
          prefix_boosts(iprefix)%particle_recorders(irec)%&
              particle_lists(ispecies)%count = 0
          CALL create_empty_partlist(&
              prefix_boosts(iprefix)%particle_recorders(irec)%&
              particle_lists(ispecies)%attached_list, .TRUE.)
        END DO
      END DO
    END DO

    IF (total_recorders == max_n_recorders) THEN
      IF (max_n_recorders == 1) THEN
        IF (rank == 0) THEN
          PRINT *, '*** CAUTION ***'
          PRINT *, 'Output in the lab frame has been requested. This requires'
          PRINT *, 'that EPOCH store copies of particles and fields at constant'
          PRINT *, 'lab time. EPOCH''s memory footprint will be above normal'
        END IF
      ELSE IF (max_n_recorders > 1) THEN
        IF (rank == 0) THEN
          PRINT *, '*** WARNING ***'
          PRINT *, 'Output in the lab frame has been requested with an interval'
          PRINT *, 'lower than the lab frame time difference across the'
          PRINT '(A,I3,A)', ' simulation domain. This will require', &
              max_n_recorders ,' copies of'
          PRINT *, 'fields and particles to be kept and will lead to very high'
          PRINT *, 'memory usage. To avoid this, reduce "dt_snapshot" in '
          PRINT *, 'lab frame diagnostics to be less then', dt_domain 
        END IF
      END IF
    ELSE
      PRINT *, '*** WARNING ***'
      PRINT *, 'Multiple lab frame outputs have been requested.'
      PRINT '(A,I3,A)', ' This will require', &
            total_recorders ,'copies of'
      PRINT *, 'fields and particles to be kept and will lead to very high'
      PRINT *, 'memory usage.'
      IF (max_n_recorders > 1) THEN
        PRINT *, 'Also some requested lab frame outputs have an interval'
        PRINT *, 'lower than the lab frame time difference across the'
        PRINT '(A,I3,A)', ' simulation domain. The most restrictive requires', &
            max_n_recorders ,' copies of'
        PRINT *, 'fields and particles to be kept. To avoid this, reduce '
        PRINT *, '"dt_snapshot" in lab frame diagnostics to be less then', &
            dt_domain 
      END IF
    END IF

    DEALLOCATE(dt_dump)

  END FUNCTION transform_io



  SUBROUTINE record_particle(part, init_x, final_x, ispecies)

    TYPE(particle), INTENT(IN) :: part
    REAL(num), INTENT(IN) :: init_x, final_x
    INTEGER, INTENT(IN) :: ispecies
    REAL(num) :: init_t, final_t, mass
    TYPE(particle), POINTER :: new_part
    INTEGER :: iprefix, irec

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
      DO irec = 1, prefix_boosts(iprefix)%n_recorders
        IF (init_t <= prefix_boosts(iprefix)%next_dump(irec) .AND. final_t &
            >= prefix_boosts(iprefix)%next_dump(irec)) THEN
          ALLOCATE(new_part, source = part)
          !Transform the particle's momentum and energy to lab frame now but
          !keep position in boosted frame. 
#ifdef PHOTONS
          IF (species_list(ispecies)%species_type /= c_species_id_photon) THEN
#endif
            new_part%part_p = transform_momentum(global_boost_info, &
                new_part%part_p, mass = mass, inverse = .TRUE.)
#ifdef PHOTONS
          ELSE
            new_part%part_p = transform_momentum(global_boost_info, &
                new_part%part_p, energy = new_part%energy, inverse = .TRUE.)
          END IF
#endif
          CALL add_particle_to_partlist(prefix_boosts(iprefix)%&
              particle_recorders(irec)%particle_lists(ispecies)%attached_list, &
              new_part)
        END IF
      END DO
    END DO

  END SUBROUTINE record_particle



  SUBROUTINE record_fields

    REAL(num) :: t_local0, t_local1
    REAL(num), DIMENSION(6) :: fields_in, fields_out
    INTEGER :: iprefix, ix, iy, irec

    DO iprefix = 1, SIZE(file_prefixes)
      IF (prefix_boosts(iprefix)%frame == c_frame_boost) CYCLE
      DO ix = 1, nx
        t_local0 = transform_time(global_boost_info, time - dt, x(ix), &
            inverse = .TRUE.)
        t_local1 = transform_time(global_boost_info, time, x(ix), &
            inverse = .TRUE.)
        DO irec = 1, prefix_boosts(iprefix)%n_recorders
          IF (t_local0 <= prefix_boosts(iprefix)%next_dump(irec) &
              .AND. t_local1 >= prefix_boosts(iprefix)%next_dump(irec)) THEN
            DO iy = 1, ny
              fields_in = [ex(ix, iy), ey(ix, iy), ez(ix, iy), &
                  bx(ix, iy), by(ix, iy), bz(ix, iy)]
              CALL transform_em_fields(global_boost_info, fields_in, &
                  fields_out, inverse = .TRUE.)
              prefix_boosts(iprefix)%field_lists(ix, iy, :, irec) = fields_out
              prefix_boosts(iprefix)%field_lists(ix, iy, 6, irec) = 1.0_num
            END DO
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE record_fields

#endif

END MODULE boosted_frame
