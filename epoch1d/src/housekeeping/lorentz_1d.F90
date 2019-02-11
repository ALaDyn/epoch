MODULE lorentz

  USE shared_data
  USE utilities
  IMPLICIT NONE

  INTERFACE
    FUNCTION transform_position_fn(boost_info, position_val, inverse, &
        time_val, time_boost)
      IMPORT boost_info_object, num
      IMPLICIT NONE
      TYPE(boost_info_object), INTENT(IN) :: boost_info
      REAL(num), INTENT(IN) :: position_val
      LOGICAL, INTENT(IN), OPTIONAL :: inverse
      REAL(num), INTENT(IN), OPTIONAL :: time_val
      REAL(num), INTENT(OUT), OPTIONAL :: time_boost
      REAL(num) :: transform_position_fn
    END FUNCTION transform_position_fn
  END INTERFACE


  CONTAINS

  !>Transform a four_vector subject to a given frame transform
  !>since this is only ever a boost in 1D, only have two components
  !>time like and x component of space like
  FUNCTION transform_four_vector(boost_info, four_vector, inverse)
    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), DIMENSION(2), INTENT(IN) :: four_vector
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num), DIMENSION(2) :: transform_four_vector
    REAL(num) :: beta

    beta = boost_info%beta
    IF (get_optional(default_val = .FALSE., optional_val = inverse)) THEN
      beta = - boost_info%beta
    ELSE
      beta = boost_info%beta
    END IF

    transform_four_vector(1) = boost_info%lorentz_gamma * four_vector(1) &
        - beta * boost_info%lorentz_gamma * four_vector(2)
    transform_four_vector(2) = boost_info%lorentz_gamma * four_vector(2) &
        - beta * boost_info%lorentz_gamma * four_vector(1)

  END FUNCTION transform_four_vector



  FUNCTION transform_interval(boost_info, interval_length, inverse)
    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), INTENT(IN) :: interval_length
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num) :: transform_interval

    IF (get_optional(default_val = .FALSE., optional_val = inverse)) THEN
      !Inverse case
      transform_interval = interval_length / boost_info%lorentz_gamma
    ELSE
      transform_interval = interval_length * boost_info%lorentz_gamma
    END IF

  END FUNCTION transform_interval



  FUNCTION transform_mass(boost_info, mass, inverse)
    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), INTENT(IN) :: mass
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num) :: transform_mass

    IF (get_optional(default_val = .FALSE., optional_val = inverse)) THEN
      !Inverse case
      transform_mass = mass / boost_info%lorentz_gamma
    ELSE
      transform_mass = mass * boost_info%lorentz_gamma
    END IF

  END FUNCTION transform_mass



  FUNCTION transform_length(boost_info, length, inverse)

    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), INTENT(IN) :: length
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num) :: transform_length

    IF (get_optional(default_val = .FALSE., optional_val = inverse)) THEN
      !Inverse case
      transform_length = length * boost_info%lorentz_gamma
    ELSE
      transform_length = length / boost_info%lorentz_gamma
    END IF

  END FUNCTION transform_length



  FUNCTION transform_volume(boost_info, volume, inverse)
    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), INTENT(IN) :: volume
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num) :: transform_volume

    IF (get_optional(default_val = .FALSE., optional_val = inverse)) THEN
      !Inverse case
      transform_volume = volume * boost_info%lorentz_gamma
    ELSE
      transform_volume = volume / boost_info%lorentz_gamma
    END IF

  END FUNCTION transform_volume



  FUNCTION transform_density(boost_info, density, inverse)
    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), INTENT(IN) :: density
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num) :: transform_density

    IF (get_optional(default_val = .FALSE., optional_val = inverse)) THEN
      !Inverse case
      transform_density = density / boost_info%lorentz_gamma
    ELSE
      transform_density = density * boost_info%lorentz_gamma
    END IF

  END FUNCTION transform_density



  !>Transform a point in space at the current simulation time into a position
  !>in a transformed frame
  FUNCTION transform_position(boost_info, position_val, inverse, time_val, &
      time_boost)
    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), INTENT(IN) :: position_val
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num), INTENT(IN), OPTIONAL :: time_val
    REAL(num), INTENT(OUT), OPTIONAL :: time_boost
    REAL(num), DIMENSION(2) :: four_position
    REAL(num) :: transform_position

    four_position(1) = get_optional(optional_val = time_val, &
        default_val = time) * c
    four_position(2) = position_val

    four_position = transform_four_vector(boost_info, &
        four_position, inverse)

    transform_position = four_position(2)

    IF (PRESENT(time_boost)) time_boost = four_position(1) / c

  END FUNCTION transform_position



  !>Transform a point in space at the current simulation time into a position
  !>in a transformed frame
  FUNCTION transform_position_at_prime(boost_info, position_val, inverse, &
      time_val, time_boost)
    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), INTENT(IN) :: position_val
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num), INTENT(IN), OPTIONAL :: time_val
    REAL(num), INTENT(OUT), OPTIONAL :: time_boost
    REAL(num) :: transform_position_at_prime, timeval, beta
    LOGICAL :: do_inverse

    do_inverse = get_optional(optional_val = inverse, default_val = .FALSE.)
    timeval = get_optional(optional_val = time_val, &
        default_val = time) * c

    beta = boost_info%beta
    IF (do_inverse) beta = - beta

    transform_position_at_prime = SQRT(1.0_num - beta**2) * position_val &
        - beta * timeval

    IF (PRESENT(time_boost)) time_boost = timeval

  END FUNCTION transform_position_at_prime



  !>Transform a time at a given position into a time in another frame
  FUNCTION transform_time(boost_info, time_val, position_val, inverse)
    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), INTENT(IN) :: time_val
    REAL(num), INTENT(IN) :: position_val
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num), DIMENSION(2) :: four_position
    REAL(num) :: transform_time

    four_position(1) = time_val * c
    four_position(2) = position_val

    four_position = transform_four_vector(boost_info, &
        four_position, inverse)

    transform_time = four_position(1) / c

  END FUNCTION transform_time



  !>Transform a point in space at the current simulation time into a position
  !>in a transformed frame
  FUNCTION transform_frequency(boost_info, omega, kx, inverse, k_boost)

    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), INTENT(IN) :: omega
    REAL(num), INTENT(IN) :: kx
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num), INTENT(OUT), OPTIONAL :: k_boost
    REAL(num), DIMENSION(2) :: four_frequency
    REAL(num) :: transform_frequency

    four_frequency(1) = omega/c
    four_frequency(2) = kx

    four_frequency = transform_four_vector(boost_info, &
        four_frequency, inverse)

    transform_frequency = four_frequency(1) * c
    IF (PRESENT(k_boost)) k_boost = four_frequency(2)

  END FUNCTION transform_frequency



  !>Lorentz transform for momentum. Must specify either total energy or mass
  !>If both are supplied then total energy is used. If neither is supplied
  !>then the routine will assume a massless, zero energy particle.
  FUNCTION transform_momentum(boost_info, p, inverse, energy, mass)
    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), DIMENSION(3), INTENT(IN) :: p
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num), INTENT(IN), OPTIONAL :: energy, mass
    REAL(num), DIMENSION(3) :: transform_momentum
    REAL(num), DIMENSION(2) :: four_p
    REAL(num) :: energy_c

    IF (PRESENT(energy)) THEN
      energy_c = energy / c
    ELSE IF (PRESENT(mass)) THEN
      energy_c = mass * c * SQRT(1.0_num + (p(1)**2 + p(2)**2 + p(3)**2) &
          / (mass * c)**2)
    ELSE
      energy_c = 0.0_num
    END IF
    four_p(1) = energy_c
    four_p(2) = p(1)

    four_p = transform_four_vector(boost_info, four_p, inverse)

    transform_momentum = [four_p(2), p(2), p(3)]

  END FUNCTION transform_momentum



  !>Transform E and B fields intro a transformed frame
  !>Fields are packed into a 6D array as [ex,ey,ez,bx,by,bz]
  SUBROUTINE transform_em_fields(boost_info, fields_in, fields_out, inverse)
    TYPE(boost_info_object), INTENT(IN) :: boost_info
    REAL(num), DIMENSION(6), INTENT(IN) :: fields_in
    REAL(num), DIMENSION(6), INTENT(OUT) :: fields_out
    LOGICAL, INTENT(IN), OPTIONAL :: inverse
    REAL(num) :: beta, v, lgamma

    beta = boost_info%beta
    IF (get_optional(default_val = .FALSE., optional_val = inverse)) &
        beta = - beta
    v = beta * c
    lgamma = boost_info%lorentz_gamma

    !Transform E
    fields_out(1) = fields_in(1)
    fields_out(2) = lgamma * (fields_in(2) - v * fields_in(6))
    fields_out(3) = lgamma * (fields_in(3) + v * fields_in(5))

    !Transform B
    fields_out(4) = fields_in(4)
    fields_out(5) = lgamma * (fields_in(5) + beta/c * fields_in(3))
    fields_out(6) = lgamma * (fields_in(6) + beta/c * fields_in(2))

  END SUBROUTINE transform_em_fields

END MODULE lorentz
