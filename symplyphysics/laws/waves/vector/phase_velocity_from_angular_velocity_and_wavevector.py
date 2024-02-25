from symplyphysics import (
    units,
    angle_type,
    Symbol,
    Quantity,
    Vector,
    QuantityVector,
    validate_input,
    validate_output,
    scale_vector,
    vector_magnitude,
)

# Description
## The phase velocity of a wave is the rate at which the wave propagates in a medium.
## It is the velocity at which the phase of one frequency component of the wave travels.
## The phase velocity is collinear with the wavevector.

# Law: v = w * k / |k|**2
## v - vector of phase velocity of wave
## w - angular frequency of wave
## k - wavevector of wave
## |k| - wavenumber of wave, where |v| is the Euclidean norm of vector v

angular_frequency = Symbol("angular_frequency", angle_type / units.time, positive=True)


def phase_velocity_law(wavevector_: Vector) -> Vector:
    wavenumber_ = vector_magnitude(wavevector_)
    return scale_vector(
        angular_frequency / wavenumber_**2,
        wavevector_,
    )


# k = |k| * e_k = (w / |v|) * (v / |v|) = w * v / |v|**2
def wavevector_law(phase_velocity_: Vector) -> Vector:
    phase_speed_ = vector_magnitude(phase_velocity_)
    return scale_vector(
        angular_frequency / phase_speed_**2,
        phase_velocity_,
    )


@validate_input(
    angular_frequency_=angular_frequency,
    wavevector_=angle_type / units.length,
)
@validate_output(units.velocity)
def calculate_phase_velocity(
    angular_frequency_: Quantity,
    wavevector_: QuantityVector,
) -> QuantityVector:
    result = phase_velocity_law(wavevector_.to_base_vector())
    return QuantityVector.from_base_vector(
        result,
        subs={angular_frequency: angular_frequency_},
    )


@validate_input(
    angular_frequency_=angular_frequency,
    phase_velocity_=units.velocity
)
@validate_output(angle_type / units.length)
def calculate_wavevector(
    angular_frequency_: Quantity,
    phase_velocity_: QuantityVector,
) -> QuantityVector:
    result = wavevector_law(phase_velocity_.to_base_vector())
    return QuantityVector.from_base_vector(
        result,
        subs={angular_frequency: angular_frequency_},
    )
