from sympy import symbols as sym_symbols
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    Vector,
    QuantityVector,
    validate_input,
    validate_output,
    scale_vector,
    vector_magnitude,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals

# Description
## The phase velocity of a wave is the rate at which the wave propagates in a medium.
## It is the velocity at which the phase of one frequency component of the wave travels.
## The phase velocity is collinear with the wavevector.

# Note
## Angular wavevector is a vector used in describing a wave. Its magnitude is the angular
## wavenumber of the wave. Its direction is perpendicular to the wavefront, and in isotropic
## media it is also the direction of wave propagation.

# Law: v = (w / |k|) * (k / |k|)
## v - [phase velocity of wave](../phase_velocity_from_angular_velocity_and_wavevector.py)
## w - angular frequency of wave
## k - angular wavevector of wave
## |k| - [angular wavenumber of wave](../../../definitions/angular_wavenumber_is_inverse_wavelength.py)

# Links: Wikipedia, derivable from here <https://en.wikipedia.org/wiki/Phase_velocity>

angular_frequency = clone_as_symbol(symbols.angular_frequency, positive=True)


def phase_velocity_law(wavevector_: Vector) -> Vector:
    wavenumber_ = vector_magnitude(wavevector_)
    return scale_vector(
        angular_frequency / wavenumber_**2,
        wavevector_,
    )


# k   =   |k| * e_k
#   =(1)= (w / |v|) * e_k
#   =(2)= (w / |v|) * (v / |v|)
#     =   w * v / |v|**2
# where
# (1) |v| = (w / |k|) * (|k| / |k|) = w / |k| => |k| = w / |v|
# (2) v / |v| =(1)= (w / |k|) * (k / |k|) / (w / |k|) = k / |k| = e_k
def wavevector_law(phase_velocity_: Vector) -> Vector:
    phase_speed_ = vector_magnitude(phase_velocity_)
    return scale_vector(
        angular_frequency / phase_speed_**2,
        phase_velocity_,
    )


# Prove the wavevector law for arbitrary v and w
_phase_velocity = Vector(sym_symbols("phase_velocity_x:z"))
_wavevector_derived = wavevector_law(_phase_velocity)
_phase_velocity_derived = phase_velocity_law(_wavevector_derived)
for _component, _derived_component in zip(_phase_velocity.components,
    _phase_velocity_derived.components):
    assert expr_equals(_component, _derived_component)


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


@validate_input(angular_frequency_=angular_frequency, phase_velocity_=units.velocity)
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
