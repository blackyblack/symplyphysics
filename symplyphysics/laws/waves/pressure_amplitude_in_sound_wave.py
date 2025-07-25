"""
Pressure amplitude in sound wave
================================

Sound waves cause a pressure change of the medium from the equilibrium pressure.
This change is proportional to the speed of sound and the density of the medium,
the angular frequency of the wave and the displacement of particles in the medium.

**Links:**

#. Equation 17-14 on p. 484 of "Fundamentals of Physics" by David Halladay et al., 10th Ed.
"""

from sympy import Eq, cos, sin, Wild, Q, solve
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol,
    clone_as_function)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.laws.hydro import bulk_stress_is_bulk_modulus_times_strain as _bulk_law
from symplyphysics.laws.quantities import (
    fractional_change_is_change_over_initial_value as _fractional_law,)
from symplyphysics.laws.waves import (
    phase_of_traveling_wave as _phase_law,
    wave_equation_general_solution_in_one_dimension as _solution_law,
    speed_of_sound_via_bulk_modulus_and_density as _sound_speed_law,
    phase_velocity_from_angular_frequency_and_wavenumber as _phase_speed_law,
)

pressure_amplitude = clone_as_symbol(
    symbols.pressure,
    display_symbol="Delta(p)_max",
    display_latex="(\\Delta p)_\\text{max}",
)
"""
Amplitude of :symbols:`pressure` change.
"""

speed_of_sound = symbols.speed
"""
:symbols:`speed` of sound in the medium.
"""

medium_density = symbols.density
"""
:symbols:`density` of the medium.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the sound wave.
"""

displacement_amplitude = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="s_max",
    display_latex="s_\\text{max}",
)
"""
Displacement amplitude of particles in the medium. See :symbols:`euclidean_distance`.
"""

law = Eq(
    pressure_amplitude,
    speed_of_sound * medium_density * angular_frequency * displacement_amplitude,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

_position = _phase_law.position
_time = _phase_law.time
_displacement = clone_as_function(symbols.distance, (_position, _time), real=True)

# 1. Find the expression for the displacement of air particles as a function of position and time.

_phase_expr = _phase_law.law.rhs.subs({
    _phase_law.angular_frequency: angular_frequency,
})

_phase_shift = symbols.phase_shift

_displacement_expr = _solution_law.law.rhs.replace(
    _solution_law.solution,
    lambda arg: displacement_amplitude * cos(arg + _phase_shift),
).subs({
    _solution_law.wave_phase: _phase_expr,
})

# 2. Consider a small element of air of thickness `dx` that gets displaced by `ds` units in the
# direction of sound wave propagation. Also see [figure](https://www.researchgate.net/publication/336555466/figure/fig1/AS:814253246271488@1571144510196/Sound-wave-propagation.jpg).

_small_element_length = clone_as_symbol(symbols.length, display_symbol="dx", positive=True)

# NOTE: `ds = (∂s/∂x) * dx`
_small_element_displacement_expr = (_displacement(_position, _time).diff(_position) *
    _small_element_length)

_area = symbols.area

_small_element_volume_expr = _area * _small_element_length

# NOTE: The area stays the same; the air element only gets displaced in the direction of wave
# propagation.
_small_element_volume_change_expr = _area * _small_element_displacement_expr

_fractional_volume_change_expr = _fractional_law.law.rhs.subs({
    _fractional_law.change: _small_element_volume_change_expr,
    _fractional_law.initial_value: _small_element_volume_expr,
})

_pressure_expr = _bulk_law.law.rhs.subs({
    _bulk_law.fractional_volume_change: _fractional_volume_change_expr,
}).subs({
    _displacement(_position, _time): _displacement_expr,
}).doit()

# 3. Find the amplitude of the pressure wave

_wild_amplitude = Wild("a")
_wild_phase = Wild("b")

# We "guess" the correct form of the expression, and if our guess is incorrect the `.match` method
# would return a `None`.
_matched_expr = _pressure_expr.match(_wild_amplitude * sin(_wild_phase))[_wild_amplitude]

_pressure_amplitude_expr = abs(_matched_expr).refine(
    Q.positive(_bulk_law.bulk_modulus) & Q.positive(displacement_amplitude) &
    Q.positive(_phase_law.angular_wavenumber))

# 4. Perform final substitutions.

_pressure_amplitude_eqn = Eq(pressure_amplitude, _pressure_amplitude_expr)

_sound_speed_eqn = _sound_speed_law.law.subs({
    _sound_speed_law.phase_velocity: speed_of_sound,
    _sound_speed_law.bulk_modulus: _bulk_law.bulk_modulus,
    _sound_speed_law.density: medium_density,
})

_phase_speed_eqn = _phase_speed_law.law.subs({
    _phase_speed_law.phase_speed: speed_of_sound,
    _phase_speed_law.angular_frequency: angular_frequency,
    _phase_speed_law.angular_wavenumber: _phase_law.angular_wavenumber,
})

_pressure_amplitude_derived = solve(
    (_pressure_amplitude_eqn, _sound_speed_eqn, _phase_speed_eqn),
    (pressure_amplitude, _bulk_law.bulk_modulus, _phase_law.angular_wavenumber),
    dict=True,
)[0][pressure_amplitude]

assert expr_equals(_pressure_amplitude_derived, law.rhs)


@validate_input(
    speed_of_sound_=speed_of_sound,
    density_of_medium_=medium_density,
    angular_frequency_=angular_frequency,
    displacement_amplitude_=displacement_amplitude,
)
@validate_output(pressure_amplitude)
def calculate_pressure_amplitude(
    speed_of_sound_: Quantity,
    density_of_medium_: Quantity,
    angular_frequency_: Quantity,
    displacement_amplitude_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        speed_of_sound: speed_of_sound_,
        medium_density: density_of_medium_,
        angular_frequency: angular_frequency_,
        displacement_amplitude: displacement_amplitude_,
    })
    return Quantity(result)
