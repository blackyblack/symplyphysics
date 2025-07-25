"""
Intensity of sound wave via displacement amplitude
==================================================

The intensity of a sound wave is the rate per unit area of energy transfer
through or onto a surface. It depends on the density of the medium, the phase
speed and the angular frequency of the wave and the amplitude of particles
in the medium.

**Conditions:**

#. The potential energy and the kinetic energy have the same average rate of propagation in the
   medium.

**Links:**

#. Equation 17-27 on p. 489 of "Fundamentals of Physics" by David Halladay et al., 10th Ed.
"""

from sympy import Eq, solve, cos, Expr
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol,
    clone_as_function)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.experimental.solvers import apply

from symplyphysics.definitions import (
    intensity_of_sound_wave_is_rate_of_energy_transfer_over_area as _intensity_def,
    density_from_mass_volume as _density_def,
    speed_is_distance_derivative as _speed_def,
    power_is_energy_derivative as _power_def,
    period_from_angular_frequency as _period_def,
    mechanical_energy_is_kinetic_and_potential_energy as _mechanical_energy_def,
)
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_speed as _kinetic_energy_law
from symplyphysics.laws.waves import (
    phase_of_traveling_wave as _phase_law,
    wave_equation_general_solution_in_one_dimension as _solution_law,
)

wave_intensity = symbols.intensity
"""
:symbols:`intensity` of the sound wave.
"""

medium_density = symbols.density
"""
:symbols:`density` of the medium in which the sound wave is being propagated.
"""

phase_speed = symbols.phase_speed
"""
:symbols:`phase_speed` of the wave.
"""

angular_frequency = clone_as_symbol(symbols.angular_frequency, real=True, zero=False)
"""
:symbols:`angular_frequency` of the wave.
"""

displacement_amplitude = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="s_max",
    display_latex="s_\\text{max}",
)
"""
Displacement amplitude of the particles in the medium. See :symbols:`euclidean_distance`.
"""

law = Eq(
    wave_intensity,
    medium_density * phase_speed * angular_frequency**2 * displacement_amplitude**2 / 2,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

# Consider a thin slice of medium, oscillating as the sound wave passes though it

# 1. Find the mass of the thin slice

_small_thickness = clone_as_symbol(symbols.thickness, display_symbol="dx", positive=True)
_area = _intensity_def.area

_small_volume_expr = _small_thickness * _area

_small_mass_expr = solve(_density_def.definition, _density_def.mass)[0].subs({
    _density_def.density: medium_density,
    _density_def.volume: _small_volume_expr,
})

# 2. Construct the expression for the displacement of the medium's particles

_time = _phase_law.time

_phase_expr = _phase_law.law.rhs.subs({
    _phase_law.angular_frequency: angular_frequency,
})

_phase_shift = clone_as_symbol(symbols.phase_shift)

# `_solution_law.solution` is an arbitrary single-valued function, so we can replace it with a
# sinusoidal function as long as it depends on `_solution_law.wave_phase`. Then we can substitute
# the wave phase symbol with the expression for it, which we've already found.
_displacement_expr = _solution_law.law.rhs.replace(
    _solution_law.solution,
    lambda arg: displacement_amplitude * cos(arg + _phase_shift),
).subs({
    _solution_law.wave_phase: _phase_expr,
})

# 3. Find the speed of an oscillating element of the medium

_speed_expr = _speed_def.definition.rhs.subs(_speed_def.time, _time).replace(
    _speed_def.distance,
    lambda _: _displacement_expr,
).doit()

# 4. Find the kinetic energy of the thin slice. Note that since the slice is thin, the speed can be
#    considered constant throughout it.

_small_kinetic_energy_expr = _kinetic_energy_law.law.rhs.subs({
    _kinetic_energy_law.mass: _small_mass_expr,
    _kinetic_energy_law.speed: _speed_expr,
})

# 4. Find the (local) power of the sound wave due to kinetic energy

_kinetic_energy_change = clone_as_symbol(symbols.kinetic_energy, display_symbol="dK", positive=True)
_time_change = clone_as_symbol(symbols.time, display_symbol="dt", real=True)

_original_power_eqn = _power_def.definition.subs(_power_def.time, _time)

_power_eqn = _original_power_eqn.subs({
    _power_def.energy(_time).diff(_time): _kinetic_energy_change / _time_change,
})

_speed_eqn = _speed_def.definition.subs(_speed_def.time, _time).subs({
    _speed_def.speed(_time): phase_speed,
    _speed_def.distance(_time).diff(_time): _small_thickness / _time_change,
})

_power = _power_def.power

_kinetic_power_expr = solve(
    (_power_eqn, _speed_eqn, Eq(_kinetic_energy_change, _small_kinetic_energy_expr)),
    (_power(_time), _small_thickness, _kinetic_energy_change),
    dict=True,
)[0][_power(_time)]

# 5. Average the (kinetic) power over the wave period

_period = _period_def.period


def _period_average(expr: Expr) -> Expr:
    return expr.integrate((_time, 0, _period)) / _period


_period_expr = _period_def.law.rhs.subs({
    _period_def.angular_frequency: angular_frequency,
})

_average_kinetic_power_expr = _period_average(_kinetic_power_expr).subs(_period, _period_expr)

# 6. Find the average mechanical energy

# Show that `avg(P_mechanical) = avg(P_kinetic) + avg(P_potential)`

_mechanical_energy = clone_as_function(symbols.mechanical_energy, (_time,), positive=True)
_potential_energy = clone_as_function(symbols.potential_energy, (_time,), positive=True)
_kinetic_energy = clone_as_function(symbols.kinetic_energy, (_time,), positive=True)

_mechanical_energy_eqn = _mechanical_energy_def.definition.subs({
    _mechanical_energy_def.mechanical_energy: _mechanical_energy(_time),
    _mechanical_energy_def.potential_energy: _potential_energy(_time),
    _mechanical_energy_def.kinetic_energy: _kinetic_energy(_time),
})

_mechanical_energy_applied_eqn = apply(
    _mechanical_energy_eqn,
    lambda x: _period_average(x.diff(_time)),
).doit()

_average_power_expr = _period_average(_original_power_eqn.rhs)

_average_mechanical_power = clone_as_symbol(symbols.power, display_symbol="avg(P)", real=True)
_average_kinetic_power = clone_as_symbol(symbols.power, display_symbol="avg(P_k)", real=True)
_average_potential_power = clone_as_symbol(symbols.power, display_symbol="avg(P_p)", real=True)
_energy = _power_def.energy

_average_mechanical_power_expr = solve(
    (
    _mechanical_energy_applied_eqn,
    Eq(_average_mechanical_power, _average_power_expr.replace(_energy, _mechanical_energy)),
    Eq(_average_potential_power, _average_power_expr.replace(_energy, _potential_energy)),
    Eq(_average_kinetic_power, _average_power_expr.replace(_energy, _kinetic_energy)),
    ),
    (
    _average_mechanical_power,
    _mechanical_energy(0),
    _potential_energy(0),
    _kinetic_energy(0),
    ),
    dict=True,
)[0][_average_mechanical_power]

# NOTE: we assume that the potential energy is carried along the wave at the same *average* rate
# as the kinetic energy.
_average_mechanical_power_expr = _average_mechanical_power_expr.subs({
    _average_potential_power: _average_kinetic_power,
}).subs({
    _average_kinetic_power: _average_kinetic_power_expr,
})

# 7. Use the definition of intensity

_intensity_expr = _intensity_def.definition.rhs.subs({
    _intensity_def.power: _average_mechanical_power_expr,
})

assert expr_equals(_intensity_expr, law.rhs), _intensity_expr


@validate_input(
    density_=medium_density,
    phase_speed_=phase_speed,
    angular_frequency_=angular_frequency,
    displacement_amplitude_=displacement_amplitude,
)
@validate_output(wave_intensity)
def calculate_intensity(
    density_: Quantity,
    phase_speed_: Quantity,
    angular_frequency_: Quantity,
    displacement_amplitude_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        medium_density: density_,
        phase_speed: phase_speed_,
        angular_frequency: angular_frequency_,
        displacement_amplitude: displacement_amplitude_,
    })
    return Quantity(result)
