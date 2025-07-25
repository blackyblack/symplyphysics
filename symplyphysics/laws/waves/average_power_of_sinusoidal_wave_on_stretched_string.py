"""
Average power of sinusoidal wave on stretched string
====================================================

The average power of a wave of any type is proportional to the square of its amplitude and to the
square of its angular frequency.

**Conditions:**

#. The wave is sinusoidal.

#. The wave is linear, i.e., assuming :math:`s(x, t)` is the transverse displacement:

   #. :math:`\\frac{|s(x, t)|}{L} \\ll 1`, where :math:`L` is the string length

   #. :math:`\\left| \\frac{\\partial s}{\\partial x} \\right| \\ll 1`

   #. the tension is constant thoughout the string

   #. the linear density is constant throughout the string

   #. the displacement must be orthogonal to the string's equilibrium shape

   #. no damping or energy dissipation within the string

**Links:**

#. `Physics LibreTexts, formula 16.5.1 <https://phys.libretexts.org/Bookshelves/University_Physics/University_Physics_(OpenStax)/Book%3A_University_Physics_I_-_Mechanics_Sound_Oscillations_and_Waves_(OpenStax)/16%3A_Waves/16.05%3A_Energy_and_Power_of_a_Wave>`__.
"""

from sympy import Eq, cos, solve, sqrt, Symbol as SymSymbol
from symplyphysics import (Quantity, validate_input, validate_output, symbols, clone_as_symbol,
    clone_as_function)
from symplyphysics.core.expr_comparisons import expr_equals

from symplyphysics.definitions import (
    speed_is_distance_derivative as _speed_def,
    power_is_energy_derivative as _power_def,
    period_from_angular_frequency as _period_law,
    mechanical_energy_is_kinetic_and_potential_energy as _mechanical_energy_def,
)
from symplyphysics.laws.dynamics import (
    kinetic_energy_from_mass_and_speed as _kinetic_energy_def,
    mechanical_work_from_force_and_distance as _work_law,
)
from symplyphysics.laws.quantities import (
    quantity_is_linear_density_times_length as _linear_density_law,)
from symplyphysics.laws.waves import (
    phase_of_traveling_wave as _phase_law,
    wave_equation_general_solution_in_one_dimension as _solution_law,
    phase_speed_of_wave_on_stretched_string as _string_speed_law,
    phase_velocity_from_angular_frequency_and_wavenumber as _phase_speed_law,
)

wave_average_power = symbols.power
"""
Average :symbols:`power`, or rate of energy transfer, of the wave.
"""

string_linear_density = symbols.linear_density
"""
:symbols:`linear_density` of the string.
"""

phase_speed = symbols.phase_speed
"""
:symbols:`phase_speed` of the wave.
"""

wave_angular_frequency = clone_as_symbol(symbols.angular_frequency, positive=True)
"""
:symbols:`angular_frequency` of the wave.
"""

wave_amplitude = clone_as_symbol(
    symbols.euclidean_distance,
    display_symbol="u_max",
    display_latex="u_\\text{max}",
)
"""
Amplitude of the wave. See :symbols:`euclidean_distance`.
"""

law = Eq(
    wave_average_power,
    string_linear_density * phase_speed * wave_angular_frequency**2 * wave_amplitude**2 / 2,
)
"""
:laws:symbol::

:laws:latex::
"""

# Derive law

_position = _phase_law.position
_time = _phase_law.time

# 1. Find the expression for the displacement function of a point on a string

_phase_expr = _phase_law.law.rhs.subs({
    _phase_law.angular_frequency: wave_angular_frequency,
})

_transverse_displacement_expr = _solution_law.law.rhs.replace(
    _solution_law.solution,
    lambda _: wave_amplitude * cos(_phase_expr),
)

# 2. Find the kinetic energy change of a small string element

# This is the speed at which a string element at a specific location moves in the transverse
# direction.
_transverse_speed_expr = _speed_def.definition.rhs.subs(_speed_def.time, _time).replace(
    _speed_def.distance,
    lambda _: _transverse_displacement_expr,
).doit()

# NOTE: we choose a small string element, which is initially horizontal (due to the equilibrium
#  shape)
_small_element_length = clone_as_symbol(symbols.distance, display_symbol="dx", positive=True)

_small_element_mass_expr = _linear_density_law.law.rhs.subs({
    _linear_density_law.linear_density: string_linear_density,
    _linear_density_law.length: _small_element_length,
})

_small_element_kinetic_energy_expr = _kinetic_energy_def.law.rhs.subs({
    _kinetic_energy_def.mass: _small_element_mass_expr,
    _kinetic_energy_def.speed: _transverse_speed_expr,
})

# 2. Find the potential energy change of that small string element

_transverse_displacement = clone_as_function(symbols.distance, (_position, _time), real=True)

# NOTE: `ds = (∂s/∂x) * dx`
_small_transverse_displacement_expr = (_transverse_displacement(_position, _time).diff(_position) *
    _small_element_length)

# A symbol we use to replace `∂s/∂x` with in the expression
_slope = SymSymbol("k")

# The string element was initially horizontal, and after a `dt` time it moves up a `ds` distance.
# Since the lengths are small we can find the length of the arc using the Pythagoras theorem.
_small_arc_length_expr = sqrt(_small_element_length**2 + _small_transverse_displacement_expr**2)

# Per condition 2 we can use the Taylor series expansion w.r.t. `k = ∂s/∂x`
_small_arc_length_expr = _small_arc_length_expr.subs(
    _transverse_displacement(_position, _time).diff(_position),
    _slope,
).series(_slope, 0, 4).removeO().simplify().subs(
    _slope,
    _transverse_displacement(_position, _time).diff(_position),
).replace(
    _transverse_displacement,
    lambda *_: _transverse_displacement_expr,
)

_arc_length_change = (_small_arc_length_expr - _small_element_length).simplify()

# The string element gets stretched out and thus, potential energy is stored in the string.
_tension = _string_speed_law.tension

# The tension force vector and the displacement of the string element are almost parallel.
_potential_energy_change_expr = _work_law.law.rhs.subs({
    _work_law.force: _tension,
    _work_law.distance: _arc_length_change,
})

# 3. Find mechanical energy change

_mechanical_energy_change_expr = _mechanical_energy_def.definition.rhs.subs({
    _mechanical_energy_def.kinetic_energy: _small_element_kinetic_energy_expr,
    _mechanical_energy_def.potential_energy: _potential_energy_change_expr,
}).simplify()

_mechanical_energy_change = clone_as_symbol(
    symbols.mechanical_energy,
    display_symbol="dE",
    real=True,
)

_mechanical_energy_change_eqn = Eq(_mechanical_energy_change, _mechanical_energy_change_expr)

_time_change = clone_as_symbol(symbols.time, display_symbol="dt", real=True)

_power_eqn = _power_def.definition.subs(_power_def.time, _time).subs({
    _power_def.energy(_time).diff(_time): _mechanical_energy_change / _time_change,
})

_speed_eqn = _speed_def.definition.subs(_speed_def.time, _time).subs({
    _speed_def.speed(_time): phase_speed,
    _speed_def.distance(_time).diff(_time): _small_element_length / _time_change,
})

_tension_eqn = _string_speed_law.law.subs({
    _string_speed_law.phase_speed: phase_speed,
    _string_speed_law.linear_density: string_linear_density,
})

_wavenumber_eqn = _phase_speed_law.law.subs({
    _phase_speed_law.phase_speed: phase_speed,
    _phase_speed_law.angular_frequency: wave_angular_frequency,
    _phase_speed_law.angular_wavenumber: _phase_law.angular_wavenumber,
})

_power_expr = solve(
    (_mechanical_energy_change_eqn, _power_eqn, _speed_eqn, _tension_eqn, _wavenumber_eqn),
    (_mechanical_energy_change, _power_def.power(_time), _small_element_length, _tension,
    _phase_law.angular_wavenumber),
    dict=True,
)[0][_power_def.power(_time)]

# 4. Average power over the wave period

_period_expr = _period_law.law.rhs.subs({
    _period_law.angular_frequency: wave_angular_frequency,
})

# Since the string performs a simple harmonic motion (i.e. the motion repeats itself every period)
# we can average it as time ranges from `0` to the oscillation period.
_average_power_expr = _power_expr.integrate((_time, 0, _period_expr)).simplify() / _period_expr

assert expr_equals(_average_power_expr, law.rhs)


@validate_input(
    string_linear_density_=string_linear_density,
    phase_speed_=phase_speed,
    wave_angular_frequency_=wave_angular_frequency,
    wave_amplitude_=wave_amplitude,
)
@validate_output(wave_average_power)
def calculate_average_power(
    string_linear_density_: Quantity,
    phase_speed_: Quantity,
    wave_angular_frequency_: Quantity,
    wave_amplitude_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        string_linear_density: string_linear_density_,
        phase_speed: phase_speed_,
        wave_angular_frequency: wave_angular_frequency_,
        wave_amplitude: wave_amplitude_,
    })
    return Quantity(result)
