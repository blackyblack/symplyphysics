"""
Speed of sound in ideal gas
===========================

Also known as the Laplace's formula of the speed of sound, it provides a correct expression
for the speed of sound for ideal gases compared to the Newton's formula. See
:doc:`laws.thermodynamics.isentropic_speed_of_sound_via_pressure_derivative`.

**Notation:**

#. :quantity_notation:`molar_gas_constant`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Speed_of_sound#Speed_of_sound_in_ideal_gases_and_air>`__.
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    symbols,
    validate_input,
    validate_output,
    quantities,
    clone_as_symbol,
    clone_as_function,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import (
    isentropic_speed_of_sound_via_pressure_derivative as sound_law,
    pressure_and_volume_in_adiabatic_process as adiabate_law,
)
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation
from symplyphysics.laws.quantities import (
    quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law,
    quantity_is_volumetric_density_times_volume as density_law,
)

speed_of_sound = clone_as_symbol(
    symbols.speed,
    display_symbol="v_s",
    display_latex="v_\\text{s}",
)
"""
:symbols:`speed` of sound in gas.
"""

temperature = symbols.temperature
"""
:symbols:`temperature` of the gas.
"""

heat_capacity_ratio = symbols.adiabatic_index
"""
Heat capacity ratio, or :symbols:`adiabatic_index`, of the gas.
"""

molar_mass = symbols.molar_mass
"""
:symbols:`molar_mass` of the gas.
"""

law = Eq(speed_of_sound,
    sqrt(heat_capacity_ratio * quantities.molar_gas_constant * temperature / molar_mass))
"""
:laws:symbol::

:laws:latex::
"""

# Derive from law of speed of sound and adiabate equation.

_gas_mass = density_law.extensive_quantity
_density = density_law.volumetric_density
_pressure = clone_as_function(symbols.pressure)(_density)

_volume_expr = solve(density_law.law, density_law.volume)[0]

_adiabate_eqn = adiabate_law.adiabatic_condition.subs({
    adiabate_law.adiabatic_index: heat_capacity_ratio,
    adiabate_law.final_volume: _volume_expr,
    adiabate_law.final_pressure: _pressure,
})

_pressure_derivative_eqn = Eq(
    _adiabate_eqn.lhs.diff(_density),
    _adiabate_eqn.rhs.diff(_density),
)

_ideal_gas_pressure = solve(ideal_gas_equation.law, ideal_gas_equation.pressure)[0].subs({
    ideal_gas_equation.volume: _volume_expr,
})

_pressure_derivative_expr = solve(_pressure_derivative_eqn, _pressure.diff(_density))[0].subs({
    _pressure: _ideal_gas_pressure,
})

_molar_qty_eqn = molar_qty_law.law.subs({
    molar_qty_law.extensive_quantity: _gas_mass,
    molar_qty_law.molar_quantity: molar_mass,
    molar_qty_law.amount_of_substance: ideal_gas_equation.amount_of_substance,
})

_pressure_derivative_expr = solve(
    (Eq(_pressure.diff(_density), _pressure_derivative_expr), _molar_qty_eqn),
    (_pressure.diff(_density), _gas_mass),
    dict=True,
)[0][_pressure.diff(_density)]

_speed_of_sound_expr = sound_law.law.rhs.subs(
    sound_law.pressure(sound_law.density, sound_law.entropy).diff(sound_law.density),
    _pressure_derivative_expr,
)

assert expr_equals(_speed_of_sound_expr, law.rhs)


@validate_input(temperature_=temperature,
    heat_capacity_ratio_=heat_capacity_ratio,
    mole_mass_=molar_mass)
@validate_output(speed_of_sound)
def calculate_speed_of_sound(
    temperature_: Quantity,
    heat_capacity_ratio_: float,
    mole_mass_: Quantity,
) -> Quantity:
    result_expr = solve(law, speed_of_sound, dict=True)[0][speed_of_sound]
    result_applied = result_expr.subs({
        temperature: temperature_,
        heat_capacity_ratio: heat_capacity_ratio_,
        molar_mass: mole_mass_,
    })
    return Quantity(result_applied)
