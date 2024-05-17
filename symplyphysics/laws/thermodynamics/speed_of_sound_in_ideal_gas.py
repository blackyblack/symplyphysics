from sympy import Eq, solve, sqrt, Function as SymFunction, symbols as sym_symbols
from symplyphysics import Quantity, Symbol, dimensionless, symbols, units, validate_input, validate_output
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import (
    isentropic_speed_of_sound as sound_law,
    zero_heat_transfer as adiabate_law,
)
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation
from symplyphysics.laws.quantities import (
    quantity_is_molar_quantity_times_amount_of_substance as molar_qty_law,
    quantity_is_volumetric_density_times_volume as density_law,
)

# Speed of sound for ideal gases
# c = sqrt( gamma * R * T / M ), where
# gamma is heat capacity ratio (adiabatic index),
# T is temperature,
# M is mass of one mole of this gas.
# R is ideal gas constant,
# c is speed of sound.

temperature = symbols.thermodynamics.temperature
heat_capacity_ratio = Symbol("heat_capacity_ratio", dimensionless)
mole_mass = Symbol("mole_mass", units.mass / units.amount_of_substance)
speed_of_sound = Symbol("speed_of_sound", units.velocity)

law = Eq(speed_of_sound,
    sqrt(heat_capacity_ratio * units.molar_gas_constant * temperature / mole_mass))

# Derive from law of speed of sound and adiabate equation.

_gas_mass = density_law.extensive_quantity
_density = density_law.volumetric_density
_pressure = sym_symbols("pressure", cls=SymFunction)(_density)

_volume_expr = solve(density_law.law, density_law.volume)[0]

_adiabate_eqn = adiabate_law.adiabatic_condition.subs({
    adiabate_law.specific_heats_ratio: heat_capacity_ratio,
    adiabate_law.volume_end: _volume_expr,
    adiabate_law.pressure_end: _pressure,
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
    molar_qty_law.molar_quantity: mole_mass,
    molar_qty_law.amount_of_substance: ideal_gas_equation.mole_count,
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
    mole_mass_=mole_mass)
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
        mole_mass: mole_mass_,
    })
    return Quantity(result_applied)
