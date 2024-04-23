from sympy import Eq, solve, Symbol as SymSymbol, Function as SymFunction, dsolve, symbols as sym_symbols
from sympy.abc import t
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol, print_expression,
    dimensionless, validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import heat_capacity_ratio
from symplyphysics.laws.thermodynamics import (
    internal_energy_change_via_amount_of_heat_and_work_done as first_law,
    internal_energy_of_ideal_gas_is_proportional_to_temperature as internal_energy_law,
    infinitesimal_work_in_quasistatic_process as work_law,
    isochoric_heat_capacity_via_isobaric_heat_capacity_for_ideal_gas as mayers_relation,
)
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

# Description
## Adiabatic process: Q = 0, P * V^y = const
## Where:
## Q is amount of transferred heat between system and environment
## P is pressure,
## V is volume,
## y is the ratio of specific heats (also known as heat capacity
##   ratio) (https://en.wikipedia.org/wiki/Heat_capacity_ratio)

specific_heats_ratio = Symbol("specific_heats_ratio", dimensionless)
# Some of these parameters depend on each other. It is up to user, which of these parameters to choose
# as known.
temperature_start = clone_symbol(symbols.thermodynamics.temperature, "temperature_start")
temperature_end = clone_symbol(symbols.thermodynamics.temperature, "temperature_end")
volume_start = Symbol("volume_start", units.volume)
volume_end = Symbol("volume_end", units.volume)
pressure_start = Symbol("pressure_start", units.pressure)
pressure_end = Symbol("pressure_end", units.pressure)

adiabatic_condition = Eq(pressure_start * (volume_start**specific_heats_ratio),
    pressure_end * (volume_end**specific_heats_ratio))

eq_start = ideal_gas_equation.law.subs({
    ideal_gas_equation.temperature: temperature_start,
    ideal_gas_equation.volume: volume_start,
    ideal_gas_equation.pressure: pressure_start
})

eq_end = ideal_gas_equation.law.subs({
    ideal_gas_equation.temperature: temperature_end,
    ideal_gas_equation.volume: volume_end,
    ideal_gas_equation.pressure: pressure_end
})

law = [eq_start, eq_end, adiabatic_condition]

# Derive from first law of thermodynamics and other relations

_temperature = ideal_gas_equation.temperature
_pressure = ideal_gas_equation.pressure
_volume = ideal_gas_equation.volume

_adiabatic_index = Symbol("adiabatic_index", positive=True)
_isobaric_heat_capacity = mayers_relation.isobaric_heat_capacity
_isochoric_heat_capacity = mayers_relation.isochoric_heat_capacity

_temperature_change = SymSymbol("temperature_change")
_pressure_change = SymSymbol("pressure_change")
_volume_change = SymSymbol("volume_change")

# _internal_energy_expr = internal_energy_law.law.rhs.subs({
#     internal_energy_law.temperature_change: _temperature,
#     internal_energy_law.isochoric_heat_capacity: _isochoric_heat_capacity,
# })

# _internal_energy_change_expr = _internal_energy_expr.diff(_temperature) * _temperature_change
_internal_energy_change_expr = internal_energy_law.law.rhs.subs({
    internal_energy_law.temperature_change: _temperature_change,
    internal_energy_law.isochoric_heat_capacity: _isochoric_heat_capacity,
})

_work_change_expr = work_law.law.rhs.subs({
    work_law.pressure_inside_system: _pressure,
    work_law.infinitesimal_volume_change: _volume_change,
})

_first_law_eqn = first_law.law.subs({
    first_law.internal_energy_change: _internal_energy_change_expr,
    first_law.heat_supplied_to_system: 0,
    first_law.work_done_by_system: _work_change_expr,
})

# Expressing the temperature differential through pressure and volume differentials

_parameterized_pressure = sym_symbols("pressure", cls=SymFunction)
_parameterized_volume = sym_symbols("volume", cls=SymFunction)

_temperature_change_expr = solve(
    ideal_gas_equation.law, _temperature
)[0].subs({
    _pressure: _parameterized_pressure(t),
    _volume: _parameterized_volume(t),
}).diff(t).subs({
    _parameterized_pressure(t).diff(t): _pressure_change,
    _parameterized_volume(t).diff(t): _volume_change,
    _parameterized_pressure(t): _pressure,
    _parameterized_volume(t): _volume,
})

_mayers_relation = mayers_relation.law.subs(
    mayers_relation.amount_of_substance, ideal_gas_equation.mole_count
)

_temperature_change_expr = solve(
    (Eq(_temperature_change, _temperature_change_expr), _mayers_relation),
    (_temperature_change, ideal_gas_equation.mole_count),
    dict=True,
)[0][_temperature_change]

_adiabatic_index_eqn = heat_capacity_ratio.definition.subs({
    heat_capacity_ratio.heat_capacity_ratio: _adiabatic_index,
    heat_capacity_ratio.isochoric_heat_capacity: _isochoric_heat_capacity,
    heat_capacity_ratio.isobaric_heat_capacity: _isobaric_heat_capacity,
})

_pressure_change_expr = solve(
    (_first_law_eqn, Eq(_temperature_change, _temperature_change_expr), _adiabatic_index_eqn),
    (_pressure_change, _temperature_change, _isobaric_heat_capacity),
    dict=True,
)[0][_pressure_change]

_diff_eqn = Eq(
    _parameterized_pressure(_volume).diff(_volume),
    _pressure_change_expr.subs({
        _pressure: _parameterized_pressure(_volume),
        _volume_change: 1,
    })
)

_pressure_eqn = dsolve(
    _diff_eqn,
    _parameterized_pressure(_volume),
    ics={_parameterized_pressure(volume_start): pressure_start}
)

_pressure_derived = _pressure_eqn.rhs.subs(_volume, volume_end)

_pressure_from_law = solve(adiabatic_condition, pressure_end)[0].subs(
    specific_heats_ratio, _adiabatic_index
)

assert expr_equals(_pressure_derived, _pressure_from_law)


def print_law() -> str:
    return print_expression(law)


@validate_input(mole_count_=ideal_gas_equation.mole_count,
    temperature_start_=temperature_start,
    volume_start_=volume_start,
    volume_end_=volume_end,
    specific_heats_ratio_=specific_heats_ratio)
@validate_output(pressure_end)
def calculate_pressure(mole_count_: Quantity, temperature_start_: Quantity, volume_start_: Quantity,
    volume_end_: Quantity, specific_heats_ratio_: float) -> Quantity:

    solved = solve(law, (pressure_start, temperature_end, pressure_end), dict=True)[0][pressure_end]
    result_pressure = solved.subs({
        ideal_gas_equation.mole_count: mole_count_,
        temperature_start: temperature_start_,
        volume_start: volume_start_,
        volume_end: volume_end_,
        specific_heats_ratio: specific_heats_ratio_
    })
    return Quantity(result_pressure)
