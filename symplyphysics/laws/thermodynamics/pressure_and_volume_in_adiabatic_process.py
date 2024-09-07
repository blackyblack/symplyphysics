"""
Adiabatic process equation via pressure and volume
==================================================

An *adiabatic process* is a type of thermodynamic process that occurs without transferring
heat or mass between the thermodynamic system and its environment.
"""

from sympy import Eq, Rational, solve, Symbol as SymSymbol, Function as SymFunction, dsolve, symbols as sym_symbols
from sympy.abc import t
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol, dimensionless,
    validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import heat_capacity_ratio
from symplyphysics.laws.thermodynamics import (
    internal_energy_change_of_ideal_gas_via_temperature as internal_energy_law,
    internal_energy_change_via_heat_and_work as first_law,
    infinitesimal_work_in_quasistatic_process as work_law,
    isochoric_and_isobaric_heat_capacity_of_ideal_gas as mayers_relation,
)
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

adiabatic_index = Symbol("adiabatic_index", dimensionless)
r"""
Adiabatic index, also known as :doc:`heat capacity ratio <definitions.heat_capacity_ratio>`, of the gas.

Symbol:
    :code:`gamma`

Latex:
    :math:`\gamma`
"""

# Some of the following parameters depend on each other. It is up to user which to choose as known.

initial_temperature = clone_symbol(symbols.thermodynamics.temperature,
    display_symbol="T0",
    display_latex="T_0")
"""
Initial :attr:`~symplyphysics.symbols.thermodynamics.temperature` of the system.
"""

final_temperature = clone_symbol(symbols.thermodynamics.temperature,
    display_symbol="T1",
    display_latex="T_1")
"""
Final :attr:`~symplyphysics.symbols.thermodynamics.temperature` of the system.
"""

initial_volume = Symbol("initial_volume", units.volume)
"""
Initial volume of the system.

Symbol:
    :code:`V0`

Latex:
    :math:`V_0`
"""

final_volume = Symbol("final_volume", units.volume)
"""
Final volume of the system.

Symbol:
    :code:`V1`

Latex:
    :math:`V_1`
"""

initial_pressure = Symbol("initial_pressure", units.pressure)
"""
Initial pressure inside the system.

Symbol:
    :code:`p0`

Latex:
    :math:`p_0`
"""

final_pressure = Symbol("final_pressure", units.pressure)
"""
Final pressure inside the system.

Symbol:
    :code:`p1`

Latex:
    :math:`p_1`
"""

adiabatic_condition = Eq(
    initial_pressure * (initial_volume**adiabatic_index),
    final_pressure * (final_volume**adiabatic_index),
)
r"""
:code:`p0 * V0^gamma = p1 * V1^gamma`

Latex:
    .. math::
        p_0 V_0^\gamma = p_1 V_1^\gamma
"""

eq_start = ideal_gas_equation.law.subs({
    ideal_gas_equation.temperature: initial_temperature,
    ideal_gas_equation.volume: initial_volume,
    ideal_gas_equation.pressure: initial_pressure
})

eq_end = ideal_gas_equation.law.subs({
    ideal_gas_equation.temperature: final_temperature,
    ideal_gas_equation.volume: final_volume,
    ideal_gas_equation.pressure: final_pressure
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

_internal_energy_change_expr = internal_energy_law.law.rhs.subs({
    internal_energy_law.temperature_change: _temperature_change,
    internal_energy_law.isochoric_heat_capacity: _isochoric_heat_capacity,
})

_work_change_expr = work_law.law.rhs.subs({
    work_law.pressure: _pressure,
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

_temperature_change_expr = solve(ideal_gas_equation.law, _temperature)[0].subs({
    _pressure: _parameterized_pressure(t),
    _volume: _parameterized_volume(t),
}).diff(t).subs({
    _parameterized_pressure(t).diff(t): _pressure_change,
    _parameterized_volume(t).diff(t): _volume_change,
    _parameterized_pressure(t): _pressure,
    _parameterized_volume(t): _volume,
})

_mayers_relation = mayers_relation.law.subs(mayers_relation.amount_of_substance,
    ideal_gas_equation.amount_of_substance)

_temperature_change_expr = solve(
    (Eq(_temperature_change, _temperature_change_expr), _mayers_relation),
    (_temperature_change, ideal_gas_equation.amount_of_substance),
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
    }))

_pressure_eqn = dsolve(_diff_eqn,
    _parameterized_pressure(_volume),
    ics={_parameterized_pressure(initial_volume): initial_pressure})

_pressure_derived = _pressure_eqn.rhs.subs(_volume, final_volume)

_pressure_from_law = solve(adiabatic_condition,
    final_pressure)[0].subs(adiabatic_index, _adiabatic_index)

assert expr_equals(_pressure_derived, _pressure_from_law)


@validate_input(mole_count_=ideal_gas_equation.amount_of_substance,
    temperature_start_=initial_temperature,
    volume_start_=initial_volume,
    volume_end_=final_volume,
    specific_heats_ratio_=adiabatic_index)
@validate_output(final_pressure)
def calculate_pressure(mole_count_: Quantity, temperature_start_: Quantity, volume_start_: Quantity,
    volume_end_: Quantity, specific_heats_ratio_: Rational) -> Quantity:

    solved = solve(law, (initial_pressure, final_temperature, final_pressure),
        dict=True)[0][final_pressure]
    result_pressure = solved.subs({
        ideal_gas_equation.amount_of_substance: mole_count_,
        initial_temperature: temperature_start_,
        initial_volume: volume_start_,
        final_volume: volume_end_,
        adiabatic_index: specific_heats_ratio_
    })
    return Quantity(result_pressure)
