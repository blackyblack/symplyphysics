"""
Adiabatic process equation via pressure and volume
==================================================

An *adiabatic process* is a type of thermodynamic process that occurs without transferring
heat or mass between the thermodynamic system and its environment.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Isentropic_process#Derivation_of_the_isentropic_relations>`__.

..
    TODO rename file
"""

from sympy import Eq, Rational, solve, dsolve
from sympy.abc import t as _t
from symplyphysics import (
    clone_as_symbol,
    clone_as_function,
    symbols,
    Quantity,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import heat_capacity_ratio
from symplyphysics.laws.thermodynamics import (
    internal_energy_change_of_ideal_gas_via_temperature as internal_energy_law,
    internal_energy_change_via_heat_and_work as first_law,
    infinitesimal_work_in_quasistatic_process as work_law,
    isochoric_and_isobaric_heat_capacity_of_ideal_gas as mayers_relation,
)
from symplyphysics.laws.thermodynamics.equations_of_state import ideal_gas_equation

adiabatic_index = symbols.adiabatic_index
"""
:symbols:`adiabatic_index`, also known as :doc:`heat capacity ratio <definitions.heat_capacity_ratio>`, of the gas.
"""

# Some of the following parameters depend on each other. It is up to user which to choose as known.

initial_temperature = clone_as_symbol(symbols.temperature, subscript="0")
"""
Initial :symbols:`temperature` of the system.
"""

final_temperature = clone_as_symbol(symbols.temperature, subscript="1")
"""
Final :symbols:`temperature` of the system.
"""

initial_volume = clone_as_symbol(symbols.volume, subscript="0")
"""
Initial :symbols:`volume` of the system.
"""

final_volume = clone_as_symbol(symbols.volume, subscript="1")
"""
Final :symbols:`volume` of the system.
"""

initial_pressure = clone_as_symbol(symbols.pressure, subscript="0")
"""
Initial :symbols:`pressure` inside the system.
"""

final_pressure = clone_as_symbol(symbols.pressure, subscript="1")
"""
Final :symbols:`pressure` inside the system.
"""

adiabatic_condition = Eq(
    initial_pressure * (initial_volume**adiabatic_index),
    final_pressure * (final_volume**adiabatic_index),
)
"""
:laws:symbol::

:laws:latex::
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

_adiabatic_index = clone_as_symbol(symbols.adiabatic_index, positive=True)
_isobaric_heat_capacity = mayers_relation.isobaric_heat_capacity
_isochoric_heat_capacity = mayers_relation.isochoric_heat_capacity

_temperature_change = clone_as_symbol(symbols.temperature)
_pressure_change = clone_as_symbol(symbols.pressure)
_volume_change = clone_as_symbol(symbols.volume)

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

_parameterized_pressure = clone_as_function(symbols.pressure)
_parameterized_volume = clone_as_function(symbols.volume)

_temperature_change_expr = solve(ideal_gas_equation.law, _temperature)[0].subs({
    _pressure: _parameterized_pressure(_t),
    _volume: _parameterized_volume(_t),
}).diff(_t).subs({
    _parameterized_pressure(_t).diff(_t): _pressure_change,
    _parameterized_volume(_t).diff(_t): _volume_change,
    _parameterized_pressure(_t): _pressure,
    _parameterized_volume(_t): _volume,
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
