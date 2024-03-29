from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import (
    internal_energy_change_via_amount_of_heat_and_work_done as first_law,
    entropy_increment_in_reversible_process as second_law,
    infinitesimal_work_in_quasistatic_process as work_law,
)

# Description
## The fundamental thermodynamic relation are fundamental equations which demonstate how important
## thermodynamic quantities depend on variables that are measurable experimentally.

# Law: dU = T * dS - P * dV
## U - internal energy
## T - absolute temperature
## S - entropy
## P - pressure
## V - volume
## Notation: d(x) - full differential of `x`

# Note
## - Entropy `S` and volume `V` are so called natural variables of internal energy as a
##   thermodynamic potential.

# Conditions
## - The system is closed
## - The system is in thermal equilibrium with its surroundings

internal_energy_change = Symbol("internal_energy_change", units.energy)
temperature = symbols.thermodynamics.temperature
entropy_change = Symbol("entropy_change", units.energy / units.temperature)
pressure = Symbol("pressure", units.pressure)
volume_change = Symbol("volume_change", units.volume)

law = Eq(
    internal_energy_change,
    temperature * entropy_change - pressure * volume_change,
)

# Derive from the first and second laws of thermodynamics

_heat_supplied_to_system = solve(
    second_law.law, second_law.infinitesimal_transfer_of_heat
)[0].subs({
    second_law.infinitesimal_entropy_change: entropy_change,
    second_law.common_temperature: temperature,
})

_work_done_by_system = work_law.law.rhs.subs({
    work_law.pressure_inside_system: pressure,
    work_law.infinitesimal_volume_change: volume_change,
})

_internal_energy_change = first_law.law.rhs.subs({
    first_law.heat_supplied_to_system: _heat_supplied_to_system,
    first_law.work_done_by_system: _work_done_by_system,
})

assert expr_equals(law.rhs, _internal_energy_change)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    temperature_=temperature,
    entropy_change_=entropy_change,
    pressure_=pressure,
    volume_change_=volume_change,
)
@validate_output(internal_energy_change)
def calculate_internal_energy_change(
    temperature_: Quantity,
    entropy_change_: Quantity,
    pressure_: Quantity,
    volume_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        entropy_change: entropy_change_,
        pressure: pressure_,
        volume_change: volume_change_,
    })
    return Quantity(result)
