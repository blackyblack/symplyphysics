from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
    symbols,
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

internal_energy_change = Function("internal_energy_change", units.energy)
temperature = symbols.thermodynamics.temperature
entropy_change = Symbol("entropy_change", units.energy / units.temperature)
pressure = Symbol("pressure", units.pressure)
volume_change = Symbol("volume_change", units.volume)

law = Eq(
    internal_energy_change(entropy_change, volume_change),
    temperature * entropy_change - pressure * volume_change,
)

# TODO: Derive from the first and second laws of thermodynamics


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
