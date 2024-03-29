from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The fundamental thermodynamic relation are fundamental equations which demonstate how important
## thermodynamic quantities depend on variables that are measurable experimentally.

# Law: dF = -S * dT - p * dV
## F - Helmholtz free energy
## S - entropy
## T - absolute temperature
## p - pressure
## V - volume
## Notation: d(x) - full differential of `x`

# Note
## Temperature `T` and volume `V` are so called natural variables of Helmholtz free energy
## as a thermodynamic potential.

# Conditions
## - The system is closed
## - The system is in thermal equilibrium with its surroundings

free_energy_change = Function("free_energy_change", units.energy)
entropy = Symbol("entropy", units.energy / units.temperature)
temperature_change = Symbol("temperature_change", units.temperature)
pressure = Symbol("pressure", units.pressure)
volume_change = Symbol("volume_change", units.volume)

law = Eq(
    free_energy_change(temperature_change, volume_change),
    -1 * entropy * temperature_change - pressure * volume_change,
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    entropy_=entropy,
    temperature_change_=temperature_change,
    pressure_=pressure,
    volume_change_=volume_change,
)
@validate_output(free_energy_change)
def calculate_free_energy_change(
    entropy_: Quantity,
    temperature_change_: Quantity,
    pressure_: Quantity,
    volume_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        entropy: entropy_,
        temperature_change: temperature_change_,
        pressure: pressure_,
        volume_change: volume_change_,
    })
    return Quantity(result)
