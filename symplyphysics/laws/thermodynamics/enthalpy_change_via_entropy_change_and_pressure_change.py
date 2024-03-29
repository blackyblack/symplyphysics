from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    symbols,
)

# Description
## The fundamental thermodynamic relation are fundamental equations which demonstate how important
## thermodynamic quantities depend on variables that are measurable experimentally.

# Law: dH = T * dS + V * dp
## H - [enthalpy](./enthalpy_is_internal_energy_plus_pressure_energy.py)
## T - absolute temperature
## S - entropy
## p - pressure
## V - volume
## Notation: `d(x)` is an exact differential of `x`

# Note
## - Entropy `S` and pressure `p` are so called natural variables of enthalpy as a
##   thermodynamic potential

# Conditions
## - The system is closed and homogeneous
## - The system is in thermal equilibrium with the environment
## - Only reversible prossesses or pure heat transfer are considered

enthalpy_change = Symbol("enthalpy_change", units.energy)
temperature = symbols.thermodynamics.temperature
entropy_change = Symbol("entropy_change", units.energy / units.temperature)
volume = Symbol("volume", units.volume)
pressure_change = Symbol("pressure_change", units.pressure)

law = Eq(
    enthalpy_change,
    temperature * entropy_change + volume * pressure_change,
)

# TODO: Derive from definition of enthalpy and internal energy differential


def print_law() -> str:
    return print_expression(law)


@validate_input(
    temperature_=temperature,
    entropy_change_=entropy_change,
    volume_=volume,
    pressure_change_=pressure_change,
)
@validate_output(enthalpy_change)
def calculate_enthalpy_change(
    temperature_: Quantity,
    entropy_change_: Quantity,
    volume_: Quantity,
    pressure_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        temperature: temperature_,
        entropy_change: entropy_change_,
        volume: volume_,
        pressure_change: pressure_change_,
    })
    return Quantity(result)
