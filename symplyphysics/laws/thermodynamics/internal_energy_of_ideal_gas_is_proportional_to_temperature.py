from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The internal energy of an ideal gas depends solely on its temperature and the number of gas particles
## and is independent of other thermodynamic quantities such as pressure or density.

# Law: dU = C_V * dT
## U - internal energy
## C_V - isochoric heat capacity, which only depends on the temperature of the ideal gas
## T - temperature
## Notation: d(x) - exact differential of `x`

internal_energy_change = Symbol("internal_energy_change", units.energy)
isochoric_heat_capacity = Symbol("isochoric_heat_capacity", units.energy / units.temperature)
temperature_change = Symbol("temperature_change", units.temperature)

law = Eq(internal_energy_change, isochoric_heat_capacity * temperature_change)

# TODO: derive law from Joule-Thompson effect


def print_law() -> str:
    return print_expression(law)


@validate_input(
    isochoric_heat_capacity_=isochoric_heat_capacity,
    temperature_change_=temperature_change,
)
@validate_output(internal_energy_change)
def calculate_internal_energy(
    isochoric_heat_capacity_: Quantity,
    temperature_change_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        isochoric_heat_capacity: isochoric_heat_capacity_,
        temperature_change: temperature_change_,
    })
    return Quantity(result)
