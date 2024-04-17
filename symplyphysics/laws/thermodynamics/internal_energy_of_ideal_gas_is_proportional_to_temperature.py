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
## The internal energy of an ideal gas depends solely on its temperature and the number of gas particles
## and is independent of other thermodynamic quantities such as pressure or density.

# Law: U = C_V * T
## U - internal energy
## C_V - isochoric heat capacity
## T - temperature

internal_energy = Symbol("internal_energy", units.energy)
isochoric_heat_capacity = Symbol("isochoric_heat_capacity", units.energy / units.temperature)
temperature = symbols.thermodynamics.temperature

law = Eq(internal_energy, isochoric_heat_capacity * temperature)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    isochoric_heat_capacity_=isochoric_heat_capacity,
    temperature_=temperature,
)
@validate_output(internal_energy)
def calculate_internal_energy(
    isochoric_heat_capacity_: Quantity,
    temperature_: Quantity,
) -> Quantity:
    result = law.rhs.subs({
        isochoric_heat_capacity: isochoric_heat_capacity_,
        temperature: temperature_,
    })
    return Quantity(result)
