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
## Critical parameters of the van der Waals equation of state are such value of volume, pressure, and
## temperature at which the isotherm has an inflection point whose tangent at that point is zero, i.e.
## the first and second derivatives of pressure with respect to volume at constant temperature are zero.

# Law: p_c = a / (27 * b**2)
## p_c - critical pressure
## a, b - parameters of van der Waals equation of state

critical_pressure = Symbol("critical_pressure", units.pressure)

bonding_forces_parameter = Symbol(
    "bonding_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance) ** 2,
)

molecules_volume_parameter = Symbol(
    "molecules_volume_parameter",
    units.volume / units.amount_of_substance,
)

law = Eq(
    critical_pressure,
    bonding_forces_parameter / (27 * molecules_volume_parameter**2),
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    bonding_forces_parameter_=bonding_forces_parameter,
    molecules_volume_parameter_=molecules_volume_parameter,
)
@validate_output(critical_pressure)
def calculate_critical_pressure(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
) -> Quantity:
    pressure_ = law.rhs.subs({
        bonding_forces_parameter: bonding_forces_parameter_,
        molecules_volume_parameter: molecules_volume_parameter_,
    })
    return Quantity(pressure_)
