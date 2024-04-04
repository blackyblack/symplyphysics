from sympy import Eq
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
)
from symplyphysics.core.quantity_decorator import validate_multiple_output

# Description
## Critical parameters of the van der Waals equation of state are such value of volume, pressure, and
## temperature at which the isotherm has an inflection point whose tangent at that point is zero, i.e.
## the first and second derivatives of pressure with respect to volume at constant temperature are zero.

# Law: V_c = 3 * b
#      p_c = a / (27 * b**2)
#      T_c = (8 * a) / (27 * R * b)
## V_c, p_c, T_c - critical volume, pressure, and temperature respectively
## a, b - parameters of van der Waals equation of state
## R - molar gas constant

critical_molar_volume = Symbol("critical_molar_volume", units.volume / units.amount_of_substance)
critical_pressure = Symbol("critical_pressure", units.pressure)
critical_temperature = Symbol("critical_temperature", units.temperature)
bonding_forces_parameter = Symbol(
    "bonding_forces_parameter",
    units.pressure * (units.volume / units.amount_of_substance) ** 2,
)
molecules_volume_parameter = Symbol(
    "molecules_volume_parameter",
    units.volume / units.amount_of_substance,
)

molar_volume_law = Eq(critical_molar_volume, 3 * molecules_volume_parameter)

pressure_law = Eq(
    critical_pressure,
    bonding_forces_parameter / (27 * molecules_volume_parameter**2),
)

temperature_law = Eq(
    critical_temperature,
    (8 * bonding_forces_parameter) / (27 * units.molar_gas_constant * molecules_volume_parameter),
)

law = molar_volume_law, pressure_law, temperature_law


def print_law() -> str:
    return print_expression(law)


@validate_input(
    bonding_forces_parameter_=bonding_forces_parameter,
    molecules_volume_parameter_=molecules_volume_parameter,
)
@validate_multiple_output(
    critical_molar_volume, critical_pressure, critical_temperature
)
def calculate_critical_parameters(
    bonding_forces_parameter_: Quantity,
    molecules_volume_parameter_: Quantity,
) -> tuple[Quantity, Quantity, Quantity]:
    subs_ = {
        bonding_forces_parameter: bonding_forces_parameter_,
        molecules_volume_parameter: molecules_volume_parameter_,
    }
    molar_volume_ = molar_volume_law.rhs.subs(subs_)
    pressure_ = pressure_law.rhs.subs(subs_)
    temperature_ = temperature_law.rhs.subs(subs_)
    return molar_volume_, pressure_, temperature_
