from sympy import Eq
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    convert_to_float,
)

# Description
## Reduced units are used in the dimensionless van der Waals equation of state.

# Law: p* = p / p_c
## p* - reduced pressure
## p - pressure
## p_c - [critical pressure](./critical_pressure.py)

reduced_pressure = Symbol("reduced_pressure", dimensionless)
pressure = Symbol("pressure", units.pressure)
critical_pressure = Symbol("critical_pressure", units.pressure)

law = Eq(reduced_pressure, pressure / critical_pressure)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    pressure_=pressure,
    critical_pressure_=critical_pressure,
)
@validate_output(reduced_pressure)
def calculate_reduced_pressure(
    pressure_: Quantity,
    critical_pressure_: Quantity,
) -> float:
    result = law.rhs.subs({
        pressure: pressure_,
        critical_pressure: critical_pressure_,
    })
    return convert_to_float(result)
