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

# Law: T* = T / T_c
## T* - reduced temperature
## T - temperature
## T_c - critical temperature

reduced_temperature = Symbol("reduced_temperature", dimensionless)
temperature = Symbol("temperature", units.temperature)
critical_temperature = Symbol("critical_temperature", units.temperature)

law = Eq(reduced_temperature, temperature / critical_temperature)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    temperature_=temperature,
    critical_temperature_=critical_temperature,
)
@validate_output(reduced_temperature)
def calculate_reduced_temperature(
    temperature_: Quantity,
    critical_temperature_: Quantity,
) -> float:
    result = law.rhs.subs({
        temperature: temperature_,
        critical_temperature: critical_temperature_,
    })
    return convert_to_float(result)
