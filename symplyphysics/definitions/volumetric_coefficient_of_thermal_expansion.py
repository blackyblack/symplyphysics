from sympy import Eq, Derivative
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
## The coefficient of thermal expansion describes how the size of an object changes with a change in temperature
## at constant pressure.

# Law: alpha_V = (dV(T)/dT)_p / V(T)
## alpha_V - volumetric coefficient of thermal expansion
## V - volume of body (gas, liquid or solid)
## T - temperature
## (d/dT)_p - temperature derivative with pressure held constant

# Conditions
## - Pressure must be constant during the expansion process

volumetric_expansion_coefficient = Symbol("volumetric_expansion_coefficient", 1 / units.temperature)
volume = Function("volume", units.volume)
temperature = symbols.thermodynamics.temperature

law = Eq(
    volumetric_expansion_coefficient,
    Derivative(volume(temperature), temperature) / volume(temperature),
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    volume_before_=volume,
    volume_after_=volume,
    temperature_before_=temperature,
    temperature_after_=temperature,
)
@validate_output(volumetric_expansion_coefficient)
def calculate_volumetric_expansion_coefficient(
    volume_before_: Quantity,
    volume_after_: Quantity,
    temperature_before_: Quantity,
    temperature_after_: Quantity,
) -> Quantity:
    # The RHS of the equation is calculated at the temperature point after expansion (`temperature_after_`)

    volume_function = (
        volume_before_
        + (volume_after_ - volume_before_)
        * (temperature - temperature_before_)
        / (temperature_after_ - temperature_before_)
    )
    result = (
        (law.rhs)
        .subs(volume(temperature), volume_function)
        .doit()
        .subs(temperature, temperature_after_)
    )
    return Quantity(result)
