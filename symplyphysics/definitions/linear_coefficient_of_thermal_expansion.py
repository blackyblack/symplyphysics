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
from symplyphysics.core.geometry.line import two_point_function, Point2D

# Description
## The coefficient of thermal expansion describes how the size of an object changes with a change in temperature
## at constant pressure.

# Law: alpha_L = (dL(T)/dT)_p / L(T)
## alpha_L - linear coefficient of thermal expansion
## L - body length
## T - temperature
## (d/dT)_p - temperature derivative with pressure held constant

# Conditions
## - Pressure must be constant during the expansion process

linear_expansion_coefficient = Symbol("linear_expansion_coefficient", 1 / units.temperature)
length = Function("length", units.length)
temperature = symbols.thermodynamics.temperature

definition = Eq(
    linear_expansion_coefficient,
    Derivative(length(temperature), temperature) / length(temperature),
)


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    length_before_=length,
    length_after_=length,
    temperature_before_=temperature,
    temperature_after_=temperature,
)
@validate_output(linear_expansion_coefficient)
def calculate_linear_expansion_coefficient(
    length_before_: Quantity,
    length_after_: Quantity,
    temperature_before_: Quantity,
    temperature_after_: Quantity,
) -> Quantity:
    # The RHS of the equation is calculated at the temperature point after expansion (`temperature_after_`)

    length_function = two_point_function(
        Point2D(temperature_before_, length_before_),
        Point2D(temperature_after_, length_after_),
        temperature,
    )
    result = ((definition.rhs).subs(length(temperature),
        length_function).doit().subs(temperature, temperature_after_))
    return Quantity(result)
