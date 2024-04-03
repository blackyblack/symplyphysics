from sympy import Eq, Derivative
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
## The compressibility, or the coefficient of compressibility, is a measure of the instantaneous
## relative volume change of a fluid or solid as a response to pressure or mean stress change.

# Definition: beta = -1/V(p) * (dV/dp)
## beta - compressibility
## V - volume (as a function of pressure and other parameters)
## p - pressure
## d/dp - partial derivative with respect to pressure

# Note
## - This definition is incomplete in the sense that the value of the compressibility coefficient
##   depends on whether the process is isentropic or isothermal, hence the partial derivative should
##   be taken at either constant entropy or constant temperature.

compressibility = Symbol("compressibility", 1 / units.pressure)
volume = Function("volume", units.volume)
pressure = Symbol("pressure", units.pressure)

definition = Eq(
    compressibility,
    -1 * Derivative(volume(pressure), pressure) / volume(pressure)
)


def print_law() -> str:
    return print_expression(definition)


@validate_input(
    volume_before_=volume,
    volume_after_=volume,
    pressure_before_=pressure,
    pressure_after_=pressure,
)
@validate_output(compressibility)
def calculate_compressibility(
    volume_before_: Quantity,
    volume_after_: Quantity,
    pressure_before_: Quantity,
    pressure_after_: Quantity,
) -> Quantity:
    # The value of the volume is calculated in the `pressure_after_` point

    volume_function = (
        volume_before_ +
        (volume_after_ - volume_before_) * (pressure - pressure_before_) / (pressure_after_ - pressure_before_)
    )
    expr = definition.rhs.subs(volume(pressure), volume_function).doit()
    result = expr.subs(pressure, pressure_after_)
    return Quantity(result)
