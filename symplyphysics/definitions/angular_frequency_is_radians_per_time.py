import numbers
from sympy import (Eq, solve)
from symplyphysics import (angle_type, units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input_symbols, validate_output_symbol)

# Description
## Angular frequency "ω" (also referred to by the terms angular speed and angular rate) is a scalar measure of the angular displacement per unit time.

# Law: w = N / t
# Where:
## N is radians per time,
## t is time,
## w is angular frequency.

radians = Symbol("radians", angle_type)
time = Symbol("period", units.time)
# This is equivalent to frequency, but radians per second is a common notation for angular frequency
angular_frequency = Symbol("angular_frequency", angle_type / units.time)

definition = Eq(angular_frequency, radians / time)

definition_units_SI = units.radian / units.second


def print() -> str:
    return print_expression(definition)


@validate_input_symbols(time_=time, radians_=radians)
@validate_output_symbol(angular_frequency)
def calculate_frequency(radians_: float | Quantity, time_: Quantity) -> Quantity:
    #HACK: SymPy angles are always in radians
    angle_radians = radians_ if isinstance(radians_, numbers.Number) else radians_.scale_factor
    solved = solve(definition, angular_frequency, dict=True)[0][angular_frequency]
    result_expr = solved.subs({time: time_, radians: angle_radians})
    return expr_to_quantity(result_expr)
