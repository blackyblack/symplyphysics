import numbers
from sympy import (Eq, solve)
from symplyphysics import (angle_type, units, expr_to_quantity, Quantity, Symbol, print_expression,
    validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.definitions import temporal_frequency_is_events_per_time as frequency_def

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

law = Eq(angular_frequency, radians / time)

# Derive the same law from temporal frequency definition

frequency_of_radian = frequency_def.definition.subs({
    frequency_def.events: radians,
    frequency_def.time: time
}).rhs
assert expr_equals(frequency_of_radian, law.rhs)


def print_law() -> str:
    return print_expression(law)


@validate_input(time_=time, radians_=radians)
@validate_output(angular_frequency)
def calculate_frequency(radians_: float | Quantity, time_: Quantity) -> Quantity:
    #HACK: SymPy angles are always in radians
    angle_radians = radians_ if isinstance(radians_, numbers.Number) else radians_.scale_factor
    solved = solve(law, angular_frequency, dict=True)[0][angular_frequency]
    result_expr = solved.subs({time: time_, radians: angle_radians})
    return expr_to_quantity(result_expr)
