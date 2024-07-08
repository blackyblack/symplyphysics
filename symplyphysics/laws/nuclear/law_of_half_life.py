from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    convert_to_float,
    dimensionless,
)

# Description
## Half-life is the time required for a quantity (of substance) to reduce to half of its initial value.
## The term is commonly used in nuclear physics to describe how quickly unstable atoms undergo radioactive
## decay or how long stable atoms survive. The term is also used more generally to characterize any type
## of exponential (or, rarely, non-exponential) decay.

## Law is: N = N0 * 2^(-t / T), where
## N - number of cores after a period of time,
## N0 - initial number of cores,
## t - an arbitrary time interval,
## T - half-life of the nuclei.

number_of_cores = Symbol("number_of_cores", dimensionless)

number_of_cores_initial = Symbol("number_of_cores_initial", dimensionless)
half_life = Symbol("half_life", units.time)
decay_time = Symbol("decay_time", units.time)

law = Eq(number_of_cores, number_of_cores_initial * 2**(-decay_time / half_life))


def print_law() -> str:
    return print_expression(law)


@validate_input(number_of_cores_initial_=number_of_cores_initial,
    half_life_=half_life,
    decay_time_=decay_time)
@validate_output(number_of_cores)
def calculate_number_of_cores(number_of_cores_initial_: int, half_life_: Quantity,
    decay_time_: Quantity) -> int:
    if number_of_cores_initial_ < 0:
        raise ValueError("Number of cores cannot be negative")
    result_expr = solve(law, number_of_cores, dict=True)[0][number_of_cores]
    result_expr = result_expr.subs({
        number_of_cores_initial: number_of_cores_initial_,
        half_life: half_life_,
        decay_time: decay_time_
    })
    return int(convert_to_float(result_expr))
