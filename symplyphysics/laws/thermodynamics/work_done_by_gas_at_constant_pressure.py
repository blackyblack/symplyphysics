from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output)

# Description
# Common formula for work done by gas is A = ∫p * dV, where
# A is work done,
# p is pressure,
# dV is change in volume.
# For process with constant pressure the equation can be transformed to
# A = p * (V2 - V1), where
# V1 is initial volume,
# V2 is resulting volume,
# p is pressure,
# A is work done.


pressure = Symbol("pressure", units.pressure)
init_volume = Symbol("init_volume", units.volume)
result_volume = Symbol("result_volume", units.volume)
work = Symbol("work", units.energy)

law = Eq(work, pressure * (result_volume - init_volume))


def print_law():
    print_expression(law)


@validate_input(
    pressure_=pressure,
    init_volume_=init_volume,
    result_volume_=result_volume
)
@validate_output(work)
def calculate_work(pressure_: Quantity, init_volume_: Quantity, result_volume_: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    result_work = result_work_expr.subs(
        {
            pressure: pressure_,
            init_volume: init_volume_,
            result_volume: result_volume_,
        }
    )
    return Quantity(result_work)
