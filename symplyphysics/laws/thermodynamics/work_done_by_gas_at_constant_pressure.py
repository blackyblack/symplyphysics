from sympy import (Eq, solve, integrate)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output)

# Description
# Common formula for work done by gas is A = ∫p * dV
# A is work done
# p is pressure
# dV is change in volume


pressure = Symbol("pressure", units.pressure)
volume = Symbol("volume", units.volume)
init_volume = Symbol("init_volume", units.volume)
result_volume = Symbol("result_volume", units.volume)
work = Symbol("work", units.energy)

common_law = integrate(pressure, (volume, init_volume, result_volume))
law = Eq(work, common_law)


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
