from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
# Work done by gas: A = p * (V_2 - V_1), where
# A is work done
# p is pressure
# V_2 is resulting volume
# V_1 is initial volume

pressure = Symbol("pressure", units.pressure)
volume_1 = Symbol("volume_1", units.volume)
volume_2 = Symbol("volume_2", units.volume)
work = Symbol("work", units.energy)

law = Eq(work, pressure * (volume_2 - volume_1))

def print_law():
    print_expression(law)


@validate_input(pressure=Quantity, volume_change=Quantity)
@validate_output(work)
def calculate_work(pressure_: Quantity, volume_1_: Quantity, volume_2_: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    result_work = result_work_expr.subs(
        {
            pressure: pressure_,
            volume_1: volume_1_,
            volume_2: volume_2_,
        }
    )
    return Quantity(result_work)
