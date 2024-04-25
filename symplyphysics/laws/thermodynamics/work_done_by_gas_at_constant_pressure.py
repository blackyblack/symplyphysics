from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import work_is_volume_integral_of_pressure as work_law

# Description
# When energy is added to gas molecules and increases their kinetic energy,
# the gas expands and does work on its surroundings. The work done by the gas
# with constant pressure can be found by: A = p * ΔV, where
# ΔV is the change in the volume of the gas,
# p is pressure,
# A is work done.
# When the volume of a gas changes from to V1 to V2, the change
# in the volume of the gas is ΔV = V2 - V1.
# Law: A = p * (V2 - V1), where
# V1 is initial volume,
# V2 is resulting volume,
# p is pressure,
# A is work done.

pressure = Symbol("pressure", units.pressure)
init_volume = Symbol("init_volume", units.volume)
result_volume = Symbol("result_volume", units.volume)
work = Symbol("work", units.energy)

law = Eq(work, pressure * (result_volume - init_volume))

# Derive law from work integral

_work_expr = work_law.law.rhs.xreplace({
    work_law.pressure(work_law.volume): pressure,
}).subs({
    work_law.volume_before: init_volume,
    work_law.volume_after: result_volume,
}).doit()

assert expr_equals(_work_expr, law.rhs)


@validate_input(pressure_=pressure, init_volume_=init_volume, result_volume_=result_volume)
@validate_output(work)
def calculate_work(pressure_: Quantity, init_volume_: Quantity,
    result_volume_: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    result_work = result_work_expr.subs({
        pressure: pressure_,
        init_volume: init_volume_,
        result_volume: result_volume_,
    })
    return Quantity(result_work)
