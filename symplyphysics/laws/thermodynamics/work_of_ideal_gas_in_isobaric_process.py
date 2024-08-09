"""
Work of ideal gas in isobaric process
=====================================

An isobaric process is a process that occurs at constant pressure. The work done
by the gas is proportional to the pressure of the gas and its volume change.

**Conditions:**

#. The gas is ideal.
"""

from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import work_is_integral_of_pressure_over_volume as work_law

work = Symbol("work", units.energy)
"""
Work done by the gas during the isobaric process.

Symbol:
    :code:`W`
"""

pressure = Symbol("pressure", units.pressure)
"""
Pressure of the gas, constant during the process.

Symbol:
    :code:`p`
"""

initial_volume = Symbol("initial_volume", units.volume)
r"""
Initial volume of the gas.

Symbol:
    :code:`V_0`

Latex:
    :math:`V_0`
"""

final_volume = Symbol("final_volume", units.volume)
r"""
Final volume of the gas

Symbol:
    :code:`V_1`

Latex:
    :math:`V_1`
"""

law = Eq(work, pressure * (final_volume - initial_volume))
r"""
:code:`W = p * (V_1 - V_0)`

Latex:
    .. math::
        W = p (V_1 - V_0)
"""

# Derive law from work integral

_work_expr = work_law.law.rhs.subs({
    work_law.pressure(work_law.volume): pressure,
    work_law.initial_volume: initial_volume,
    work_law.final_volume: final_volume,
}).doit()

assert expr_equals(_work_expr, law.rhs)


@validate_input(pressure_=pressure, init_volume_=initial_volume, result_volume_=final_volume)
@validate_output(work)
def calculate_work(pressure_: Quantity, init_volume_: Quantity,
    result_volume_: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    result_work = result_work_expr.subs({
        pressure: pressure_,
        initial_volume: init_volume_,
        final_volume: result_volume_,
    })
    return Quantity(result_work)
