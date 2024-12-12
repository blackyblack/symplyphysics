"""
Work of ideal gas in isobaric process
=====================================

An isobaric process is a process that occurs at constant pressure. The work done
by the gas is proportional to the pressure of the gas and its volume change.

**Conditions:**

#. The gas is ideal.
#. The process is isobaric.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Isobaric_process>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import (
    work_is_integral_of_pressure_over_volume as work_law,
)

work = symbols.work
"""
:symbols:`work` done by the gas during the isobaric process.
"""

pressure = symbols.pressure
"""
:symbols:`pressure` of the gas, constant during the process.
"""

initial_volume = clone_as_symbol(symbols.volume, subscript="0")
"""
Initial :symbols:`volume` of the gas.
"""

final_volume = clone_as_symbol(symbols.volume, subscript="1")
"""
Final :symbols:`volume` of the gas
"""

law = Eq(work, pressure * (final_volume - initial_volume))
"""
:laws:symbol::

:laws:latex::
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
