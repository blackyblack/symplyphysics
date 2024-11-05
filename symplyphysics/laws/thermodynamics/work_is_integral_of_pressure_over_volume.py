"""
Work is integral of pressure over volume
========================================

The pressure-volume work is the work done by a thermodynamic system when the volume of the
system changes.

**Conditions:**

#. The process is reversible or quasi-static.
#. The system is closed.
"""

from sympy import Eq, Integral, Point2D, solve
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
    clone_as_symbol,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.core.geometry.line import two_point_function
from symplyphysics.laws.thermodynamics import infinitesimal_work_in_quasistatic_process as work_law

work = symbols.work
"""
Pressure-volume :symbols:`work` done by the system.
"""

volume = symbols.volume
"""
:symbols:`volume` of the system.
"""

pressure = clone_as_function(symbols.pressure, [volume])
"""
:symbols:`pressure` inside the system as a function of :attr:`~volume`.
"""

initial_volume = clone_as_symbol(symbols.volume, subscript="0")
"""
Initial :symbols:`volume` of the system.
"""

final_volume = clone_as_symbol(symbols.volume, subscript="1")
"""
Final :symbols:`volume` of the system.
"""

law = Eq(work, Integral(pressure(volume), (volume, initial_volume, final_volume)))
"""
:laws:symbol::

:laws:latex::
"""

# Derive law from infinitesimal counterpart

_work_function = clone_as_function(symbols.work)

# Prepare law for integration
# Note that `dW = (dW/dV) * dV` and we set `dV` to 1
_infinitesimal_law = work_law.law.subs({
    work_law.infinitesimal_work_done: _work_function(volume).diff(volume),
    work_law.infinitesimal_volume_change: 1,
    work_law.pressure: pressure(volume),
})

_integrated_law = Eq(_infinitesimal_law.lhs.integrate((volume, initial_volume, final_volume)),
    _infinitesimal_law.rhs.integrate((volume, initial_volume, final_volume))).subs({
    _work_function(initial_volume): 0,  # there was no work done at the start of the process
    _work_function(final_volume):
            work,  # the work at the end of the process is the work done during the process
    })

_work_derived = solve(_integrated_law, work)[0]

assert expr_equals(_work_derived, law.rhs)


@validate_input(
    volume_before_=initial_volume,
    volume_after_=final_volume,
    pressure_before_=pressure,
    pressure_after_=pressure,
)
@validate_output(work)
def calculate_work(
    volume_before_: Quantity,
    volume_after_: Quantity,
    pressure_before_: Quantity,
    pressure_after_: Quantity,
) -> Quantity:
    pressure_ = two_point_function(
        Point2D(volume_before_, pressure_before_),
        Point2D(volume_after_, pressure_after_),
        volume,
    )
    result = law.rhs.subs({
        pressure(volume): pressure_,
        initial_volume: volume_before_,
        final_volume: volume_after_,
    }).doit()
    return Quantity(result)
