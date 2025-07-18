"""
Volume flux is constant
=======================

The product of the area and the fluid speed, which is called the *volume flux*, is constant
at all points along the tube of flow of an incompressible liquid. This equation is also
known as the *equation of continuity*.

**Conditions:**

#. The fluid is :ref:`ideal <ideal_fluid_def>`.

**Links:**

#. `Engineering LibreTexts, derivable from here <https://eng.libretexts.org/Bookshelves/Aerospace_Engineering/Fundamentals_of_Aerospace_Engineering_(Arnedo)/03%3A_Aerodynamics/3.01%3A_Fundamentals_of_fluid_mechanics/3.1.02%3A_Continuity_equation>`__.
"""

from sympy import Eq, solve, dsolve, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

time = symbols.time
"""
:symbols:`time`.
"""

tube_area = clone_as_function(symbols.area, [time])
"""
Cross-sectional :symbols:`area` of the tube of flow as a function of :attr:`~time`.
"""

flow_speed = clone_as_function(symbols.flow_speed, [time])
"""
:symbols:`flow_speed` of the fluid as a function of :attr:`~time`.
"""

law = Eq(Derivative(tube_area(time) * flow_speed(time), time), 0)
"""
:laws:symbol::

:laws:latex::
"""

# Derivable from a [more general continuity equation](https://en.wikipedia.org/wiki/Continuity_equation#Fluid_dynamics)
# for ideal fluids.


@validate_input(tube_area_before_=tube_area,
    fluid_speed_before_=flow_speed,
    tube_area_after_=tube_area)
@validate_output(flow_speed)
def calculate_fluid_speed(tube_area_before_: Quantity, fluid_speed_before_: Quantity,
    tube_area_after_: Quantity) -> Quantity:
    dsolved = dsolve(law, flow_speed(time))
    c1_value = solve(dsolved, "C1")[0].subs({
        tube_area(time): tube_area_before_,
        flow_speed(time): fluid_speed_before_,
    })
    result_expr = dsolved.subs({
        "C1": c1_value,
        tube_area(time): tube_area_after_,
    }).rhs
    return Quantity(result_expr)
