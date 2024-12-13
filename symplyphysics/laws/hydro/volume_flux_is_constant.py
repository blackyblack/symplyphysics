r"""
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
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)

tube_area = Function("tube_area", units.area)
"""
Cross-sectional area of the tube of flow as a function of time.

Symbol:
    :code:`A(t)`
"""

fluid_speed = Function("fluid_speed", units.velocity)
"""
Fluid speed as a function of time.

Symbol:
    :code:`u(t)`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

law = Eq(Derivative(tube_area(time) * fluid_speed(time), time), 0)
r"""
:code:`Derivative(A(t) * u(t), t) = 0`

Latex:
    .. math::
        \frac{d (A u)}{d t} = 0
"""


@validate_input(tube_area_before_=tube_area,
    fluid_speed_before_=fluid_speed,
    tube_area_after_=tube_area)
@validate_output(fluid_speed)
def calculate_fluid_speed(tube_area_before_: Quantity, fluid_speed_before_: Quantity,
    tube_area_after_: Quantity) -> Quantity:
    dsolved = dsolve(law, fluid_speed(time))
    c1_value = solve(dsolved, "C1")[0].subs({
        tube_area(time): tube_area_before_,
        fluid_speed(time): fluid_speed_before_,
    })
    result_expr = dsolved.subs({
        "C1": c1_value,
        tube_area(time): tube_area_after_,
    }).rhs
    return Quantity(result_expr)
