r"""
Equation in homogeneous medium in one dimension
===============================================

Heat equation governs heat diffusion, as well as other diffusive processes. It describes
the evolution of heat transferred from hotter to colder environments in time and space.

**Notes:**

#. There is no straghtforward solution to this equation, and it depends on initial
   conditions as well.

**Conditions:**

#. There are no heat sources in the system, i.e. the heat distribution only depends on
   the initial conditions.
#. Thermal diffusivity :math:`\chi` does not depend on position.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Heat_equation#Heat_flow_in_a_uniform_rod>`__.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    units,
    Symbol,
    Function,
)

temperature = Function("temperature", units.temperature)
"""
Temperature as a function of position and time.

Symbol:
    :code:`T(x, t)`
"""

position = Symbol("position", units.length)
"""
Position, or spatial variable.

Symbol:
    :code:`x`
"""

time = Symbol("time", units.time)
"""
Time.

Symbol:
    :code:`t`
"""

# TODO: Create law of thermal diffusivity via thermal conductivity

thermal_diffusivity = Symbol("thermal_diffusivity", units.area / units.time)
r"""
`Thermal diffusivity <https://en.wikipedia.org/wiki/Thermal_diffusivity>`_.

Symbol:
    :code:`chi`

Latex:
    :math:`\chi`
"""

law = Eq(Derivative(temperature(position, time), time),
    thermal_diffusivity * Derivative(temperature(position, time), position, 2))
r"""
:code:`Derivative(T(x, t), t) = chi * Derivative(T(x, t), (x, 2))`

Latex:
    .. math::
        \frac{\partial T}{\partial t} = \chi \frac{\partial^2 T}{\partial x^2}
"""
