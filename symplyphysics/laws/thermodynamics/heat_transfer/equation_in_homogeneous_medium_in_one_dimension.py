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
from symplyphysics import symbols, clone_as_function

position = symbols.position
"""
:symbols:`position`, or spatial variable.
"""

time = symbols.time
"""
:symbols:`time`.
"""

temperature = clone_as_function(symbols.temperature, [position, time])
"""
Temperature as a function of :attr:`~position` and :attr:`~time`.
"""

# TODO: Create law of thermal diffusivity via thermal conductivity

thermal_diffusivity = symbols.thermal_diffusivity
"""
:symbols:`thermal_diffusivity`.
"""

law = Eq(Derivative(temperature(position, time), time),
    thermal_diffusivity * Derivative(temperature(position, time), position, 2))
"""
:laws:symbol::

:laws:latex::
"""
