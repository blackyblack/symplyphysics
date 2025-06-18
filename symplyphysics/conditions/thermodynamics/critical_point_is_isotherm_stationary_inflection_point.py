"""
Critical point is stationary inflection point of isotherm
=========================================================

Critical point (in the thermodynamic sense) is such values of volume, pressure, and temperature at
which only one phase exists and at the vicinity of which the physical properties of the phases of
the substance change dramatically. Algebraically, the critical point is the stationary inflection
point of the isothermal pressure-volume dependency line.

**Note:**

#. These equations need to be solved together with the :ref:`equation of state <**Equations of state**>`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Critical_point_(thermodynamics)#Overview>`__.
"""

from sympy import Eq, Derivative
from symplyphysics import symbols, clone_as_function

volume = symbols.volume
"""
:symbols:`volume`
"""

pressure = clone_as_function(symbols.pressure, [volume])
"""
:symbols:`pressure` as a function of :attr:`~volume`.
"""

inflection_point_condition = Eq(Derivative(pressure(volume), volume, 1), 0)
"""
:laws:symbol::

:laws:latex::
"""

flat_tangent_condition = Eq(Derivative(pressure(volume), volume, 2), 0)
"""
:laws:symbol::

:laws:latex::
"""

conditions = inflection_point_condition, flat_tangent_condition
