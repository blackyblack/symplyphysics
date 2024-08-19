r"""
Speed of light via vacuum permittivity and permeability
=======================================================

Speed of light can be expressed using the permittivity and permeability
of the vacuum.

**Notation:**

#. :math:`c` is the speed of light.
#. :math:`\varepsilon_0` (:code:`epsilon_0`) is vacuum permittivity, or electric constant.
#. :math:`\mu_0` (:code:`mu_0`) is vacuum permeability, or magnetic constant.
"""

from sympy import (Eq, sqrt)
from sympy.physics.units import speed_of_light, magnetic_constant, electric_constant
from symplyphysics import units, convert_to

law = Eq(speed_of_light, 1 / sqrt(electric_constant * magnetic_constant))
r"""
:code:`c = 1 / sqrt(epsilon_0 * mu_0)`

Latex:
    .. math::
        c = \frac{1}{\sqrt{\varepsilon_0 \mu_0}}
"""

assert convert_to(law.lhs, units.meter / units.second) == convert_to(law.rhs,
    units.meter / units.second)
