r"""
Speed of light via vacuum permittivity and permeability
=======================================================

Speed of light can be expressed using the permittivity and permeability
of the vacuum.

**Notation:**

#. :quantity_notation:`speed_of_light`.
#. :quantity_notation:`vacuum_permittivity`.
#. :quantity_notation:`vacuum_permeability`.
"""

from sympy import (Eq, sqrt)
from symplyphysics import units, convert_to, quantities

law = Eq(quantities.speed_of_light, 1 / sqrt(quantities.vacuum_permittivity * quantities.vacuum_permeability))
r"""
:code:`c = 1 / sqrt(epsilon_0 * mu_0)`

Latex:
    .. math::
        c = \frac{1}{\sqrt{\varepsilon_0 \mu_0}}
"""

assert convert_to(law.lhs, units.meter / units.second) == convert_to(law.rhs,
    units.meter / units.second)
