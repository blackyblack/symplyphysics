"""
Speed of light via vacuum permittivity and permeability
=======================================================

Speed of light can be expressed using the permittivity and permeability
of the vacuum.

**Notation:**

#. :quantity_notation:`speed_of_light`.
#. :quantity_notation:`vacuum_permittivity`.
#. :quantity_notation:`vacuum_permeability`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Speed_of_light#Propagation_of_light>`__.
"""

from sympy import (Eq, sqrt)
from symplyphysics import units, convert_to, quantities

law = Eq(quantities.speed_of_light,
    1 / sqrt(quantities.vacuum_permittivity * quantities.vacuum_permeability))
"""
:laws:symbol::

:laws:latex::
"""

assert convert_to(law.lhs, units.meter / units.second) == convert_to(law.rhs,
    units.meter / units.second)
