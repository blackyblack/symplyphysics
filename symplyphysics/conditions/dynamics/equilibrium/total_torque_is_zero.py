"""
Total torque is zero
====================

For an object to be in rotational equilibrium, the total torque on the object must be zero.

**Links:**

#. `Physics LibreTexts <https://phys.libretexts.org/Bookshelves/University_Physics/Physics_(Boundless)/8%3A_Static_Equilibrium_Elasticity_and_Torque/8.2%3A_Conditions_for_Equilibrium>`__.
"""

from sympy import Eq
from symplyphysics import symbols

total_torque = symbols.torque
"""
Total :symbols:`torque` is the magnitue of the vector sum of all torques acting on the body.
"""

law = Eq(total_torque, 0)
"""
:laws:symbol::

:laws:latex::
"""
