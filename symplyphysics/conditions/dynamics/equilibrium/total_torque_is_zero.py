from sympy import Eq
from symplyphysics import units, Symbol

# Description
## For an object to be in rotational equilibrium, the total torque on the object must be zero.

# Law: total_torque = 0
## total_torque - magnitude of the vector sum of all torques acting on a body

# Links:
## Physics LibreTexts <https://phys.libretexts.org/Bookshelves/University_Physics/Physics_(Boundless)/8%3A_Static_Equilibrium_Elasticity_and_Torque/8.2%3A_Conditions_for_Equilibrium>

total_torque = Symbol("total_torque", units.length * units.force)

law = Eq(total_torque, 0)
