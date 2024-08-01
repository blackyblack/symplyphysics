from sympy import Eq
from symplyphysics import units, Symbol

# Description
## For an object to be in rotational equilibrium, the total torque on the object must be zero.

# Law: total_torque = 0
## total_torque - magnitude of the vector sum of all torques acting on a body

total_torque = Symbol("total_torque", units.length * units.force)

law = Eq(total_torque, 0)
