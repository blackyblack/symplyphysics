from sympy import Eq
from symplyphysics import (
    units,
    Symbol,
    Quantity,
    validate_input,
    validate_output,
    print_expression,
)

# Description
## Also called Hooke's law, it is an empirical law which states that the force needed to extend
## or compress a spring by some distance is proportional with respect to the deformation of the
## spring.

# Law: F = -k * x
## F - restoring force exerted by the spring due to the deformation
## k - proportionality coefficient, characteristic of the spring, also known as stiffness
## x - deformation of the spring (positive or negative, see Note)

# Conditions
## - `x` is small compared to the total possible deformation of the spring
## - Only applies to elastic deformations of the body, i.e. the body reverts to its initial
##   state after the removal of force or load applied onto it.

# Note
## - The spring is aligned in the positive direction of the x-axis, thus the deformation can
##   be positive, in which case the spring is stretched, or negative, in which case the spring
##   is compressed. The sign of the force indicates its direction along the x-axis.

# Also see its [vector counterpart](../vector/spring_reaction_from_deformation.py)

restoring_force = Symbol("restoring_force", units.force)
stiffness = Symbol("stiffness", units.force / units.length)
deformation = Symbol("deformation", units.length)

law = Eq(restoring_force, -1 * stiffness * deformation)


def print_law() -> str:
    return print_expression(law)


@validate_input(stiffness_=stiffness, deformation_=deformation)
@validate_output(restoring_force)
def calculate_restoring_force(stiffness_: Quantity, deformation_: Quantity) -> Quantity:
    result = law.rhs.subs({stiffness: stiffness_, deformation: deformation_})
    return Quantity(result)
