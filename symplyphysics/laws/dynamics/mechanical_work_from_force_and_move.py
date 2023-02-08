from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)
from symplyphysics.core.vectors.vectors import extended_express
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from sympy.vector import CoordSys3D, Dot

# Description
## Work is measured result of force applied. Mechanical work is the only reason for the object energy to be changed.
## Work is scalar value equal to force multiplied by movement and by cosine of angle between force and movement vectors. 
## Law: A = F * S, where
## A is mechanical work
## F is force vector applied to object
## S is movement vector caused by this force
## * is a scalar multiplication of vectors (dot product)


work, force, distance = symbols('work force distance')
law = Eq(work, Dot(force, distance))

def print():
    return pretty(law, use_unicode=False)

@validate_input(force_=units.force, distance_=units.length, force_angle=angle_type, distance_angle=angle_type)
@validate_output(units.energy)
def calculate_work(force_: Quantity, distance_: Quantity, force_angle: Quantity, distance_angle: Quantity) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    Ca = CoordSys3D("Ca")
    Cy = Ca.create_new("Cy", transformation="cylindrical")
    #HACK: sympy angles are always in radians
    force_angle_radians = force_angle.scale_factor
    distance_angle_radians = distance_angle.scale_factor

    # Convert to Cartesian coordinates as Dot product does not work properly for Cylindrical coordinates
    force_vector = force_ * Cy.i + force_angle_radians * Cy.j
    distance_vector = distance_ * Cy.i + distance_angle_radians * Cy.j
    force_vector_cartesian = extended_express(force_vector, Ca)
    distance_vector_cartesian = extended_express(distance_vector, Ca)
    result_expr = result_work_expr.subs({force: force_vector_cartesian, distance: distance_vector_cartesian}).doit()
    return expr_to_quantity(result_expr, "vector_work")