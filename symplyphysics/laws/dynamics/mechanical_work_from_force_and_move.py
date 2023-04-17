import numbers
from sympy import Expr
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from sympy.vector import Dot
from symplyphysics import (
    Eq, pretty, solve, units, validate_input, expr_to_quantity,
    CoordinateSystem, coordinates_transform, Vector, vector_rebase, sympy_vector_from_vector
)
from symplyphysics.core.quantity_decorator import validate_input_symbols, validate_output_symbol
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.core.symbols.symbols import Symbol, to_printable

# Description
## Work is measured result of force applied. Mechanical work is the only reason for the object energy to be changed.
## Work is scalar value equal to force multiplied by movement and by cosine of angle between force and movement vectors. 
## Law: A = F * S, where
## A is mechanical work
## F is force vector applied to object
## S is movement vector caused by this force
## * is a scalar multiplication of vectors (dot product)

work = Symbol("work", units.energy)
force = Symbol("force", units.force)
distance = Symbol("distance", units.length)

law = Eq(work, Dot(force, distance))

def print(expr: Expr) -> str:
    symbols = [work, force, distance]
    return pretty(to_printable(expr, symbols), use_unicode=False)

#TODO: consider removing 'validate_input' and pass angles as float numbers
@validate_input_symbols(force_=force, distance_=distance)
@validate_input(force_angle=angle_type, distance_angle=angle_type)
@validate_output_symbol(work)
def calculate_work(force_: Quantity, distance_: Quantity, force_angle: Quantity | float, distance_angle: Quantity | float) -> Quantity:
    result_work_expr = solve(law, work, dict=True)[0][work]
    coordinates_polar = CoordinateSystem(CoordinateSystem.System.CYLINDRICAL)
    coordinates_cartesian = coordinates_transform(coordinates_polar, CoordinateSystem.System.CARTESIAN)
    #HACK: sympy angles are always in radians
    force_angle_radians = force_angle if isinstance(force_angle, numbers.Number) else force_angle.scale_factor
    distance_angle_radians = distance_angle if isinstance(distance_angle, numbers.Number) else distance_angle.scale_factor

    # Convert to Cartesian coordinates as Dot product does not work properly for Cylindrical coordinates
    force_vector = Vector([force_, force_angle_radians], coordinates_polar)
    distance_vector = Vector([distance_, distance_angle_radians], coordinates_polar)
    force_vector_cartesian = vector_rebase(force_vector, coordinates_cartesian)
    distance_vector_cartesian = vector_rebase(distance_vector, coordinates_cartesian)
    force_vector_cartesian_sympy = sympy_vector_from_vector(force_vector_cartesian)
    distance_vector_cartesian_sympy = sympy_vector_from_vector(distance_vector_cartesian)
    result_expr = result_work_expr.subs({force: force_vector_cartesian_sympy, distance: distance_vector_cartesian_sympy}).doit()
    return expr_to_quantity(result_expr)
