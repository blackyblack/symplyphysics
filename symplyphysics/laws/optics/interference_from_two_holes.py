from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression,
                           validate_input, validate_output)

# Description
## The coherence of vibrations that occur at the same
## time at different points in the plane perpendicular to the direction
## of wave propagation is called spatial coherence.
## Spatial coherence depends on the conditions of radiation and
## the formation of light waves in extended sources.

# Law: delta = x * d / l
# Where:
## delta - optical travel difference for two coherent waves
## x - coordinate in interference picture
## d - distance between holes (or source parameter)
## l - distance from holes to picture


travel_difference = Symbol("travel_difference", units.length)
coordinate = Symbol("coordinate", units.length)
distance_holes = Symbol("distance_holes", units.length)
distance_to_picture = Symbol("distance_to_picture", units.length)

law = Eq(travel_difference, coordinate * distance_holes / distance_to_picture)


def print_law():
    return print_expression(law)


@validate_input(coordinate_=coordinate, distance_holes_=distance_holes, distance_to_picture_=distance_to_picture)
@validate_output(travel_difference)
def calculate_travel_difference(coordinate_: Quantity, distance_holes_: Quantity, distance_to_picture_: Quantity) -> Quantity:
    solved = solve(law, travel_difference, dict=True)[0][travel_difference]
    result_expr = solved.subs({
        coordinate: coordinate_,
        distance_holes: distance_holes_,
        distance_to_picture: distance_to_picture_,
    })
    result = Quantity(result_expr)
    return result
