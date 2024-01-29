from sympy import (Eq, solve)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
                           validate_output, dimensionless)

# Description
## A diffraction grating is an optical device whose operation is based on the use of the phenomenon of light diffraction. It is a collection of a large number of regularly spaced strokes (slits, protrusions) applied to a certain surface.
## The distance through which the strokes on the grating are repeated is called the period of the diffraction grating.


## Definition: d = 1 / N
## Where:
## d is period of the diffraction grating
## N is number of strokes


diffraction_grating_period = Symbol("diffraction_grating_period", units.length)
number_of_strokes = Symbol("number_of_strokes", dimensionless)

definition = Eq(diffraction_grating_period, 1 / number_of_strokes)

definition_units_SI = units.millimeter

def print_definition() -> str:
    return print_expression(definition)


@validate_input(number_of_strokes_=number_of_strokes)
@validate_output(diffraction_grating_period)
def calculate_diffraction_grating_period(number_of_strokes_: int) -> Quantity:
    result_expr = definition.rhs
    result_diffraction_grating_period = result_expr.subs(number_of_strokes, number_of_strokes_)
    return Quantity(result_diffraction_grating_period)
