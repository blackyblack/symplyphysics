from sympy import Eq
from sympy.core import Equality
from symplyphysics.core.fields.operators import divergence_operator
from symplyphysics.core.fields.vector_field import VectorField

## Description
## The divergence of magnetic induction is equal to zero. A simpler formulation of the law
## is that magnetic charges not exist.

## Law is: div(B) = 0, where
## B - magnetic induction (vector field),
## div - divergence (the sum of partial derivatives in coordinates).


def magnetic_field_divergence_condition(magnetic_induction_: VectorField) -> Equality:
    divergence_magnetic_induction = divergence_operator(magnetic_induction_)
    return Eq(divergence_magnetic_induction, 0)
