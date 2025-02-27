from symplyphysics import Quantity
from symplyphysics.core.fields.operators import divergence_operator
from symplyphysics.core.fields.vector_field import VectorField

## Description
## The divergence of magnetic induction is equal to zero. A simpler formulation of the law
## is that magnetic charges do not exist.

## Law is: div(B) = 0, where
## B - magnetic induction (vector field),
## div - divergence (the sum of partial derivatives in coordinates).

# Links:
## Wikipedia, second equation <https://en.wikipedia.org/wiki/Maxwell%27s_equations#>


def magnetic_field_divergence_condition(magnetic_induction_: VectorField) -> bool:
    divergence_magnetic_induction = divergence_operator(magnetic_induction_)
    divergence_magnetic_induction_quantity = Quantity(divergence_magnetic_induction)
    return divergence_magnetic_induction_quantity.scale_factor == 0  # type: ignore[no-any-return]
