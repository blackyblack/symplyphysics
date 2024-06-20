from sympy import Eq
from sympy.core import Equality
from symplyphysics import (units, Symbol)
from symplyphysics.core.fields.operators import divergence_operator
from symplyphysics.core.fields.vector_field import VectorField

## Description
## The divergence of magnetic induction at some point is equal to zero. A simpler formulation of the law
## is that magnetic charges not exist.

## Law is: div(B) = 0, where
## B - magnetic induction (vector field),
## div - divergence (the sum of partial derivatives in coordinates).

time = Symbol('time', units.time)

def gauss_law_for_magnetic_field(magnetic_induction_: VectorField) -> Equality:
    divergence_magnetic_induction = divergence_operator(magnetic_induction_)
    law = Eq(divergence_magnetic_induction, 0)
    return law
