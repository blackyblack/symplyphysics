from sympy import Eq, Derivative, Matrix
from sympy.core import Equality
from symplyphysics import (units, Symbol)
from symplyphysics.core.fields.operators import curl_operator
from symplyphysics.core.fields.vector_field import VectorField


## Description
## Faraday's law of induction states that a change in magnetic induction generates a vortex electric field.
## This law is valid for any magnetic field that changes over time.

## Law is: curl(E(r, t)) = -d(B(r, t)) / dt, where
## E - the field of electric intensity,
## B - the field of magnetic induction,
## r - the vector of a point in space,
## t - time,
## curl - curl of the vector of magnetic induction,
## d / dt - time partial derivative.

time = Symbol('time', units.time)


def faradays_law_of_induction(electric_intensity_: VectorField, 
                              magnetic_induction_: VectorField) -> tuple[Equality, Equality, Equality]:
    curl_electric_intensity = curl_operator(electric_intensity_)

    magnetic_induction_space = magnetic_induction_.apply_to_basis()
    magnetic_induction_time_derivative = [-Derivative(c, time) for c in magnetic_induction_space.components]

    return Eq(Matrix(curl_electric_intensity.apply_to_basis().components), Matrix(magnetic_induction_time_derivative))
