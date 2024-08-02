from symplyphysics import (units, Quantity, validate_input, validate_output)
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.fields.operators import divergence_operator
from symplyphysics.core.fields.scalar_field import ScalarField
from symplyphysics.core.fields.vector_field import VectorField

## Description
## The divergence of electric induction (or electric displacement field) at some point
## is equal to the bulk charge density at that point. A simpler formulation of the law
## is that electric charges exist.

## Law is: div(D) = p, where
## D - electric induction (vector field),
## p - electric charge volumetric density (scalar field),
## div - divergence (the sum of partial derivatives in coordinates).


def charge_volumetric_density_law(electric_induction: VectorField) -> ScalarField:
    return ScalarField.from_expression(divergence_operator(electric_induction),
        electric_induction.coordinate_system)


@validate_input(cartesian_point_=units.length)
@validate_output(units.charge / units.volume)
def calculate_charge_volumetric_density_at_point(
        electric_induction_: VectorField, cartesian_point_: tuple[Quantity, Quantity,
    Quantity]) -> Quantity:
    electric_induction_vector = electric_induction_.apply(cartesian_point_)
    for i, c in enumerate(electric_induction_vector.components):
        assert_equivalent_dimension(c, f"electric_induction_vector[{i}]",
            "calculate_charge_volumetric_density_at_point", units.charge / units.area)
    charge_density_scalar_field_ = charge_volumetric_density_law(electric_induction_)
    result = charge_density_scalar_field_.apply(cartesian_point_)
    return Quantity(result)
