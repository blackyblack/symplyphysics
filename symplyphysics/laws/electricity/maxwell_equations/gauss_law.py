from sympy import Eq, solve
from symplyphysics import (units, Quantity, Symbol, validate_input, validate_output, SI)
from symplyphysics.core.fields.operators import divergence_operator as div
from symplyphysics.core.fields.vector_field import VectorField

## Description
## The divergence of electrical induction at some point is equal to the bulk charge density at that
## point. A simpler formulation of the law is that electric charges exist.

## Law is: div(D(x, y, z)) = p, where
## D - the vector of electrical induction,
## x, y, z - coordinates in a rectangular coordinate system,
## p - the bulk density of the electric charge,
## div - divergence (the sum of partial derivatives in coordinates).

x_coor = Symbol('x_coor', units.length)
y_coor = Symbol('y_coor', units.length)
z_coor = Symbol('z_coor', units.length)
divergence_electrical_induction_vector = Symbol('divergence_electrical_induction_vector', units.charge / units.volume)
bulk_density_of_electric_charge = Symbol('bulk_density_of_electric_charge', units.charge / units.volume)

law = Eq(divergence_electrical_induction_vector, bulk_density_of_electric_charge)


@validate_input(x_=x_coor, y_=y_coor, z_=z_coor)
@validate_output(bulk_density_of_electric_charge)
def calculate_bulk_density_of_electric_charge(electrical_induction_: VectorField, x_: Quantity, y_: Quantity, z_: Quantity) -> Quantity:
    x = electrical_induction_.coordinate_system.coord_system.base_scalars()[0]
    y = electrical_induction_.coordinate_system.coord_system.base_scalars()[1]
    z = electrical_induction_.coordinate_system.coord_system.base_scalars()[2]
    dict_xyz = {
        x: x_,
        y: y_,
        z: z_,
    }

    field_components = list(electrical_induction_.apply_to_basis().components) + [0] * (3 - len(electrical_induction_.apply_to_basis().components))
    assert SI.get_dimension_system().equivalent_dims(Quantity(field_components[0].subs(dict_xyz)).dimension, units.charge / units.area)
    assert SI.get_dimension_system().equivalent_dims(Quantity(field_components[1].subs(dict_xyz)).dimension, units.charge / units.area)
    assert SI.get_dimension_system().equivalent_dims(Quantity(field_components[2].subs(dict_xyz)).dimension, units.charge / units.area)

    divergence_electrical_induction_vector_ = div(electrical_induction_)

    applied_law = law.subs(divergence_electrical_induction_vector, divergence_electrical_induction_vector_)
    result = applied_law.subs(dict_xyz)
    result = solve(result, bulk_density_of_electric_charge, dict=True)[0][bulk_density_of_electric_charge]
    return Quantity(result)
