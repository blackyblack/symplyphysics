from sympy import Eq
from sympy.core import Equality
from symplyphysics import (units, Quantity, Symbol, validate_input)
from symplyphysics.core.dimensions import assert_equivalent_dimension
from symplyphysics.core.fields.operators import divergence_operator
from symplyphysics.core.fields.vector_field import VectorField

## Description
## The divergence of magnetic induction at some point is equal to zero. A simpler formulation of the law
## is that magnetic charges not exist.

## Law is: div(B) = 0, where
## B - magnetic induction (vector field),
## div - divergence (the sum of partial derivatives in coordinates).

# x_coor = Symbol('x_coor', units.length)
# y_coor = Symbol('y_coor', units.length)
# z_coor = Symbol('z_coor', units.length)
time = Symbol('time', units.time)


# def gauss_law_for_magnetic_field(magnetic_induction_: VectorField, time_: Quantity) -> Equality:
#     divergence_magnetic_induction = divergence_operator(magnetic_induction_)
#     x = magnetic_induction_.coordinate_system.coord_system.base_scalars()[0]
#     y = magnetic_induction_.coordinate_system.coord_system.base_scalars()[1]
#     z = magnetic_induction_.coordinate_system.coord_system.base_scalars()[2]
#     dict_xyz = {
#         x: x_coor,
#         y: y_coor,
#         z: z_coor,
#         time_: time,
#     }
#     law = Eq(divergence_magnetic_induction.subs(dict_xyz), 0)
#     return law
def gauss_law_for_magnetic_field(magnetic_induction_: VectorField) -> Equality:
    divergence_magnetic_induction = divergence_operator(magnetic_induction_)
    law = Eq(divergence_magnetic_induction, 0)
    return law


@validate_input(cartesian_point_=units.length,
                time_=units.time)
def get_gauss_law_for_magnetic_field(magnetic_induction_: VectorField, time_: Quantity, 
                                     cartesian_point_: tuple[Quantity, Quantity, Quantity]) -> Equality:
    magnetic_induction_vector = magnetic_induction_.apply(cartesian_point_)
    for i, c in enumerate(magnetic_induction_vector.components):
        assert_equivalent_dimension(c, f"magnetic_induction_vector[{i}]",
            "calculate_charge_volumetric_density_at_point", units.magnetic_density)
    result = gauss_law_for_magnetic_field(magnetic_induction_, time_)
    
    return result
