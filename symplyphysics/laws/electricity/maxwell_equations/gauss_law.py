from sympy import Eq, solve, sqrt, log, pi, sin, Expr, Derivative
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, Function, Vector, SI)
from symplyphysics.core.fields.operators import divergence_operator as div
from symplyphysics.core.coordinate_systems.coordinate_systems import CoordinateSystem
from symplyphysics.core.fields.vector_field import VectorField
from sympy.physics.units import electric_constant

## Description
## The divergence of electrical induction at some point is equal to the bulk charge density at that
## point. A simpler formulation of the law is that electric charges exist.

## Law is: div(D(x, y, z)) = p, where
## D - the vector of electrical induction,
## x, y, z - coordinates in a rectangular coordinate system,
## p - the bulk density of the electric charge,
## div - divergence (the sum of partial derivatives in coordinates).

electrical_induction_x = Function("D_x", units.charge / units.area)
electrical_induction_y = Function("D_y", units.charge / units.area)
electrical_induction_z = Function("D_z", units.charge / units.area)
x = Symbol('x', units.length)
y = Symbol('y', units.length)
z = Symbol('z', units.length)
bulk_density_of_electric_charge = Symbol('bulk_density_of_electric_charge', units.charge / units.volume)

law = Eq(Derivative(electrical_induction_x(x, y, z), (x, 1)) + Derivative(electrical_induction_y(x, y, z), (y, 1)) + Derivative(electrical_induction_z(x, y, z), (z, 1)), bulk_density_of_electric_charge)


def apply_electrical_induction_function(electrical_induction_x_: Expr, electrical_induction_y_: Expr, electrical_induction_z_: Expr) -> Expr:
    applied_law = law.subs(electrical_induction_x(x, y, z), electrical_induction_x_)
    applied_law = applied_law.subs(electrical_induction_y(x, y, z), electrical_induction_y_)
    applied_law = applied_law.subs(electrical_induction_z(x, y, z), electrical_induction_z_)
    return applied_law


@validate_input(x_=x, y_=y, z_=z)
@validate_output(bulk_density_of_electric_charge)
def calculate_bulk_density_of_electric_charge(electrical_induction_: VectorField, x_: Quantity, y_: Quantity, z_: Quantity) -> Quantity:
    electrical_induction_x_ = electrical_induction_.field_function[0]
    electrical_induction_y_ = electrical_induction_.field_function[1]
    electrical_induction_z_ = electrical_induction_.field_function[2]
    assert SI.get_dimension_system().equivalent_dims(Quantity(electrical_induction_x_).dimension,
        electrical_induction_x.dimension)
    assert SI.get_dimension_system().equivalent_dims(Quantity(electrical_induction_y_).dimension,
        electrical_induction_y.dimension)
    assert SI.get_dimension_system().equivalent_dims(Quantity(electrical_induction_z_).dimension,
        electrical_induction_z.dimension)
    applied_law = apply_electrical_induction_function(electrical_induction_x_, electrical_induction_y_, electrical_induction_z_)
    result_expr = applied_law.subs({
        x: x_,
        y: y_,
        z: z_,
    })
    result = solve(result_expr, bulk_density_of_electric_charge, dict=True)[0][bulk_density_of_electric_charge]
    return Quantity(result)
