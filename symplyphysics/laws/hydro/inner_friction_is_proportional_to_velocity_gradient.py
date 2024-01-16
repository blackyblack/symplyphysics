from sympy import Eq, Derivative
from symplyphysics import units, Quantity, Symbol, Function, print_expression, validate_input, validate_output

# Description
## An element of flowing fluid will endure forces from the surrounding fluid (stress forces), that
## will gradually deform the element over time. In a steady (laminar) flow, these stress forces
## occur between fluid layers. For an incompressible and isotropic fluid, the shear stress exerted
## to an element of fluid is proportional to the spatial derivative of the velocity component
## that is perpendicular to the velocity vector (and as such, parallel to the direction of shear).
## This law is also called the Newton's law of viscosity, and the fluids that follow it are said to
## be Newtonian.

## As an example of this law, consider two solid flat plates that contain water in between.
## The bottom plate is fixed in place, while the top plate moves parallel to the bottom one with
## a small speed v (such that the flow of the water is steady). If we were to measure the force that
## needed to make the top plate continue to move, we would find that it is proportional to the area
## of the plate and the ratio v/d where d is the distance between the plates.

# Law: tau = eta * dv/dy
## tau - shear stress, i.e. the component of stress coplanar with a material cross section
## eta - dynamic viscosity
## v - speed of fluid flow
## d/dy - spatial derivative in the direction perpendicular to velocity vector


shear_stress = Symbol("shear_stress", units.pressure)
dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)
fluid_speed = Function("fluid_speed", units.velocity)
layer_position = Symbol("layer_position", units.length)

law = Eq(shear_stress, dynamic_viscosity * Derivative(fluid_speed(layer_position), layer_position))


def print_law() -> str:
    return print_expression(law)


@validate_input(
    dynamic_viscosity_=dynamic_viscosity,
    fluid_speed_before_=fluid_speed,
    fluid_speed_after_=fluid_speed,
    layer_separation_=layer_position,
)
@validate_output(shear_stress)
def calculate_inner_friction(
    dynamic_viscosity_: Quantity,
    fluid_speed_before_: Quantity,
    fluid_speed_after_: Quantity,
    layer_separation_: Quantity,
) -> Quantity:
    fluid_speed_function = layer_position * (fluid_speed_after_ - fluid_speed_before_) / layer_separation_
    applied_law = law.subs({
        dynamic_viscosity: dynamic_viscosity_,
        fluid_speed(layer_position): fluid_speed_function,
    })
    dsolved = applied_law.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
