from sympy import Eq, Derivative
from symplyphysics import units, Quantity, Symbol, Function, print_expression, validate_input, validate_output

# Description
## In a steady (laminar) flow in a tube, the speed of fluid particles varies from its maximum value
## (along the tube axis) to zero (near tube walls). Each layer decelerates the layer closer to the axis,
## and accelerates the layer further from the axis. In between the fluid layers there exist forces of
## inner friction.

# Law: F = eta * dv/dy * S
## F - forces of inner friction between fluid layers
## eta - dynamic viscosity
## v - speed of fluid flow
## d/dy - spatial derivative in the direction perpendicular to velocity vector (i.e., a gradient)
## S - area of fluid layers

# Also called the Newton's law of viscosity

# If a fluid follows this law, it is called a Newtonian one. If not, it's called non-Newtonian.
# In such fluids, the viscosity is a function of pressure and the velocity gradient.

inner_friction = Symbol("inner_friction", units.force)
dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)
fluid_speed = Function("fluid_speed", units.length / units.time)
layer_position = Symbol("layer_position", units.length)
layer_area = Symbol("layer_area", units.area)

law = Eq(inner_friction, dynamic_viscosity * Derivative(fluid_speed(layer_position), layer_position) * layer_area)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    dynamic_viscosity_=dynamic_viscosity,
    fluid_speed_before_=fluid_speed,
    fluid_speed_after_=fluid_speed,
    layer_separation_=layer_position,
    layer_area_=layer_area,
)
@validate_output(inner_friction)
def calculate_inner_friction(
    dynamic_viscosity_: Quantity,
    fluid_speed_before_: Quantity,
    fluid_speed_after_: Quantity,
    layer_separation_: Quantity,
    layer_area_: Quantity,
) -> Quantity:
    fluid_speed_function = layer_position * (fluid_speed_after_ - fluid_speed_before_) / layer_separation_
    applied_law = law.subs({
        dynamic_viscosity: dynamic_viscosity_,
        fluid_speed(layer_position): fluid_speed_function,
        layer_area: layer_area_,
    })
    dsolved = applied_law.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
