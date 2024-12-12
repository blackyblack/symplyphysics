r"""
Shear stress is proportional to speed gradient
==============================================

An element of flowing fluid will endure forces from the surrounding fluid (stress forces) that
will gradually deform the element over time. In a steady (laminar) flow, these stress forces
occur between fluid layers. For an incompressible and isotropic fluid, the shear stress exerted
to an element of fluid is proportional to the spatial derivative of the velocity component
that is perpendicular to the velocity vector (and as such, parallel to the direction of shear).
This law is also called the *Newton's law of viscosity*, and the fluids that follow it are said to
be *Newtonian*.

As an example of this law, consider two solid flat plates that contain water in between.
The bottom plate is fixed in place, while the top plate moves parallel to the bottom one with
a small speed :code:`u` such that the flow of the water is steady. If we were to measure the force that
needed to make the top plate continue to move, we would find that it is proportional to the area
of the plate and the ratio :math:`\frac{u}{d}` where :math:`d` is the distance between the plates.

**Conditions:**

#. The flow is one-dimensional. For a two-dimensional flow, replace the derivative with a sum of partial
   derivatives with respect to both perpendicular directions.
#. The fluid is incompressible and isotropic.
#. The fluid flow is steady (laminar).

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Viscosity#Dynamic_viscosity>`__.
"""

from sympy import Eq, Derivative, solve
from symplyphysics import units, Quantity, Symbol, Function, validate_input, validate_output

shear_stress = Symbol("shear_stress", units.pressure)
r"""
Shear stress. See :doc:`laws.dynamics.pressure_from_force_and_area`.

Symbol:
    :code:`tau`

Latex:
    :math:`\tau`
"""

dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)
r"""
Dynamic viscosity of the fluid.

Symbol:
    :code:`eta`

Latex:
    :math:`\eta`
"""

fluid_speed = Function("fluid_speed", units.velocity)
"""
Fluid speed as a function of position perpendicular to the fluid flow, or layer position.

Symbol:
    :code:`u(y)`
"""

layer_position = Symbol("layer_position", units.length)
"""
Layer position, or position in the direction perpendicular to fluid velocity.

Symbol:
    :code:`y`
"""

law = Eq(shear_stress, dynamic_viscosity * Derivative(fluid_speed(layer_position), layer_position))
r"""
:code:`tau = eta * Derivative(u(y), y)`

Latex:
    .. math::
        \tau = \eta \frac{d u}{d y}
"""


@validate_input(
    dynamic_viscosity_=dynamic_viscosity,
    fluid_speed_before_=fluid_speed,
    fluid_speed_after_=fluid_speed,
    layer_separation_=layer_position,
)
@validate_output(shear_stress)
def calculate_shear_stress(
    dynamic_viscosity_: Quantity,
    fluid_speed_before_: Quantity,
    fluid_speed_after_: Quantity,
    layer_separation_: Quantity,
) -> Quantity:
    fluid_speed_function = layer_position * (fluid_speed_after_ -
        fluid_speed_before_) / layer_separation_
    applied_law = law.subs({
        dynamic_viscosity: dynamic_viscosity_,
        fluid_speed(layer_position): fluid_speed_function,
    })
    result_expr = solve(applied_law, shear_stress)[0]
    return Quantity(result_expr)
