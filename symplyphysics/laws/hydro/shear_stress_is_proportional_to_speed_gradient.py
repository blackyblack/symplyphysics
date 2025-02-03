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
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

shear_stress = symbols.shear_stress
"""
:symbols:`shear_stress`.
"""

dynamic_viscosity = symbols.dynamic_viscosity
"""
:symbols:`dynamic_viscosity` of the fluid.
"""

layer_position = symbols.position
"""
Layer :symbols:`position`, or position in the direction perpendicular to fluid velocity.
"""

flow_speed = clone_as_function(symbols.flow_speed, [layer_position])
"""
:symbols:`flow_speed` as a function of position perpendicular to the fluid flow, or
:attr:`~layer_position`.
"""

law = Eq(shear_stress, dynamic_viscosity * Derivative(flow_speed(layer_position), layer_position))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(
    dynamic_viscosity_=dynamic_viscosity,
    fluid_speed_before_=flow_speed,
    fluid_speed_after_=flow_speed,
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
        flow_speed(layer_position): fluid_speed_function,
    })
    result_expr = solve(applied_law, shear_stress)[0]
    return Quantity(result_expr)
