"""
Buoyant force from density and volume
=====================================

Any object, totally or partially immersed in a fluid (i.e. liquid or gas), is buoyed up by a force equal to the
weight of the fluid displaced by the object. Also known as the Archimedes principle. The *buoyant force*
vector is directed opposite to the gravity vector.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_as_symbol, symbols, quantities, Quantity, validate_input,
    validate_output, mul)

buoyant_force = clone_as_symbol(symbols.force,
    display_symbol="F_A",
    display_latex="F_\\text{A}")
"""
The buoyant (Archimedes) :symbols:`force`.
"""

fluid_density = symbols.density
"""
The :symbols:`density` of the fluid.
"""

displaced_volume = symbols.volume
"""
The :symbols:`volume` of the displaced fluid. Equivalently, the volume of the part of the body
immersed in the fluid.
"""

law = Eq(buoyant_force, mul(-1, fluid_density, quantities.acceleration_due_to_gravity, displaced_volume))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(fluid_density_=fluid_density, displaced_volume_=displaced_volume)
@validate_output(buoyant_force)
def calculate_force_buoyant(fluid_density_: Quantity, displaced_volume_: Quantity) -> Quantity:
    result_force_expr = solve(law, buoyant_force, dict=True)[0][buoyant_force]
    result_expr = result_force_expr.subs({
        fluid_density: fluid_density_,
        displaced_volume: displaced_volume_
    })
    return Quantity(abs(result_expr))
