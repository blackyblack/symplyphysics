"""
Buoyant force from density and volume
=====================================

Any object, totally or partially immersed in a fluid (i.e. liquid or gas), is buoyed up by a force equal to the
weight of the fluid displaced by the object. Also known as the Archimedes principle. The *buoyant force*
vector is directed opposite to the gravity vector.
"""

from sympy import (Eq, solve)
from symplyphysics import (clone_symbol, symbols, units, Quantity, Symbol,
    validate_input, validate_output)

force_buoyant = clone_symbol(symbols.dynamics.force, "force_buoyant")
r"""
The buoyant (Archimedes) :attr:`~symplyphysics.symbols.dynamics.force`.

Symbol:
    Fa

Latex:
    :math:`F_\text{A}`
"""

fluid_density = Symbol("fluid_density", units.mass / units.volume)
r"""
The density of the fluid.

Symbol:
    rho

Latex:
    :math:`\rho`
"""

displaced_volume = Symbol("displaced_volume", units.volume)
"""
The volume of the displaced fluid. Equivalently, the volume of the part of the body immersed in the fluid.

Symbols:
    V
"""

law = Eq(force_buoyant, -1 * fluid_density * units.acceleration_due_to_gravity * displaced_volume)
r"""
Fa = -1 * rho * g * V

Latex:
    .. math::
        F_\text{A} = - \rho g V
"""

@validate_input(fluid_density_=fluid_density, displaced_volume_=displaced_volume)
@validate_output(force_buoyant)
def calculate_force_buoyant(fluid_density_: Quantity, displaced_volume_: Quantity) -> Quantity:
    result_force_expr = solve(law, force_buoyant, dict=True)[0][force_buoyant]
    result_expr = result_force_expr.subs({
        fluid_density: fluid_density_,
        displaced_volume: displaced_volume_
    })
    return Quantity(abs(result_expr))
