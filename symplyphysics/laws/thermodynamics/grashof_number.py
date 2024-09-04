r"""
Grashof number
==============

The Grashof number (Gr) is a dimensionless number which approximates the ratio of the
buoyancy to viscous forces acting on a fluid. It arises in situations involving natural convection
and is analogous to the Reynolds number.

**Notation:**

#. :math:`g` is the acceleration due to gravity.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_symbol,
    symbols,
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    dimensionless,
    convert_to_float,
)

grashof_number = Symbol("grashof_number", dimensionless)
r"""
Grashof number of the fluid.

Symbol:
    :code:`Gr`

Latex:
    :math:`\text{Gr}`
"""

volumetric_expansion_coefficient = Symbol("volumetric_expansion_coefficient", 1 / units.temperature)
r"""
Volumetric coefficient of thermal expansion of the body.

Symbol:
    :code:`beta`

Latex:
    :math:`\beta`
"""

surface_temperature = clone_symbol(symbols.thermodynamics.temperature,
    display_symbol="T_s",
    display_latex="T_\\text{s}")
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the surface of the fluid.
"""

bulk_temperature = clone_symbol(symbols.thermodynamics.temperature,
    display_symbol="T_b",
    display_latex="T_\\text{b}")
"""
Average :attr:`~symplyphysics.symbols.thermodynamics.temperature` of the inside of the fluid.
"""

characteristic_length = Symbol("characteristic_length", units.length)
"""
An important dimension of the given fluid, usually defined as the volume of the body
divided by its area.

Symbol:
    :code:`L`
"""

kinematic_viscosity = Symbol("kinematic_viscosity", units.area / units.time)
r"""
Kinematic viscosity of the fluid.

Symbol:
    :code:`nu`

Latex:
    :math:`\nu`
"""

law = Eq(
    grashof_number,
    units.acceleration_due_to_gravity * volumetric_expansion_coefficient *
    (surface_temperature - bulk_temperature) * characteristic_length**3 / (kinematic_viscosity**2))
r"""
:code:`Gr = g * beta * (T_s - T_b) * L^3 / nu^2`

Latex:
    .. math::
        \text{Gr} = g \beta (T_\text{s} - T_\text{b}) \frac{L^3}{\nu^2}
"""


@validate_input(coefficient_of_volume_expansion_=volumetric_expansion_coefficient,
    surface_temperature_=surface_temperature,
    fluid_temperature_=bulk_temperature,
    characteristic_length_=characteristic_length,
    viscosity_=kinematic_viscosity)
@validate_output(grashof_number)
def calculate_grashof_number(coefficient_of_volume_expansion_: Quantity,
    surface_temperature_: Quantity, fluid_temperature_: Quantity, characteristic_length_: Quantity,
    viscosity_: Quantity) -> float:
    result_expr = solve(law, grashof_number, dict=True)[0][grashof_number]
    result_applied = result_expr.subs({
        volumetric_expansion_coefficient: coefficient_of_volume_expansion_,
        surface_temperature: surface_temperature_,
        bulk_temperature: fluid_temperature_,
        characteristic_length: characteristic_length_,
        kinematic_viscosity: viscosity_
    })
    return convert_to_float(result_applied)
