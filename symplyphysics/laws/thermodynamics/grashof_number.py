r"""
Grashof number
==============

The Grashof number (Gr) is a dimensionless number which approximates the ratio of the
buoyancy to viscous forces acting on a fluid. It arises in situations involving natural convection
and is analogous to the Reynolds number.

**Notation:**

#. :quantity_notation:`acceleration_due_to_gravity`.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Grashof_number#Definition>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    clone_as_symbol,
    symbols,
    Quantity,
    validate_input,
    validate_output,
    convert_to_float,
    quantities,
)

grashof_number = symbols.grashof_number
"""
:symbols:`grashof_number` of the fluid.
"""

volumetric_expansion_coefficient = clone_as_symbol(symbols.thermal_expansion_coefficient, subscript="V")
r"""
Volumetric (see :symbols:`volume`) :symbols:`thermal_expansion_coefficient` of the body.
"""

surface_temperature = clone_as_symbol(symbols.temperature,
    display_symbol="T_s",
    display_latex="T_\\text{s}")
"""
:symbols:`temperature` of the surface of the fluid.
"""

bulk_temperature = clone_as_symbol(symbols.temperature,
    display_symbol="T_b",
    display_latex="T_\\text{b}")
"""
Average :symbols:`temperature` of the inside of the fluid.
"""

characteristic_length = symbols.characteristic_length
"""
:symbols:`characteristic_length` of the fluid body.
"""

kinematic_viscosity = symbols.kinematic_viscosity
"""
:symbols:`kinematic_viscosity` of the fluid.
"""

law = Eq(
    grashof_number,
    quantities.acceleration_due_to_gravity * volumetric_expansion_coefficient *
    (surface_temperature - bulk_temperature) * characteristic_length**3 / (kinematic_viscosity**2))
"""
:laws:symbol::

:laws:latex::
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
