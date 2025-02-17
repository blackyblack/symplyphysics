r"""
Prandtl number via dynamic viscosity and thermal conductivity
=============================================================

*Prandtl number* is a dimensionless quantity defined as the ratio of kinetic viscosity
(momentum diffusivity) to thermal diffusivity. It can also be expressed using dynamic viscosity
and thermal conductivity.

**Links:**

#. `Wikipedia, last formula within the box <https://en.wikipedia.org/wiki/Prandtl_number>`__.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    Symbol,
    units,
    validate_input,
    validate_output,
    convert_to_float,
    symbols,
)

# TODO Create "definition" law: Pr = momentum diffusivity / thermal diffusivity

prandtl_number = symbols.prandtl_number
"""
:symbols:`prandtl_number` of the fluid.
"""

isobaric_specific_heat_capacity = Symbol("c_p", units.energy / units.mass / units.temperature)
"""
:symbols:`heat_capacity` at constant :symbols:`pressure` per unit :symbols:`mass`.
"""

dynamic_viscosity = symbols.dynamic_viscosity
"""
:symbols:`dynamic_viscosity` of the fluid.
"""

thermal_conductivity = symbols.thermal_conductivity
"""
:symbols:`thermal_conductivity` of the fluid.
"""

law = Eq(prandtl_number, isobaric_specific_heat_capacity * dynamic_viscosity / thermal_conductivity)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(heat_capacity_=isobaric_specific_heat_capacity,
    dynamic_viscosity_=dynamic_viscosity,
    thermal_conductivity_=thermal_conductivity)
@validate_output(prandtl_number)
def calculate_prandtl_number(heat_capacity_: Quantity, dynamic_viscosity_: Quantity,
    thermal_conductivity_: Quantity) -> float:
    result_expr = solve(law, prandtl_number, dict=True)[0][prandtl_number]
    result_applied = result_expr.subs({
        isobaric_specific_heat_capacity: heat_capacity_,
        dynamic_viscosity: dynamic_viscosity_,
        thermal_conductivity: thermal_conductivity_
    })
    return convert_to_float(result_applied)
