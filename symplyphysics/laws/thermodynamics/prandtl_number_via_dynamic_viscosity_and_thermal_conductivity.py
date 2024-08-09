r"""
Prandtl number via dynamic viscosity and thermal conductivity
=============================================================

*Prandtl number* is a dimensionless quantity defined as the ratio of kinetic viscosity
(momentum diffusivity) to thermal diffusivity. It can also be expressed using dynamic viscosity
and thermal conductivity.
"""

from sympy import Eq, solve
from symplyphysics import (
    Quantity,
    Symbol,
    dimensionless,
    units,
    validate_input,
    validate_output,
    convert_to_float,
)

# TODO Create "definition" law: Pr = momentum diffusivity / thermal diffusivity

isobaric_specific_heat_capacity = Symbol("isobaric_specific_heat_capacity", units.energy / units.mass / units.temperature)
r"""
Heat capacity at constant pressure per unit mass.

Symbol:
    :code:`c_p`

Latex:
    :math:`c_p`
"""

dynamic_viscosity = Symbol("dynamic_viscosity", units.pressure * units.time)
r"""
Dynamic viscosity of the fluid.

Symbol:
    :code:`mu`

Latex:
    :math:`\mu`
"""

thermal_conductivity = Symbol("thermal_conductivity",
    units.power / units.length / units.temperature)
"""
Thermal conductivity of the fluid.

Symbol:
    :code:`k`
"""

prandtl_number = Symbol("prandtl_number", dimensionless)
r"""
Prandtl number of the fluid.

Symbol:
    :code:`Pr`

Latex:
    :math:`\text{Pr}`
"""

law = Eq(prandtl_number, isobaric_specific_heat_capacity * dynamic_viscosity / thermal_conductivity)
r"""
:code:`Pr = c_p * mu / k`

Latex:
    .. math::
        \text{Pr} = \frac{c_p \mu}{k}
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
    result = Quantity(result_applied)
    return convert_to_float(result)
