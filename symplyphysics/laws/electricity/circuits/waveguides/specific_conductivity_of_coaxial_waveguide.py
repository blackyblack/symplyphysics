"""
Specific conductivity of coaxial waveguide
==========================================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The specific conductivity of a
coaxial waveguide depends on the frequency of signal and the specific capacitance of
coaxial waveguide, as well as on the tangent of the dielectric loss angle of the
insulator material.

..
    TODO: find link
"""

from sympy import Eq, solve
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
)

specific_conductance = SymbolNew("G", units.conductance / units.length)
"""
:symbols:`electrical_conductance` of the waveguide per unit :symbols:`length`.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the signal.
"""

specific_capacitance = SymbolNew("C", units.capacitance / units.length)
"""
:symbols:`capacitance` per unit :symbols:`length`.
"""

loss_tangent = symbols.dielectric_loss_tangent
"""
:symbols:`dielectric_loss_tangent`.
"""

law = Eq(specific_conductance,
    angular_frequency * specific_capacitance * loss_tangent)
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(angular_frequency_=angular_frequency,
    specific_capacitance_=specific_capacitance,
    tangent_dielectric_loss_angle_=loss_tangent)
@validate_output(specific_conductance)
def calculate_specific_conductivity(angular_frequency_: Quantity, specific_capacitance_: Quantity,
    tangent_dielectric_loss_angle_: float) -> Quantity:
    result_velocity_expr = solve(law, specific_conductance, dict=True)[0][specific_conductance]
    result_expr = result_velocity_expr.subs({
        angular_frequency: angular_frequency_,
        specific_capacitance: specific_capacitance_,
        loss_tangent: tangent_dielectric_loss_angle_
    })
    return Quantity(result_expr)
