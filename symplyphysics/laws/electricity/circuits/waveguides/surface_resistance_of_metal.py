"""
Surface resistance of metal
===========================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The resistance formed by the
surface of the metal sheath of the cable can be calculated by knowing the signal
frequency, magnetic permeability and specific conductivity of the metal.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    validate_input,
    validate_output,
    symbols,
)

resistance = symbols.electrical_resistance
"""
Surface :symbols:`electrical_resistance`.
"""

absolute_permeability = symbols.absolute_permeability
"""
:symbols:`absolute_permeability` of the insulator.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the signal.
"""

specific_conductance = Symbol("G", units.conductance / units.length)
"""
Specific conductance of the conductor, i.e. :symbols:`electrical_conductance` per unit
:symbols:`length`.
"""

law = Eq(
    resistance,
    sqrt(angular_frequency * absolute_permeability /
    (2 * specific_conductance)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(absolute_permeability_=absolute_permeability,
    angular_frequency_=angular_frequency,
    specific_conductivity_=specific_conductance)
@validate_output(resistance)
def calculate_surface_resistance(absolute_permeability_: Quantity, angular_frequency_: Quantity,
    specific_conductivity_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_velocity_expr.subs({
        absolute_permeability: absolute_permeability_,
        angular_frequency: angular_frequency_,
        specific_conductance: specific_conductivity_
    })
    return Quantity(result_expr)
