"""
Surface resistance of metal
===========================

A coaxial waveguide is an electrical cable consisting of a central conductor and a
shield arranged coaxially and separated by an insulating material or an air gap. It is
used to transmit radio frequency electrical signals. The resistance formed by the
surface of the metal sheath of the cable can be calculated by knowing the signal
frequency, magnetic permeability and specific conductivity of the metal.

**Notation:**

#. :quantity_notation:`vacuum_permeability`.

..
    TODO: replace `mu_0 * mu_r` with `mu`
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    units,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    quantities,
    symbols,
)

resistance = symbols.electrical_resistance
"""
Surface :symbols:`electrical_resistance`.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the insulator.
"""

angular_frequency = symbols.angular_frequency
"""
:symbols:`angular_frequency` of the signal.
"""

specific_conductivity = SymbolNew("sigma", units.conductance / units.length, display_latex="\\sigma")
"""
Specific conductivity of the conductor, i.e. :symbols:`electrical_conductivity` per unit
:symbols:`length`.
"""

law = Eq(
    resistance,
    sqrt(angular_frequency * quantities.vacuum_permeability * relative_permeability /
    (2 * specific_conductivity)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permeability_=relative_permeability,
    angular_frequency_=angular_frequency,
    specific_conductivity_=specific_conductivity)
@validate_output(resistance)
def calculate_surface_resistance(relative_permeability_: float, angular_frequency_: Quantity,
    specific_conductivity_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, resistance, dict=True)[0][resistance]
    result_expr = result_velocity_expr.subs({
        relative_permeability: relative_permeability_,
        angular_frequency: angular_frequency_,
        specific_conductivity: specific_conductivity_
    })
    return Quantity(result_expr)
