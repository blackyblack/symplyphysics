"""
Group speed of wave in rectangular waveguide
============================================

A rectangular waveguide is a rectangular metal waveguide capable of supporting waves
propagating along it. There is a critical wavelength. Signals with a wavelength greater
than the critical one are attenuated and do not propagate in the waveguide.

**Notation:**

#. :quantity_notation:`speed_of_light`.

..
    TODO: find link
"""

from sympy import Eq, solve, sqrt
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    quantities,
    symbols,
    clone_as_symbol,
)

group_speed = symbols.group_speed
"""
:symbols:`group_speed` of the wave in the waveguide.
"""

relative_permittivity = symbols.relative_permittivity
"""
:symbols:`relative_permittivity` of the insulator.
"""

relative_permeability = symbols.relative_permeability
"""
:symbols:`relative_permeability` of the insulator.
"""

wavelength = symbols.wavelength
"""
:symbols:`wavelength` of the signal.
"""

critical_wavelength = clone_as_symbol(symbols.wavelength, subscript="\\text{c}")
"""
Critical :symbols:`wavelength` of the system.
"""

law = Eq(
    group_speed,
    quantities.speed_of_light * sqrt(1 - (wavelength / critical_wavelength)**2) /
    sqrt(relative_permittivity * relative_permeability))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability,
    wavelength_=wavelength,
    critical_wavelength_=critical_wavelength)
@validate_output(group_speed)
def calculate_group_velocity(relative_permittivity_: float, relative_permeability_: float,
    wavelength_: Quantity, critical_wavelength_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, group_speed, dict=True)[0][group_speed]
    result_expr = result_velocity_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        wavelength: wavelength_,
        critical_wavelength: critical_wavelength_
    })
    return Quantity(result_expr)
