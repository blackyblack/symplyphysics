"""
Phase speed of wave in rectangular waveguide
============================================

A rectangular waveguide is a rectangular metal waveguide capable of supporting waves
propagating along it. There is a critical wavelength. Signals with a wavelength greater
than the critical one are attenuated and do not propagate in the waveguide.

**Notation:**

#. :quantity_notation:`speed_of_light`.

..
    TODO: fix file name
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

## Description

## Law is: Vph = c / sqrt(er * mur * (1 - (L / L1)^2)), where
## Vph - phase velocity of wave in rectangular waveguide,
## c - speed of light,
## er - relative permittivity of insulating material,
## mur - relative permeability of the insulating material,
## L - wavelength,
## L1 - critical wavelength.

phase_speed = symbols.phase_speed
"""
:symbols:`phase_speed` of the wave in the waveguide.
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
    phase_speed, quantities.speed_of_light / sqrt(relative_permittivity * relative_permeability * (1 -
    (wavelength / critical_wavelength)**2)))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(relative_permittivity_=relative_permittivity,
    relative_permeability_=relative_permeability,
    wavelength_=wavelength,
    critical_wavelength_=critical_wavelength)
@validate_output(phase_speed)
def calculate_phase_velocity(relative_permittivity_: float, relative_permeability_: float,
    wavelength_: Quantity, critical_wavelength_: Quantity) -> Quantity:
    result_velocity_expr = solve(law, phase_speed, dict=True)[0][phase_speed]
    result_expr = result_velocity_expr.subs({
        relative_permittivity: relative_permittivity_,
        relative_permeability: relative_permeability_,
        wavelength: wavelength_,
        critical_wavelength: critical_wavelength_
    })
    return Quantity(result_expr)
