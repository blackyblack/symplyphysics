"""
Peak wavelength via temperature
===============================

The wavelength of a blackbody's radiation peak is inversely proportional to its
temperature. This law is known as the *Wien's displacement law*.
"""

from sympy import Eq, solve
from sympy.physics.units import speed_of_light, planck, boltzmann
from symplyphysics import (symbols, units, Quantity, Symbol, validate_input, validate_output)

peak_wavelength = Symbol("peak_wavelength", units.length)
r"""
Wavelength at which the peak of the spectral radiance of a blackbody occurs.

Symbol:
    :code:`lambda_peak`

Latex:
    :math:`\lambda_\text{peak}`
"""

temperature = symbols.thermodynamics.temperature
"""
:attr:`~symplyphysics.symbols.thermodynamics.temperature` of the blackbody.
"""

wien_displacement_constant = speed_of_light * planck / (boltzmann * 4.965114)
"""
A constant of proportionality.

Symbol:
    :code:`b`
"""

law = Eq(peak_wavelength, wien_displacement_constant / temperature)
r"""
:code:`lambda_peak = b / T`

Latex:
    .. math::
        \lambda_\text{peak} = \frac{b}{T}
"""


@validate_input(object_temperature_=temperature)
@validate_output(peak_wavelength)
def calculate_intensive_wavelength(object_temperature_: Quantity) -> Quantity:
    result_wavelength_expr = solve(law, peak_wavelength, dict=True)[0][peak_wavelength]
    result_expr = result_wavelength_expr.subs({temperature: object_temperature_})
    return Quantity(result_expr)
