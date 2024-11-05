"""
Diffusion flux from diffusion coefficient and concentration gradient
====================================================================

*Fick's first law* relates the diffusion flux to the gradient of the concentration.
It postulates that the flux goes from regions of high concentration to regions of
low concentration, with a magnitude that is proportional to the concentration gradient
(spatial derivative), or, in simplistic terms, the concept that a solute will move from
a region of high concentration to a region of low concentration across a concentration
gradient. The minus sign in the law indicates that the directions of the diffusion
flow and the gradient are opposite: the gas diffuses towards a lower concentration,
and the gradient is directed towards a higher concentration.

**Conditions:**

#. The mixture is ideal.
"""

from sympy import Eq, Derivative
from symplyphysics import (
    Quantity,
    validate_input,
    validate_output,
    symbols,
    clone_as_function,
)

position = symbols.position
"""
:symbols:`position` of particles.
"""

diffusion_flux = clone_as_function(symbols.diffusion_flux, [position])
"""
:symbols:`diffusion_flux` as a function of :attr:`~position`.
"""

diffusion_coefficient = symbols.diffusion_coefficient
"""
:symbols:`diffusion_coefficient`.
"""

concentration = clone_as_function(symbols.molar_concentration, [position])
"""
:symbols:`molar_concentration` of particles as a function of :attr:`~position`.
"""

law = Eq(diffusion_flux(position),
    -diffusion_coefficient * Derivative(concentration(position), position))
"""
:laws:symbol::

:laws:latex::
"""


@validate_input(diffusion_coefficient_=diffusion_coefficient,
    concentration_start_=concentration,
    concentration_end_=concentration,
    position_=position)
@validate_output(diffusion_flux)
def calculate_diffusion_flux(diffusion_coefficient_: Quantity, concentration_start_: Quantity,
    concentration_end_: Quantity, position_: Quantity) -> Quantity:
    concentration_function_ = position * (concentration_end_ - concentration_start_) / position_
    applied_definition = law.subs({
        concentration(position): concentration_function_,
        diffusion_coefficient: diffusion_coefficient_
    })
    dsolved = applied_definition.doit()
    result_expr = dsolved.rhs
    return Quantity(result_expr)
