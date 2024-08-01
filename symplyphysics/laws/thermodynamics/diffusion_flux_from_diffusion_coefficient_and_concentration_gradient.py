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

from sympy import (Eq, Derivative)
from symplyphysics import (units, Function, Quantity, Symbol, validate_input,
    validate_output)

diffusion_flux = Function("diffusion_flux", units.amount_of_substance / (units.area * units.time))
"""
Diffusion flux.

Symbol:
    :code:`J`
"""

diffusion_coefficient = Symbol("diffusion_coefficient", units.area / units.time)
"""
Diffusion coefficient.

Symbol:
    :code:`D`
"""

concentration = Function("concentration", units.amount_of_substance / units.volume)
"""
Concentration of particles as a function of position.

Symbol:
    :code:`n(x)`
"""

position = Symbol("position", units.length)
"""
Position of particles.

Symbol:
    :code:`x`
"""

law = Eq(diffusion_flux(position),
    -diffusion_coefficient * Derivative(concentration(position), position))
r"""
:code:`J = -1 * D * Derivative(n(x), x)`

Latex:
    .. math::
        J = - D \frac{d n}{d x}
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
