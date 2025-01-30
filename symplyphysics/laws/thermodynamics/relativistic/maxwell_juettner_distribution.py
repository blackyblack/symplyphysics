r"""
Maxwell—Jüttner distribution
============================

The Maxwell—Jüttner distribution is the distribution of speeds of particles in a hypothetical
gas of relativistic particles. It is similar to the Maxwell—Boltzmann distribution in that it considers 
an ideal gas where particles are dilute and do not interact with each other, but the effects of special
relativity are taken into account. In the limit of low temperatures, this distribution becomes identical
to the Maxwell—Boltzmann distribution.

**Notation:**

#. :math:`K_2` is the `modified Bessel function of the second kind <https://docs.sympy.org/latest/modules/functions/special.html#sympy.functions.special.bessel.besselk>`_.

**Conditions:**

#. The system is in thermal equilibrium with the environment.
#. Particle interactions are not taken into account.
#. No quantum effects occur in the system.
#. Antiparticles cannot occur in the system.
#. Temperature must be isotropic, i.e. each degree of freedom has to have the same translational kinetic energy.

**Links:**

#. `Wikipedia <https://en.wikipedia.org/wiki/Maxwell%E2%80%93J%C3%BCttner_distribution#Definition>`__.

..
    TODO: add `besselk` to code printers
"""

from sympy import Eq, exp, sqrt
from sympy.functions.special.bessel import besselk
from symplyphysics import (
    convert_to_float,
    dimensionless,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
)

distribution_function = SymbolNew("f(gamma)", dimensionless, display_latex="f(\\gamma)")
"""
:attr:`~lorentz_factor` distribution function.
"""

lorentz_factor = symbols.lorentz_factor
"""
:symbols:`lorentz_factor` of relativistic particles.
"""

reduced_temperature = SymbolNew("theta", dimensionless, display_latex="\\theta")
"""
:doc:`Reduced temperature <laws.thermodynamics.relativistic.reduced_temperature_in_maxwell_juettner_statistics>` of the system.
"""

law = Eq(
    distribution_function,
    lorentz_factor * sqrt(lorentz_factor**2 - 1) /
    (reduced_temperature * besselk(2, 1 / reduced_temperature)) *
    exp(-1 * lorentz_factor / reduced_temperature))
r"""
..
    The printers do not yet recognize the `besselk` function.

:code:`f(gamma) = (gamma * sqrt(gamma^2 - 1)) / (theta * K_2(1 / theta)) * exp(-1 * gamma / theta)`

Latex:
    .. math::
        f(\gamma) = \frac{\gamma \sqrt{\gamma^2 - 1}}
                         {\theta K_2 \left( \frac{1}{\theta} \right)} 
                    \exp{\left( -\frac{\gamma}{\theta} \right)}
"""


@validate_input(
    lorentz_factor_=lorentz_factor,
    dimensionless_temperature_=reduced_temperature,
)
@validate_output(distribution_function)
def calculate_distribution_function(
    lorentz_factor_: float,
    dimensionless_temperature_: float,
) -> float:
    result = law.rhs.subs({
        lorentz_factor: lorentz_factor_,
        reduced_temperature: dimensionless_temperature_,
    })
    return convert_to_float(result)
