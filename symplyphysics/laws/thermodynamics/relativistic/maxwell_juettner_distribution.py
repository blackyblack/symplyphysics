from sympy import Eq, exp, sqrt
from sympy.functions.special.bessel import besselk
from symplyphysics import (
    dimensionless,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## The Maxwell-Juettner distribution is the distribution of speeds of particles in a hypothetical
## gas of relativistic particles. It is similar to Maxwell-Boltzmann distribution in that it consideres
## an ideal gas where particles are dilute and do not interact with each other, but the effects of special
## relativity are taken into account. In the limit of low temperatures, this distribution becomes identical
## to the Maxwell-Boltzmann distribution.

# Law: f(gamma) = (gamma * sqrt(gamma**2 - 1) / (theta * K_2(1 / theta)) * exp(-gamma/theta)
## f(gamma) - probability distribution function
## gamma - [Lorentz factor](../../../definitions/lorentz_factor.py)
## theta = (k * T)/(m * c**2) - dimensionless temperature
## K_2 - modified Bessel function of the second kind

# Conditions
## - The system is in thermal equilibrium with the environment.

# Limitations
## - Particle interactions are not taken into account
## - No quantum effects in the system
## - Antiparticles cannot occur in the system
## - Temperature must be isotropic (each degree of freedom has to have the same translational kinetic energy)

distribution_function = Function("distribution_function", dimensionless)
lorentz_factor = Symbol("lorentz_factor", dimensionless)
dimensionless_temperature = Symbol("dimensionless_temperature", dimensionless)

law = Eq(
    distribution_function(lorentz_factor),
    lorentz_factor
    * sqrt(lorentz_factor**2 - 1)
    / (dimensionless_temperature * besselk(2, 1 / dimensionless_temperature))
    * exp(-1 * lorentz_factor / dimensionless_temperature)
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    lorentz_factor_=lorentz_factor,
    dimensionless_temperature_=dimensionless_temperature,
)
@validate_output(distribution_function)
def calculate_distribution_function(
    lorentz_factor_: float,
    dimensionless_temperature_: float,
) -> float:
    result = law.rhs.subs({
        lorentz_factor: lorentz_factor_,
        dimensionless_temperature: dimensionless_temperature_,
    })
    return float(result)
