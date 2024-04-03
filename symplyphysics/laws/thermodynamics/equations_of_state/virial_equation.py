from sympy import Eq, S
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output,
    convert_to,
)

# Description
## Also called the virial expansion, the virial equation of state expresses the compressibility factor
## (and therefore the pressure) of a real gas in local equilibrium as a power series of the density.

# Law: Z = 1 + B * rho + C * rho**2 + ...
## Z - compressibility factor; TODO: add link to definition
## B, C, ... - virial coefficients
## rho - molar density (amount of substance divided by volume)

# Note
## - The first virial coefficient A is defined as the constant of value 1 in order to enforce that the equation
##   reduces to the ideal gas equation as the gas density approaches zero.
## - N-th virial coefficient represents the non-additive n-body interactions of particles plus all
##   mutual interactions of 2 up to `n - 1` particles.
## - 4- and more body interactions are quite rare to happen, so the expansion is truncated to contain only
##   the second and third virial coefficients. Moreover, the latter have been extensively studied and tabulated
##   for many fluids.

compressibility_factor = Symbol("compressibility_factor", dimensionless, positive=True)
second_virial_coefficient = Symbol("second_virial_coefficient", units.volume / units.amount_of_substance, real=True)
third_virial_coefficient = Symbol("third_virial_coefficient", (units.volume / units.amount_of_substance)**2, real=True)
molar_density = Symbol("molar_density", units.amount_of_substance / units.volume, positive=True)

law = Eq(
    compressibility_factor,
    1 + second_virial_coefficient * molar_density + third_virial_coefficient * molar_density**2
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    second_virial_coefficient_=second_virial_coefficient,
    third_virial_coefficient_=third_virial_coefficient,
    molar_density_=molar_density,
)
@validate_output(compressibility_factor)
def calculate_compressibility_factor(
    second_virial_coefficient_: Quantity,
    third_virial_coefficient_: Quantity,
    molar_density_: Quantity,
) -> float:
    result = law.rhs.subs({
        second_virial_coefficient: second_virial_coefficient_,
        third_virial_coefficient: third_virial_coefficient_,
        molar_density: molar_density_,
    })
    return float(convert_to(Quantity(result), S.One))
