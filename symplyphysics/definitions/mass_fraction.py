from sympy import (Eq, solve, S)
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output, dimensionless, convert_to)
from symplyphysics.core.symbols.fraction import Fraction

# Description
## The mass fraction is the ratio of the mass of the mixture to the total mass of the mixture
## omega = m_i / m
## Where:
## m_i - mass of component
## m - mass of mixture
## omega - mass fraction of the mixture

mass_fraction = Symbol("mass_fraction", dimensionless)
mass_of_component = Symbol("mass_of_component", units.mass)
mass_of_mixture = Symbol("mass_of_mixture", units.mass)

definition = Eq(mass_fraction, mass_of_component / mass_of_mixture)

definition_units_SI = dimensionless


def print_law() -> str:
    return print_expression(definition)


@validate_input(mass_of_component_=mass_of_component, mass_of_mixture_=mass_of_mixture)
@validate_output(mass_fraction)
def calculate_mass_fraction(mass_of_component_: Quantity,
    mass_of_mixture_: Quantity) -> Fraction:
    result_mass_fraction_expr = solve(definition, mass_fraction, dict=True)[0][mass_fraction]
    result_expr = result_mass_fraction_expr.subs({
        mass_of_component: mass_of_component_,
        mass_of_mixture: mass_of_mixture_
    })
    result = Quantity(result_expr)
    return Fraction(convert_to(result, S.One).evalf())
