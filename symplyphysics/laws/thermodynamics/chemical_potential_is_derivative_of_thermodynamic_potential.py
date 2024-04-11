from sympy import Eq, Symbol as SymSymbol, Derivative
from symplyphysics import (
    units,
    dimensionless,
    Quantity,
    Symbol,
    Function,
    print_expression,
    validate_input,
    validate_output,
)

# Description
## Chemical potential can be described as the amount of energy the system absorbs or releases
## due to the change of particle number. Thermodynamically it can be found by differentiating
## thermodynamic potentials with respect to the particle number keeping the other two natural
## variables of the potential constant.

# Law: mu = (dX/dN)_(a, b)
## mu - chemical potential
## X - thermodynamic potential
## N - particle count
## a, b - other two natural variables of X

# Note
## - mu = (dU/dN)_(S, V) = (dF/dN)_(T, V) = (dH/dN)_(S, p) = (dG/dN)_(T, p) = (dΩ/dN)_(T, V)
##   where U is internal energy, F is Helmholtz free energy, H is enthalpy, G is Gibbs energy,
##   Ω is grand potential.

chemical_potential = Symbol("chemical_potential", units.energy)
thermodynamic_potential = Function("thermodynamic_potential", units.energy)
particle_count = Symbol("particle_count", dimensionless)
first_natural_variable = SymSymbol("first_natural_variable")
second_natural_variable = SymSymbol("second_natural_variable")

law = Eq(
    chemical_potential,
    Derivative(thermodynamic_potential(first_natural_variable, second_natural_variable, particle_count), particle_count)
)


def print_law() -> str:
    return print_expression(law)


@validate_input(
    thermodynamic_potential_before_=thermodynamic_potential,
    thermodynamic_potential_after_=thermodynamic_potential,
    particle_count_before_=particle_count,
    particle_count_after_=particle_count,
)
@validate_output(chemical_potential)
def calculate_chemical_potential(
    thermodynamic_potential_before_: Quantity,
    thermodynamic_potential_after_: Quantity,
    particle_count_before_: int,
    particle_count_after_: int,
) -> Quantity:
    thermodynamic_potential_ = (
        thermodynamic_potential_before_ +
        (thermodynamic_potential_after_ - thermodynamic_potential_before_)
        * (particle_count - particle_count_before_)
        / (particle_count_after_ - particle_count_before_)
    )
    result = law.rhs.subs(
        thermodynamic_potential(first_natural_variable, second_natural_variable, particle_count),
        thermodynamic_potential_,
    ).doit()
    return Quantity(result)
