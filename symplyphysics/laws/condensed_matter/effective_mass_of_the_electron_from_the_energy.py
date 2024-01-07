from sympy import (Derivative, Expr, Eq, solve, pi,)
from sympy.physics.units import planck as planck_constant
from symplyphysics import (SI, units, Quantity, Symbol, Function, print_expression, validate_input,
    validate_output,)

# Description
## In solid state physics, a particle's effective mass is the mass that it seems to have when responding to forces,
## or the mass that it seems to have when interacting with other identical particles in a thermal distribution.
## The effective mass is a quantity that is used to simplify band structures by modeling the behavior
## of a free particle with that mass.

# Law: m_effective = (h/2pi)^2 / (d^2(E)/dk^2)
## Where:
## E - energy function of the electron,
## h is Planck constant,
## k is number of wave cycles per length,
## (d^2(E)/dk^2) - derivative of the second order of energy by the component of the wave vector.
## m_effective is effective mass of the electron.

wavenumber = Symbol("wavenumber", 1/units.length)
energy_function = Function("energy_function", units.energy)
mass = Symbol("mass", units.mass)

law = Eq(mass, ((planck_constant/2/pi)**2)
         / Derivative(energy_function(wavenumber), (wavenumber, 2)))


def print_law() -> str:
    return print_expression(law)


def apply_energy_function(energy_function_: Expr) -> Expr:
    applied_law = law.subs(energy_function(wavenumber), energy_function_)
    return applied_law


@validate_input(wavenumber_=wavenumber)
@validate_output(mass)
def calculate_mass(energy_function_: Expr, wavenumber_: Quantity) -> Quantity:
    energy_function_quantity = Quantity(energy_function_)
    assert SI.get_dimension_system().equivalent_dims(energy_function_quantity.dimension,
        energy_function.dimension)
    applied_law = apply_energy_function(energy_function_)
    result_expr = applied_law.subs(wavenumber, wavenumber_)
    result = solve(result_expr, mass, dict=True)[0][mass]
    return Quantity(result)
