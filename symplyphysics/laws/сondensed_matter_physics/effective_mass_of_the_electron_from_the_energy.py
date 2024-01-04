from sympy import (Derivative, Expr, Eq, solve, pi,)
from sympy.physics.units import planck as planck_constant
from symplyphysics import (SI, units, Quantity, Symbol, Function, dimensionless, print_expression, validate_input,
    validate_output,)

# Description
## In solid state physics, a particle's effective mass is the mass that it seems to have when responding to forces,
## or the mass that it seems to have when interacting with other identical particles in a thermal distribution.
## The effective mass is a quantity that is used to simplify band structures by modeling the behavior
## of a free particle with that mass.

# Law: (2pi/h)*d^2(E)/dk^2 = 1/m_effective
## Where:
## E - energy of the electron,
## h is Planck constant,
## k is component of the propagation vector,
## m_effective is effective mass of the electron.

propogation_vec_axis = Symbol("propogation_vec_axis", dimensionless/units.length)

energy = Function("energy", units.energy)
mass = Symbol("mass", units.mass)

law = Eq(mass, ((planck_constant/2/pi)**2)
         / Derivative(Derivative(energy(propogation_vec_axis), propogation_vec_axis), propogation_vec_axis))

def print_law() -> str:
    return print_expression(law)


def apply_energy_func(energy_func_: Expr) -> Expr:
    applied_law = law.subs(energy(propogation_vec_axis), energy_func_)
    return applied_law


@validate_input(propogation_vec_axis_=propogation_vec_axis)
@validate_output(mass)
def calculate_mass(energy_func_: Expr, propogation_vec_axis_: Quantity) -> Quantity:
    energy_func_quantity = Quantity(energy_func_)
    assert SI.get_dimension_system().equivalent_dims(energy_func_quantity.dimension,
        energy.dimension)

    applied_law = apply_energy_func(energy_func_)

    result_expr = applied_law.subs({
        propogation_vec_axis: propogation_vec_axis_
    })
    result = solve(result_expr, mass,
        dict=True)[0][mass]

    return Quantity(result)
