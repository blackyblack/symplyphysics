from sympy import (Eq, solve)
from sympy.physics.units import planck as planck_constant
from symplyphysics import (units, Quantity, Symbol, print_expression, validate_input,
    validate_output)

# Description
## The photoelectric effect is a theory proposed by Einstein.
## It states that in the range of kinetic energies of the electrons removed
## from their respective atomic bindings by the absorption of a photon, the
## highest kinetic energy is equal to the difference between the energy of
## the photon and the work function.
## Law is: K_max = h * nu - W, where
## K_max is the highest kinetic energy of the electrons,
## h is Planck constant,
## nu is frequency of photon,
## W is work function of the surface.

max_kinetic_energy = Symbol("max_kinetic_energy", units.energy)
photon_frequency = Symbol("frequency", units.frequency)
work_function = Symbol("work_function", units.energy)

law = Eq(max_kinetic_energy, planck_constant * photon_frequency - work_function)


def print_law() -> str:
    return print_expression(law)


@validate_input(photon_frequency_=photon_frequency, work_function_=work_function)
@validate_output(max_kinetic_energy)
def calculate_max_kinetic_energy(photon_frequency_: Quantity, work_function_: Quantity) -> Quantity:
    result_energy_expr = solve(law, max_kinetic_energy, dict=True)[0][max_kinetic_energy]
    result_expr = result_energy_expr.subs({
        photon_frequency: photon_frequency_,
        work_function: work_function_
    })
    return Quantity(result_expr)
