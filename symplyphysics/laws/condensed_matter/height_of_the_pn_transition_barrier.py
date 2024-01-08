from sympy import (Eq, solve, log)
from sympy.physics.units import boltzmann
from symplyphysics import (
    units,
    Quantity,
    Symbol,
    print_expression,
    validate_input,
    validate_output
)

# Description
## The p-n junction has a potential barrier preventing the movement of charge carriers. If the concentration
## of donor impurities is significantly higher than the concentration of acceptor impurities, then the law for
## the height of the barrier takes the form below.

## Law is: F = (kT/q) * ln(Nd * Na / n), where
## F - height of the potential transition,
## n is concentration intrinsic charge carriers,
## q - charge of the electron
## Nd, Na, n - concentrations of donors, acceptors and intrinsic charge carriers, respectively,
## k - Boltzmann constant,
## T - temperature.

height_barrier = Symbol("height_barrier", units.voltage)

concentration_donors = Symbol("concentration_donors", 1 / units.length**3)
concentration_acceptors = Symbol("concentration_acceptors", 1 / units.length**3)
concentration_intrinsic = Symbol("concentration_intrinsic", 1 / units.length**3)
temperature = Symbol("temperature", units.temperature)

charge_electron = Quantity(1.6e-19 * units.coulomb)

law = Eq(height_barrier, (boltzmann * temperature / charge_electron) * log(concentration_donors * concentration_acceptors / concentration_intrinsic**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(concentration_donors_=concentration_donors, concentration_acceptors_=concentration_acceptors, concentration_intrinsic_=concentration_intrinsic, temperature_=temperature)
@validate_output(height_barrier)
def calculate_height_barrier(concentration_donors_: Quantity, concentration_acceptors_: Quantity, concentration_intrinsic_: Quantity, temperature_: Quantity) -> Quantity:
    result_momentum_expr = solve(law, height_barrier, dict=True)[0][height_barrier]
    result_expr = result_momentum_expr.subs({
        concentration_donors: concentration_donors_,
        concentration_acceptors: concentration_acceptors_,
        concentration_intrinsic: concentration_intrinsic_,
        temperature: temperature_,
    })
    return Quantity(result_expr)
