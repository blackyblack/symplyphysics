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

## Law is: F = (k * T / q) * ln(Nd * Na / n), where
## F - height of the potential transition (barrier),
## n is concentration of intrinsic charge carriers,
## q - charge of electron or hole,
## Nd - concentrations of donors,
## Nd - concentrations of acceptors,
## k - Boltzmann constant,
## T - temperature.

height_barrier = Symbol("height_barrier", units.voltage)

donors_concentration = Symbol("donors_concentration", 1 / units.volume)
acceptors_concentration = Symbol("acceptors_concentration", 1 / units.volume)
charge_carriers_concentration = Symbol("charge_carriers_concentration", 1 / units.volume)
temperature = Symbol("temperature", units.temperature)

charge_electron = Quantity(1.6e-19 * units.coulomb)

law = Eq(height_barrier, (boltzmann * temperature / charge_electron) * log(donors_concentration * acceptors_concentration / charge_carriers_concentration**2))


def print_law() -> str:
    return print_expression(law)


@validate_input(donors_concentration_=donors_concentration, acceptors_concentration_=acceptors_concentration, charge_carriers_concentration_=charge_carriers_concentration, temperature_=temperature)
@validate_output(height_barrier)
def calculate_height_barrier(donors_concentration_: Quantity, acceptors_concentration_: Quantity, charge_carriers_concentration_: Quantity, temperature_: Quantity) -> Quantity:
    result_expr = solve(law, height_barrier, dict=True)[0][height_barrier]
    result_expr = result_expr.subs({
        donors_concentration: donors_concentration_,
        acceptors_concentration: acceptors_concentration_,
        charge_carriers_concentration: charge_carriers_concentration_,
        temperature: temperature_,
    })
    return Quantity(result_expr)
